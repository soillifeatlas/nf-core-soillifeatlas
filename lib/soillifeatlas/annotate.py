"""Stage 4: merge lipid annotations from several sources into one table.

Every source is turned into the same small table (one row per candidate):

    feature_id, lipid_name, lipid_class, sum_composition, level, rank, score,
    source, evidence

and one written priority rule chooses, per feature, the name and the
confidence level that is reported:

* ``level`` is how much of the structure a source claims: ``species``
  (each chain resolved), ``sum`` (class plus total carbons and double bonds)
  or ``class``.
* ``rank`` is the kind of evidence, lower is stronger:
  1 = species with diagnostic fragments (MZmine lipid rules at species level,
  LipidSearch grade A/B, MS-DIAL with matched MS/MS); 2 = spectral library
  match; 3 = rule-based class or sum assignment; 4 = in-silico prediction
  (SIRIUS / CANOPUS), last because it confuses lipid isomers
  (Zulfiqar et al. 2023).

Confidence (the labels the Soil Life Atlas uses):

* **Gold**: species level from a rank-1 source, or the same species from two
  independent sources.
* **Silver**: sum composition (reported at sum level even if a weaker source
  named a species).
* **Bronze**: lipid class only.
* **Unidentified**: no lipid annotation.

Sources that come from another processing run (for example the LipidSearch,
GNPS or SIRIUS results of an earlier study) are matched to this run's
features by m/z and retention time, after removing a constant retention-time
offset estimated from the data.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

PROTON = 1.007276
ADDUCT_SHIFT = {
    "M+H": PROTON,
    "M+NH4": 18.033823,
    "M+Na": 22.989218,
    "M+K": 38.963158,
    "M+H-H2O": PROTON - 18.010565,
    "M+2H": PROTON,  # singly treated; rare for lipids
    "M-H": -PROTON,
    "M+HCOO": 44.998201,
    "M+CH3COO": 59.013851,
    "M-CH3": -15.022926,
    "M+Cl": 34.969402,
}

# Classes whose names carry a single chain: one chain = complete species.
SINGLE_CHAIN = {
    "LPC", "LPE", "LPG", "LPI", "LPS", "LPA", "MG", "FA", "CE", "SPH", "SPB",
    "CAR", "NAE", "LPC-O", "LPE-O", "LPC-P", "LPE-P", "MGMG", "MGDG-L", "LDGTS",
    "HexCer-L", "So", "Sa",
}
CLASS_ALIAS = {
    "LYSOPC": "LPC", "LYSO-PC": "LPC", "LPC": "LPC", "LYSOPE": "LPE", "LYSO-PE": "LPE",
    "TAG": "TG", "DAG": "DG", "MAG": "MG", "CHE": "CE", "CHOLESTERYLESTER": "CE",
    "SPHINGOMYELIN": "SM", "CERAMIDE": "Cer", "CER": "Cer", "HEXCER": "HexCer",
    "GLCCER": "HexCer", "GALCER": "HexCer",
}
LEVEL_ORDER = {"species": 0, "sum": 1, "class": 2}

# Two-chain phospholipid classes. A name with one chain of at most
# LYSO_MAX_CARBONS carbons ("PC 16:0", "PC(16:0/0:0)") is the lyso form.
DIACYL_PL = {"PC", "PE", "PG", "PI", "PS", "PA"}
LYSO_MAX_CARBONS = 24

# Class family used to judge agreement between sources: a SIRIUS
# "glycerophosphocholine" (PC) agrees with an LPC species.
CLASS_FAMILY = {f"L{c}": c for c in DIACYL_PL} | {
    "Hex1Cer": "HexCer", "Hex2Cer": "HexCer", "Hex3Cer": "HexCer",
    "LPC-O": "PC", "LPE-O": "PE", "LPC-P": "PC", "LPE-P": "PE",
}


def class_family(cls: str) -> str:
    return CLASS_FAMILY.get(cls, cls)

# CANOPUS / ClassyFire class names that map unambiguously to a lipid class.
CLASSYFIRE_MAP = {
    "glycerophosphocholines": "PC", "glycerophosphoethanolamines": "PE",
    "glycerophosphoglycerols": "PG", "glycerophosphoinositols": "PI",
    "glycerophosphoserines": "PS", "glycerophosphates": "PA",
    "lysophosphatidylcholines": "LPC", "lysophosphatidylethanolamines": "LPE",
    "triradylcglycerols": "TG", "triacylglycerols": "TG", "diradylglycerols": "DG",
    "diacylglycerols": "DG", "monoradylglycerols": "MG", "monoacylglycerols": "MG",
    "ceramides": "Cer", "phosphosphingolipids": "SM", "sphingomyelins": "SM",
    "glycosphingolipids": "HexCer", "cholesteryl esters": "CE", "sterol esters": "CE",
    "glycosylglycerols": "MGDG", "betaine lipids": "DGTS",
}


# ---------------------------------------------------------------------------
# Lipid names
# ---------------------------------------------------------------------------
@dataclass
class Lipid:
    lipid_class: str
    chains: list[str]
    modification: str = ""

    @property
    def level(self) -> str:
        if not self.chains:
            return "class"
        if len(self.chains) > 1 or self.lipid_class in SINGLE_CHAIN:
            return "species"
        return "sum"

    @property
    def sum_composition(self) -> str:
        if not self.chains:
            return ""
        carbons = dbs = 0
        prefix = ""
        extra = ""
        for i, ch in enumerate(self.chains):
            m = re.match(r"^([OPdtm]-?|[OP]-)?(\d+):(\d+)(.*)$", ch)
            if not m:
                return ""
            if i == 0 and m.group(1):
                prefix = m.group(1) if m.group(1).endswith("-") or m.group(1) in "dtm" else m.group(1) + "-"
            carbons += int(m.group(2))
            dbs += int(m.group(3))
            if m.group(4) and not extra:
                extra = m.group(4)
        return f"{prefix}{carbons}:{dbs}{extra}"

    def name(self, level: str | None = None) -> str:
        level = level or self.level
        mod = self.modification
        if level == "class" or not self.chains:
            return self.lipid_class
        if level == "sum":
            return f"{self.lipid_class} {self.sum_composition}{mod}"
        return f"{self.lipid_class} {'_'.join(self.chains)}{mod}"


_CHAIN = r"(?:[OoPp]-|[dtm])?\d{1,2}:\d{1,2}(?:;O\d?|;\d?O)?"
_NAME = re.compile(
    rf"^(?P<cls>[A-Za-z][A-Za-z0-9]*(?:-[A-Za-z][A-Za-z]?)?)\s*\(?\s*"
    rf"(?P<chains>{_CHAIN}(?:\s*[_/]\s*{_CHAIN})*)\s*\)?\s*(?P<mod>\+.*)?$"
)


def normalise_class(cls: str) -> str:
    key = cls.replace(" ", "").upper()
    return CLASS_ALIAS.get(key, cls)


def parse_lipid(text: str) -> Lipid | None:
    """Parse common lipid shorthand into class, chains and modification.

    Accepts LIPID MAPS style ("PC 16:0_18:1", "PC 34:1"), LipidSearch style
    ("TG(16:0_16:0_16:1)", "TG(20:4_14:0)+O:(s)"), library names with
    prefixes or suffixes ("Massbank:LQB00208 PE 36:2",
    "PC(16:0/18:1); [M+H]+ C42H83N1O8P1", "Spectral Match to Lyso-PC(16:0)
    from NIST14") and ether or sphingoid prefixes ("PE O-16:0_18:1",
    "Cer d18:1/15:0"). Returns None for anything that is not a lipid name.
    """
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return None
    s = str(text).strip()
    s = re.sub(r"^Spectral Match to\s+", "", s, flags=re.I)
    s = re.sub(r"\s+from\s+\S+$", "", s)
    s = re.sub(r";\s*\[M.*$", "", s)  # "; [M+H]+ C42H..."
    s = re.sub(r"^[A-Za-z]+:[A-Za-z]+\d*\s+", "", s)  # library prefix such as Massbank:LQB00208
    s = re.sub(r"\s+CollisionEnergy.*$", "", s)
    s = s.replace("(s)", "s").strip()
    m = _NAME.match(s)
    if not m:
        bare = re.match(r"^([A-Za-z][A-Za-z0-9]{0,8})$", s)
        if bare and normalise_class(bare.group(1)) in set(CLASS_ALIAS.values()) | SINGLE_CHAIN:
            return Lipid(normalise_class(bare.group(1)), [])
        return None
    cls = normalise_class(m.group("cls"))
    chains = [re.sub(r"^([op])-", lambda x: x.group(1).upper() + "-", c.strip())
              for c in re.split(r"\s*[_/]\s*", m.group("chains"))]
    chains = [c for c in chains if c and not re.fullmatch(r"0:0", c)]
    mod = (m.group("mod") or "").replace(":", "").strip()
    if cls in DIACYL_PL and len(chains) == 1:
        c = re.search(r"(\d+):", chains[0])
        if c and int(c.group(1)) <= LYSO_MAX_CARBONS:
            cls = "L" + cls
    return Lipid(cls, chains, mod)


def normalise_adduct(adduct: str) -> str:
    a = str(adduct).strip().replace(" ", "")
    a = re.sub(r"^\[(.*)\]\d*[+-]$", r"\1", a)
    return a


def ion_mz(neutral_mass: float, adduct: str) -> float | None:
    shift = ADDUCT_SHIFT.get(normalise_adduct(adduct))
    return None if shift is None else neutral_mass + shift


# ---------------------------------------------------------------------------
# Matching external results to this run's features
# ---------------------------------------------------------------------------
def match_by_mz_rt(features: pd.DataFrame, ext: pd.DataFrame, ppm: float = 5.0, rt_window: float = 0.2,
                   coarse_window: float = 1.0) -> pd.DataFrame:
    """Attach the nearest feature (``feature_id``) to each row of ``ext``.

    ``features`` has columns id, mz, rt; ``ext`` has mz and rt (minutes).
    A constant retention-time offset between the two runs is estimated from
    unambiguous coarse matches and removed before the final match.
    Adds ``feature_id``, ``match_ppm`` and ``match_drt``; unmatched rows get
    ``feature_id`` = NaN.
    """
    f = features[["id", "mz", "rt"]].dropna().sort_values("mz").reset_index(drop=True)
    fmz, frt, fid = f["mz"].to_numpy(), f["rt"].to_numpy(), f["id"].to_numpy()

    def nearest(mz, rt, window, offset):
        tol = mz * ppm * 1e-6
        lo, hi = np.searchsorted(fmz, mz - tol), np.searchsorted(fmz, mz + tol, side="right")
        if hi <= lo or not np.isfinite(rt):
            return None
        drt = frt[lo:hi] - (rt + offset)
        ok = np.abs(drt) <= window
        if not ok.any():
            return None
        dppm = (fmz[lo:hi] - mz) / mz * 1e6
        cost = np.where(ok, np.abs(dppm) / ppm + np.abs(drt) / window, np.inf)
        j = int(np.argmin(cost))
        return lo + j, float(dppm[j]), float(drt[j])

    offsets = []
    for mz, rt in zip(ext["mz"], ext["rt"]):
        hit = nearest(mz, rt, coarse_window, 0.0)
        if hit:
            offsets.append(hit[2])
    offset = float(np.median(offsets)) if len(offsets) >= 10 else 0.0

    out = ext.copy()
    ids, ppms, drts = [], [], []
    for mz, rt in zip(out["mz"], out["rt"]):
        hit = nearest(mz, rt, rt_window, offset)
        if hit:
            ids.append(fid[hit[0]])
            ppms.append(hit[1])
            drts.append(hit[2])
        else:
            ids.append(np.nan)
            ppms.append(np.nan)
            drts.append(np.nan)
    out["feature_id"] = ids
    out["match_ppm"] = ppms
    out["match_drt"] = drts
    out.attrs["rt_offset"] = offset
    return out


# ---------------------------------------------------------------------------
# Source adapters: each returns the common candidate table
# ---------------------------------------------------------------------------
COLUMNS = ["feature_id", "lipid_name", "lipid_class", "sum_composition", "level", "rank", "score",
           "source", "evidence"]


def _candidates(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=COLUMNS)


def _row(fid, lip: Lipid, level: str, rank: int, score: float, source: str, evidence: str) -> dict:
    return {"feature_id": int(fid), "lipid_name": lip.name(level), "lipid_class": lip.lipid_class,
            "sum_composition": lip.sum_composition, "level": level, "rank": rank,
            "score": score, "source": source, "evidence": evidence}


def from_mzmine(full: pd.DataFrame) -> pd.DataFrame:
    """MZmine lipid annotation module (rule-based, MS/MS fragments).

    Species-level names are rank 1 (chains resolved from fragments); sum or
    class names are rank 3 (rule-based assignment).
    """
    col = "lipid_annotations:lipid_annotations"
    if col not in full.columns:
        return _candidates([])
    rows = []
    sub = full[full[col].notna()]
    for _, r in sub.iterrows():
        for name in str(r[col]).split(";"):
            lip = parse_lipid(name)
            if not lip:
                continue
            lvl = lip.level
            rank = 1 if lvl == "species" and len(lip.chains) > 1 else 3
            if lip.lipid_class in SINGLE_CHAIN and lvl == "species":
                rank = 1
            score = r.get("lipid_annotations:explained_intensity_percent", np.nan)
            ev = f"{r.get('lipid_annotations:ion_adduct', '')}; explained MS/MS intensity {score}"
            rows.append(_row(r["id"], lip, lvl, rank, float(score) if pd.notna(score) else np.nan,
                             "MZmine lipid rules", ev))
    return _candidates(rows)


def read_table(path: str | Path, polarity: str, sheet_hint: tuple[str, ...] = ()) -> pd.DataFrame:
    """Read CSV/TSV, or the sheet of an Excel file that matches the polarity."""
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        sheets = pd.read_excel(p, sheet_name=None)
        key = "POS" if polarity == "pos" else "NEG"
        for name, df in sheets.items():
            if key in name.upper() and not any(h in name.upper() for h in ("COMBINED",)):
                return df
        return next(iter(sheets.values()))
    sep = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    return pd.read_csv(p, sep=sep, low_memory=False)


def from_lipidsearch(table: pd.DataFrame, features: pd.DataFrame, ppm: float, rt_window: float) -> pd.DataFrame:
    """LipidSearch export (LipidMolec, CalcMass, BaseRt, MainIon, TotalGrade).

    Grade A or B with resolved chains = rank 1; grade C or D, or no resolved
    chains, = rank 3, and grade C/D names are reported at sum level at most.
    """
    t = table.copy()
    if "Rej" in t.columns:
        t = t[~t["Rej"].astype(str).str.lower().isin(["true", "1"])]
    t["mz"] = [ion_mz(m, a) for m, a in zip(t["CalcMass"], t["MainIon"])]
    t["rt"] = pd.to_numeric(t["BaseRt"], errors="coerce")
    t = t[t["mz"].notna()]
    t = match_by_mz_rt(features, t, ppm, rt_window)
    rows = []
    for _, r in t[t["feature_id"].notna()].iterrows():
        lip = parse_lipid(r["LipidMolec"])
        if not lip:
            continue
        grade = str(r.get("TotalGrade", "")).strip().upper()
        lvl = lip.level
        if grade in ("A", "B") and lvl == "species":
            rank = 1
        else:
            rank = 3
            if lvl == "species" and grade not in ("A", "B") and lip.lipid_class not in SINGLE_CHAIN:
                lvl = "sum"
        ev = f"grade {grade}; {r['MainIon']}; {r['match_ppm']:.1f} ppm, {r['match_drt']:+.2f} min"
        rows.append(_row(r["feature_id"], lip, lvl, rank, {"A": 1.0, "B": 0.8, "C": 0.5}.get(grade, 0.3),
                         "LipidSearch", ev))
    out = _candidates(rows)
    out.attrs["rt_offset"] = t.attrs.get("rt_offset", 0.0)
    out.attrs["n_input"] = len(table)
    return out


def from_msdial(table: pd.DataFrame, features: pd.DataFrame, ppm: float, rt_window: float) -> pd.DataFrame:
    """MS-DIAL alignment export (Metabolite name, Average Mz, Average Rt(min),
    Total score, MS/MS matched)."""
    t = table.rename(columns={"Average Mz": "mz", "Average Rt(min)": "rt"}).copy()
    t = t[t["Metabolite name"].notna() & ~t["Metabolite name"].astype(str).str.lower().isin(["unknown", "null", ""])]
    t = match_by_mz_rt(features, t, ppm, rt_window)
    rows = []
    for _, r in t[t["feature_id"].notna()].iterrows():
        lip = parse_lipid(str(r["Metabolite name"]).split("|")[-1])
        if not lip:
            continue
        msms = str(r.get("MS/MS matched", "")).lower() in ("true", "1")
        lvl = lip.level
        rank = 1 if (msms and lvl == "species") else 3
        score = pd.to_numeric(r.get("Total score", np.nan), errors="coerce")
        rows.append(_row(r["feature_id"], lip, lvl, rank, float(score) if pd.notna(score) else np.nan, "MS-DIAL",
                         f"MS/MS matched {msms}; {r['match_ppm']:.1f} ppm"))
    return _candidates(rows)


def from_library(table: pd.DataFrame, features: pd.DataFrame, ppm: float, rt_window: float,
                 min_score: float = 0.7, min_peaks: int = 4, same_run: bool = False) -> pd.DataFrame:
    """Spectral library matches, e.g. GNPS library results (Compound_Name,
    MQScore, SharedPeaks, SpecMZ, RT_Query in seconds, #Scan#).

    With ``same_run`` the ``#Scan#`` column is taken as this run's feature id;
    otherwise matches are placed by precursor m/z and retention time.
    """
    t = table.copy()
    t = t[(pd.to_numeric(t["MQScore"], errors="coerce") >= min_score)
          & (pd.to_numeric(t["SharedPeaks"], errors="coerce") >= min_peaks)]
    if same_run:
        t["feature_id"] = t["#Scan#"]
        t["match_ppm"] = np.nan
    else:
        t["mz"] = pd.to_numeric(t["SpecMZ"], errors="coerce")
        rtq = pd.to_numeric(t["RT_Query"], errors="coerce")
        t["rt"] = rtq / 60.0 if rtq.median() > 60 else rtq
        t = match_by_mz_rt(features, t, ppm, rt_window)
    rows = []
    for _, r in t[t["feature_id"].notna()].iterrows():
        lip = parse_lipid(r["Compound_Name"])
        if not lip:
            continue
        rows.append(_row(r["feature_id"], lip, lip.level, 2, float(r["MQScore"]), "spectral library",
                         f"cosine {float(r['MQScore']):.2f}, {int(r['SharedPeaks'])} peaks"))
    return _candidates(rows)


def from_sirius(table: pd.DataFrame, features: pd.DataFrame, ppm: float, rt_window: float,
                min_probability: float = 0.7, same_run: bool = False) -> pd.DataFrame:
    """SIRIUS/CANOPUS summary: class-level only (rank 4).

    Uses the most specific ClassyFire class when it maps to a lipid class and
    its probability is at least ``min_probability``.
    """
    t = table.copy()
    if same_run and "featureId" in t.columns:
        t["feature_id"] = t["featureId"]
    else:
        t = t.rename(columns={"row m/z": "mz", "row retention time": "rt"})
        t = match_by_mz_rt(features, t, ppm, rt_window)
    rows = []
    for _, r in t[t["feature_id"].notna()].iterrows():
        cls = None
        for col in ("ClassyFire#level 5", "ClassyFire#subclass", "ClassyFire#most specific class", "ClassyFire#class"):
            name = str(r.get(col, "")).strip().lower()
            prob = pd.to_numeric(r.get(f"{col} Probability", np.nan), errors="coerce")
            if name in CLASSYFIRE_MAP and (pd.isna(prob) or prob >= min_probability):
                cls = CLASSYFIRE_MAP[name]
                break
        if not cls:
            continue
        lip = Lipid(cls, [])
        rows.append(_row(r["feature_id"], lip, "class", 4, float(prob) if pd.notna(prob) else np.nan, "SIRIUS/CANOPUS",
                         f"ClassyFire {name}"))
    return _candidates(rows)


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------
def merge(cands: pd.DataFrame, features: pd.DataFrame, is_classes: set[str] | None = None) -> pd.DataFrame:
    """One row per annotated feature with name, confidence and quantification route.

    The best candidate is the lowest rank, then the deepest level, then the
    highest score. Confidence follows the module docstring. ``quant_route``
    is the highest level of the Q1-Q4 ladder the annotation allows: Q1 needs
    a standard of the same class in the run (``is_classes``).
    """
    is_classes = is_classes or set()
    if cands.empty:
        return pd.DataFrame(columns=["feature_id"])
    c = cands.copy()
    c["lvl"] = c["level"].map(LEVEL_ORDER)
    c["score_f"] = pd.to_numeric(c["score"], errors="coerce").fillna(0)
    c = c.sort_values(["feature_id", "rank", "lvl", "score_f"], ascending=[True, True, True, False])
    fmeta = features.set_index("id")[["mz", "rt"]]
    rows = []
    for fid, g in c.groupby("feature_id", sort=False):
        best = g.iloc[0]
        sources = sorted(set(g["source"]))
        classes = {class_family(x) for x in g["lipid_class"]}
        species = g[g["level"] == "species"]
        sp_counts = species.groupby("lipid_name")["source"].nunique() if len(species) else pd.Series(dtype=int)
        agreed_species = sp_counts[sp_counts >= 2].index.tolist()
        if best["level"] == "species" and best["rank"] == 1:
            conf, name = "Gold", best["lipid_name"]
        elif agreed_species:
            conf, name = "Gold", agreed_species[0]
        elif best["level"] in ("species", "sum") and best["sum_composition"]:
            conf, name = "Silver", f"{best['lipid_class']} {best['sum_composition']}"
        else:
            conf, name = "Bronze", best["lipid_class"]
        cls = best["lipid_class"]
        if conf in ("Gold", "Silver"):
            route = "Q1" if cls in is_classes else "Q2"
        else:
            route = "Q3"
        alts = [n for n in dict.fromkeys(g["lipid_name"]) if n != name][:5]
        rows.append({
            "feature_id": fid,
            "mz": fmeta["mz"].get(fid, np.nan),
            "rt": fmeta["rt"].get(fid, np.nan),
            "lipid_name": name,
            "lipid_class": cls,
            "confidence": conf,
            "quant_route": route,
            "best_source": best["source"],
            "best_rank": int(best["rank"]),
            "n_sources": len(sources),
            "sources": "; ".join(sources),
            "class_agreement": len(classes) == 1,
            "alternatives": "; ".join(alts),
        })
    return pd.DataFrame(rows)


def _carbons_dbs(name: str) -> tuple[float, float]:
    lip = parse_lipid(name)
    if lip is None or not lip.chains:
        return np.nan, np.nan
    m = re.search(r"(\d+):(\d+)", lip.sum_composition)
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)


# Number of fatty-acyl or sphingoid chains per class, for the composition check.
N_CHAINS = {
    "TG": 3, "DG": 2, "MG": 1, "CL": 4, "FA": 1, "CE": 1, "SPH": 1, "SPB": 1,
    "PC": 2, "PE": 2, "PG": 2, "PI": 2, "PS": 2, "PA": 2,
    "LPC": 1, "LPE": 1, "LPG": 1, "LPI": 1, "LPS": 1, "LPA": 1,
    "Cer": 2, "SM": 2, "HexCer": 2, "Hex1Cer": 2, "Hex2Cer": 2, "Hex3Cer": 2,
    "MGDG": 2, "DGDG": 2, "SQDG": 2, "DGTS": 2, "DGTA": 2, "DGGA": 2, "LDGTS": 1,
}
MIN_CHAIN_CARBONS = 10   # shortest acyl or sphingoid chain accepted as a membrane or storage lipid
MIN_CHAIN_CARBONS_MODIFIED = 5  # oxidation can truncate a chain (e.g. 9:0 aldehyde or acid cores)
MAX_CHAIN_DBS = 6        # most double bonds in one chain (22:6)


def check_composition(name: str, lipid_class: str) -> str:
    """``ok``, ``implausible`` or ``untested`` for a lipid name.

    Implausible: any named chain shorter than MIN_CHAIN_CARBONS
    (MIN_CHAIN_CARBONS_MODIFIED for oxidised or otherwise modified species),
    or a sum composition with fewer carbons or more than MAX_CHAIN_DBS double
    bonds per chain on average, given the class's number of chains. Classes
    not in N_CHAINS and class-only names are untested.
    """
    lip = parse_lipid(name)
    n = N_CHAINS.get(lipid_class)
    if lip is None or not lip.chains or n is None:
        return "untested"
    min_c = MIN_CHAIN_CARBONS_MODIFIED if lip.modification else MIN_CHAIN_CARBONS
    for ch in lip.chains:
        m = re.search(r"(\d+):(\d+)", ch)
        if len(lip.chains) > 1 and m and int(m.group(1)) < min_c:
            return "implausible"
    c, d = _carbons_dbs(name)
    if np.isnan(c):
        return "untested"
    if c < min_c * n or d > MAX_CHAIN_DBS * n:
        return "implausible"
    return "ok"


def _rt_group(name: str, lipid_class: str) -> str:
    """Class, split into diacyl, ether and modified (oxidised) species, which elute apart."""
    lip = parse_lipid(name)
    if lip is None:
        return lipid_class
    key = lipid_class
    if lip.sum_composition.startswith(("O-", "P-")):
        key += " ether"
    if lip.modification:
        key += " modified"
    return key


def check_rt(merged: pd.DataFrame, min_anchors: int = 6, min_tol: float = 1.0,
             max_tol: float = 2.0, n_sd: float = 4.0, max_extrapolation: int = 4) -> pd.DataFrame:
    """Flag annotations whose composition or retention time does not fit their class.

    ``composition_check`` (see check_composition) marks chain lengths and
    double-bond numbers that membrane and storage lipids do not have.

    ``rt_check`` uses the equivalent carbon number rule: within a lipid class
    on reversed phase, retention time rises with carbon number and falls with
    double bonds. Diacyl, ether and modified species of a class are fitted
    separately. For each group with at least ``min_anchors`` Gold
    annotations of plausible composition, fit rt ~ carbons + double bonds
    (refit once without anchors beyond 3 robust SD), then mark every row of
    the group:

    - ``ok``: within the tolerance, n_sd x robust SD bounded to
      [min_tol, max_tol] minutes;
    - ``outlier``: further away;
    - ``untested``: no fit for the group, class-only name, or carbon number
      more than ``max_extrapolation`` beyond the anchors' range.

    Adds ``composition_check``, ``rt_expected`` and ``rt_check``; confidence
    is not changed.
    """
    out = merged.copy()
    out["composition_check"] = "untested"
    out["rt_expected"] = np.nan
    out["rt_check"] = "untested"
    if out.empty:
        return out
    out["composition_check"] = [check_composition(n, c) for n, c in zip(out["lipid_name"], out["lipid_class"])]
    cd = out["lipid_name"].map(_carbons_dbs)
    out["_c"] = [x[0] for x in cd]
    out["_d"] = [x[1] for x in cd]
    out["_g"] = [_rt_group(n, c) for n, c in zip(out["lipid_name"], out["lipid_class"])]
    for _grp, g in out.groupby("_g"):
        anchors = g[(g["confidence"] == "Gold") & g["_c"].notna() & (g["composition_check"] != "implausible")]
        if len(anchors) < min_anchors or anchors["_c"].nunique() < 2:
            continue

        def fit(a):
            X = np.column_stack([np.ones(len(a)), a["_c"], a["_d"]])
            coef, *_ = np.linalg.lstsq(X, a["rt"].to_numpy(float), rcond=None)
            res = a["rt"].to_numpy(float) - X @ coef
            sd = 1.4826 * np.median(np.abs(res - np.median(res)))
            return coef, res, sd

        coef, res, sd = fit(anchors)
        keep = np.abs(res) <= max(3 * sd, min_tol)
        if keep.sum() >= min_anchors and not keep.all():
            anchors = anchors[keep]
            coef, res, sd = fit(anchors)
        tol = min(max_tol, max(min_tol, n_sd * sd))
        lo, hi = anchors["_c"].min() - max_extrapolation, anchors["_c"].max() + max_extrapolation
        rows = g[g["_c"].notna()]
        exp = coef[0] + coef[1] * rows["_c"] + coef[2] * rows["_d"]
        inside = rows["_c"].between(lo, hi)
        out.loc[rows.index, "rt_expected"] = exp.round(2)
        out.loc[rows.index[inside], "rt_check"] = np.where(
            (rows["rt"] - exp).abs()[inside] <= tol, "ok", "outlier")
    return out.drop(columns=["_c", "_d", "_g"])


def summarise(merged: pd.DataFrame, cands: pd.DataFrame, features: pd.DataFrame,
              source_inputs: dict[str, int] | None = None) -> dict:
    n_feat = int(len(features))
    n_ms2 = int((pd.to_numeric(features.get("fragment_scans", pd.Series(dtype=float)), errors="coerce") > 0).sum())
    conf = merged["confidence"].value_counts().to_dict() if len(merged) else {}
    by_source = cands.groupby("source")["feature_id"].nunique().to_dict() if len(cands) else {}
    only_source = {}
    if len(merged):
        for s in by_source:
            only_source[s] = int((merged["sources"] == s).sum())
    return {
        "features": n_feat,
        "features_with_ms2": n_ms2,
        "annotated_features": int(len(merged)),
        "confidence": {k: int(conf.get(k, 0)) for k in ("Gold", "Silver", "Bronze")},
        "quant_route": merged["quant_route"].value_counts().to_dict() if len(merged) else {},
        "features_per_source": {k: int(v) for k, v in by_source.items()},
        "features_only_from_source": only_source,
        "source_inputs": source_inputs or {},
        "class_conflicts": int((~merged["class_agreement"]).sum()) if len(merged) else 0,
        "classes": merged["lipid_class"].value_counts().head(25).to_dict() if len(merged) else {},
        "rt_check": (pd.crosstab(merged["confidence"], merged["rt_check"]).to_dict(orient="index")
                     if "rt_check" in merged and len(merged) else {}),
        "composition_check": (pd.crosstab(merged["confidence"], merged["composition_check"]).to_dict(orient="index")
                              if "composition_check" in merged and len(merged) else {}),
    }

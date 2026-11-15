"""Stage 3 quality checks for raw-to-feature processing (process mode).

Reads an MZmine "full" feature table (one row per aligned feature, with
per-file columns ``datafile:<file>:area``, ``:mz`` and ``:rt``) together with
the sample sheet, and computes run-level quality checks:

* internal standard (EquiSPLASH) retention-time stability, mass accuracy and
  peak-area variation across injections;
* pooled-QC feature reproducibility (relative standard deviation);
* blank contribution to sample features;
* per-injection summaries (features detected, total area) in injection order.

Default thresholds follow Samrat et al. (2025, Soil Biol. Biochem. 209,
109892): retention time within ±0.1 min, internal-standard peak area RSD
below 20 %, mass accuracy within ±5 ppm. The pooled-QC RSD limit (30 %) and the
blank ratio follow common practice (Broadhurst et al. 2018, Metabolomics).
All thresholds are parameters.

The functions are pure (DataFrames in, DataFrames out) so they can be unit
tested without MZmine output; ``bin/qc_report.py`` is the command-line wrapper.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROLES = {
    "sample",
    "pooled_qc",
    "blank",
    "extraction_blank",
    "nist",
    "standard_cal",
}



@dataclass
class Thresholds:
    rt_tolerance_min: float = 0.1
    is_area_rsd_pct: float = 20.0
    mass_accuracy_ppm: float = 5.0
    qc_feature_rsd_pct: float = 30.0
    qc_feature_pass_fraction: float = 0.7
    blank_ratio: float = 0.2
    match_ppm: float = 10.0
    match_rt_window_min: float | None = None


@dataclass
class Check:
    check: str
    threshold: str
    value: str
    status: str  # pass | warn | fail | not_run
    detail: str = ""


@dataclass
class QCResult:
    checks: list[Check] = field(default_factory=list)
    internal_standards: pd.DataFrame | None = None
    is_per_injection: pd.DataFrame | None = None
    injections: pd.DataFrame | None = None
    features: pd.DataFrame | None = None

    def checks_frame(self) -> pd.DataFrame:
        return pd.DataFrame([c.__dict__ for c in self.checks])


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def normalise_role(role: str) -> str:
    r = str(role).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {"qc": "pooled_qc", "pooledqc": "pooled_qc", "pool_qc": "pooled_qc",
               "srm": "nist", "nist_srm_1950": "nist", "srm1950": "nist",
               "solvent_blank": "blank", "process_blank": "extraction_blank",
               "standard": "standard_cal", "calibration": "standard_cal",
               "is_calibration": "standard_cal"}
    return aliases.get(r, r)


def read_samplesheet(path: str | Path) -> pd.DataFrame:
    """Read and validate the sample sheet.

    Required columns: ``file``, ``polarity``, ``role``. Optional columns are
    kept as they are. ``file`` may be a path; only its base name (without
    extension) is used to match MZmine columns.
    """
    sheet = pd.read_csv(path)
    sheet.columns = [c.strip().lower() for c in sheet.columns]
    missing = {"file", "polarity", "role"} - set(sheet.columns)
    if missing:
        raise ValueError(f"Sample sheet is missing column(s): {', '.join(sorted(missing))}")
    sheet["polarity"] = sheet["polarity"].str.strip().str.lower()
    bad_pol = sorted(set(sheet["polarity"]) - {"pos", "neg"})
    if bad_pol:
        raise ValueError(f"polarity must be 'pos' or 'neg', found: {bad_pol}")
    sheet["role"] = sheet["role"].map(normalise_role)
    bad_roles = sorted(set(sheet["role"]) - ROLES)
    if bad_roles:
        raise ValueError(f"Unknown role(s) {bad_roles}; use one of {sorted(ROLES)}")
    sheet["file_id"] = sheet["file"].map(file_id)
    dup = sheet["file_id"][sheet["file_id"].duplicated()]
    if not dup.empty:
        raise ValueError(f"Duplicate file names in sample sheet: {sorted(set(dup))}")
    return sheet


TARGET_REQUIRED = {"compound", "polarity", "adduct", "mz"}


def read_targets(path: str | Path, compounds: list[str] | None = None) -> pd.DataFrame:
    """Read an internal-standard target table.

    The default table (``assets/equisplash_targets.csv``) lists all 13
    EquiSPLASH components with several adducts each. A lab can use its own
    table instead. Required columns: ``compound``, ``polarity`` (pos/neg),
    ``adduct``, ``mz``. Optional: ``lipid_class`` (defaults to the compound
    name), ``priority`` (primary/secondary; default primary), ``spike_ng``,
    and ``rt``: the expected retention time on the lab's own method. When
    ``rt`` is given, a standard is only accepted within the RT window.

    ``compounds`` keeps only the named compounds (for example the subset of
    EquiSPLASH a lab actually spiked).
    """
    tg = pd.read_csv(path)
    tg.columns = [c.strip().lower() for c in tg.columns]
    missing = TARGET_REQUIRED - set(tg.columns)
    if missing:
        raise ValueError(f"Internal-standard table is missing column(s): {', '.join(sorted(missing))}")
    tg["polarity"] = tg["polarity"].str.strip().str.lower()
    tg["mz"] = pd.to_numeric(tg["mz"], errors="raise")
    if "lipid_class" not in tg.columns:
        tg["lipid_class"] = tg["compound"]
    if "priority" not in tg.columns:
        tg["priority"] = "primary"
    tg["priority"] = tg["priority"].fillna("primary").str.strip().str.lower()
    if "spike_ng" not in tg.columns:
        tg["spike_ng"] = np.nan
    if "rt" in tg.columns:
        tg["rt"] = pd.to_numeric(tg["rt"], errors="coerce")
    if compounds:
        wanted = [c.strip() for c in compounds if c.strip()]
        unknown = sorted(set(wanted) - set(tg["compound"]))
        if unknown:
            raise ValueError(f"Internal standard(s) not in the table: {unknown}. "
                             f"Available: {sorted(tg['compound'].unique())}")
        tg = tg[tg["compound"].isin(wanted)]
    return tg.reset_index(drop=True)


def file_id(path_or_name: str) -> str:
    name = Path(str(path_or_name).strip()).name
    return re.sub(r"\.(raw|mzml|mzxml)$", "", name, flags=re.IGNORECASE)


def read_mzml_start_time(path: str | Path, max_bytes: int = 200_000) -> datetime | None:
    """Return the run start time stored in an mzML header, if present."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(max_bytes).decode("utf-8", errors="ignore")
    except OSError:
        return None
    m = re.search(r'startTimeStamp="([^"]+)"', head)
    if not m:
        return None
    try:
        return datetime.fromisoformat(m.group(1).replace("Z", "+00:00"))
    except ValueError:
        return None


def per_file_matrix(full: pd.DataFrame, value: str) -> pd.DataFrame:
    """Return a features x files matrix for one per-file value (area, mz, rt).

    Columns are renamed to the file id (file name without extension).
    """
    pattern = re.compile(rf"^datafile:(.+?):{re.escape(value)}$")
    cols = {c: file_id(pattern.match(c).group(1)) for c in full.columns if pattern.match(c)}
    if not cols:
        raise ValueError(f"No per-file '{value}' columns found in feature table")
    mat = full[list(cols)].rename(columns=cols)
    mat.index = full["id"] if "id" in full.columns else full.index
    return mat.apply(pd.to_numeric, errors="coerce")


# ---------------------------------------------------------------------------
# Internal standards
# ---------------------------------------------------------------------------
def match_internal_standards(
    full: pd.DataFrame,
    targets: pd.DataFrame,
    polarity: str,
    area: pd.DataFrame,
    ppm: float = 10.0,
    rt_window: float | None = None,
    reference_files: list[str] | None = None,
) -> pd.DataFrame:
    """Find the feature that best represents each internal standard.

    For every compound, all target adducts of the given polarity are matched
    by m/z within ``ppm``. Soil contains endogenous lipids within a few ppm
    of the deuterated standards, so m/z alone is not enough:

    * If ``reference_files`` are given (injections of the standard on its own,
      role ``standard_cal``), the candidate with the largest area in those
      injections wins, whatever its adduct. This fixes the standard's feature
      and retention time from data where it dominates.
    * Otherwise a primary adduct beats a secondary one, then the highest
      median area across all injections wins.

    Returns one row per compound (matched or not), with ``selected_by``.
    """
    tg = targets[targets["polarity"].str.lower() == polarity].copy()
    mz = pd.to_numeric(full["mz"], errors="coerce").to_numpy()
    rt = pd.to_numeric(full["rt"], errors="coerce").to_numpy()
    ids = full["id"].to_numpy() if "id" in full.columns else full.index.to_numpy()
    median_area = area.median(axis=1, skipna=True).reindex(ids).to_numpy()
    ref_cols = [f for f in (reference_files or []) if f in area.columns]
    ref_area = (
        area[ref_cols].max(axis=1, skipna=True).reindex(ids).to_numpy() if ref_cols else None
    )

    rows = []
    for compound, grp in tg.groupby("compound", sort=True):
        best = None
        for _, t in grp.iterrows():
            err = (mz - t["mz"]) / t["mz"] * 1e6
            hit = np.abs(err) <= ppm
            has_rt = "rt" in t and pd.notna(t.get("rt"))
            if has_rt:
                hit &= np.abs(rt - float(t["rt"])) <= (rt_window if rt_window is not None else 0.3)
            for i in np.flatnonzero(hit):
                if has_rt:
                    key = (abs(rt[i] - float(t["rt"])), 0 if t["priority"] == "primary" else 1)
                    basis = "expected RT"
                elif ref_area is not None:
                    ra = np.nan_to_num(ref_area[i], nan=0.0)
                    if ra <= 0:
                        continue
                    key = (-ra, 0 if t["priority"] == "primary" else 1)
                    basis = "standard injection"
                else:
                    key = (0 if t["priority"] == "primary" else 1, -np.nan_to_num(median_area[i], nan=0.0))
                    basis = "m/z only (unconfirmed)"
                cand = (key, i, t, err[i], basis)
                if best is None or cand[0] < best[0]:
                    best = cand
        first = grp.iloc[0]
        if best is None:
            rows.append({"compound": compound, "lipid_class": first["lipid_class"], "matched": False})
            continue
        _, i, t, e, basis = best
        rows.append(
            {
                "compound": compound,
                "lipid_class": t["lipid_class"],
                "matched": True,
                "adduct": t["adduct"],
                "priority": t["priority"],
                "selected_by": basis,
                "target_mz": float(t["mz"]),
                "feature_id": ids[i],
                "feature_mz": float(mz[i]),
                "feature_rt": float(rt[i]),
                "feature_ppm_error": float(e),
                "spike_ng": t.get("spike_ng"),
            }
        )
    return pd.DataFrame(rows)


def internal_standard_metrics(
    matches: pd.DataFrame,
    area: pd.DataFrame,
    mzs: pd.DataFrame,
    rts: pd.DataFrame,
    sheet: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-standard summary and per-injection detail for matched standards.

    Following Samrat et al. (2025), the pass/fail criteria apply to QC
    injections and standards, not to soils, whose standard signal varies with
    matrix effects and extraction recovery:

    * retention time and mass accuracy: pooled QCs, NIST and standard
      calibration injections;
    * peak-area RSD: pooled QCs (a repeated identical sample); if there are
      none, NIST injections; if neither, samples, flagged in ``qc_basis``.

    Soil values are reported alongside (``sample_*`` columns) for information.
    """
    roles = sheet.set_index("file_id")["role"]
    files = [f for f in area.columns if f in roles.index]
    by_role = {r: [f for f in files if roles[f] == r] for r in ROLES}
    rt_ppm_set = by_role["pooled_qc"] + by_role["nist"] + by_role["standard_cal"]
    if by_role["pooled_qc"]:
        rsd_set, basis = by_role["pooled_qc"], "pooled QC"
    elif by_role["nist"]:
        rsd_set, basis = by_role["nist"], "NIST (no pooled QC)"
    else:
        rsd_set, basis = by_role["sample"], "samples (no pooled QC or NIST)"
    if not rt_ppm_set:
        rt_ppm_set = by_role["sample"]
    samples = by_role["sample"]

    def rsd(values: pd.Series) -> float:
        v = values.replace(0, np.nan).dropna()
        return float(v.std(ddof=1) / v.mean() * 100) if len(v) >= 2 else np.nan

    def max_rt_dev(values: pd.Series) -> float:
        v = values.dropna()
        return float((v - v.median()).abs().max()) if len(v) else np.nan

    detail_rows = []
    summary_rows = []
    for _, m in matches.iterrows():
        if not m.get("matched", False):
            summary_rows.append({"compound": m["compound"], "lipid_class": m["lipid_class"], "matched": False})
            continue
        fid = m["feature_id"]
        a = area.loc[fid, files]
        x = mzs.loc[fid, files]
        r = rts.loc[fid, files]
        ppm = (x - m["target_mz"]) / m["target_mz"] * 1e6
        for f in files:
            detail_rows.append(
                {
                    "compound": m["compound"],
                    "file_id": f,
                    "role": roles[f],
                    "area": a[f],
                    "rt": r[f],
                    "ppm_error": ppm[f],
                }
            )
        present = a > 0
        rc = r[rt_ppm_set][present[rt_ppm_set]]
        pc = ppm[rt_ppm_set][present[rt_ppm_set]]
        rs = r[samples][present[samples]]
        ps = ppm[samples][present[samples]]
        summary_rows.append(
            {
                "compound": m["compound"],
                "lipid_class": m["lipid_class"],
                "matched": True,
                "adduct": m["adduct"],
                "priority": m["priority"],
                "selected_by": m.get("selected_by"),
                "feature_id": fid,
                "qc_basis": basis,
                "median_rt": float(rc.median()) if len(rc) else np.nan,
                "max_rt_deviation_min": max_rt_dev(rc),
                "median_ppm_error": float(pc.median()) if len(pc) else np.nan,
                "max_abs_ppm_error": float(pc.abs().max()) if len(pc) else np.nan,
                "area_rsd_pct": rsd(a[rsd_set]),
                "detected_in_qc": int(present[rt_ppm_set].sum()),
                "qc_injections": len(rt_ppm_set),
                "sample_area_rsd_pct": rsd(a[samples]),
                "sample_max_rt_deviation_min": max_rt_dev(rs),
                "sample_max_abs_ppm_error": float(ps.abs().max()) if len(ps) else np.nan,
                "detected_in_samples": int(present[samples].sum()),
                "samples": len(samples),
            }
        )
    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


# ---------------------------------------------------------------------------
# Features
# ---------------------------------------------------------------------------
def feature_quality(area: pd.DataFrame, sheet: pd.DataFrame) -> pd.DataFrame:
    """Per-feature pooled-QC RSD and blank-to-sample ratio."""
    roles = sheet.set_index("file_id")["role"]
    by_role = {r: [f for f in area.columns if roles.get(f) == r] for r in ROLES}
    out = pd.DataFrame(index=area.index)
    qc = area[by_role["pooled_qc"]].replace(0, np.nan)
    if qc.shape[1] >= 2:
        out["qc_detected"] = qc.notna().sum(axis=1)
        out["qc_rsd_pct"] = qc.std(axis=1, ddof=1) / qc.mean(axis=1) * 100
    samples = area[by_role["sample"]].replace(0, np.nan)
    out["sample_detected"] = samples.notna().sum(axis=1)
    out["sample_median_area"] = samples.median(axis=1)
    blanks = area[by_role["blank"] + by_role["extraction_blank"]].replace(0, np.nan)
    if blanks.shape[1] >= 1:
        out["blank_max_area"] = blanks.max(axis=1)
        out["blank_ratio"] = out["blank_max_area"] / out["sample_median_area"]
    return out


def injection_summary(area: pd.DataFrame, sheet: pd.DataFrame, start_times: dict[str, datetime | None]) -> pd.DataFrame:
    """One row per injection: role, features detected, total area, run order."""
    rows = []
    roles = sheet.set_index("file_id")["role"]
    for f in area.columns:
        col = area[f]
        rows.append(
            {
                "file_id": f,
                "role": roles.get(f, "not_in_sheet"),
                "start_time": start_times.get(f),
                "features_detected": int((col > 0).sum()),
                "total_area": float(col.fillna(0).sum()),
            }
        )
    inj = pd.DataFrame(rows)
    if inj["start_time"].notna().all():
        inj = inj.sort_values("start_time")
        inj["order_source"] = "mzML start time"
    else:
        order = {fid: i for i, fid in enumerate(sheet["file_id"])}
        inj["_o"] = inj["file_id"].map(order)
        inj = inj.sort_values("_o").drop(columns="_o")
        inj["order_source"] = "sample sheet row order"
    inj["injection_order"] = range(1, len(inj) + 1)
    return inj.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
def _status(ok: bool, warn: bool = False) -> str:
    return "pass" if ok else ("warn" if warn else "fail")


def _threshold_check(name: str, threshold: str, values: pd.Series, limit: float, compounds: pd.Series,
                     note: str = "") -> Check:
    """Judge one metric per standard; standards without a value are 'not measured', never 'within'."""
    measured = values.notna()
    bad = compounds[measured & (values > limit)].tolist()
    n_meas = int(measured.sum())
    missing = compounds[~measured].tolist()
    if n_meas == 0:
        return Check(name, threshold, "not measured", "not_run", note or "no QC injection contains these standards")
    value = f"{n_meas - len(bad)} of {n_meas} measured within"
    detail = ", ".join(bad)
    if missing:
        value += f"; {len(missing)} not measured"
        detail = (detail + "; " if detail else "") + "not measured: " + ", ".join(missing)
    status = _status(not bad, warn=len(bad) < n_meas)
    if missing and status == "pass":
        status = "warn"
    return Check(name, threshold, value, status, (detail + (" " + note if note else "")).strip())


def build_checks(
    is_summary: pd.DataFrame,
    feats: pd.DataFrame,
    inj: pd.DataFrame,
    t: Thresholds,
) -> list[Check]:
    checks: list[Check] = []

    roles = inj["role"].value_counts().to_dict()
    for role, need, why in [
        ("pooled_qc", 3, "feature reproducibility and drift"),
        ("blank", 1, "background and carryover"),
        ("extraction_blank", 1, "contaminants from extraction"),
        ("nist", 1, "accuracy against reference values (stage 5)"),
        ("standard_cal", 1, "standard calibration and recovery (stage 5)"),
    ]:
        n = roles.get(role, 0)
        checks.append(Check(f"Injections with role '{role}'", f">= {need}", str(n),
                            "pass" if n >= need else "warn", why))
    not_in_sheet = int((inj["role"] == "not_in_sheet").sum())
    if not_in_sheet:
        checks.append(Check("Files missing from the sample sheet", "0", str(not_in_sheet), "warn",
                            "these injections were processed but have no role"))

    m = is_summary[is_summary.get("matched", False) == True] if not is_summary.empty else is_summary  # noqa: E712
    n_total = len(is_summary)
    checks.append(Check("Internal standards found", f"{n_total} of {n_total}", f"{len(m)} of {n_total}",
                        _status(len(m) == n_total, warn=len(m) > 0),
                        "matched by m/z; a secondary adduct is used if the primary is absent"))
    if len(m) and "selected_by" in m:
        unconfirmed = m[m["selected_by"] == "m/z only (unconfirmed)"]["compound"].tolist()
        checks.append(Check("Internal standard identity confirmed", "by expected RT or a standard injection",
                            f"{len(m) - len(unconfirmed)} of {len(m)} confirmed",
                            _status(not unconfirmed, warn=True),
                            ("matched by m/z only; soil lipids can sit within a few ppm. Add an 'rt' column "
                             "to the standards table or a standard-only injection: " + ", ".join(unconfirmed))
                            if unconfirmed else ""))
    if len(m):
        comp = m["compound"].reset_index(drop=True)
        checks.append(_threshold_check("Internal standard retention time (QC injections and standards)",
                                       f"±{t.rt_tolerance_min} min from median",
                                       m["max_rt_deviation_min"].reset_index(drop=True), t.rt_tolerance_min, comp))
        checks.append(_threshold_check("Internal standard mass accuracy (QC injections and standards)",
                                       f"±{t.mass_accuracy_ppm} ppm",
                                       m["max_abs_ppm_error"].reset_index(drop=True), t.mass_accuracy_ppm, comp))
        basis = str(m["qc_basis"].iloc[0]) if "qc_basis" in m else "pooled QC"
        c = _threshold_check("Internal standard peak area RSD", f"< {t.is_area_rsd_pct:g}% over {basis}",
                             m["area_rsd_pct"].reset_index(drop=True), t.is_area_rsd_pct, comp)
        if basis != "pooled QC" and c.status in ("pass", "fail"):
            c.status = "warn"
            c.detail = (c.detail + " (no pooled QC: fallback set, so a warning only)").strip()
        checks.append(c)
        if "sample_area_rsd_pct" in m:
            med = float(m["sample_area_rsd_pct"].median())
            checks.append(Check("Internal standard variation across samples (information)", "reported, not judged",
                                f"median RSD {med:.0f}%", "pass",
                                "soils differ in matrix effects and recovery; stage 5 corrects for this"))

    if "qc_rsd_pct" in feats:
        detected = feats[feats["qc_detected"] >= 2]
        frac = float((detected["qc_rsd_pct"] < t.qc_feature_rsd_pct).mean()) if len(detected) else np.nan
        checks.append(Check(f"Pooled-QC features with RSD < {t.qc_feature_rsd_pct:g}%",
                            f">= {t.qc_feature_pass_fraction:.0%} of features",
                            f"{frac:.0%} of {len(detected)}" if len(detected) else "no features",
                            _status(bool(frac >= t.qc_feature_pass_fraction)) if len(detected) else "not_run",
                            "features detected in at least two pooled QCs"))
    else:
        checks.append(Check("Pooled-QC feature reproducibility", "needs >= 2 pooled QCs", "-", "not_run", ""))

    if "blank_ratio" in feats:
        in_samples = feats[feats["sample_detected"] > 0]
        flagged = int((in_samples["blank_ratio"] > t.blank_ratio).sum())
        frac = flagged / len(in_samples) if len(in_samples) else np.nan
        checks.append(Check(f"Sample features also high in blanks (blank > {t.blank_ratio:.0%} of sample median)",
                            "reported, not removed", f"{flagged} of {len(in_samples)} ({frac:.0%})", "pass",
                            "listed in qc_features.tsv for filtering in later stages"))
    return checks


def run_qc(
    full: pd.DataFrame,
    sheet: pd.DataFrame,
    targets: pd.DataFrame,
    polarity: str,
    start_times: dict[str, datetime | None] | None = None,
    thresholds: Thresholds | None = None,
) -> QCResult:
    t = thresholds or Thresholds()
    sheet = sheet[sheet["polarity"] == polarity].copy()
    area = per_file_matrix(full, "area")
    mzs = per_file_matrix(full, "mz")
    rts = per_file_matrix(full, "rt")
    ref_files = sheet.loc[sheet["role"] == "standard_cal", "file_id"].tolist()
    matches = match_internal_standards(full, targets, polarity, area, ppm=t.match_ppm,
                                       rt_window=t.match_rt_window_min, reference_files=ref_files)
    is_summary, is_detail = internal_standard_metrics(matches, area, mzs, rts, sheet)
    feats = feature_quality(area, sheet)
    inj = injection_summary(area, sheet, start_times or {})
    checks = build_checks(is_summary, feats, inj, t)
    return QCResult(checks=checks, internal_standards=is_summary, is_per_injection=is_detail,
                    injections=inj, features=feats)

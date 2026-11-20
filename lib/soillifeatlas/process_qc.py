"""Stage 3 quality checks for raw-to-feature processing (process mode).

Reads an MZmine "full" feature table (one row per aligned feature, with
per-file columns ``datafile:<file>:area``, ``:mz`` and ``:rt``) together with
the sample sheet, and computes run-level quality checks:

* internal standard (EquiSPLASH) retention-time stability, mass accuracy and
  peak-area variation across injections;
* pooled-QC feature reproducibility (relative standard deviation), signal
  drift with injection order, and the percentage-difference (Bland–Altman)
  view used in Samrat et al. (2025);
* blank contribution to sample features;
* per-injection summaries (features detected, total area) in injection order.

Only peaks that MZmine actually detected (``feature_state`` DETECTED) count
for the QC metrics. Gap-filled values (ESTIMATED) fill holes for statistics in
later stages, but they are noise-level estimates and would make every
reproducibility measure look worse than the data are.

Default thresholds follow Samrat et al. (2025, Soil Biol. Biochem. 209,
109892): retention time within ±0.1 min, internal-standard peak area RSD
below 20 %, mass accuracy within ±5 ppm. The pooled-QC RSD limit (30 %), the
70 % detection rule and the blank ratio follow Broadhurst et al. (2018,
Metabolomics 14:72). All thresholds are parameters.

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
    "conditioning",  # system-conditioning QCs at the start: run, but left out of every QC statistic
}


@dataclass
class Thresholds:
    rt_tolerance_min: float = 0.1
    is_area_rsd_pct: float = 20.0
    mass_accuracy_ppm: float = 5.0
    qc_feature_rsd_pct: float = 30.0
    qc_feature_pass_fraction: float = 0.7
    qc_min_detect_fraction: float = 0.7  # feature must be peak-picked in this share of pooled QCs
    is_min_qc_detected: int = 3  # a standard RSD needs at least this many detected pooled QCs ...
    is_min_qc_fraction: float = 0.5  # ... and at least this share of them
    qc_drift_pct: float = 20.0  # change in summed pooled-QC signal over the run
    is_qc_shift_window_min: float = 0.2  # how far the standard may move between calibration and QCs
    block_gap_hours: float = 12.0  # a longer pause starts a new sequence block
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
    drift: pd.DataFrame | None = None
    agreement: dict | None = None

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
               "is_calibration": "standard_cal", "cond_qc": "conditioning",
               "conditioning_qc": "conditioning", "condqc": "conditioning"}
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


def detected_matrix(full: pd.DataFrame, area: pd.DataFrame) -> pd.DataFrame:
    """True where MZmine peak-picked the feature in that file.

    Uses the per-file ``feature_state`` columns (DETECTED vs ESTIMATED from gap
    filling). Tables without them (older exports, other tools) fall back to
    "area above zero".
    """
    pattern = re.compile(r"^datafile:(.+?):feature_state$")
    cols = {c: file_id(pattern.match(c).group(1)) for c in full.columns if pattern.match(c)}
    if not cols:
        return area > 0
    state = full[list(cols)].rename(columns=cols)
    state.index = area.index
    det = state.apply(lambda col: col.astype(str).str.upper().eq("DETECTED"))
    return det.reindex(columns=area.columns, fill_value=False)


def sequence_blocks(start_times: dict[str, datetime | None], files: list[str], gap_hours: float) -> dict[str, int]:
    """Group injections into sequence blocks: a pause longer than ``gap_hours`` starts a new block.

    Retention time and mass accuracy drift between blocks run days apart (column
    re-equilibration, recalibration), so they are judged within a block.
    Without start times every file is in block 1.
    """
    known = sorted((t, f) for f, t in start_times.items() if t is not None and f in files)
    if len(known) < len(files):
        return {f: 1 for f in files}
    blocks, block, prev = {}, 1, None
    for t, f in known:
        if prev is not None and (t - prev).total_seconds() > gap_hours * 3600:
            block += 1
        blocks[f] = block
        prev = t
    return blocks


# ---------------------------------------------------------------------------
# Internal standards
# ---------------------------------------------------------------------------
def _prio(t) -> int:
    return 0 if t["priority"] == "primary" else 1


def match_internal_standards(
    full: pd.DataFrame,
    targets: pd.DataFrame,
    polarity: str,
    area: pd.DataFrame,
    ppm: float = 10.0,
    rt_window: float | None = None,
    reference_files: list[str] | None = None,
) -> pd.DataFrame:
    """Find the feature that identifies each internal standard.

    For every compound, all target adducts of the given polarity are matched
    by m/z within ``ppm``. Soil contains endogenous lipids within a few ppm
    of the deuterated standards, so m/z alone is not enough:

    * With an ``rt`` column in the targets, the candidate closest to the
      expected retention time wins (within ``rt_window``).
    * Otherwise, if ``reference_files`` are given (injections of the standard
      itself, role ``standard_cal``), only candidates seen in those injections
      are kept; a primary adduct beats a secondary one, then the largest area
      in the reference injections wins.
    * Otherwise a primary adduct beats a secondary one, then the highest
      median area across all injections wins (flagged as unconfirmed).

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
                    key = (abs(rt[i] - float(t["rt"])), _prio(t))
                    basis = "expected RT"
                elif ref_area is not None:
                    ra = np.nan_to_num(ref_area[i], nan=0.0)
                    if ra <= 0:
                        continue
                    key = (_prio(t), -ra)
                    basis = "standard injection"
                else:
                    key = (_prio(t), -np.nan_to_num(median_area[i], nan=0.0))
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


def qc_features_for_standards(
    matches: pd.DataFrame,
    full: pd.DataFrame,
    detected: pd.DataFrame,
    area: pd.DataFrame,
    qc_files: list[str],
    ppm: float = 10.0,
    window: float = 0.2,
) -> pd.DataFrame:
    """Pick the feature row on which each standard is measured in the pooled QCs.

    The identity feature comes from the calibration injections. When those ran
    on another day, the standard can elute a little earlier or later in the QCs
    and soils, and the aligner may put it on a separate row. Among rows of the
    same adduct within ``window`` minutes and ``ppm`` of the identity feature,
    the one peak-picked in the most pooled QCs (then the largest QC area) is
    used for the QC metrics. Columns ``qc_feature_id`` and ``qc_feature_rt``
    are added; they equal the identity feature when it is already the best.
    """
    out = matches.copy()
    if out.empty or "feature_id" not in out:
        return out
    mz = pd.to_numeric(full["mz"], errors="coerce").to_numpy()
    rt = pd.to_numeric(full["rt"], errors="coerce").to_numpy()
    ids = full["id"].to_numpy() if "id" in full.columns else full.index.to_numpy()
    qcs = [f for f in qc_files if f in detected.columns]
    n_det = detected[qcs].sum(axis=1).reindex(ids).to_numpy() if qcs else np.zeros(len(ids))
    qc_area = area[qcs].where(detected[qcs]).median(axis=1).reindex(ids).fillna(0).to_numpy() if qcs else np.zeros(len(ids))
    qc_id, qc_rt = [], []
    for _, m in out.iterrows():
        if not m.get("matched", False):
            qc_id.append(np.nan)
            qc_rt.append(np.nan)
            continue
        hit = (np.abs((mz - m["target_mz"]) / m["target_mz"] * 1e6) <= ppm) & (np.abs(rt - m["feature_rt"]) <= window)
        cand = np.flatnonzero(hit)
        i_id = int(np.flatnonzero(ids == m["feature_id"])[0])
        best = max(cand, key=lambda i: (n_det[i], qc_area[i], i == i_id)) if len(cand) else i_id
        if n_det[best] <= n_det[i_id]:
            best = i_id
        qc_id.append(ids[best])
        qc_rt.append(float(rt[best]))
    out["qc_feature_id"] = qc_id
    out["qc_feature_rt"] = qc_rt
    return out


def _rsd(values: pd.Series, min_n: int = 2) -> float:
    v = values.dropna()
    v = v[v > 0]
    return float(v.std(ddof=1) / v.mean() * 100) if len(v) >= max(min_n, 2) else np.nan


def internal_standard_metrics(
    matches: pd.DataFrame,
    area: pd.DataFrame,
    mzs: pd.DataFrame,
    rts: pd.DataFrame,
    sheet: pd.DataFrame,
    detected: pd.DataFrame | None = None,
    blocks: dict[str, int] | None = None,
    min_qc_detected: int = 3,
    min_qc_fraction: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-standard summary and per-injection detail for matched standards.

    Following Samrat et al. (2025), the pass/fail criteria apply to QC
    injections and standards, not to soils, whose standard signal varies with
    matrix effects and extraction recovery:

    * retention time and mass accuracy: pooled QCs, NIST and standard
      calibration injections, judged within each sequence block (``blocks``);
    * peak-area RSD: pooled QCs (a repeated identical sample); if there are
      none, NIST injections; if neither, samples, flagged in ``qc_basis``.
      It needs at least ``min_qc_detected`` peak-picked QCs and
      ``min_qc_fraction`` of them, otherwise it is "not measured".

    Only peak-picked values count (``detected``). Calibration injections are
    read on the identity feature, all other injections on the QC feature (see
    :func:`qc_features_for_standards`). Soil values are reported alongside
    (``sample_*`` columns) for information.
    """
    roles = sheet.set_index("file_id")["role"]
    files = [f for f in area.columns if f in roles.index]
    by_role = {r: [f for f in files if roles[f] == r] for r in ROLES}
    if detected is None:
        detected = area > 0
    blocks = blocks or {}
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
    cal = set(by_role["standard_cal"])
    need = max(min_qc_detected, int(np.ceil(min_qc_fraction * len(rsd_set))))

    def max_rt_dev(values: pd.Series) -> float:
        v = values.dropna()
        if not len(v):
            return np.nan
        dev = [(g - g.median()).abs().max() for _, g in v.groupby(lambda f: blocks.get(f, 1))]
        return float(max(dev))

    detail_rows = []
    summary_rows = []
    for _, m in matches.iterrows():
        if not m.get("matched", False):
            summary_rows.append({"compound": m["compound"], "lipid_class": m["lipid_class"], "matched": False})
            continue
        fid = m["feature_id"]
        qid = m.get("qc_feature_id", fid)
        qid = fid if pd.isna(qid) else qid
        row_for = {f: (fid if f in cal else qid) for f in files}
        a = pd.Series({f: area.at[row_for[f], f] for f in files}, dtype=float)
        x = pd.Series({f: mzs.at[row_for[f], f] for f in files}, dtype=float)
        r = pd.Series({f: rts.at[row_for[f], f] for f in files}, dtype=float)
        det = pd.Series({f: bool(detected.at[row_for[f], f]) for f in files})
        ppm = (x - m["target_mz"]) / m["target_mz"] * 1e6
        for f in files:
            detail_rows.append(
                {
                    "compound": m["compound"],
                    "file_id": f,
                    "role": roles[f],
                    "feature_id": row_for[f],
                    "block": blocks.get(f, 1),
                    "detected": bool(det[f]),
                    "area": a[f],
                    "rt": r[f],
                    "ppm_error": ppm[f],
                }
            )
        a_det, r_det, p_det = a.where(det), r.where(det), ppm.where(det)
        rc, pc = r_det[rt_ppm_set], p_det[rt_ppm_set]
        n_rsd = int(det[rsd_set].sum())
        cal_det = a_det[list(cal)].dropna() if cal else pd.Series(dtype=float)
        qc_med = float(a_det[rsd_set].median()) if n_rsd else np.nan
        summary_rows.append(
            {
                "compound": m["compound"],
                "lipid_class": m["lipid_class"],
                "matched": True,
                "adduct": m["adduct"],
                "priority": m["priority"],
                "selected_by": m.get("selected_by"),
                "feature_id": fid,
                "qc_feature_id": qid,
                "rt_identity": m.get("feature_rt"),
                "rt_qc": m.get("qc_feature_rt", m.get("feature_rt")),
                "qc_basis": basis,
                "median_rt": float(rc.median()) if rc.notna().any() else np.nan,
                "max_rt_deviation_min": max_rt_dev(rc),
                "median_ppm_error": float(pc.median()) if pc.notna().any() else np.nan,
                "max_abs_ppm_error": float(pc.abs().max()) if pc.notna().any() else np.nan,
                "area_rsd_pct": _rsd(a_det[rsd_set], need) if n_rsd >= need else np.nan,
                "detected_in_qc": n_rsd,
                "qc_injections": len(rsd_set),
                "qc_median_area": qc_med,
                "lowest_cal_area": float(cal_det.min()) if len(cal_det) else np.nan,
                "qc_to_lowest_cal": qc_med / float(cal_det.min()) if len(cal_det) and n_rsd else np.nan,
                "sample_area_rsd_pct": _rsd(a_det[samples]),
                "sample_max_abs_ppm_error": float(p_det[samples].abs().max()) if p_det[samples].notna().any() else np.nan,
                "detected_in_samples": int(det[samples].sum()),
                "samples": len(samples),
            }
        )
    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


# ---------------------------------------------------------------------------
# Features
# ---------------------------------------------------------------------------
def feature_quality(area: pd.DataFrame, sheet: pd.DataFrame, detected: pd.DataFrame | None = None,
                    min_detect_fraction: float = 0.7) -> pd.DataFrame:
    """Per-feature pooled-QC detection and RSD, and blank-to-sample ratio.

    ``qc_detected`` counts pooled QCs where the feature was peak-picked.
    ``qc_rsd_pct`` uses those values only, and only for features peak-picked
    in at least ``min_detect_fraction`` of the pooled QCs (``qc_eligible``).
    The blank ratio uses all values, gap-filled included, to stay cautious.
    """
    roles = sheet.set_index("file_id")["role"]
    by_role = {r: [f for f in area.columns if roles.get(f) == r] for r in ROLES}
    if detected is None:
        detected = area > 0
    out = pd.DataFrame(index=area.index)
    qcs = by_role["pooled_qc"]
    if len(qcs) >= 2:
        qc = area[qcs].where(detected[qcs] & (area[qcs] > 0))
        need = max(2, int(np.ceil(min_detect_fraction * len(qcs))))
        out["qc_detected"] = qc.notna().sum(axis=1)
        out["qc_eligible"] = out["qc_detected"] >= need
        rsd = qc.std(axis=1, ddof=1) / qc.mean(axis=1) * 100
        out["qc_rsd_pct"] = rsd.where(out["qc_eligible"])
    samples = area[by_role["sample"]].replace(0, np.nan)
    # Peak-picked only, like qc_detected; the median and the blank ratio keep all values (cautious).
    out["sample_detected"] = (samples.notna() & detected[by_role["sample"]]).sum(axis=1)
    out["sample_median_area"] = samples.median(axis=1)
    blanks = area[by_role["blank"] + by_role["extraction_blank"]].replace(0, np.nan)
    if blanks.shape[1] >= 1:
        out["blank_max_area"] = blanks.max(axis=1)
        out["blank_ratio"] = out["blank_max_area"] / out["sample_median_area"]
    return out


def qc_agreement(area: pd.DataFrame, detected: pd.DataFrame, qc_files: list[str], eligible: pd.Series) -> dict:
    """Percentage difference of each pooled QC from the feature's QC mean (Bland–Altman view).

    This is the pooled-QC figure of Samrat et al. (2025): limits of agreement
    are the mean difference ± 1.96 SD over all feature × QC points.
    """
    qcs = [f for f in qc_files if f in area.columns]
    if len(qcs) < 2 or not eligible.any():
        return {}
    qc = area.loc[eligible, qcs].where(detected.loc[eligible, qcs])
    pct = qc.sub(qc.mean(axis=1), axis=0).div(qc.mean(axis=1), axis=0) * 100
    log_mean = np.log(qc.mean(axis=1))
    vals = pct.stack()
    mean, sd = float(vals.mean()), float(vals.std(ddof=1))
    lo, hi = mean - 1.96 * sd, mean + 1.96 * sd
    within = float(((vals >= lo) & (vals <= hi)).mean())
    pts = pct.stack().rename("pct_diff").reset_index()
    pts.columns = ["feature_id", "file_id", "pct_diff"]
    pts["log_mean_area"] = pts["feature_id"].map(log_mean)
    return {"loa_low": lo, "loa_high": hi, "mean_diff": mean, "within_fraction": within,
            "n_points": int(len(vals)), "n_features": int(eligible.sum()), "points": pts}


def qc_drift(area: pd.DataFrame, detected: pd.DataFrame, inj: pd.DataFrame, qc_files: list[str]) -> pd.DataFrame:
    """Summed pooled-QC signal of features peak-picked in every pooled QC, in injection order.

    Adds a straight-line fit; the change from first to last fitted value is the
    run drift (``drift_pct`` on every row).
    """
    qcs = [f for f in qc_files if f in area.columns]
    if len(qcs) < 3:
        return pd.DataFrame()
    allqc = detected[qcs].all(axis=1)
    if not allqc.any():
        return pd.DataFrame()
    total = area.loc[allqc, qcs].sum(axis=0)
    order = inj.set_index("file_id")["injection_order"]
    d = pd.DataFrame({"file_id": qcs, "summed_area": total[qcs].to_numpy(),
                      "injection_order": [order.get(f, np.nan) for f in qcs]}).sort_values("injection_order")
    slope, icpt = np.polyfit(d["injection_order"], d["summed_area"], 1)
    d["fitted"] = slope * d["injection_order"] + icpt
    first, last = d["fitted"].iloc[0], d["fitted"].iloc[-1]
    d["drift_pct"] = (last - first) / first * 100 if first else np.nan
    d["features_used"] = int(allqc.sum())
    return d.reset_index(drop=True)


def injection_summary(area: pd.DataFrame, sheet: pd.DataFrame, start_times: dict[str, datetime | None],
                      detected: pd.DataFrame | None = None, blocks: dict[str, int] | None = None) -> pd.DataFrame:
    """One row per injection: role, features detected (peak-picked), total area, run order, block."""
    if detected is None:
        detected = area > 0
    rows = []
    roles = sheet.set_index("file_id")["role"]
    for f in area.columns:
        col = area[f].where(detected[f])
        rows.append(
            {
                "file_id": f,
                "role": roles.get(f, "not_in_sheet"),
                "start_time": start_times.get(f),
                "block": (blocks or {}).get(f, 1),
                "features_detected": int(detected[f].sum()),
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
    drift: pd.DataFrame | None = None,
    agreement: dict | None = None,
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
    if roles.get("conditioning", 0):
        checks.append(Check("Conditioning injections", "left out of QC statistics", str(roles["conditioning"]),
                            "pass", "system-conditioning QCs at the start of the sequence"))
    not_in_sheet = int((inj["role"] == "not_in_sheet").sum())
    if not_in_sheet:
        checks.append(Check("Files missing from the sample sheet", "0", str(not_in_sheet), "warn",
                            "these injections were processed but have no role"))
    n_blocks = int(inj["block"].nunique()) if "block" in inj else 1
    if n_blocks > 1:
        checks.append(Check("Sequence blocks", "information", str(n_blocks), "pass",
                            f"pauses over {t.block_gap_hours:g} h split the run; RT and mass accuracy are judged within each block"))

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
                                       f"±{t.rt_tolerance_min} min from block median",
                                       m["max_rt_deviation_min"].reset_index(drop=True), t.rt_tolerance_min, comp))
        checks.append(_threshold_check("Internal standard mass accuracy (QC injections and standards)",
                                       f"±{t.mass_accuracy_ppm} ppm",
                                       m["max_abs_ppm_error"].reset_index(drop=True), t.mass_accuracy_ppm, comp))
        basis = str(m["qc_basis"].iloc[0]) if "qc_basis" in m else "pooled QC"
        c = _threshold_check("Internal standard peak area RSD",
                             f"< {t.is_area_rsd_pct:g}% over {basis} (needs ≥ {t.is_min_qc_detected} detected)",
                             m["area_rsd_pct"].reset_index(drop=True), t.is_area_rsd_pct, comp)
        if basis != "pooled QC" and c.status in ("pass", "fail"):
            c.status = "warn"
            c.detail = (c.detail + " (no pooled QC: fallback set, so a warning only)").strip()
        checks.append(c)
        if "qc_to_lowest_cal" in m and m["lowest_cal_area"].notna().any():
            ratio = m["qc_to_lowest_cal"]
            low = m.loc[ratio.notna() & (ratio < 1), "compound"].tolist()
            absent = m.loc[ratio.isna(), "compound"].tolist()
            n_in = int((ratio >= 1).sum())
            detail = []
            if low:
                detail.append("below the lowest calibration level: " + ", ".join(low))
            if absent:
                detail.append("not detected in pooled QCs: " + ", ".join(absent))
            if low or absent:
                detail.append("consider a larger spike or a smaller reconstitution volume")
            checks.append(Check("Internal standard signal in pooled QCs vs calibration range",
                                "at or above the lowest calibration level", f"{n_in} of {len(m)}",
                                "pass" if n_in == len(m) else "warn", "; ".join(detail)))
        if "sample_area_rsd_pct" in m:
            med = float(m["sample_area_rsd_pct"].median())
            checks.append(Check("Internal standard variation across samples (information)", "reported, not judged",
                                f"median RSD {med:.0f}%" if med == med else "not measured", "pass",
                                "soils differ in matrix effects and recovery; stage 5 corrects for this"))

    if "qc_rsd_pct" in feats:
        eligible = feats[feats["qc_eligible"]]
        frac = float((eligible["qc_rsd_pct"] < t.qc_feature_rsd_pct).mean()) if len(eligible) else np.nan
        checks.append(Check(f"Pooled-QC features with RSD < {t.qc_feature_rsd_pct:g}%",
                            f">= {t.qc_feature_pass_fraction:.0%} of features",
                            f"{frac:.1%} of {len(eligible)}" if len(eligible) else "no features",
                            _status(bool(frac >= t.qc_feature_pass_fraction)) if len(eligible) else "not_run",
                            f"features peak-picked in at least {t.qc_min_detect_fraction:.0%} of pooled QCs "
                            f"(gap-filled values not counted); {int((feats['qc_detected'] > 0).sum())} features "
                            "peak-picked in at least one pooled QC"))
    else:
        checks.append(Check("Pooled-QC feature reproducibility", "needs >= 2 pooled QCs", "-", "not_run", ""))

    if agreement:
        checks.append(Check("Pooled-QC difference from mean (Bland–Altman, as in Samrat et al. 2025)",
                            "information",
                            f"{agreement['within_fraction']:.1%} within {agreement['loa_low']:+.1f}% to "
                            f"{agreement['loa_high']:+.1f}%",
                            "pass", f"{agreement['n_features']} features × pooled QCs"))

    if drift is not None and not drift.empty:
        dp = float(drift["drift_pct"].iloc[0])
        checks.append(Check("Pooled-QC signal drift over the run", f"within ±{t.qc_drift_pct:g}%", f"{dp:+.1f}%",
                            "pass" if abs(dp) <= t.qc_drift_pct else "warn",
                            f"straight-line fit of summed signal of {int(drift['features_used'].iloc[0])} features "
                            "found in every pooled QC; drift can be corrected with the QCs in a later stage"))

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
    start_times = start_times or {}
    area = per_file_matrix(full, "area")
    mzs = per_file_matrix(full, "mz")
    rts = per_file_matrix(full, "rt")
    detected = detected_matrix(full, area) & (area.fillna(0) > 0)
    blocks = sequence_blocks(start_times, list(area.columns), t.block_gap_hours)
    ref_files = sheet.loc[sheet["role"] == "standard_cal", "file_id"].tolist()
    qc_files = sheet.loc[sheet["role"] == "pooled_qc", "file_id"].tolist()
    area_det = area.where(detected)
    matches = match_internal_standards(full, targets, polarity, area_det, ppm=t.match_ppm,
                                       rt_window=t.match_rt_window_min, reference_files=ref_files)
    matches = qc_features_for_standards(matches, full, detected, area, qc_files, ppm=t.match_ppm,
                                        window=t.is_qc_shift_window_min)
    is_summary, is_detail = internal_standard_metrics(
        matches, area, mzs, rts, sheet, detected=detected, blocks=blocks,
        min_qc_detected=t.is_min_qc_detected, min_qc_fraction=t.is_min_qc_fraction)
    feats = feature_quality(area, sheet, detected, t.qc_min_detect_fraction)
    inj = injection_summary(area, sheet, start_times, detected, blocks)
    drift = qc_drift(area, detected, inj, qc_files)
    agreement = qc_agreement(area, detected, qc_files, feats["qc_eligible"]) if "qc_eligible" in feats else {}
    checks = build_checks(is_summary, feats, inj, t, drift, agreement)
    return QCResult(checks=checks, internal_standards=is_summary, is_per_injection=is_detail,
                    injections=inj, features=feats, drift=drift, agreement=agreement)

#!/usr/bin/env python
"""Figure captions for the manuscript, with every number read from paper/ (same sources as the
figures). Writes paper/figures/captions.md and captions.json.

    python paper/figures/captions.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import panels_atlas as atlas  # noqa: E402
from common import FIG_DIR, read_table, read_yaml  # noqa: E402


def _rsd(pol):
    s = read_table(f"01_qc/{pol}/rsd_summary.tsv").set_index("mode")
    thr = read_yaml(f"01_qc/{pol}/provenance.yaml")["params"]["thresholds_from_run"]
    return s, thr


def build() -> dict[str, str]:
    c = {}
    sp, thr = _rsd("pos")
    sn, _ = _rsd("neg")
    pp, gf = sp.loc["peak_picked_only"], sp.loc["including_gap_filled"]
    elig = (f"≥ {int(pp.min_qc_values_required)} of {int(pp.n_pooled_qc)} (POS) and "
            f"≥ {int(sn.loc['peak_picked_only', 'min_qc_values_required'])} of "
            f"{int(sn.loc['peak_picked_only', 'n_pooled_qc'])} (NEG)")
    rsd_thr = thr["qc_feature_rsd_pct"]
    c["Figure 2"] = (
        "Integrated workflow for soil lipidomics. (A) Seven stages from Thermo raw files to "
        "Atlas-linked lipid profiles; stages 1–4 run as Nextflow processes, stages 5–7 (dashed) "
        "as library code and scripts in this version. The Soil Life Atlas enters at three "
        "stages: reference spectra and names during annotation (when available), class-level "
        "response during concentration estimation, and identity and occurrence during "
        "interpretation. (B) Pooled-QC relative standard deviation (RSD) per feature in the "
        "ClimGrass run, positive (POS) and negative (NEG) mode, computed from detected peaks only "
        "(dark) and including gap-filled values (grey); a feature is eligible when it has a "
        f"value in {elig} pooled QCs; dotted line, {rsd_thr:g} % RSD. Including gap-filled values "
        f"admits {gf.n_eligible / pp.n_eligible:.0f}× more POS features, of which "
        f"{gf.pct_of_eligible_le_30:g} % have an RSD ≤ {rsd_thr:g} %, compared with "
        f"{pp.pct_of_eligible_le_30:g} % of detected peaks. (C) Summary of the one-page report "
        "generated for the ClimGrass positive-mode run.")

    cov = read_table("06_atlas_match/pos/coverage.tsv").iloc[0]
    conf = read_table("02_annotation/pos/confidence.tsv")
    m = atlas._matched()
    n_match, n_ms2 = int(cov.soil_features_matched), int(cov.soil_features_with_ms2)
    n_unnamed = int((~m.named).sum())
    one = (m.n_groups == 1).mean() * 100
    det = read_table("06_atlas_match/pos/phylum_detection.tsv")
    units = det.groupby("phylum").n_units_phylum.first()
    c["Figure 3"] = (
        "What the Soil Life Atlas adds to a soil lipid dataset (ClimGrass, positive mode). "
        f"(A) Of {n_ms2:,} soil features with an MS/MS spectrum ({int(conf.n_features.sum()):,} "
        f"detected), {n_match:,} ({cov.pct_matched:.0f} %) match an Atlas reference spectrum "
        f"({atlas.match_settings()}); {n_unnamed:,} ({n_unnamed / n_match * 100:.0f} %) of the "
        "matched features carry no lipid annotation. (B) Number of organism groups in which the "
        f"matched reference feature was detected; {one:.0f} % of matches point to a single group. "
        "(C) Share of matched soil features whose reference occurs in each organism group (pale) "
        "and only in that group (solid); reference samples per group in brackets. Counts are not "
        "normalised for the number of reference samples. (D) For the eight most frequent classes "
        "among annotated matches, share of features whose reference occurs in each group (dot "
        f"area). Organism groups follow the Atlas release; occurrence uses {len(units)} analysed "
        f"phyla and {int(units.sum())} reference samples with detected features.")

    ex = read_table("06_atlas_match/examples/examples.tsv").set_index("role")
    b, r = ex.loc["broad"], ex.loc["restricted"]
    c["Figure 4"] = (
        "Reference-based interpretation of individual soil lipid features. (A, B) Mirror plots "
        "of a soil MS/MS spectrum (top) and its matched Atlas reference spectrum (bottom) for a "
        f"widely distributed lipid, {b.lipid_name} (cosine {b.cosine:.3f}, "
        f"{int(b.matched_peaks)} shared fragments; reference detected in "
        f"{int(b.n_ref_phyla)} phyla), and a restricted one, {r.lipid_name} (cosine "
        f"{r.cosine:.3f}, {int(r.matched_peaks)} shared fragments; "
        f"{int(r.n_ref_phyla)} phylum). The strongest shared fragments are labelled with their "
        "m/z; peaks beyond the precursor range are counted in the panel but not drawn. "
        f"(C) Detection of the two reference features across the {len(units)} analysed phyla, "
        "grouped by organism group; dot area is the share of reference samples of a phylum in "
        "which the lipid was detected (detected/analysed samples below each dot); open circles, "
        "not detected. Occurrence within the reference collection indicates potential sources, "
        "not proof of origin.")

    ts = read_table("08_application/pos/test_summary.tsv").iloc[0]
    c["Figure 5"] = (
        "Application to the ClimGrass climate and drought experiment. (A) Two-by-two design: "
        "ambient climate or +3 °C and +300 ppm CO₂, each with or without drought, three plots per "
        "treatment. The independent drying–rewetting validation is pre-registered and blinded and "
        "is described here only in general terms. (B) Lipids that change with drought "
        "(two-way ANOVA on log10 amounts per gram dry soil, Benjamini–Hochberg q < 0.05; "
        f"{int(ts.n_lipids_tested)} POS lipids tested), shown as log2 fold change, drought vs "
        "ambient, for both polarities. (C) Number of phyla in which the matched Atlas reference "
        "of each drought-responsive POS lipid was detected. (D) One lipid traced from its measured "
        "response to its Atlas occurrence. [Figure 5 is drawn once the stage-5 quantification "
        "is confirmed.]")

    s1, _ = _rsd("pos")
    ipi = {p: read_table(f"01_qc/{p}/is_per_injection.tsv") for p in ("pos", "neg")}
    inj = {p: read_table(f"01_qc/{p}/injections.tsv") for p in ("pos", "neg")}
    rt_tol, ppm_tol = thr["rt_tolerance_min"], thr["mass_accuracy_ppm"]
    within = []
    for p in ipi:
        d = ipi[p][ipi[p].detected.astype(bool)]
        within.append(((d.ppm_error.abs() <= ppm_tol).sum(), len(d)))
    pct_ppm = sum(a for a, _ in within) / sum(b_ for _, b_ in within) * 100
    c["Figure S1"] = (
        "Analytical quality of the ClimGrass run. (A) Injection sequences in positive "
        f"(n = {len(inj['pos'])}) and negative (n = {len(inj['neg'])}) mode, in two blocks: neat "
        "calibration standards first, then soils, pooled QCs, blanks and NIST SRM 1950. "
        "(B) Pooled-QC RSD per feature, detected peaks only vs including gap-filled values "
        "(as Fig. 2B). (C) Internal-standard retention time relative to the median of its "
        "sequence block, and mass error, for every injection; dashed lines, acceptance limits "
        f"(±{rt_tol:g} min, ±{ppm_tol:g} ppm); shaded, calibration block. {pct_ppm:.0f} % of mass "
        "errors lie within the limit; the pipeline applies the limits to pooled-QC and standard "
        "injections.")

    inter = read_table("02_annotation/pos/upset_intersections.tsv")
    single = inter[inter.n_sources == 1].n_features.sum() / inter.n_features.sum() * 100
    rt = read_table("02_annotation/pos/rt_model_PC.tsv")
    c["Figure S2"] = (
        "Lipid annotation of the ClimGrass soils. (A) Annotation sources and their exclusive "
        f"intersections (UpSet), POS and NEG; {single:.0f} % of annotated POS features are named "
        "by a single source. (B) Confidence levels: Gold (molecular species), Silver (sum "
        "composition) and Bronze (lipid class), with the number of annotated of all detected "
        "features. (C) The twelve largest lipid classes per polarity, as features and as unique "
        "lipids after collapsing adducts. (D) PC retention-time model: observed against fitted "
        f"retention time; crosses mark outliers ({int((rt.rt_check == 'outlier').sum())} of "
        f"{int(rt.rt_check.isin(['ok', 'outlier']).sum())} tested POS PC names); untested "
        "lipids, including ether PC, are counted in the panel but not drawn.")

    resp = read_table("03_response/response_by_class_adduct.tsv")
    st = read_table("03_response/ie_model_stats.tsv")
    best = st[st.best.astype(bool)].iloc[0]
    nist = read_table("05_nist/pos/nist_comparison.tsv")
    q1 = nist[nist.routes == "Q1"]
    c["Figure S3"] = (
        "Response and quantification routes. (A) Ionisation response of "
        f"{len(resp)} standard measurements as relative ionisation efficiency (RIE) relative to "
        "LPE, by class and adduct (Samrat et al., 2025, Table S10); log scale. (B) Measured "
        "against predicted RIE, leave-one-class-out; best model "
        f"{best.model.replace('_', ' ')} (R² = {best.r2:.2f}, RMSE = {best.rmse_log10:.2f} "
        f"log10 units, {int(best.n_standards)} POS standards in {int(best.n_classes)} classes). "
        "(C) Route assigned to each feature: Q1, class-matched internal standard; Q3, "
        "class-level response; Q4, relative signal only. Q2 (lipid-specific predicted response) "
        "is not used because cross-class prediction is weak (B). (D) Q1 estimates for NIST SRM "
        "1950 against consensus values (Bowden et al., 2017); solid line 1:1, dashed lines "
        f"2-fold; POS median estimate/consensus ratio {np.median(q1.ratio):.2f} (n = {len(q1)}). "
        "Q3 estimates are added once the stage-5 quantification is confirmed.")

    c["Supplementary note"] = (
        "[Figures S4 (internal-standard recovery against the neat calibration) and S5 (class "
        "amounts per gram dry soil) are added once the reconstitution volume and the stage-5 "
        "quantification are confirmed.]")
    return c


def main() -> int:
    c = build()
    (FIG_DIR / "captions.json").write_text(json.dumps(c, ensure_ascii=False, indent=1))
    md = ["# Figure captions (generated by captions.py from paper/ data)", ""]
    md += [f"**{k}.** {v}\n" if not k.startswith("Supp") else f"{v}\n" for k, v in c.items()]
    (FIG_DIR / "captions.md").write_text("\n".join(md))
    print("\n\n".join(f"{k}. {v}" for k, v in c.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

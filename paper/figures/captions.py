#!/usr/bin/env python
"""Figure captions for the manuscript. Every number comes from analysis.py's outputs
(figures/derived/numbers.tsv and derived tables); nothing is computed here.

    python paper/figures/captions.py     -> paper/figures/captions.md and captions.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import FIG_DIR, PAPER  # noqa: E402

DERIVED = FIG_DIR / "derived"


def _numbers() -> dict:
    t = pd.read_csv(DERIVED / "numbers.tsv", sep="\t", dtype=str)
    out = {}
    for k, v in zip(t.key, t.value):
        try:
            f = float(v)
            out[k] = int(f) if f.is_integer() and "." not in v else f
        except ValueError:
            out[k] = v
    return out


def _table(name: str) -> pd.DataFrame:
    return pd.read_csv(DERIVED / f"{name}.tsv", sep="\t")


def build() -> dict[str, str]:
    n = _numbers()
    c = {}
    s = _table("qc_rsd_summary").set_index(["pol", "mode"])
    pp, pn = s.loc[("pos", "peak_picked_only")], s.loc[("neg", "peak_picked_only")]
    elig = (f"≥ {int(pp.min_qc)} of {int(pp.n_qc)} (POS) and ≥ {int(pn.min_qc)} of "
            f"{int(pn.n_qc)} (NEG)")
    thr = n["qc_rsd_threshold_pct_pos"]
    c["Figure 1"] = (
        f"The Soil Life Atlas, {n['release_label']}. (A) Reference collection: samples "
        "measured in positive mode per phylum, grouped by organism group; hatched phyla were "
        f"collected but are not yet analysed ({n['release_phyla_analysed']} of "
        f"{n['release_phyla_collected']} phyla analysed; {n['release_genera']} genera; "
        f"{n['release_reference_samples_neg']} samples in negative mode). (B) From everything "
        "measured to the curated Atlas: consensus features aligned across batches, the "
        "quality-filtered matrix, phylum-enriched lipids (shaded; the Atlas), and of these the "
        "lipids with a usable MS/MS spectrum and with an accepted annotation (Gold, Silver or "
        "Bronze), in positive (dark) and negative (light) mode; log scale. Of "
        f"{n['release_features_considered']:,} features detected in at least two samples and "
        f"two batches, only {n['release_exclusive_features']} "
        f"({n['release_exclusive_pct']:g} %) occur in a single phylum: phyla are told apart "
        "mainly by shared lipids present in different amounts. (C) Each Atlas lipid carries an "
        "identity (m/z, retention time, MS/MS spectrum and name), its occurrence across phyla "
        "and organism groups, and a class-level ionisation response; the workflow uses these for "
        "annotation, interpretation and amounts (Fig. 2). Counts are those of the draft release "
        "and may change while the reference data are reprocessed.")

    c["Figure 2"] = (
        "Integrated workflow for soil lipidomics. (A) Seven stages from Thermo raw files to "
        "Atlas-linked lipid profiles; stages 1–4 run as Nextflow processes, stages 5–7 (dashed) "
        "as library code and scripts in this version (main tools: ThermoRawFileParser, MZmine 4.9, "
        "LipidSearch, MS-DIAL, GNPS, SIRIUS). The Soil Life Atlas enters at three "
        "stages: reference spectra and names during annotation (when available), class-level "
        "response during concentration estimation, and identity and occurrence during "
        "interpretation (solid arrows, used by default; dashed, when available). (B) Pooled-QC relative standard deviation (RSD) per feature in the "
        "ClimGrass run, positive (POS) and negative (NEG) mode, computed from detected peaks only "
        "(dark) and including gap-filled values (grey); a feature is eligible when it has a "
        f"value in {elig} pooled QCs; dotted line, {thr:g} % RSD. Including gap-filled values "
        f"admits {n['qc_gapfill_feature_ratio_pos']:.0f}× more POS features, of which "
        f"{n['qc_rsd_pct_le_threshold_including_gap_filled_pos']:g} % have an RSD ≤ {thr:g} %, "
        f"compared with {n['qc_rsd_pct_le_threshold_peak_picked_only_pos']:g} % of detected "
        "peaks. (C) Summary of the one-page report generated for the ClimGrass positive-mode run.")

    settings = (f"modified cosine ≥ {n['atlas_match_min_cosine']}, ≥ "
                f"{n['atlas_match_min_peaks']} shared fragments, precursor ± "
                f"{n['atlas_match_ppm']:g} ppm")
    c["Figure 3"] = (
        "What the Soil Life Atlas adds to a soil lipid dataset (ClimGrass, positive mode). "
        f"(A) Of {n['atlas_n_ms2']:,} soil features with an MS/MS spectrum "
        f"({n['atlas_n_features_all']:,} detected), {n['atlas_n_matched']:,} "
        f"({n['atlas_pct_matched']:.0f} %) match an Atlas reference spectrum ({settings}); "
        f"{n['atlas_n_matched_unnamed']:,} ({n['atlas_pct_unnamed_of_matched']:.0f} %) of the "
        "matched features carry no lipid annotation. (B) Number of organism groups in which the "
        f"matched reference feature was detected; {n['atlas_pct_single_group']:.0f} % of matches "
        "point to a single group. (C) Share of matched soil features whose reference occurs in "
        "each organism group (pale) and only in that group (solid); reference samples per group "
        "in brackets. Counts are not normalised for the number of reference samples. (D) For "
        "the eight most frequent classes among annotated matches, share of features whose "
        "reference occurs in each group (dot area). Organism groups follow the Atlas release; "
        f"occurrence uses {n['atlas_n_phyla_matching']} analysed phyla and "
        f"{n['atlas_n_ref_samples_matching']} reference samples with detected features.")

    ex = _table("interp_examples").set_index("role")
    b, r = ex.loc["broad"], ex.loc["restricted"]
    c["Figure 4"] = (
        "Reference-based interpretation of individual soil lipid features. (A, B) Mirror plots "
        "of a soil MS/MS spectrum (top) and its matched Atlas reference spectrum (bottom) for a "
        f"widely distributed lipid, {b.lipid_name} (cosine {b.cosine:.3f}, "
        f"{int(b.matched_peaks)} shared fragments; reference detected in "
        f"{int(b.n_ref_phyla)} phyla), and a restricted one, {r.lipid_name} (cosine "
        f"{r.cosine:.3f}, {int(r.matched_peaks)} shared fragments; "
        f"{int(r.n_ref_phyla)} phylum). The strongest shared fragments are labelled with their "
        "m/z; reference peaks beyond the precursor range are not drawn "
        f"({n['interp_hidden_peaks_broad']} for {b.lipid_name}, "
        f"{n['interp_hidden_peaks_restricted']} for {r.lipid_name}). "
        f"(C) Detection of the two reference features across the {n['atlas_n_phyla_matching']} "
        "analysed phyla, grouped by organism group; dot area is the share of reference samples "
        "of a phylum in which the lipid was detected; "
        "open circles, not detected. Occurrence within the reference collection indicates "
        "potential sources, not proof of origin.")

    c["Figure 5"] = (
        "Application to the ClimGrass climate and drought experiment. (A) Two-by-two design: "
        "ambient climate or +3 °C and +300 ppm CO₂, each with or without drought, three plots per "
        f"treatment and polarity ({n['app_soils_pos']} soils per polarity; in the future + drought "
        "treatment one plot was measured in positive and another in negative mode). An "
        "independent drying–rewetting validation is pre-registered and blinded and is not shown. "
        "(B) Lipids that change with drought (two-way ANOVA on log10 amounts per gram dry soil, "
        "Benjamini–Hochberg q < 0.05; "
        f"{n['app_n_responsive_pos']} of {n['app_n_lipids_tested_pos']} POS and "
        f"{n['app_n_responsive_neg']} of {n['app_n_lipids_tested_neg']} NEG lipids), shown as "
        "log2 fold change, drought vs ambient. (C) Number of phyla in which the matched Atlas reference of each "
        "drought-responsive POS lipid occurs; "
        f"{n['app_n_responsive_pos_matched']} of {n['app_n_responsive_pos']} have a match."
        + ("" if n["app_trace_lipid"] == "none" else
           f" (D) {n['app_trace_lipid']} traced from its measured amount to its Atlas occurrence: "
           "mean amount per treatment (ng per g dry soil, mean of plots; q for drought above) "
           "and the analysed phyla of each organism group in which its matched reference occurs "
           "(filled)."))

    c["Figure S1"] = (
        "Analytical quality of the ClimGrass run. (A) Injection sequences in positive "
        f"(n = {n['qc_n_injections_pos']}) and negative (n = {n['qc_n_injections_neg']}) mode, "
        "in two blocks: neat calibration standards first, then soils, pooled QCs, blanks and "
        "NIST SRM 1950. (B) Pooled-QC RSD per feature, detected peaks only vs including "
        "gap-filled values (as Fig. 2B). (C) Internal-standard retention time relative to the "
        "median of its sequence block, and mass error, for every injection; dashed lines, "
        f"acceptance limits (±{n['qc_rt_tolerance_min_pos']:g} min, "
        f"±{n['qc_mass_tolerance_ppm_pos']:g} ppm); shaded, calibration block. "
        f"{n['qc_is_pct_ppm_within']:.0f} % of mass errors lie within the limit; the pipeline "
        "applies the limits to pooled-QC and standard injections.")

    c["Figure S2"] = (
        "Lipid annotation of the ClimGrass soils. (A) Annotation sources and their exclusive "
        f"intersections (UpSet), POS and NEG; {n['ann_single_source_pct_pos']:.0f} % of "
        "annotated POS features are named by a single source. (B) Confidence levels: Gold "
        "(molecular species), Silver (sum composition) and Bronze (lipid class), with the "
        "number of annotated of all detected features. (C) The twelve largest lipid classes per "
        "polarity, as features and as unique lipids after collapsing adducts. (D) PC "
        "retention-time model: observed against fitted retention time; crosses mark outliers "
        f"({n['ann_rt_outliers_pos']} of {n['ann_rt_tested_pos']} tested POS PC names); untested "
        "lipids, including ether PC, are counted in the panel but not drawn.")

    c["Figure S3"] = (
        f"Response and quantification routes. (A) Ionisation response of "
        f"{n['resp_n_standards']} standard measurements as relative ionisation efficiency (RIE) "
        "relative to LPE, by class and adduct (Samrat et al., 2025, Table S10); log scale. "
        "(B) Measured against predicted RIE, leave-one-class-out; best model "
        f"{n['resp_ie_model'].replace('_', ' ')} (R² = {n['resp_ie_r2']:.2f}, RMSE = "
        f"{n['resp_ie_rmse_log10']:.2f} log10 units, {n['resp_ie_n_standards']} POS standards "
        f"in {n['resp_ie_n_classes']} classes). (C) Route assigned to each feature: Q1, "
        "class-matched internal standard; Q3, class-level response; Q4, relative signal only. Q2 "
        "(lipid-specific predicted response) is not used because cross-class prediction is weak "
        "(B). (D) Q1 estimates for NIST SRM 1950 against consensus values (Bowden et al., 2017); "
        "solid line 1:1, dashed lines 2-fold; POS median estimate/consensus ratio "
        f"{n['nist_q1_median_ratio_pos']:.2f} (n = {n['nist_q1_n_pos']}). Open circles, Q3 "
        f"estimates (NEG, n = {n['nist_q3_n_neg']}; median ratio "
        f"{n['nist_q3_median_ratio_neg']:.2f}, {n['nist_q3_pct_within_2fold_neg']:.0f} % within "
        f"2-fold; {n['nist_q3_n_pos']} in POS).")

    c["Figure S5"] = (
        "Lipid amounts in the ClimGrass soils. (A) Amount per lipid class in µg per g dry soil, "
        "summed over the lipids quantified in every soil (common lipid set), mean ± SD of soils, "
        "positive and negative mode; log scale. (B) The same per treatment (positive mode), "
        "mean ± SD of the three plots.")

    if n["figS4_shown"] == "True":
        c["Figure S4"] = (
            "Internal standards in soil extracts against the neat calibration (ClimGrass). (A) Median "
            "concentration in the soil extracts of each internal standard detected in soils "
            "(points; detected/analysed soils in brackets), read from a log–log calibration of the "
            "neat EquiSPLASH dilution series, against the calibration range (grey); "
            f"{n['qc_is_n_within_cal']} of {n['qc_is_n_detected']} "
            f"{'lies' if n['qc_is_n_within_cal'] == 1 else 'lie'} within the range. "
            "(B) Apparent recovery, the amount found (extract concentration × "
            f"{n['qc_recon_volume_ul']:g} µL reconstitution volume) relative to the 10 ng spiked "
            "before extraction, median and range over soils; dashed line, 100 %. Median values range "
            f"from {n['qc_recovery_median_min_pct']:.0f} % to {n['qc_recovery_median_max_pct']:.0f} %. "
            "Without post-extraction spikes, losses during extraction and matrix effects on the "
            "signal cannot be separated.")
    return c


def main() -> int:
    c = build()
    (FIG_DIR / "captions.json").write_text(json.dumps(c, ensure_ascii=False, indent=1))
    md = ["# Figure captions (generated by captions.py from analysis.py outputs)", ""]
    md += [f"**{k}.** {v}\n" if not k.startswith("Supp") else f"{v}\n" for k, v in c.items()]
    (FIG_DIR / "captions.md").write_text("\n".join(md))
    print(f"wrote {FIG_DIR / 'captions.md'} ({len(c)} captions; inputs under "
          f"{DERIVED.relative_to(PAPER)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

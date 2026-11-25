# Results for the manuscript (first draft, 2026-11-25)

Every value below comes from a recorded run. Pipeline run: ClimGrass stages 1–4, both
polarities, commit `2ed5d70`, LiSC `sla_workflow_dev/runs/climgrass_both_2ed5d70`
(finished 2026-11-25 12:10, exit 0). Each export folder has a `provenance.yaml`
(export commit, parameters, input sha256, command). Paths are relative to `paper/`.
POS = positive mode, NEG = negative mode.

**Not to be copied:** old draft numbers (90.6 %, 2,689, 743 …) were recomputed; where a
recomputed value happens to equal an old one (POS 90.6 %), the source is now the file named.

## Methods: bracketed items

| Manuscript item | Value / sentence | Source |
|---|---|---|
| §2.1 reference collection, replicates, alignment, consensus spectra | Not settled here (issue #7). Release 0.9 draft counts: 168 POS / 195 NEG reference samples, 19 phyla collected (16 analysed), 106 genera, 6 groups, 6 batches; 11,371 POS / 5,697 NEG Atlas lipids. Which count goes where (collected vs analysed vs matching denominators): `06_atlas_match/ATLAS_COUNTS.md`. | `~/Desktop/SoilLifeAtlas/release/release.json` (draft, SLA `fd7e820`) |
| §2.2 sample mass, spike, volumes | ~2 g freeze-dried soil (1.965–2.082 g; plot 50: 0.178 g); 10 µL of 1 ng µL⁻¹ EquiSPLASH = 10 ng per standard before extraction; NIST SRM 1950: 50 µL plasma per extraction batch. **Reconstitution volume, injection volume, column, gradient, MS settings: still to fill** (not found in Notion E5). | `00_inputs/climgrass_samples.csv`, `00_inputs/SOURCES.md` |
| §2.6 thresholds, ranking, flagged annotations | Unchanged from the pipeline at `2ed5d70` (see `lib/soillifeatlas/annotate.py`); flagged annotations are kept and marked (`flags` column), confidence is not lowered. | `02_annotation/*/annotations.tsv` |
| §2.7 recovery, volume factors, calibration fitting | R_app = m_found / m_spiked; m_found from a log–log fit of area on nominal ESTD-1 concentration (EquiSPLASH 10 ng µL⁻¹ at A1, 2-fold to A7), × reconstitution volume (**assumed 100 µL, to confirm**). No post-extraction spikes → matrix effect not separable. | `01_qc/*/is_vs_calibration.tsv` |
| §2.7 response units, adducts, redundant features, extrapolation | f = RIE(standard) / RIE(analyte), RIE = class mean per adduct from Samrat et al. (2025) Table S10, relative to LPE; [M+NH4]⁺ rows bridged to the LPE scale through PS, [M+HCOO]⁻ through HexCer; positive mode: standard at its measured adduct, analyte at its class's usual adduct ([M+NH4]⁺ for TG, DG, CE, PG, else [M+H]⁺); MG and PI without a standard of their class → Q4; Q1 may use a head-group family standard (e.g. LPC for PC, flagged `is_match = family`); class sums over lipids with an amount in every soil that has a standard; f limited to 0.1–10 and flagged when limited. Q3 uses the standard nearest in RT among standards detected in ≥ half of the soils. Adducts of one lipid (same name, RT within 0.1 min) are collapsed to one unique lipid; class sums use one feature per lipid. Only peak-picked areas are quantified. | `lib/soillifeatlas/quantify.py`, `04_quant/*/provenance.yaml` |
| §2.8 similarity estimator, weighting, unexplained share | Library: release reprocessing reference spectra (one per reference feature), limited to features detected in ≥ 1 reference sample. Spectral matching implemented (modified cosine, matchms greedy, top 50 peaks). The sample-level similarity/unexplained-share calculation is **not implemented**; keep as planned analysis. | `scripts/atlas_match.py` |
| §2.10 sample numbers, contrasts, criteria | 12 soils = 2 × 2 (drought × +3 °C/+300 ppm CO₂), 3 plots each (ambient 19, 22, 36; future 35, 40, 52; ambient+drought 28, 44, 49; future+drought 27, 43, 48). NEG: 48 failed (no IS), replaced by 50 (same treatment). Test: two-way ANOVA (type II) on log10 ng g⁻¹ per unique lipid, Benjamini–Hochberg per effect, q < 0.05; lipids quantified in ≥ 2 plots of every treatment. | `00_inputs/climgrass_samples.csv`, `08_application/*/provenance.yaml` |
| Data and code availability | Repository `soillifeatlas/nf-core-soillifeatlas`; results branch `paper/results-2026-11`. Zenodo DOI, release tag, licences: **open** (see PR). | — |

## Results to complete

### §3.2 Workflow
Stages 1–4 ran end to end in Nextflow on ClimGrass (both polarities). Stage 5 (amounts),
stage 6 (one-page report) and stage 7 (spectral matching) ran as library code and scripts,
not yet as Nextflow processes, so the text should not claim a seven-stage Nextflow run.
Report screenshot: `09_report/report_pos.html`. Real sample sheet for Fig. 2B:
`09_report/samplesheet_climgrass_pos.csv` (40 injections, in run order).

### §3.3 Analytical quality (ClimGrass)
- Sequence: POS 40 injections (12 soils, 9 pooled QCs, 1 solvent blank, 1 extraction blank, 3 NIST SRM 1950, 14 calibration) in two blocks: calibration (sequence OE12-1, 27–28 Mar 2023) and soils/QCs/SRM (OE12-3, 30–31 Mar 2023). NEG 39 (8 pooled QCs). → `01_qc/*/injections.tsv`
- Pooled-QC RSD, peak-picked only: POS median 13.8 %, 90.6 % of 6,513 eligible features ≤ 30 % (eligible = peak-picked in ≥ 7 of 9 QCs; 96,729 features in total). NEG median 9.9 %, 92.9 % of 2,128 (≥ 6 of 8; 20,280 features). → `01_qc/*/rsd_summary.tsv`
- Including gap-filled values: POS median 47.4 %, 16.0 % of 66,961 eligible ≤ 30 %; NEG 38.2 %, 37.1 % of 8,095. On the same features the effect is small (POS 13.8 vs 15.9 %), so the main effect of gap filling is admitting many more, poorly measured features. → same file, row `paired_same_features`
- Internal standards: identity, RT and mass accuracy confirmed for 12/12 (POS) and 9/9 (NEG) found standards; |mass error| median 0.70 ppm (POS, max 7.9) and 1.9 ppm (NEG, max 3.8). RT shifted between the calibration block and the soil block by a median −0.14 min (POS, n = 7 standards) and −0.09 min (NEG, n = 4). → `01_qc/*/is_by_block.tsv`, `is_per_injection.tsv`
- Standards in soil: detected in ≥ half of the soils for 4 of 12 POS standards (Cer, LPC, PE, SM) and 1 of 9 NEG standards (Cer); MG, PG, PI, PS, TG (POS) and all but Cer (NEG) not detected in soil. Against the neat calibration (POS): Cer within range; LPC, PC, SM, DG partly below; LPE, PE below the lowest level. NEG Cer partly below. → `01_qc/*/is_vs_calibration.tsv`, `04_quant/*/standards_used.tsv`
- Background: sample features with blank > 20 % of the sample median: POS 10,408 of 64,923 (16.0 %), NEG 468 of 15,193 (3.1 %). → `01_qc/*/blank_summary.tsv`
- Drift of summed pooled-QC signal over the run: POS +14.2 %, NEG +8.3 %. → `01_qc/*/drift.tsv`
- Couvillion (drying–rewetting): use only the existing blinded stages 1–4 report (Drive `06_results_and_reports/couvillion_blinded_2026-11-20`); not re-opened here (#11).

### §3.4 Chemical characterisation
- Annotated features: POS 2,692 of 96,729 (Gold 595, Silver 1,578, Bronze 519; 2,079 unique Gold/Silver lipids after collapsing adducts); NEG 1,013 of 20,280 (Gold 309, Silver 510, Bronze 194; 805 unique). → `02_annotation/*/confidence.tsv`, `summary.tsv`
- Sources (POS set sizes): LipidSearch 1,741, SIRIUS/CANOPUS 801, MZmine lipid rules 487, spectral library 174; largest exclusive intersections LipidSearch only 1,388, SIRIUS only 519, MZmine only 328, LipidSearch + SIRIUS 179. NEG: 698 / 286 / 262 / 27. → `02_annotation/*/upset_sets.tsv`, `upset_intersections.tsv`
- Classes (POS features / unique lipids): TG 720/696, DG 519/516, PC 478/440, Cer 168/165, PE 122/120, PI 115/111, SM 102/96, MG 97/96. NEG: FA 252/245, PE 143/143, Cer 106/105, PC 96/96, PI 88/87. → `02_annotation/*/classes.tsv`
- Chemical consistency: RT outliers POS 591 of 1,468 tested (TG alone 258; the TG fit itself is adequate — accepted TG R² 0.96, median |residual| 0.18 min — and 243 of the 246 diacyl/triacyl TG outliers elute earlier than predicted, 245 named by LipidSearch, so they mark misassigned names: `02_annotation/pos/rt_model_check.md`), NEG 7 of 367; implausible compositions POS 191, NEG 15. PC RT model (Fig. 4D): POS 177 within tolerance, 92 outliers, 60 untested (ether PC 149 untested); NEG in `rt_model_PC.tsv`. → `02_annotation/*/rt_model_PC.tsv`
- Atlas matching contribution: see §3.6 (separate from the annotation counts).

### §3.5 Response, recovery, routes
- Response by class and adduct (Fig. 5A): 109 standard measurements, Samrat et al. (2025) Table S10. → `03_response/response_by_class_adduct.tsv`
- Measured vs predicted ionisation efficiency (Fig. 5B), leave-one-class-out over 57 POS standards in 24 classes: best model random forest, R² = 0.16, RMSE = 1.27 log10 units (~18-fold). **The placeholder R² = 0.89 is not supported**; cross-class prediction is weak, so Q2 is not used in this draft. → `03_response/ie_model_stats.tsv`, `ie_loco_predictions.tsv`
- Apparent recovery R_app (Fig. 5C), median across soils, **with the assumed 100 µL reconstitution volume**: POS Cer 375 %, LPC 228 %, SM 166 %, PC 153 %, PE 95 %, LPE 75 % (1 soil); NEG Cer 109 %. Values > 100 % point to the volume assumption or to the calibration being run 3 days earlier, not to real gains; do not report until the volume is confirmed. → `01_qc/*/is_vs_calibration.tsv`
- Matrix effect: not computed (no post-extraction spikes).
- Routes: POS Q1 916 features (865 unique lipids), Q2 0, Q3 1,561 (1,523), Q4 94,252 (97.4 % of features; includes 97 MG and 118 PI features without a standard of their class, decisions 2026-11-25); NEG Q1 106 (105), Q2 0, Q3 907 (894), Q4 19,267. No Q3 response factor hits the 0.1–10 limit any more (median f TG 0.14, DG 1.08, PG 0.25). → `04_quant/*/route_counts.tsv`, `per_lipid.tsv`
- How Q1 standards were matched (`is_match`): POS Q1 class 428 features (415 unique lipids) + family 488 (450: 478 PC on the LPC standard, 10 LPE on the PE standard; the PC standard was detected in only 4 soils); NEG Q1 all class (106). Q3 always nearest RT. Family matches stay Q1 and are reported as such (decision 2026-11-25). → `04_quant/*/route_match_counts.tsv`, `per_lipid.tsv` (`is_match`)
- Review 2026-11-25, applied: (1) positive-mode Q3 analytes are scored at their class's usual adduct ([M+NH4]⁺ for TG, DG, CE, PG) instead of the standard's [M+H]⁺, which had produced f ≈ 50–100,000 and 89–98 % capped values; MG (no usable response) and PI (NIST SRM 1950 contradicts both adducts) are Q4 in positive mode (`04_quant/CAP_DIAGNOSIS.md`). (2) Class sums and totals use the common lipid set: lipids with an amount in every soil that has a detected standard — POS 279 of 2,388 unique lipids, 12 soils; NEG 177 of 999, 11 soils (plot 49 NEG has no standard) — because the PE standard was detected only in the drought soils and plot 19 (`04_quant/PLOT19_CHECK.md`, `common_set.tsv`).
- Class amounts per g dry soil (Fig. 5F), common lipid set, mean of soils: POS PC 28.9, DG 1.3, Cer 1.1, TG 1.0, LPC 0.5 µg g⁻¹; NEG FA 31.5, PC 8.4, Cer 1.2, PE 0.45 µg g⁻¹. These are sums over the common set, i.e. lower bounds of class totals. → `04_quant/*/class_by_sample.tsv`, `08_application/*/class_profiles.tsv`
- NIST SRM 1950 (Fig. 5E), three extracts of 50 µL plasma in the ClimGrass sequence, against the consensus MEDM values (nmol/mL, ≥ 5 laboratories; Bowden et al. 2017 / NISTIR 8185): POS 76 lipids, all Q1, Spearman ρ = 0.86 (log scale), median estimate/consensus ratio 1.02, 37 % within 2-fold. The 8 PIs that were Q3 are no longer quantified: scored at [M+NH4]⁺ they were 4–180-fold low, at [M+H]⁺ ~4-fold high, so PI is Q4 (decision 2026-11-25). TG/DG cannot be checked: none of their NIST estimates matches a consensus composition and the PE standard is not detected in the SRM extracts. NEG 69 lipids, ρ = 0.69, median ratio 0.34 (67 of 69 by Q3). The most abundant PCs (34:1, 36:2, 36:4, 38:4, 38:6) are 4–5-fold high; SM and LPC within 0.8–1.5-fold. → `05_nist/*/nist_summary.tsv`, `nist_comparison.tsv`

### §3.6 Interpretation and application
- Soil-to-Atlas matching (POS; 10 ppm, modified cosine ≥ 0.7, ≥ 4 fragments, 0.02 Da): 9,618 of 21,982 soil features with MS/MS (43.8 %) match a reference spectrum, hitting 7,597 reference features; library = 149,580 reference spectra detected in ≥ 1 of 161 reference samples (16 phyla) of the release reprocessing. Matched features occur in a median of 3 phyla. 1,773 matched features also carry an annotation (Gold 395, Silver 1,026, Bronze 352). NEG not run. → `06_atlas_match/pos/coverage.tsv`, `best_match.tsv`
- Examples (proposed, to confirm): broad, **LPC 18:0** (Gold, RT 5.44 min) – cosine 0.999, 20 fragments, reference feature in 13 phyla (e.g. Evosea 5/5, Arthropoda 9/11, Ascomycota 13/27 samples); restricted, **SM 34:3** (Silver, RT 11.92 min) – cosine 0.966, 17 fragments, Arthropoda only (7/11 samples). Mirror spectra and detection per phylum → `06_atlas_match/examples/`
- Atlas occurrence of the drought-responsive POS lipids: 15 of 18 match a reference spectrum (e.g. PC 30:1, 10 phyla; Cer d18:1_17:0, 5 phyla; Cer d17:0_15:0+Os, 2 phyla). → `06_atlas_match/examples/responsive_atlas_occurrence.tsv`
- Total lipid amounts (common lipid set; POS Q1–Q3, µg g⁻¹, mean ± SD of 3 plots): ambient 31.0 ± 3.8, future 33.6 ± 2.2, ambient+drought 35.3 ± 5.7, future+drought 35.0 ± 9.7; NEG 38.3 ± 6.9, 48.3 ± 10.1, 39.9 ± 5.3 (n = 2, plot 49 without standard), 42.3 ± 13.1 (includes plot 50, biased high). No clear treatment difference in totals. Plot 19 (34.0) is an ordinary ambient soil. → `08_application/*/total_by_treatment.tsv`
- Responsive lipids: of 393 POS and 514 NEG unique lipids tested (≥ 2 plots with an amount in every treatment), 18 POS and 15 NEG change with drought (q < 0.05; POS 13 down, 5 up; NEG 13 down, 2 up); none with future climate or the interaction. → `08_application/*/test_summary.tsv`, `responsive.tsv`
- Kept separate from the companion paper's kingdom/phylum decomposition (not used here).

## Placeholder register: status after this work

| Figure / panel | Now | Source |
|---|---|---|
| 1A | REAL as before (168 POS samples / 19 phyla collected); counts reconciled — occurrence statements use 16 analysed phyla | `06_atlas_match/ATLAS_COUNTS.md` |
| 1B, 1D | unchanged (issue #7, not in this pass) | — |
| 2A | still SHOWN AS BUILT for stages 5–7 in Nextflow; scripts exist | `scripts/` |
| 2B | **REAL** | `09_report/samplesheet_climgrass_pos.csv` |
| 2C | unchanged (#10) | — |
| 2D | **REAL** (minimal report) | `09_report/report_pos.html` |
| 3A ClimGrass | **REAL** | `01_qc/*/injections.tsv` |
| 3B | **REAL** | `01_qc/*/rsd_per_feature.tsv`, `rsd_summary.tsv` |
| 3C | **REAL** | `01_qc/*/is_per_injection.tsv`, `is_by_block.tsv` |
| 3D | **REAL areas**; amounts depend on the assumed reconstitution volume | `01_qc/*/is_vs_calibration.tsv` |
| 4A | **REAL** | `02_annotation/*/upset_*.tsv` |
| 4B ClimGrass | **REAL** | `02_annotation/*/confidence.tsv` |
| 4C ClimGrass | **REAL** | `02_annotation/*/classes.tsv` |
| 4D | **REAL** (PC) | `02_annotation/*/rt_model_PC.tsv` |
| 5A | **REAL** (RIE relative to LPE, not response per pmol) | `03_response/response_by_class_adduct.tsv` |
| 5B | **REAL** (R² 0.16, not 0.89) | `03_response/ie_model_stats.tsv` |
| 5C | PARTLY: depends on reconstitution volume | `01_qc/*/is_vs_calibration.tsv` |
| 5D | SCHEMATIC, now matches the code (option A) | `lib/soillifeatlas/quantify.py` |
| 5E | **REAL** (keep the panel) | `05_nist/*/nist_comparison.tsv` |
| 5F | **REAL** (semi-quantitative; common lipid set; Q3 adduct fixed) | `04_quant/*/class_by_sample.tsv` |
| 6A–C | **REAL** (examples proposed, to confirm) | `06_atlas_match/examples/` |
| 6D, 1D | unchanged (Explorer screenshots, #7) | — |
| 7A ClimGrass | REAL (unchanged) | — |
| 7A drying–rewetting | unchanged (#11, blinded) | — |
| 7B | **REAL** (common lipid set) | `08_application/*/class_profiles.tsv` |
| 7C | **REAL** | `08_application/*/responsive.tsv` |
| 7D | **REAL** | `06_atlas_match/examples/responsive_atlas_occurrence.tsv` |
| 7E | REAL pieces exist (route, amount, match, occurrence); pick one lipid | `08_application/pos/responsive.tsv` + above |
| Graphical abstract bars | unchanged | — |

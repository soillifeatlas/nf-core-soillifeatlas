# Figure status (2026-11-25)

Built by `python paper/figures/build_all.py` on branch `paper/figures-2026-11` at `cd726c1` (paper/figures has uncommitted changes).
Every plotted number is read from a file under `paper/` by the figure script; nothing is typed by hand.

Status: **REAL** = drawn from recorded data · **SCHEMATIC** = cropped from the 2026-11-24 draft, no data · **PENDING** = grey box, waits for the named item · **UNCHANGED** = not rebuilt in this pass.

Data SHA = last commit that changed the file in this repository; 'export' = `export_git_sha` in the folder's provenance.yaml.

Gates (gates.yaml): `reconstitution_volume_confirmed: False`, `stage5_fix_confirmed: False`, `fig7E_lipid: None`

| Fig | Panel | Status | Source file(s) under paper/ | Data SHA | Note |
|---|---|---|---|---|---|
| 1 | A–D | **UNCHANGED** | – | – | not rebuilt in this pass (issue #7) |
| 2 | A | **SCHEMATIC** | `figures/schematics/source/Fig2_workflow_2026-11-24.png` | cd726c1 | crop of Fig2_workflow_2026-11-24.png; stages 5–7 overlaid dashed ('script, not yet in Nextflow') |
| 2 | B | **REAL** | `01_qc/pos/injections.tsv`<br>`09_report/samplesheet_climgrass_pos.csv` | 756333b (export ac355b6)<br>33eebd2 (export f476c4c) | first two calibration rows and the first nine soil-block rows shown |
| 2 | C | **SCHEMATIC** | `figures/schematics/source/Fig2_workflow_2026-11-24.png` | cd726c1 | crop of Fig2_workflow_2026-11-24.png, unchanged (#10) |
| 2 | D | **REAL** | `09_report/report_pos.html`<br>`figures/rendered/report_pos.png` | 24f98ec (export f476c4c)<br>cd726c1 | headless-Chrome render of report_pos.html (top section); the 'Lipids with an amount' tile is stage-5 output |
| 3 | A | **REAL** | `01_qc/neg/injections.tsv`<br>`01_qc/pos/injections.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6) | POS and NEG sequences, both blocks |
| 3 | B | **REAL** | `01_qc/neg/provenance.yaml`<br>`01_qc/neg/rsd_per_feature.tsv`<br>`01_qc/neg/rsd_summary.tsv`<br>`01_qc/pos/provenance.yaml`<br>`01_qc/pos/rsd_per_feature.tsv`<br>`01_qc/pos/rsd_summary.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6) | ECDF of RSD, x-axis cut at 150 %; n and % ≤ 30 % from rsd_summary.tsv |
| 3 | C | **REAL** | `01_qc/neg/injections.tsv`<br>`01_qc/neg/is_by_block.tsv`<br>`01_qc/neg/is_per_injection.tsv`<br>`01_qc/neg/provenance.yaml`<br>`01_qc/pos/injections.tsv`<br>`01_qc/pos/is_by_block.tsv`<br>`01_qc/pos/is_per_injection.tsv`<br>`01_qc/pos/provenance.yaml` | 756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6) | detected IS in all injections, RT relative to the recorded block median (is_by_block.tsv); tolerances from provenance.yaml; dates in Europe/Vienna time |
| 3 | D | **PENDING** | `01_qc/neg/is_vs_calibration.tsv`<br>`01_qc/pos/is_vs_calibration.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6) | gated on reconstitution_volume_confirmed |
| 4 | A | **REAL** | `02_annotation/neg/upset_intersections.tsv`<br>`02_annotation/neg/upset_sets.tsv`<br>`02_annotation/pos/upset_intersections.tsv`<br>`02_annotation/pos/upset_sets.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export 72ba0e8)<br>756333b (export 72ba0e8) | exclusive intersections, POS and NEG |
| 4 | B | **REAL** | `02_annotation/neg/confidence.tsv`<br>`02_annotation/pos/confidence.tsv` | 756333b (export ac355b6)<br>756333b (export 72ba0e8) | ClimGrass only; drying–rewetting values omitted (not in a recorded paper/ file) |
| 4 | C | **REAL** | `02_annotation/neg/classes.tsv`<br>`02_annotation/pos/classes.tsv` | 756333b (export ac355b6)<br>756333b (export 72ba0e8) | 12 largest classes per polarity by feature count |
| 4 | D | **REAL** | `02_annotation/neg/rt_model_PC.tsv`<br>`02_annotation/pos/rt_model_PC.tsv` | 756333b (export ac355b6)<br>90a4383 (export 72ba0e8) | points with an RT check; untested (incl. ether PC) counted, not drawn |
| 5 | A | **REAL** | `03_response/response_by_class_adduct.tsv` | 57fa0ce (export 8012e4c) | axis: RIE relative to LPE (not response per pmol) |
| 5 | B | **REAL** | `03_response/ie_loco_predictions.tsv`<br>`03_response/ie_model_stats.tsv` | 57fa0ce (export 8012e4c)<br>57fa0ce (export 8012e4c) | R², RMSE from ie_model_stats.tsv (best model) |
| 5 | C | **PENDING** | `01_qc/neg/is_vs_calibration.tsv`<br>`01_qc/pos/is_vs_calibration.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6) | gated on reconstitution_volume_confirmed |
| 5 | D | **SCHEMATIC** | `figures/schematics/source/Fig5_quantification_2026-11-24.png` | cd726c1 | crop of Fig5_quantification_2026-11-24.png, unchanged; Q2 drawn but not used in this run (0 features) |
| 5 | E | **REAL** | `05_nist/neg/nist_comparison.tsv`<br>`05_nist/pos/nist_comparison.tsv` | e5dce35 (export f476c4c)<br>24f98ec (export f476c4c) | Q1 lipids only; Q3 points gated on stage5_fix_confirmed (caption: Q3 pending) |
| 5 | F | **PENDING** | `04_quant/neg/class_by_sample.tsv`<br>`04_quant/pos/class_by_sample.tsv` | 24f98ec<br>24f98ec | gated on stage5_fix_confirmed |
| 6 | A | **REAL** | `06_atlas_match/examples/examples.tsv`<br>`06_atlas_match/examples/mirror_21718.tsv` | a288361 (export 24f98ec)<br>a288361 (export 24f98ec) | LPC 18:0; examples proposed, to confirm |
| 6 | B | **REAL** | `06_atlas_match/examples/examples.tsv`<br>`06_atlas_match/examples/mirror_41748.tsv` | a288361 (export 24f98ec)<br>a288361 (export 24f98ec) | SM 34:3; examples proposed, to confirm |
| 6 | C | **REAL** | `06_atlas_match/examples/examples.tsv`<br>`06_atlas_match/examples/examples_phylum_detection.tsv` | a288361 (export 24f98ec)<br>a288361 (export 24f98ec) | only phyla where either reference feature is detected are listed |
| 6 | D | **PENDING** | – | – | Explorer screenshot, issue #7 |
| 7 | A | **REAL** | `00_inputs/climgrass_samples.csv` | def3188 | ClimGrass design from climgrass_samples.csv; drying–rewetting as a generic text box (blinded, #11) |
| 7 | B | **PENDING** | `08_application/neg/class_profiles.tsv`<br>`08_application/pos/class_profiles.tsv` | 24f98ec (export f476c4c)<br>24f98ec (export f476c4c) | gated on stage5_fix_confirmed |
| 7 | C | **PENDING** | `08_application/neg/responsive.tsv`<br>`08_application/neg/test_summary.tsv`<br>`08_application/pos/responsive.tsv`<br>`08_application/pos/test_summary.tsv` | b711a7f (export f476c4c)<br>b711a7f (export f476c4c)<br>24f98ec (export f476c4c)<br>24f98ec (export f476c4c) | gated on stage5_fix_confirmed |
| 7 | D | **PENDING** | `06_atlas_match/examples/responsive_atlas_occurrence.tsv` | ed8495e (export 24f98ec) | gated on stage5_fix_confirmed |
| 7 | E | **PENDING** | `06_atlas_match/examples/responsive_atlas_occurrence.tsv`<br>`08_application/pos/responsive.tsv` | ed8495e (export 24f98ec)<br>24f98ec (export f476c4c) | gated on stage5_fix_confirmed and fig7E_lipid |
| GA | – | **UNCHANGED** | – | – | graphical abstract not rebuilt (issue #7) |

Outputs:

- Fig 2: `figures/output/Fig2_workflow_2026-11-25.pdf`, `figures/output/Fig2_workflow_2026-11-25.png`
- Fig 3: `figures/output/Fig3_quality_2026-11-25.pdf`, `figures/output/Fig3_quality_2026-11-25.png`
- Fig 4: `figures/output/Fig4_annotation_2026-11-25.pdf`, `figures/output/Fig4_annotation_2026-11-25.png`
- Fig 5: `figures/output/Fig5_quantification_2026-11-25.pdf`, `figures/output/Fig5_quantification_2026-11-25.png`
- Fig 6: `figures/output/Fig6_interpretation_2026-11-25.pdf`, `figures/output/Fig6_interpretation_2026-11-25.png`
- Fig 7: `figures/output/Fig7_application_2026-11-25.pdf`, `figures/output/Fig7_application_2026-11-25.png`

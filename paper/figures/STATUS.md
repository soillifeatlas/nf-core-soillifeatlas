# Figure status (2026-11-25)

Built by `python paper/figures/build_all.py` on branch `paper/figures-2026-11` at `42a902e`.
Every plotted number is read from a file under `paper/` by the figure script; nothing is typed by hand.

Main text: Fig 1 Atlas, Fig 2 what the Atlas adds to a soil dataset, Fig 3 workflow, Fig 4 interpreting single features, Fig 5 ClimGrass demonstration. Supplement: S1 quality, S2 annotation, S3 response and routes, S4 recovery, S5 amounts.

Status: **REAL** = drawn from recorded data · **SCHEMATIC** = cropped from the 2026-11-24 draft, no data · **NOT BUILT** = figure not drawn yet, waits for the named switch in gates.yaml · **UNCHANGED** = not rebuilt in this pass.

Data SHA = last commit that changed the file in this repository; 'export' = `export_git_sha` in the folder's provenance.yaml.

Gates (gates.yaml): `reconstitution_volume_confirmed: False`, `stage5_fix_confirmed: False`, `fig7E_lipid: None`

| Fig | Panel | Status | Source file(s) under paper/ | Data SHA | Note |
|---|---|---|---|---|---|
| GA graphical_abstract | – | **REAL** | `02_annotation/pos/annotations.tsv`<br>`06_atlas_match/pos/best_match.tsv`<br>`06_atlas_match/pos/coverage.tsv`<br>`figures/inputs/atlas_groups.csv` | 756333b (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>3e25938 | Atlas group counts from the release export (figures/inputs/); match numbers from 06_atlas_match/pos |
| 1 atlas | A–D | **UNCHANGED** | – | – | Atlas overview; not rebuilt here (issue #7, release data outside paper/) |
| 2 atlas_matching | A | **REAL** | `02_annotation/pos/annotations.tsv`<br>`02_annotation/pos/confidence.tsv`<br>`06_atlas_match/pos/best_match.tsv`<br>`06_atlas_match/pos/coverage.tsv`<br>`06_atlas_match/pos/provenance.yaml`<br>`figures/inputs/atlas_groups.csv` | 756333b (export 72ba0e8)<br>756333b (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>3e25938 | POS; modified cosine ≥ 0.7, ≥ 4 shared fragments, precursor ± 10 ppm |
|  | B | **REAL** | `02_annotation/pos/annotations.tsv`<br>`06_atlas_match/pos/best_match.tsv`<br>`figures/inputs/atlas_groups.csv` | 756333b (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>3e25938 | organism groups from the release (figures/inputs/atlas_groups.csv) |
|  | C | **REAL** | `02_annotation/pos/annotations.tsv`<br>`06_atlas_match/pos/best_match.tsv`<br>`06_atlas_match/pos/phylum_detection.tsv`<br>`figures/inputs/atlas_groups.csv` | 756333b (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>3e25938 | % of matched features; reference samples per group in brackets (matching denominators) |
|  | D | **REAL** | `02_annotation/pos/annotations.tsv`<br>`06_atlas_match/pos/best_match.tsv`<br>`figures/inputs/atlas_groups.csv` | 756333b (export 72ba0e8)<br>55e0087 (export 72ba0e8)<br>3e25938 | 8 most frequent classes among named matches |
| 3 workflow | A | **SCHEMATIC** | – | – | drawn natively; stage and tool text from the 2026-11-24 draft |
|  | B | **REAL** | `01_qc/neg/provenance.yaml`<br>`01_qc/neg/rsd_per_feature.tsv`<br>`01_qc/neg/rsd_summary.tsv`<br>`01_qc/pos/provenance.yaml`<br>`01_qc/pos/rsd_per_feature.tsv`<br>`01_qc/pos/rsd_summary.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6) | pooled-QC RSD per feature; eligibility rule in rsd_summary.tsv |
|  | C | **REAL** | `09_report/report_pos.html`<br>`figures/rendered/report_pos.png` | 24f98ec (export f476c4c)<br>f12f08c | summary tiles of report_pos.html; the 'lipids with an amount' tile is stage-5 output |
| 4 interpretation | A | **REAL** | `06_atlas_match/examples/examples.tsv`<br>`06_atlas_match/examples/mirror_21718.tsv` | a288361 (export 24f98ec)<br>a288361 (export 24f98ec) | zoomed to m/z 40–precursor+25; example proposed, to confirm |
|  | B | **REAL** | `06_atlas_match/examples/examples.tsv`<br>`06_atlas_match/examples/mirror_41748.tsv` | a288361 (export 24f98ec)<br>a288361 (export 24f98ec) | zoomed; one reference peak above the range noted in the panel |
|  | C | **REAL** | `06_atlas_match/examples/examples.tsv`<br>`06_atlas_match/examples/examples_phylum_detection.tsv`<br>`06_atlas_match/pos/phylum_detection.tsv`<br>`figures/inputs/atlas_groups.csv` | a288361 (export 24f98ec)<br>a288361 (export 24f98ec)<br>55e0087 (export 72ba0e8)<br>3e25938 | all 16 analysed phyla, grouped by organism group |
| 5 application | A | **NOT BUILT** | – | – | waits for stage5_fix_confirmed |
|  | B | **NOT BUILT** | – | – | waits for stage5_fix_confirmed |
|  | C | **NOT BUILT** | – | – | waits for stage5_fix_confirmed |
| S1 quality | A | **REAL** | `01_qc/neg/injections.tsv`<br>`01_qc/pos/injections.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6) | POS and NEG sequences, both blocks |
|  | B | **REAL** | `01_qc/neg/provenance.yaml`<br>`01_qc/neg/rsd_per_feature.tsv`<br>`01_qc/neg/rsd_summary.tsv`<br>`01_qc/pos/provenance.yaml`<br>`01_qc/pos/rsd_per_feature.tsv`<br>`01_qc/pos/rsd_summary.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6) | ECDF, x cut at 150 %; n and % ≤ 30 % from rsd_summary.tsv |
|  | C | **REAL** | `01_qc/neg/injections.tsv`<br>`01_qc/neg/is_by_block.tsv`<br>`01_qc/neg/is_per_injection.tsv`<br>`01_qc/neg/provenance.yaml`<br>`01_qc/pos/injections.tsv`<br>`01_qc/pos/is_by_block.tsv`<br>`01_qc/pos/is_per_injection.tsv`<br>`01_qc/pos/provenance.yaml` | 756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export ac355b6) | RT relative to recorded block median; tolerances from provenance.yaml |
| S2 annotation | A | **REAL** | `02_annotation/neg/upset_intersections.tsv`<br>`02_annotation/neg/upset_sets.tsv`<br>`02_annotation/pos/upset_intersections.tsv`<br>`02_annotation/pos/upset_sets.tsv` | 756333b (export ac355b6)<br>756333b (export ac355b6)<br>756333b (export 72ba0e8)<br>756333b (export 72ba0e8) | exclusive intersections, POS and NEG |
|  | B | **REAL** | `02_annotation/neg/confidence.tsv`<br>`02_annotation/pos/confidence.tsv` | 756333b (export ac355b6)<br>756333b (export 72ba0e8) | ClimGrass only |
|  | C | **REAL** | `02_annotation/neg/classes.tsv`<br>`02_annotation/pos/classes.tsv` | 756333b (export ac355b6)<br>756333b (export 72ba0e8) | 12 largest classes per polarity |
|  | D | **REAL** | `02_annotation/neg/rt_model_PC.tsv`<br>`02_annotation/pos/rt_model_PC.tsv` | 756333b (export ac355b6)<br>90a4383 (export 72ba0e8) | PC; untested lipids counted, not drawn |
| S3 response_routes | A | **REAL** | `03_response/response_by_class_adduct.tsv` | 57fa0ce (export 8012e4c) | RIE relative to LPE |
|  | B | **REAL** | `03_response/ie_loco_predictions.tsv`<br>`03_response/ie_model_stats.tsv` | 57fa0ce (export 8012e4c)<br>57fa0ce (export 8012e4c) | leave-one-class-out; statistics from ie_model_stats.tsv |
|  | C | **SCHEMATIC** | – | – | drawn natively; Q2 not used in this run |
|  | D | **REAL** | `05_nist/neg/nist_comparison.tsv`<br>`05_nist/pos/nist_comparison.tsv` | e5dce35 (export f476c4c)<br>24f98ec (export f476c4c) | Q1 only; Q3 waits for stage5_fix_confirmed |
| S4 recovery | A | **NOT BUILT** | – | – | waits for reconstitution_volume_confirmed |
|  | B | **NOT BUILT** | – | – | waits for reconstitution_volume_confirmed |
| S5 amounts | A | **NOT BUILT** | – | – | waits for stage5_fix_confirmed |
|  | B | **NOT BUILT** | – | – | waits for stage5_fix_confirmed |

Outputs:

- Fig GA: `figures/output/FigGA_graphical_abstract_2026-11-25.pdf`, `figures/output/FigGA_graphical_abstract_2026-11-25.png`
- Fig 2: `figures/output/Fig2_atlas_matching_2026-11-25.pdf`, `figures/output/Fig2_atlas_matching_2026-11-25.png`
- Fig 3: `figures/output/Fig3_workflow_2026-11-25.pdf`, `figures/output/Fig3_workflow_2026-11-25.png`
- Fig 4: `figures/output/Fig4_interpretation_2026-11-25.pdf`, `figures/output/Fig4_interpretation_2026-11-25.png`
- Fig S1: `figures/output/FigS1_quality_2026-11-25.pdf`, `figures/output/FigS1_quality_2026-11-25.png`
- Fig S2: `figures/output/FigS2_annotation_2026-11-25.pdf`, `figures/output/FigS2_annotation_2026-11-25.png`
- Fig S3: `figures/output/FigS3_response_routes_2026-11-25.pdf`, `figures/output/FigS3_response_routes_2026-11-25.png`

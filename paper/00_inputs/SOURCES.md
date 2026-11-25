# Inputs for the ClimGrass paper runs

- `climgrass_samples.csv`: treatment and dry soil mass from the Notion page
  "Extraction-5-climgrass" (extraction E5, weighed soil table, "Factorial" list), read 2026-11-25.
  Treatment codes there: C0T0 = ambient, C2T2 = +300 ppm CO2 and +3 °C ("future").
- Spike: 10 µL of EquiSPLASH (EQ102, 10 ng/µL diluted 1:10 to 1 ng/µL) per sample
  = 10 ng per internal standard (same Notion page).
- SRM: NIST SRM 1950, 50 µL plasma, extracted with E5 (same Notion page).
- Soil 48 NEG failed (no IS); soil 50 replaces it in negative mode only (decision 2026-10-13,
  Paper 2 negative-mode reprocessing). Soil 50 had 0.178 g, so its per-gram values are biased high.
- Calibration (ESTD1 A1–A7 and duplicates) was acquired 27–28 Mar 2023 in sequence OE12-1;
  soils, pooled QCs and SRM on 30–31 Mar 2023 in sequence OE12-3 (raw-file creation dates).
- Raw files: LiSC `/lisc/data/work/ter/rahul/4-LC-HRMS-Orbitrap/Orbitrap-data/`.
- `rie_table_s10.csv`: Samrat et al. (2025), Soil Biology & Biochemistry, Supplementary Table S10
  (109 standard measurements; class, species, adduct, RT, RIE relative to LPE). Copied from
  `~/Desktop/SoilMass/analysis/analysis/analysis-19/00_inputs/` (sha256 056da89dfa158615…, same as the SLA repo copy).
  Stage 5 uses unfloored class means per adduct; [M+NH4]+ rows bridged to the LPE scale via PS, [M+HCOO]- via HexCer; f capped to 0.1–10 and flagged (decision 2026-11-25).
- `nist_srm1950_consensus_nistir8185.csv`: consensus estimates (nmol/mL; DSL, VR, MEDM with uncertainties, COD)
  for lipids measured by >= 5 laboratories, Tables 2–6 of NISTIR 8185 (Bowden et al. 2017, J Lipid Res 58:2275),
  parsed from https://nvlpubs.nist.gov/nistpubs/ir/2017/NIST.IR.8185.pdf (downloaded 2026-11-25) with
  `scripts/nist_compare.py --parse`. 335 rows parsed (the paper reports 339; a few wrapped rows were not parsed).

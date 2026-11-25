# Plot 19 (ambient): why its total was 477 µg g⁻¹ (review item 4)

2026-11-25. POS, one feature per unique lipid. Data: `pos/estimates_long.tsv`,
`cap_options_pos/per_lipid_f_option_A.tsv`. Stage 5 unchanged; options as in CAP_DIAGNOSIS.md.

## Cause: which standards were detected, not the soil

The PE standard (PE 15:0–18:1-d7) was peak-picked in plot 19 and in all six drought soils, and
in none of the other five non-drought soils. Q3 lipids whose nearest standard is PE (438 TG,
195 DG, …) therefore get an amount only in those seven soils; elsewhere they are missing, and
class sums treat missing as absent.

| plot | treatment | PE standard detected | lipids with an amount | TG with an amount |
|---|---|---|---|---|
| 22, 36 | ambient | no | 547, 573 | 77, 81 |
| 35, 40, 52 | future | no | 912, 896, 918 | 185, 181, 168 |
| **19** | **ambient** | **yes** | **1,168** | **319** |
| 27, 43, 48 | future + drought | yes | 1,499, 1,163, 1,523 | 426, 331, 431 |
| 28, 44, 49 | ambient + drought | yes | 1,487, 1,517, 1,150 | 401, 409, 325 |

Plot 19 by class (µg g⁻¹, current → option A): TG 344 → 6.0 (plots 22/36: 13 / 11 → 0.8 / 0.8),
DG 50 → 3.6, PE 11 → 11 (0 in 22/36). PC (Q1, LPC standard, detected everywhere) is similar
across the three plots: 33, 33, 25.

## Does it remain an outlier?

| plot 19 vs 22 / 36 (µg g⁻¹) | current | A | A + MG→Q4 | B |
|---|---|---|---|---|
| all quantified lipids | 477 vs 70 / 57 | 77 vs 49 / 40 | 67 vs 38 / 31 | 56 vs 35 / 28 |
| **lipids with an amount in all 12 soils (n = 326)** | 62 vs 60 / 49 | 43 vs 43 / 35 | 34 vs 33 / 27 | 31 vs 30 / 25 |

Options A and B shrink the gap (because TG's f drops) but plot 19 stays ~1.6–2.2× the other
ambient plots. Only summing over the same lipids in every soil removes it: then plot 19 is an
ordinary ambient soil.

## The same artefact inflates the drought effect

Treatment means (µg g⁻¹, mean ± SD of 3 plots) over the 326 lipids with an amount in all 12 soils:

| | current | A | A + MG→Q4 | B |
|---|---|---|---|---|
| ambient | 57 ± 7 | 41 ± 4 | 31 ± 4 | 29 ± 4 |
| future | 61 ± 5 | 43 ± 3 | 34 ± 2 | 31 ± 2 |
| ambient + drought | 71 ± 18 | 47 ± 8 | 36 ± 6 | 33 ± 5 |
| future + drought | 72 ± 20 | 46 ± 13 | 35 ± 10 | 32 ± 9 |

With a common lipid set, drought raises totals by ~1.1–1.2×, not ~3×. Class sums and Fig. 7B
should use a common lipid set (or impute missing standards) — a definition change, not applied.

## Applied (decision 2026-11-25: common lipid set, with A + MG→Q4)

POS sums now use the 287 lipids with an amount in all 12 soils. Plot 19 = 34.1 µg g⁻¹ vs 32.5
(plot 22) and 26.8 (plot 36): no longer an outlier. Treatment means: ambient 31.1 ± 3.8, future
33.7 ± 2.2, ambient + drought 35.5 ± 5.7, future + drought 35.1 ± 9.7. NEG uses 177 lipids over
11 soils (plot 49 has no detected standard in negative mode).

After PI → Q4 (same day): common set 279 lipids; plot 19 = 34.0 vs 32.4 / 26.7; treatment means 31.0 / 33.6 / 35.3 / 35.0 µg g⁻¹ (ambient / future / ambient + drought / future + drought).

# Why ~98 % of positive-mode TG, DG, MG, PI and PG features hit f = 10

2026-11-25, review item 1. Data: `cap_options_pos/` (script `paper/scripts/cap_options.py`).
Nothing in stage 5 has been changed yet.

## Diagnosis

The PS bridge is **not** the cause. The cause is which adduct the analyte is scored at.

`quantify._f` scores the analyte at the **standard's adduct**. Three of the four standards
used in positive mode are [M+H]⁺ (PE, SM, LPC; Cer is [M+H−H2O]⁺, which has no table row).
So a TG feature quantified against PE gets
f = RIE(PE, [M+H]⁺) / RIE(TG, [M+H]⁺) = 13.17 / 0.00605 = 2,177 (`f_source` column),
although TG, DG, MG, PI and PG ionise mainly as [M+NH4]⁺ in this method, and Table S10
shows how weak their [M+H]⁺ response is.

RIE relative to LPE (Samrat et al. 2025, Table S10, class means; [M+NH4]⁺ rows in the table are
relative to PS = 1, bridged here ×7.024 = RIE(PS,[M+H]⁺) / RIE(PS,[M+NH4]⁺)):

| Class | [M+H]⁺ | [M+NH4]⁺ as in table (PS = 1) | [M+NH4]⁺ bridged (LPE scale) | Q3 features | capped now | median f uncapped |
|---|---|---|---|---|---|---|
| TG | 0.0061 | 13.38 | 93.95 | 720 | 98.1 % | 2,177 |
| DG | 0.0323 | 4.09 | 28.71 | 519 | 97.3 % | 963 |
| MG | 0.00048 | 0.0032 | 0.022 | 97 | 97.9 % | 100,042 |
| PI | 0.264 | 11.17 | 78.45 | 115 | 90.4 % | 50 |
| PG | 0.129 | 17.95 | 126.0 | 63 | 88.9 % | 242 |
| PE (std) | 13.17 | — | — | | | |
| SM (std) | 31.13 | — | — | | | |
| LPC (std) | 976.5 | — | — | | | |

The bridge values agree with the table: they are the table's own [M+NH4]⁺ numbers × 7.02.

## Options (not applied)

**A. Fix:** score each analyte at its class's usual positive adduct ([M+NH4]⁺ for TG, DG, MG,
CE, PI, PG; [M+H]⁺ otherwise), same bridge, same 0.1–10 cap. Per-feature adducts would be
better, but only 487 MZmine-rule annotations carry one; LipidSearch's main ion is not in the
merged table yet. Median f becomes TG 0.14, DG 1.08, PI 0.40, PG 0.25, none capped.
**MG stays capped (98 %)**: its [M+NH4]⁺ response is also tiny, so under A MG should still be
demoted to Q4 (call this A + MG→Q4).

**B. Demote:** Q3 features of TG, DG, MG, PI, PG become Q4 (relative only, no mass).

## Effect (POS, µg g⁻¹ dry soil, one feature per unique lipid)

| | current | A | A + MG→Q4 | B |
|---|---|---|---|---|
| TG (mean of 12 soils) | 239.0 | 4.7 | 4.7 | — |
| DG | 65.7 | 5.6 | 5.6 | — |
| MG | 16.2 | 16.2 (capped) | — | — |
| PI | 8.5 | 0.36 | 0.36 | — |
| PG | 8.2 | 0.46 | 0.46 | — |
| PC / PE (unchanged, Q1) | 34.5 / 8.6 | same | same | same |
| Total, ambient | 202 ± 239 | 56 ± 19 | 45 ± 19 | 40 ± 15 |
| Total, future | 164 ± 37 | 68 ± 9 | 51 ± 6 | 43 ± 4 |
| Total, ambient + drought | 548 ± 138 | 94 ± 12 | 76 ± 9 | 62 ± 6 |
| Total, future + drought | 642 ± 133 | 95 ± 4 | 77 ± 3 | 60 ± 2 |

Drought still raises totals under both options (≈ 1.4–1.7×, not 3×).

Recommendation: **A + MG→Q4** (keeps TG/DG/PI/PG as estimates on the response their real adduct has; MG has no usable response). Per-feature adducts from LipidSearch would be the proper follow-up.

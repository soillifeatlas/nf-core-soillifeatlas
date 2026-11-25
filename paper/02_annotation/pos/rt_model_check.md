# Is the per-class RT model adequate for TG? (POS, review item 3)

2026-11-25. Data: `annotations.tsv` (run on 2ed5d70), refit with the same rule as
`annotate.check_rt` (rt ~ carbons + double bonds per class group, Gold anchors of plausible
composition, one refit without anchors beyond 3 robust SD; tolerance 4 × robust SD within
1–2 min). Thresholds unchanged. Plot: `rt_model_TG_residuals.png`.

| Group | Gold anchors | R² (all anchors) | robust SD | tested | outliers |
|---|---|---|---|---|---|
| TG (diacyl/triacyl) | 109 | 0.63 | 0.37 min | 404 | 246 |
| TG modified (+O…) | 78 | 0.93 | 1.33 min | 83 | 12 |
| TG ether | 2 | no fit | | | |
| DG | 70 | 0.64 | 0.46 min | 133 | 65 |
| DG ether | 6 | 0.77 | 0.94 min | 125 | 70 |
| PC | 13 | 0.62 | 0.39 min | 269 | 92 |

**Verdict: the TG fit is adequate; the outliers are mostly misassigned names, not a bad model.**

- The 158 TG features within tolerance follow the model closely: R² (observed vs expected) = 0.96,
  median |residual| 0.18 min, along ECN 30–51 (plot, green).
- The R² of 0.63 over all anchors is pulled down by Gold anchors that themselves elute far too
  early (anchor RT down to 5.3 min); 9 of the 109 anchors end up flagged as outliers.
- 243 of 246 outliers elute **earlier** than predicted (median RT 15.9 min vs 22.7 min for
  accepted TG; median residual −4.2 min), many with compositions no TG has at that RT
  (e.g. "TG 54:15" at 1.7 min, "TG 50:13" at 1.6 min). 215 are Silver, 31 Gold; 245 of 246 come
  from LipidSearch. These look like other classes, in-source fragments or adducts named as TG.
- Same pattern for DG (65 outliers, all LipidSearch) and PC (92, mostly Silver from MZmine rules
  and LipidSearch).

Consequence for the paper: report the TG outliers as flagged LipidSearch names rather than as a
model failure. Possible follow-ups (not done, they would change thresholds or ranking): require
two sources for Silver TG, or fit anchors only from the late-eluting TG band.

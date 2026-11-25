# Stage 5 (minimal): Q1–Q4 design note

Status: **decided 2026-11-25: option A (follow the manuscript).**

## The mismatch

| Route | Manuscript §2.7 | `annotate.merge()` at 2ed5d70 |
|---|---|---|
| Q1 | class-matched internal standard | Gold/Silver **and** an IS of the same class in the run |
| Q2 | lipid-specific *predicted* relative response | Gold/Silver **without** a class-matched IS |
| Q3 | class-level response | Bronze (class only) |
| Q4 | relative signal, no mass concentration | (not assigned; features without annotation) |

The code ranks routes by annotation confidence. The manuscript ranks them by the
**response information available**. They agree on Q1 and Q4 but not on Q2 and Q3.

## Proposal (option A, recommended): follow the manuscript, decide per feature by response evidence

1. **Q1**: the lipid's class has an internal standard detected in that sample → f = 1.
2. **Q2**: no class-matched IS, but a *lipid-specific* predicted response exists
   (sum composition known, i.e. Gold or Silver, **and** a predicted IE for that
   composition is available from the Samrat et al. 2025 model) → f = RF_IS / RF_pred.
3. **Q3**: class known (Gold/Silver/Bronze), no Q1/Q2 → nearest-class IS with the
   class-level response ratio from the 2025 response table; if the class has no
   table entry, f = 1 and flag `f_default`.
4. **Q4**: no class → normalised area only; excluded from mass sums.

Confidence stays a separate column; it no longer decides the route.
If the 2025 supplement has no per-composition IE predictions, **Q2 will be empty in
this draft** and the paper says so (the workflow supports it; the data for this run
do not).

## Option B: keep the code, change the manuscript

Q2 = "identified species, response from the nearest IS class", Q3 = "class-level
only". Faster, but then "lipid-specific predicted response" disappears from the
paper.

## Formula (both options)

c(i,j) = (A(i,j) / A(i,k)) × m_spiked(i,k) × f(j,k) / M_dry(i), no extra recovery
factor. Per lipid the output carries: route, IS used (k), f and its source, and a
calibration-range flag (IS area of sample i inside / below / above the neat
calibration range of standard k).

Inputs still needed per soil: dry mass M_dry, spiked amount per IS, and the
Samrat et al. (2025) response tables.

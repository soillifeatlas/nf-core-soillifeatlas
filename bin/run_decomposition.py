#!/usr/bin/env python3
"""CLI dispatcher for the four soillifeatlas decomposition methods.

One invocation per method — designed to be fanned out by Nextflow's `each`
operator so all four methods run in parallel against the same
(soil, atlas, SIMPER) triplet.

Methods (`--method`):
    nnls           scipy NNLS on L1-normalised intensity
    std_bc         Bray-Curtis similarity (1 - BC) on full feature space
    enriched_bc    BC restricted to SIMPER-enriched features per phylum
    fc_weighted_bc fold-change-weighted BC (weights from SIMPER atlas)

Inputs:
    --soil-intensity      parquet, features x samples  (unknown mixes)
    --atlas-intensity     parquet, features x samples  (phylum-labelled atlas)
    --simper-atlas        parquet SIMPER fingerprint atlas with columns
                          feature_id, phylum, kingdom, direction, fold_change
    --sample-phylum-map   CSV with columns (sample_id, phylum) linking atlas
                          sample columns to phylum labels
    --method              one of {nnls, std_bc, enriched_bc, fc_weighted_bc}
    --min-phylum-samples  phyla with fewer atlas samples than this are dropped
                          (default 2)

Output:
    --output              parquet in long (tidy) format with columns
                          (sample_id, kingdom, proportion_pct); each sample's
                          proportions sum to 100.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from soillifeatlas.decomposition import (
    build_direction_masks,
    build_phylum_reference_array,
    build_soil_matrix_array,
    decompose_enriched_only_bc,
    decompose_fc_weighted_bc,
    decompose_nnls,
    decompose_standard_bc,
    phylum_to_kingdom_array,
)


# Dispatch table: normalise a small set of short method names to the
# underlying framework call. The Nextflow `each` channel will pass these
# shorthand names on the CLI.
METHOD_DISPATCH = {
    "nnls":           lambda ref, soil, masks, weights, phyla: decompose_nnls(ref, soil),
    "std_bc":         lambda ref, soil, masks, weights, phyla: decompose_standard_bc(ref, soil),
    "enriched_bc":    lambda ref, soil, masks, weights, phyla: decompose_enriched_only_bc(ref, soil, masks, phyla),
    "fc_weighted_bc": lambda ref, soil, masks, weights, phyla: decompose_fc_weighted_bc(ref, soil, weights, phyla),
}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--soil-intensity", required=True, type=Path,
                   help="Soil intensity parquet (features x samples)")
    p.add_argument("--atlas-intensity", required=True, type=Path,
                   help="Atlas intensity parquet (features x samples)")
    p.add_argument("--simper-atlas", required=True, type=Path,
                   help="SIMPER fingerprint atlas parquet (feature_id, phylum, kingdom, direction, fold_change)")
    p.add_argument("--sample-phylum-map", required=True, type=Path,
                   help="CSV mapping atlas sample_id -> phylum")
    p.add_argument("--method", required=True, choices=sorted(METHOD_DISPATCH),
                   help="Which decomposition method to run")
    p.add_argument("--min-phylum-samples", type=int, default=2,
                   help="Drop phyla with fewer atlas samples than this (default: 2)")
    p.add_argument("--output", required=True, type=Path,
                   help="Output parquet (long-format: sample_id, kingdom, proportion_pct)")
    args = p.parse_args(argv)

    # --- Load inputs ---------------------------------------------------------
    soil_df = pd.read_parquet(args.soil_intensity)
    atlas_df = pd.read_parquet(args.atlas_intensity)
    simper = pd.read_parquet(args.simper_atlas)
    sp_df = pd.read_csv(args.sample_phylum_map)
    sample_phylum = dict(zip(sp_df["sample_id"], sp_df["phylum"]))

    # --- Step 1: phylum reference (atlas) -----------------------------------
    ref, phyla, feature_ids = build_phylum_reference_array(
        atlas_df, sample_phylum, min_samples=args.min_phylum_samples,
    )
    if ref.shape[0] == 0:
        sys.exit(
            "No phyla passed the --min-phylum-samples threshold; nothing to decompose. "
            f"(atlas samples in map: {len(sample_phylum)}, threshold: {args.min_phylum_samples})"
        )

    # --- Step 2: align soil to atlas feature order --------------------------
    soil, sample_cols = build_soil_matrix_array(soil_df, feature_ids)
    if soil.shape[0] == 0:
        sys.exit("Soil intensity matrix has zero samples; nothing to decompose.")

    # --- Step 3: direction masks + FC weights from SIMPER atlas -------------
    enriched_masks, fc_weights = build_direction_masks(simper, phyla, feature_ids)

    # --- Step 4: run the selected method ------------------------------------
    comp = METHOD_DISPATCH[args.method](ref, soil, enriched_masks, fc_weights, phyla)

    # --- Step 5: roll up phylum -> kingdom ----------------------------------
    # phylum_to_kingdom_array reads the module-level PHYLUM_KINGDOM map and
    # takes just (comp, phyla). The SIMPER atlas's `kingdom` column is not
    # consulted here — the framework treats the phylum -> kingdom mapping as
    # curation-level knowledge pinned in the decomposition module.
    kingdom_comp, kingdoms = phylum_to_kingdom_array(comp, phyla)

    # --- Step 6: emit long-format parquet -----------------------------------
    rows = [
        {
            "sample_id": sample_cols[i],
            "kingdom": kingdoms[k],
            "proportion_pct": 100.0 * float(kingdom_comp[i, k]),
        }
        for i in range(len(sample_cols))
        for k in range(len(kingdoms))
    ]
    out_df = pd.DataFrame(rows, columns=["sample_id", "kingdom", "proportion_pct"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""One-shot: assemble a tiny atlas bundle for the FETCH_ATLAS nf-test + APPLY
end-to-end nf-test.

Writes into tests/fixtures/tiny/tiny_atlas_bundle/ the same canonical artefacts
a real Zenodo-published atlas would contain, but built from the existing tiny
fixtures (renamed where necessary). This keeps the FETCH_ATLAS --atlas-path
override branch exercised in CI without hitting the network.

Artefacts laid out here:
  simper_fingerprint_atlas.parquet   <- tiny_simper_atlas.parquet
                                        (the feat_NNN-keyed one, used by the
                                        decomposition subworkflow; not the
                                        MGF-keyed variant)
  atlas_consensus.mgf                <- tiny_atlas.mgf
  consensus_aligned_table.parquet    <- tiny_atlas_intensity.parquet
  consensus_unified_annotations.csv  <- tiny_annotation.csv
                                        (kept as CSV for v0.1; downstream
                                        apply_rie_correction.py reads CSV.
                                        Renamed on the parquet/csv boundary
                                        when the atlas publishing pipeline
                                        switches format in v0.2+.)
  rie_table_s10.csv                  <- tiny_rie_table.csv
  equisplash_IS_masses_POS.csv       <- tiny_is_features.csv
  archlips_validated_features.csv    <- tiny_archlips_validated.csv
  expected_kingdom_composition.csv   <- tiny_expected_kingdom_composition.csv
  sample_metadata.tsv                <- tiny_sample_metadata.tsv
  sample_phylum_map.csv              <- tiny_sample_phylum_map.csv

Note: v0.1 FETCH_ATLAS only *requires* simper_fingerprint_atlas.parquet (see
REQUIRED_ARTEFACTS in bin/fetch_atlas.py). Extra artefacts are shipped so that
the downstream DECOMPOSE_APPLY subworkflow has everything it needs when the
APPLY workflow wires atlas_dir/<file> paths into its channels.

Commit the resulting directory so CI reproduces the bundle byte-for-byte.
Re-run only if a fixture layout needs to change.
"""
from pathlib import Path
import shutil


FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "tiny"
BUNDLE = FIXTURES / "tiny_atlas_bundle"

# (source_fixture, target_name_in_bundle) pairs that are a direct copy.
DIRECT_COPIES = [
    ("tiny_simper_atlas.parquet",               "simper_fingerprint_atlas.parquet"),
    ("tiny_atlas.mgf",                          "atlas_consensus.mgf"),
    ("tiny_atlas_intensity.parquet",            "consensus_aligned_table.parquet"),
    ("tiny_annotation.csv",                     "consensus_unified_annotations.csv"),
    ("tiny_rie_table.csv",                      "rie_table_s10.csv"),
    ("tiny_is_features.csv",                    "equisplash_IS_masses_POS.csv"),
    ("tiny_archlips_validated.csv",             "archlips_validated_features.csv"),
    ("tiny_expected_kingdom_composition.csv",   "expected_kingdom_composition.csv"),
    ("tiny_sample_metadata.tsv",                "sample_metadata.tsv"),
    ("tiny_sample_phylum_map.csv",              "sample_phylum_map.csv"),
]


def main() -> None:
    if not FIXTURES.is_dir():
        raise FileNotFoundError(f"tiny fixtures dir not found: {FIXTURES}")

    # Wipe + recreate so stale layouts from previous runs don't bleed through.
    if BUNDLE.exists():
        shutil.rmtree(BUNDLE)
    BUNDLE.mkdir(parents=True)

    for src, dst in DIRECT_COPIES:
        src_path = FIXTURES / src
        dst_path = BUNDLE / dst
        if not src_path.exists():
            raise FileNotFoundError(f"missing source fixture: {src_path}")
        shutil.copy(src_path, dst_path)

    print(f"tiny atlas bundle assembled at: {BUNDLE}")
    for f in sorted(BUNDLE.iterdir()):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()

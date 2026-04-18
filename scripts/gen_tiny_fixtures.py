"""One-shot: generate tiny synthetic fixtures for bin/ wrapper unit tests.

Run once manually; commit the generated files under tests/fixtures/tiny/.
Re-run only if the fixture schemas need to change.
"""
import numpy as np
import pandas as pd
from pathlib import Path

FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "tiny"
FIXTURES.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

# ---------- intensity matrix: 10 features x 4 samples ----------
intensity = pd.DataFrame(
    rng.random((10, 4)) * 1e6,
    columns=["S01", "S02", "S03", "S04"],
    index=[f"feat_{i:03d}" for i in range(10)],
)
intensity.index.name = "feature_id"
intensity.to_parquet(FIXTURES / "tiny_intensity.parquet")

# ---------- IS feature table (feat_000 is our IS) ----------
# Column `feature_id` links to the intensity matrix index; `compound` is the
# name consumed by --is-reference on the CLI.
is_features = pd.DataFrame({
    "feature_id": ["feat_000"],
    "compound":  ["LPE_18d7"],
    "adduct":    ["[M+H]+"],
    "spiked_pmol": [100.0],
})
is_features.to_csv(FIXTURES / "tiny_is_features.csv", index=False)

# ---------- RIE table (columns match what apply_RIE_correction expects) ----------
# Real signature: rie_lookup columns = ["class", "adduct", "RIE_LPE"].
# Very low MG RIE exercises the floor at 0.20 (wrapper CLI arg).
# NOTE on adduct format: corrections.apply_RIE_correction preprocesses the
# lookup-side adduct by stripping only [ and ] (see lines 125-130), but the
# query-side also strips the trailing charge (+/-) at line 140. To keep those
# keys aligned, we store adducts here pre-stripped as "M+H" (no brackets, no
# charge suffix). A bulk / real-world RIE table would require the same
# preprocessing before being fed to the pipeline.
rie_table = pd.DataFrame([
    {"class": "PE",      "adduct": "M+H", "RIE_LPE": 0.85},
    {"class": "PC",      "adduct": "M+H", "RIE_LPE": 1.00},
    {"class": "MG",      "adduct": "M+H", "RIE_LPE": 0.0005},  # very low — tests the floor
    {"class": "DG",      "adduct": "M+H", "RIE_LPE": 0.03},
    {"class": "UNKNOWN", "adduct": "M+H", "RIE_LPE": 0.5},
])
rie_table.to_csv(FIXTURES / "tiny_rie_table.csv", index=False)

# ---------- feature annotation (maps feature_id -> class + adduct) ----------
# First 5 features mapped; rest UNKNOWN.
annotation = pd.DataFrame({
    "feature_id": [f"feat_{i:03d}" for i in range(10)],
    "class":  ["PE", "PC", "MG", "DG", "UNKNOWN"] + ["UNKNOWN"] * 5,
    "adduct": ["[M+H]+"] * 10,
})
annotation.to_csv(FIXTURES / "tiny_annotation.csv", index=False)

# ---------- ArchLips validated features (feat_007, feat_008, feat_009) ----------
archlips = pd.DataFrame({"feature_id": ["feat_007", "feat_008", "feat_009"]})
archlips.to_csv(FIXTURES / "tiny_archlips_validated.csv", index=False)

# ---------- SIMPER fingerprint atlas (just enough for ref-quality filter test) ----------
# 10 features x 3 phyla, some Archaea, some not
simper_atlas = pd.DataFrame({
    "feature_id": [f"feat_{i:03d}" for i in range(10)] * 3,
    "phylum": ["Actinomycetota"] * 10 + ["Ascomycota"] * 10 + ["Euryarchaeota"] * 10,
    "kingdom": ["Bacteria"] * 10 + ["Fungi"] * 10 + ["Archaea"] * 10,
    "direction": ["enriched"] * 30,
    "fc_weight": rng.random(30),
    "simper_rank": list(range(10)) * 3,
})
simper_atlas.to_parquet(FIXTURES / "tiny_simper_atlas.parquet")

print(f"Fixtures written to: {FIXTURES}")

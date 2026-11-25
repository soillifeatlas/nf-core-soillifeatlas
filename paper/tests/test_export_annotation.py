import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import export_annotation as ea  # noqa: E402


def test_adducts_of_one_lipid_collapse_but_isomers_do_not():
    ann = pd.DataFrame({
        "feature_id": [1, 2, 3, 4],
        "lipid_name": ["PC 34:1", "PC 34:1", "PC 34:1", "PC"],
        "confidence": ["Gold", "Gold", "Silver", "Bronze"],
        "rt": [15.00, 15.05, 16.00, 15.0],
    })
    ids = ea.collapse_unique(ann, 0.1)
    assert ids[0] == ids[1]          # [M+H]+ and [M+Na]+ of one PC
    assert ids[2] != ids[0]          # same name 1 min later: a separate isomer
    assert ids[3] == "feature:4"     # class-only names stay separate


def test_upset_counts_each_feature_once():
    ann = pd.DataFrame({"sources": ["A; B", "A", "A; B", "B"]})
    sizes, inter = ea.upset(ann)
    assert dict(zip(sizes["source"], sizes["n_features"])) == {"A": 3, "B": 3}
    assert inter["n_features"].sum() == 4

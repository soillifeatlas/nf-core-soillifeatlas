import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import nist_compare as nc  # noqa: E402

TEXT = """Table 3. Consensus estimates
      TAG 52:2        16      nmol/mL         49            23            71            14           44             14           33
Table 4. Consensus estimates
      PC 34:1                      19     nmol/mL          33            33           1,500         1,300           120             21             17
Table 9. other
      PC 36:2      5   nmol/mL   1  1  1  1  1  1  1
"""


def test_parse_reads_medm_and_maps_names():
    t = nc.parse(TEXT).set_index("sum_name")
    assert t.loc["TG 52:2", "medm"] == 44
    assert t.loc["PC 34:1", "vr"] == 1500 and t.loc["PC 34:1", "n_labs"] == 19
    assert "PC 36:2" not in t.index          # tables outside 2-6 ignored


def test_sum_names():
    assert nc.sum_name("PC 16:0_18:1", "PC") == "PC 34:1"
    assert nc.sum_name("TG 16:0_18:1_18:2", "TG") == "TG 52:3"
    assert nc.sum_name("SM d18:1_16:0", "SM") == "SM 34:1"
    assert nc.sum_name("PC O-34:1", "PC") is None

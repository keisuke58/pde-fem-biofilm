"""compare_cal5.py on the committed CH run (the other conditions are skipped
when they are not there)."""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "ansys_usermat" / "apdl" / "results" / "2026-10-cal5"


def test_compare_cal5_ch(tmp_path):
    if not (RES / "w8_CH_cal_k0882_pm.json").exists():
        import pytest
        pytest.skip("CH run not in this checkout")
    r = subprocess.run([sys.executable, str(REPO / "ansys_usermat" / "apdl" / "compare_cal5.py"),
                        "--tag=k0882", "--results", str(RES), "--out", str(tmp_path), "--no-fig"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = (tmp_path / "cal5_compare_k0882.csv").read_text().splitlines()
    head = rows[0].split(",")
    ch = dict(zip(head, rows[1].split(",")))
    assert ch["condition"] == "CH"
    assert int(ch["n_seed"]) == 32
    assert abs(float(ch["seqv_mean_Pa"]) - 4.3192e-4) < 1e-8
    assert abs(float(ch["alpha_minus_1_mean"]) - 6.7518e-4) < 1e-8
    assert abs(float(ch["seqv_vs_CH"]) - 1.0) < 1e-12
    md = (tmp_path / "cal5_compare_k0882.md").read_text(encoding="utf-8").splitlines()
    table = [l for l in md if l.startswith("|")]
    assert len({l.count("|") for l in table}) == 1

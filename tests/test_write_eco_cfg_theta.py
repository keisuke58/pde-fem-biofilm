"""write_eco_cfg.py --theta-json: a calibrated five-species MAP file becomes a
native point-model configuration (n_active 5, c* and alpha* of
ecology_constants, eta 1, theta unchanged)."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ecology_constants as ec  # noqa: E402


def test_theta_json(tmp_path):
    theta = [0.1 * (k + 1) for k in range(20)]
    src = tmp_path / "theta_MAP.json"
    src.write_text(json.dumps({"theta_full": theta, "theta_sub": theta}))
    out = tmp_path / "eco.txt"
    subprocess.run([sys.executable, str(ROOT / "abaqus_composition" / "write_eco_cfg.py"),
                    "--theta-json", str(src), str(out)], check=True, capture_output=True)
    lines = out.read_text().split("\n")
    assert int(lines[0]) == 5
    assert [float(v) for v in lines[1].split()] == [ec.C_STAR, ec.ALPHA_STAR]
    assert [float(v) for v in lines[2].split()] == [1.0] * 5
    assert [float(v) for v in lines[3].split()] == theta


def test_theta_json_wrong_length(tmp_path):
    src = tmp_path / "bad.json"
    src.write_text(json.dumps({"theta_full": [0.0] * 14}))
    r = subprocess.run([sys.executable, str(ROOT / "abaqus_composition" / "write_eco_cfg.py"),
                        "--theta-json", str(src), str(tmp_path / "x.txt")], capture_output=True)
    assert r.returncode != 0

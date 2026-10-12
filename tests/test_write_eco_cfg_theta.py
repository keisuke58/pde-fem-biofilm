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


def test_theta_json_index_dict(tmp_path):
    """{"0": v0, ..., "19": v19}: the values, in index order, not the keys."""
    theta = [0.5 - 0.05 * k for k in range(20)]
    src = tmp_path / "theta_MAP.json"
    src.write_text(json.dumps({str(k): theta[k] for k in reversed(range(20))}))
    out = tmp_path / "eco.txt"
    subprocess.run([sys.executable, str(ROOT / "abaqus_composition" / "write_eco_cfg.py"),
                    "--theta-json", str(src), str(out)], check=True, capture_output=True)
    assert [float(v) for v in out.read_text().split("\n")[3].split()] == theta


def test_keyword_lines(tmp_path):
    """--phi-init-config prepares metadata.phi_init_exp as the paper pipeline
    does (sum 1, clip, scaled to 0.999999); --newton and --theta-tol add
    their lines; the first four lines are unchanged."""
    import numpy as np
    theta = [0.1 * (k + 1) for k in range(20)]
    src = tmp_path / "theta_MAP.json"
    src.write_text(json.dumps({"theta_full": theta}))
    cfg = tmp_path / "config.json"
    exp = [0.375, 0.025, 0.05, 0.05, 0.0025]
    cfg.write_text(json.dumps({"metadata": {"phi_init_exp": exp}}))
    out = tmp_path / "eco.txt"
    subprocess.run([sys.executable, str(ROOT / "abaqus_composition" / "write_eco_cfg.py"),
                    "--theta-json", str(src), "--phi-init-config", str(cfg),
                    "--newton", "12,1e-20", "--theta-tol", "1e-12", str(out)],
                   check=True, capture_output=True)
    lines = out.read_text().strip().split("\n")
    assert len(lines) == 7
    assert [float(v) for v in lines[3].split()] == theta
    p = np.asarray(exp)
    p = p / p.sum()          # numpy's sum, as estimate_paper_jax.py
    p = np.clip(p, 0.001, 0.99)
    p = p * (0.999999 / min(p.sum(), 0.999999))
    key, *vals = lines[4].split()
    assert key == "phi_init" and [float(v) for v in vals] == [float(x) for x in p]
    # make_initial_state scales by 0.999999 / min(sum, 0.999999): a sum of
    # 1 stays 1 (phi_0 = 0, clipped to 1e-10 at the first Newton iteration)
    assert abs(sum(float(v) for v in vals) - 1.0) < 1e-12
    assert lines[5].split() == ["newton", "12", repr(1e-20)]
    assert lines[6].split() == ["theta_tol", repr(1e-12)]


def test_theta_json_other_dict(tmp_path):
    src = tmp_path / "bad.json"
    src.write_text(json.dumps({"theta_sub": [0.0] * 20}))
    r = subprocess.run([sys.executable, str(ROOT / "abaqus_composition" / "write_eco_cfg.py"),
                        "--theta-json", str(src), str(tmp_path / "x.txt")], capture_output=True)
    assert r.returncode != 0


def test_theta_json_wrong_length(tmp_path):
    src = tmp_path / "bad.json"
    src.write_text(json.dumps({"theta_full": [0.0] * 14}))
    r = subprocess.run([sys.executable, str(ROOT / "abaqus_composition" / "write_eco_cfg.py"),
                        "--theta-json", str(src), str(tmp_path / "x.txt")], capture_output=True)
    assert r.returncode != 0


def test_phi_init_state_matches_config_path(tmp_path):
    """run_calibrated.ps1 (12 Oct) passes phi_init_state of Tmcmc202601
    final_theta_MAP/phi_init.json through --phi-init. That entry is the
    estimator's phi_init (Day-1 composition, sum 1, clipped) scaled as
    make_initial_state does; written that way it must give the same
    phi_init line, digit for digit, as --phi-init-config on the same
    Day-1 composition."""
    import numpy as np
    theta = [0.1 * (k + 1) for k in range(20)]
    src = tmp_path / "theta_MAP.json"
    src.write_text(json.dumps({"theta_full": theta}))
    exp = [0.375, 0.025, 0.05, 0.05, 0.0025]
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps({"metadata": {"phi_init_exp": exp}}))
    p = np.asarray(exp)
    p = np.clip(p / p.sum(), 0.001, 0.99)          # run_record.json "phi_init"
    state = p * (0.999999 / min(p.sum(), 0.999999))  # phi_init.json "phi_init_state"
    script = str(ROOT / "abaqus_composition" / "write_eco_cfg.py")
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    subprocess.run([sys.executable, script, "--theta-json", str(src), "--phi-init-config", str(cfg), str(a)],
                   check=True, capture_output=True)
    subprocess.run([sys.executable, script, "--theta-json", str(src),
                    "--phi-init", ",".join(repr(float(x)) for x in state), str(b)],
                   check=True, capture_output=True)
    assert a.read_text() == b.read_text()

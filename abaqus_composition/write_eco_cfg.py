"""write_eco_cfg.py -- the case constants of the point model for the Fortran
version (ecology_native.f, eco_native_init), taken from material_server.set_case
so the native run uses exactly what the server would.

    python abaqus_composition/write_eco_cfg.py CASE OUT.txt
    python abaqus_composition/write_eco_cfg.py --theta-json theta_MAP.json OUT.txt

OUT: n_active / c* alpha* / eta(5) / theta(20).

--theta-json: a calibrated five-species parameter set from TMCMC (the file
theta_MAP.json of a run, key "theta_full", 20 values in the order of the
point model: the 15 entries of A and the 5 entries of b), with c* and
alpha* from ecology_constants.py (the calibration values) and eta = 1, as
the server uses without a case. This is the step from the calibrated point
model of Chapter 3 to the element (research idea 7, RESEARCH_IDEAS.ja.md).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "ansys_usermat" / "coupling"), str(ROOT / "ansys_usermat"), str(ROOT)]
import material_server as ms  # noqa: E402


def theta_json_cfg(path):
    """(n_active, c*, alpha*, eta(5), theta(20)) for a calibrated MAP file."""
    import ecology_constants as ec
    d = json.loads(Path(path).read_text())
    if isinstance(d, dict) and "theta_full" in d:
        theta = d["theta_full"]
    elif isinstance(d, dict) and d and all(str(k).isdigit() for k in d):
        # {"0": v0, ..., "19": v19}, as in tmcmc202601 main/_runs/phibar_fix_*;
        # iterating the dict itself would take the keys 0..19 as theta
        theta = [d[k] for k in sorted(d, key=int)]
        if [int(k) for k in sorted(d, key=int)] != list(range(len(d))):
            raise ValueError(f"{path}: index keys are not 0..{len(d) - 1}")
    elif isinstance(d, list):
        theta = d
    else:
        raise ValueError(f"{path}: no 'theta_full' and not a list or an index dict")
    theta = [float(v) for v in theta]
    if len(theta) != 20:
        raise ValueError(f"{path}: expected 20 values (5 species), got {len(theta)}")
    return 5, float(ec.C_STAR), float(ec.ALPHA_STAR), [1.0] * 5, theta


def main(case, out):
    ms.set_case(case)
    c = ms.ECOLOGY_CASE
    hp = c["hp"]
    f = lambda v: " ".join(repr(float(x)) for x in v)  # noqa: E731
    Path(out).write_text(f"{ms.ECOLOGY_ACTIVE}\n{f([hp['c'], hp['alpha']])}\n{f(hp['eta'])}\n{f(c['theta'])}\n")
    print(f"wrote {out}: {case}, {ms.ECOLOGY_ACTIVE} species, c* {hp['c']}, alpha* {hp['alpha']}")


def main_theta(path, out):
    n, cs, als, eta, theta = theta_json_cfg(path)
    f = lambda v: " ".join(repr(float(x)) for x in v)  # noqa: E731
    Path(out).write_text(f"{n}\n{f([cs, als])}\n{f(eta)}\n{f(theta)}\n")
    print(f"wrote {out}: {path}, {n} species, c* {cs}, alpha* {als}")


if __name__ == "__main__":
    if sys.argv[1] == "--theta-json":
        main_theta(sys.argv[2], sys.argv[3])
    else:
        main(sys.argv[1], sys.argv[2])

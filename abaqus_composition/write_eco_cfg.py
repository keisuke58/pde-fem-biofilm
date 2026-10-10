"""write_eco_cfg.py -- the case constants of the point model for the Fortran
version (ecology_native.f, eco_native_init), taken from material_server.set_case
so the native run uses exactly what the server would.

    python abaqus_composition/write_eco_cfg.py CASE OUT.txt
    python abaqus_composition/write_eco_cfg.py --theta-json theta_MAP.json OUT.txt
        [--phi-init-config config.json | --phi-init v1,v2,v3,v4,v5]
        [--newton 12,1e-20] [--theta-tol 1e-12]

OUT: n_active / c* alpha* / eta(5) / theta(20), then optional keyword lines
(ecology_native.f reads them in any order):
  phi_init v1 .. v5   the point model's state at its first call. With
                      --phi-init-config, metadata.phi_init_exp of the
                      calibration's config.json (the Day-1 composition)
                      prepared as estimate_paper_jax.py and
                      hamilton_ode_jax_paper.make_initial_state do: divided
                      by its sum, clipped to [0.001, 0.99], scaled to a sum
                      of 0.999999. With --phi-init, the five values as given.
  newton N TOL        the paper pipeline's Newton control (12, 1e-20)
  theta_tol X         relative tolerance of the deck-vs-file theta check

--theta-json: a calibrated five-species parameter set from TMCMC (the file
theta_MAP.json of a run, key "theta_full", 20 values in the order of the
point model: the 15 entries of A and the 5 entries of b), with c* and
alpha* from ecology_constants.py (the calibration values) and eta = 1, as
the server uses without a case. This is the step from the calibrated point
model of Chapter 3 to the element (research idea 7, RESEARCH_IDEAS.ja.md).
"""
import argparse
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


def phi_init_from_config(path):
    """metadata.phi_init_exp of a calibration config.json, prepared as the
    paper pipeline prepares its initial state (estimate_paper_jax.py with
    --use-exp-init, then make_initial_state)."""
    import numpy as np
    d = json.loads(Path(path).read_text())
    p = np.asarray(d["metadata"]["phi_init_exp"], dtype=np.float64)
    if p.shape != (5,):
        raise ValueError(f"{path}: phi_init_exp has shape {p.shape}, expected (5,)")
    total = p.sum()
    if total > 0:
        p = p / total
    p = np.clip(p, 0.001, 0.99)
    phi_sum = min(p.sum(), 0.999999)
    p = p * (0.999999 / phi_sum)
    return [float(x) for x in p]


def fmt(v):
    return " ".join(repr(float(x)) for x in v)


def keyword_lines(phi_init=None, newton=None, theta_tol=None):
    out = []
    if phi_init is not None:
        if len(phi_init) != 5:
            raise ValueError(f"phi_init needs 5 values, got {len(phi_init)}")
        out.append("phi_init " + fmt(phi_init))
    if newton is not None:
        n, tol = newton
        out.append(f"newton {int(n)} {repr(float(tol))}")
    if theta_tol is not None:
        out.append(f"theta_tol {repr(float(theta_tol))}")
    return out


def main(case, out, **kw):
    ms.set_case(case)
    c = ms.ECOLOGY_CASE
    hp = c["hp"]
    lines = [str(ms.ECOLOGY_ACTIVE), fmt([hp["c"], hp["alpha"]]), fmt(hp["eta"]), fmt(c["theta"])]
    lines += keyword_lines(**kw)
    Path(out).write_text("\n".join(lines) + "\n")
    print(f"wrote {out}: {case}, {ms.ECOLOGY_ACTIVE} species, c* {hp['c']}, alpha* {hp['alpha']}")


def main_theta(path, out, **kw):
    n, cs, als, eta, theta = theta_json_cfg(path)
    lines = [str(n), fmt([cs, als]), fmt(eta), fmt(theta)] + keyword_lines(**kw)
    Path(out).write_text("\n".join(lines) + "\n")
    extra = ", ".join(k for k, v in kw.items() if v is not None)
    print(f"wrote {out}: {path}, {n} species, c* {cs}, alpha* {als}" + (f", {extra}" if extra else ""))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("args", nargs="+", help="CASE OUT, or OUT with --theta-json")
    ap.add_argument("--theta-json")
    ap.add_argument("--phi-init-config", help="calibration config.json (metadata.phi_init_exp)")
    ap.add_argument("--phi-init", help="v1,v2,v3,v4,v5 as given")
    ap.add_argument("--newton", help="N,TOL (the paper pipeline: 12,1e-20)")
    ap.add_argument("--theta-tol", type=float)
    a = ap.parse_args()
    if a.phi_init_config and a.phi_init:
        ap.error("--phi-init-config and --phi-init exclude each other")
    phi_init = None
    if a.phi_init_config:
        phi_init = phi_init_from_config(a.phi_init_config)
    elif a.phi_init:
        phi_init = [float(x) for x in a.phi_init.split(",")]
    newton = None
    if a.newton:
        n, tol = a.newton.split(",")
        newton = (int(n), float(tol))
    kw = dict(phi_init=phi_init, newton=newton, theta_tol=a.theta_tol)
    if a.theta_json:
        if len(a.args) != 1:
            ap.error("with --theta-json: one positional, OUT")
        main_theta(a.theta_json, a.args[0], **kw)
    else:
        if len(a.args) != 2:
            ap.error("CASE OUT")
        main(a.args[0], a.args[1], **kw)

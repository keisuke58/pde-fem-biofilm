"""p3_posterior.py -- step 2 of the stress-interval idea (KEIO_PLAN.ja.md section 5): the P3 element (unified UEL,
Python emulation) with point-model parameters drawn from the final posterior of chapter 3. In the UEL the growth
follows the living amount (xi_dot = L0 phi sum s_i psi_i (Y ln c - mu_b)), so the posterior reaches growth and
stress directly, and with the homeostatic-pressure law (P3) it does so through a nonlinear gate.

Particles: data/posterior_final/<C>_particles.npz (EML, 19 Oct 2026). Chosen by the integral of the living amount
(uq/posterior_growth_band.py): the 5, 50 and 95 % particles and the MAP. Case: the tooth-like sector of
p3/uel/cervical_study.py (bulge 1, bottom u_z = 0, outer face free, nutrient from the crown), law off and on.

    python3 uq/p3_posterior.py [--cond CH] [--nr 2] [--p-ref 0 0.05]    -> uq/p3_posterior_<C>_<mesh>.json
One subprocess per particle (the point model's constants are read once per process from BIOFILM_ECO_CASE).
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
UEL = ROOT / "p3" / "uel"
DATA = ROOT / "data" / "posterior_final"
SCRATCH = Path(os.environ.get("UQ_SCRATCH", "/tmp/uq_p3"))


def eco_text(theta, phi_init):
    g = np.asarray(phi_init, float); g = g * 0.9 / g.sum()          # as make_column_inp.eco_file
    return "\n".join(["5", "25.0 0.0", " ".join(["1.0"] * 5), " ".join(f"{v:.17g}" for v in theta),
                      "phi_init " + " ".join(f"{v:.17g}" for v in g), "newton 60 1e-26"]) + "\n"


def pick(cond, n_random=0):
    d = np.load(DATA / f"{cond}_particles.npz")
    T = d["t_days"] / 21.0
    I = np.trapezoid(d["sum_phipsi"], T, axis=1)
    order = np.argsort(I)
    if n_random:      # a random subset of the 200 posterior particles (equal weights), seed 0
        idx = np.random.default_rng(0).choice(len(I), n_random, replace=False)
        picks = {f"p{i:03d}": int(i) for i in idx}
    else:
        picks = {f"q{q}": int(order[int(round(q / 100 * (len(I) - 1)))]) for q in (5, 50, 95)}
    out = {k: (d["theta"][i], float(I[i])) for k, i in picks.items()}
    out["MAP"] = (d["map_theta"], float(np.trapezoid(d["map_sum_phipsi"], T)))
    return out, d["phi_init"]


WORKER = r'''
import json, sys
sys.path.insert(0, sys.argv[1])
import make_implant_uel_inp as mi
from cervical_study import measures
nr, p_ref = int(sys.argv[2]), float(sys.argv[3])
f = nr // 2
nodes, elems, layer, ijk = mi.sector(nr, 3 * nr, 4 * nr, 1.0)
r = mi.run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, seed_layers=f, p_ref=p_ref, outer="free",
           bottom="uz", nut="outer")
r.update(nr=nr, nt=3 * nr, nz=4 * nr, bulge=1.0)
s = measures(r)
s["seed_alpha_minus_1"] = r["seed_alpha_minus_1"]
print("RESULT " + json.dumps(s))
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", default="CH")
    ap.add_argument("--nr", type=int, default=2)
    ap.add_argument("--p-ref", type=float, nargs="+", default=[0.0, 0.05])
    ap.add_argument("--n-random", type=int, default=0, help="random posterior particles instead of the I quantiles")
    a = ap.parse_args()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    parts, phi_init = pick(a.cond, a.n_random)
    res = []
    for name, (theta, I) in parts.items():
        eco = SCRATCH / f"eco_{a.cond}_{name}.txt"
        eco.write_text(eco_text(theta, phi_init))
        for p_ref in a.p_ref:
            env = dict(os.environ, BIOFILM_ECO_CASE=str(eco))
            o = subprocess.run([sys.executable, "-c", WORKER, str(UEL), str(a.nr), str(p_ref)], env=env,
                               capture_output=True, text=True, cwd=UEL)
            line = [ln for ln in o.stdout.splitlines() if ln.startswith("RESULT ")]
            if not line:
                print(o.stderr[-2000:]); raise SystemExit(f"{name} p_ref {p_ref} failed")
            s = json.loads(line[0][7:])
            row = dict(particle=name, I_live=I, p_ref=p_ref, core_alpha_minus_1=s["core_alpha_minus_1"],
                       seed_alpha_minus_1=s["seed_alpha_minus_1"], core_vM_over_E=s["core_vM_over_E"],
                       band_vM_over_E=s["band_vM_over_E"], ratio=s["ratio"],
                       seed_max_p_over_E=s["seed_max_p_over_E"], cutbacks=s["cutbacks"])
            res.append(row)
            print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    out = ROOT / "uq" / (f"p3_posterior_{a.cond}_{a.nr}x{3 * a.nr}x{4 * a.nr}" + (f"_r{a.n_random}" if a.n_random else "") + ".json")
    out.write_text(json.dumps(res, indent=1))
    print("wrote", out)


if __name__ == "__main__":
    main()

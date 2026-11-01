#!/usr/bin/env python3
"""Pilot of a surrogate model for the coupled example (separate start with antibiotic), Keio continuation, not thesis.

Two parameters are varied on a Latin hypercube: the antibiotic strength alpha*_0 (log-uniform 10 to 300; the effect
switches on between 30 and 100) and the nutrient uptake g l^2/d (0.5 to 2). Every run is the 32^3 coupled example of
the thesis (Python, --fast, TVD, dc/dt with --transient-c 1e8, outer step 0.025, about 35 s on one core). The seed
means at T* = 1 go to runs/; fit_pilot.py fits a Gaussian process and checks it by leave-one-out.

    python papers/surrogate_pilot/run_pilot.py [--n 16] [--jobs 4] [--seed 0]            # first design, 32^3 (tag p)
    python papers/surrogate_pilot/run_pilot.py --tag q --n 12 --seed 1 --a0-range 40 150  # more points near the switch
    python papers/surrogate_pilot/run_pilot.py --tag f --nf 64 --reuse papers/surrogate_pilot/design_p.json
"""
import argparse
import json
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOLVER = ROOT / "ansys_usermat" / "apdl" / "coupled_living_consumption.py"
RES = ROOT / "ansys_usermat" / "apdl" / "results" / "2026-10-3d"
BASE = ["--fast", "--fortran", "--front", "1", "--gab", "1", "--tab", "0.333333", "--phi-min", "0.001",
        "--tvd", "--wall-grad", "mirror", "--v2", "--ansys-bc", "--inherit-living", "--transient-c", "1e8",
        "--only", "separate_ab"]


def lhs(n, rng):
    """Latin hypercube in [0, 1]^2"""
    u = (np.arange(n)[:, None] + rng.random((n, 2))) / n
    for j in range(2):
        u[:, j] = u[rng.permutation(n), j]
    return u


def one(i, a0, gc, nf, tag):
    out = HERE / "runs" / f"{tag}{i:03d}_a{a0:.4g}_gc{gc:.4g}.json"
    if out.exists():
        return out
    args = ["--nf", str(nf)] + BASE + ["--alpha0", f"{a0:.6g}", "--gc", f"{gc:.6g}"]
    log = subprocess.run([sys.executable, str(SOLVER), *args], cwd=ROOT, capture_output=True, text=True)
    pre = "" if nf == 32 else f"nf{nf}_"
    made = sorted(RES.glob(f"coupled_living_consumption_{pre}front1_*a{a0:.6g}_gab1_gc{gc:.6g}_*transc1e+08*separate_ab.json"),
                  key=lambda p: p.stat().st_mtime)
    if log.returncode or not made:
        print(f"run {i} failed:\n{log.stderr[-800:]}", flush=True)
        return None
    shutil.move(str(made[-1]), out)            # keep the thesis results folder free of pilot runs
    print(f"run {i}: alpha*_0 = {a0:.3g}, g l^2/d = {gc:.3g} done", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--nf", type=int, default=32)
    ap.add_argument("--tag", default="p", help="file prefix: p = first design (32^3), q = refined (32^3), f = 64^3")
    ap.add_argument("--a0-range", type=float, nargs=2, default=[10.0, 300.0])
    ap.add_argument("--reuse", default=None, help="design.json of an earlier set: run the same points")
    a = ap.parse_args()
    (HERE / "runs").mkdir(exist_ok=True)
    if a.reuse:
        d = json.loads(Path(a.reuse).read_text()); a0, gc = np.array(d["alpha0"]), np.array(d["gc"]); a.n = len(a0)
    else:
        u = lhs(a.n, np.random.default_rng(a.seed))
        lo, hi = np.log10(a.a0_range)
        a0 = 10 ** (lo + u[:, 0] * (hi - lo))         # log-uniform
        gc = 0.5 + 1.5 * u[:, 1]                      # 0.5 .. 2
    with ThreadPoolExecutor(a.jobs) as ex:
        list(ex.map(lambda k: one(k, a0[k], gc[k], a.nf, a.tag), range(a.n)))
    (HERE / f"design_{a.tag}.json").write_text(json.dumps({"alpha0": a0.tolist(), "gc": gc.tolist(), "seed": a.seed}, indent=1))


if __name__ == "__main__":
    main()

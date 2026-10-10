#!/usr/bin/env python3
"""cal5_v2.py -- verification V2 of the P3 design (tmcmc202601
docs/fem_coupling/P3_point_model_alpha_design.md): the point model inside the
element (mode 9, comp_trace.csv of one integration point, local nutrient off)
against the paper pipeline's own solver (hamilton_ode_jax_paper.simulate_0d_full,
the MAP theta, the Day-1 state, dt = 1e-4), step by step.

    python ansys_usermat/apdl/cal5_v2.py RUN --eco eco_CH_cal.txt
        --paper-module <path to hamilton_ode_jax_paper.py> [--elem 220]
        [--workdir F:\\biofilm_upf_cal5] [--tol 1e-10]

The element's coupling step is s * deltim in point-model time, split into
ceiling(s deltim / 1e-4) sub-steps of equal length; the deck is made so that
this is an integer number of 1e-4 steps (s = 0.2375, deltim = 1/95: 25). The
paper solver is run with the same number of steps and compared at every
coupling step. With the element's local nutrient off (prop(33) = 0) and
newton 12 1e-20 in the configuration file, the two are the same algorithm, so
the difference is round-off (the order of the sums), and V4 (one element,
one point) is the same comparison: a point of the element never sees another.

Prints the largest absolute difference over the 12 state components per
coupling step and overall, and exits 1 above --tol.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np


def load_eco(path):
    lines = [l.strip() for l in open(path) if l.strip()]
    theta = [float(x) for x in lines[3].split()]
    phi_init = None
    for l in lines[4:]:
        w = l.split()
        if w[0] == "phi_init":
            phi_init = [float(x) for x in w[1:6]]
    return theta, phi_init


def seed_rows(path, elem):
    rows = []
    with open(path) as f:
        for r in csv.reader(f):
            r = [x.strip() for x in r]
            if len(r) != 36 or r[0] != str(elem) or r[1] != "1":
                continue
            rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--eco", required=True)
    ap.add_argument("--paper-module", required=True)
    ap.add_argument("--elem", type=int, default=220)
    ap.add_argument("--workdir", default=os.environ.get("BIOFILM_WORKDIR") or r"F:\biofilm_upf_cal5")
    ap.add_argument("--tol", type=float, default=1e-10)
    a = ap.parse_args()
    w = Path(a.workdir)
    theta, phi_init = load_eco(w / a.eco)
    if phi_init is None:
        sys.exit("the configuration file has no phi_init line")
    rows = seed_rows(w / f"comp_trace_{a.run}.csv", a.elem)
    if not rows:
        sys.exit(f"no rows of element {a.elem}, integration point 1")
    # one row per coupling step: the last row of each (ldstep, isubst)
    last = {}
    for r in rows:
        last[(int(r[2]), int(r[3]))] = r
    keys = sorted(last)
    nsub = [int(last[k][4]) for k in keys]
    dt_pm = [float(last[k][7]) for k in keys]
    sub = nsub[0]
    if any(n != sub for n in nsub):
        sys.exit(f"sub-steps per coupling step vary: {sorted(set(nsub))}")
    dts = dt_pm[0] / sub
    print(f"{len(keys)} coupling steps, {sub} sub-steps each, dt_pm/sub = {dts:.16e}")

    spec = importlib.util.spec_from_file_location("hamilton_ode_jax_paper", a.paper_module)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    import jax.numpy as jnp
    g = mod.simulate_0d_full(jnp.array(theta), n_steps=len(keys) * sub, dt=dts,
                             phi_init=jnp.array(phi_init), K_hill=0.0, c_const=25.0,
                             alpha_const=0.0)
    g = np.asarray(g)
    g0 = np.array([float(x) for x in last[keys[0]][12:24]])
    d0 = np.abs(g0 - g[0]).max()
    print(f"start state: largest difference {d0:.3e}")
    worst = d0
    for i, k in enumerate(keys):
        gn = np.array([float(x) for x in last[k][24:36]])
        d = np.abs(gn - g[(i + 1) * sub]).max()
        worst = max(worst, d)
        if i < 3 or i == len(keys) - 1 or d > a.tol:
            print(f"  step {i + 1:4d} (ld {k[0]}, sub {k[1]}): t_pm {(i + 1) * sub * dts:.6f}  "
                  f"max |element - paper| {d:.3e}   phi {gn[:5].round(6)}")
    print(f"largest difference over all steps: {worst:.3e}  ({'OK' if worst <= a.tol else 'ABOVE'} tol {a.tol:g})")
    sys.exit(0 if worst <= a.tol else 1)


if __name__ == "__main__":
    main()

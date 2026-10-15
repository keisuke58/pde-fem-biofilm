#!/usr/bin/env python3
"""a12_check.py -- A-12 (CLOUD_TO_IKMHIWI03_2026-10-14_C.ja.md): the CH MAP replaced by
the raw-data one (Tmcmc202601 6e18e0b).

    python ansys_usermat/apdl/a12_check.py [--results DIR]

1. 0D check: the scheme (c) run w8_CH_cal_k0882_raw runs the point model free at each
   Gauss point from the Day-1 state (phi_init.json, sum 1.001, psi 0.999). Its state at
   Gauss point 1 of the first traced seed element after every step, against
   ecology_jax.ecology_substeps started from the same state (phi_0 = 1 - sum phi =
   -0.001, as the estimator starts) with the same theta, s = 0.2375 and 25 inner
   steps per coupling step. Prints the largest difference over phi_i and psi_i in
   step 1 and from Day 3 (T* = 0.1) on.
2. old vs new MAP, 8^3, scheme (a) (the Ch. 4 five-species run) and (c): seed
   alpha - 1 mean, von Mises mean, and the composition at T* = 1 (seed mean of phi_i
   and of phi_i psi_i, as compare_cal5.seed_summary).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "ansys_usermat" / "coupling"))
sys.path.insert(0, str(ROOT))
from compare_cal5 import seed_summary  # noqa: E402

SP = ("So", "An", "Vei", "Fn", "Pg")
PHI_INIT_NEW = (0.761860974108381, 0.024780725603456608, 0.11331891049687755,
                0.1000393897912847, 0.001)


def check_0d(res, job="w8_CH_cal_k0882_raw", s=0.2375, n_sub=25):
    import ecology_jax as ej  # noqa: E402  (JAX)
    pm = json.load(open(res / f"{job}_pm.json"))
    C = {k: i for i, k in enumerate(pm["columns"])}
    rows = pm["seed_trace"]
    e0 = rows[0][C["elem"]]
    by = {}
    for r in rows:
        if r[C["elem"]] == e0 and r[C["ip"]] == "1":
            by.setdefault(int(r[C["isubst"]]), r)
    steps = sorted(by)
    gold0 = np.array([float(by[steps[0]][C[f"gold_{i}"]]) for i in range(1, 13)])
    theta = np.array([float(x) for x in pm["theta"]])
    g = np.zeros(12)
    g[0:5] = PHI_INIT_NEW
    g[5] = 1.0 - sum(PHI_INIT_NEW)
    g[6:11] = 0.999
    print(f"0D check, {job}, element {e0} GP 1, {len(steps)} steps")
    print("  element's state before step 1:", " ".join(f"{x:.6g}" for x in gold0[0:6]))
    print("  python start (estimator):     ", " ".join(f"{x:.6g}" for x in g[0:6]))
    nstep = len(steps)
    d1, dlate, dmax = None, 0.0, 0.0
    for k, st in enumerate(steps):
        g = np.asarray(ej.ecology_substeps(g, theta, s / nstep, n_sub)[0], float).copy()
        ge = np.array([float(by[st][C[f"gnew_{i}"]]) for i in range(1, 13)])
        d = np.max(np.abs(ge[[0, 1, 2, 3, 4, 6, 7, 8, 9, 10]] - g[[0, 1, 2, 3, 4, 6, 7, 8, 9, 10]]))
        dmax = max(dmax, d)
        if k == 0:
            d1 = d
        if (k + 1) / nstep >= 0.1:
            dlate = max(dlate, d)
    print(f"  max |element - python| over phi_i, psi_i: step 1 {d1:.3e}, from Day 3 on {dlate:.3e}, all steps {dmax:.3e}")
    print("  at T* = 1, element phi:", " ".join(f"{x:.5f}" for x in ge[0:5]),
          "| python phi:", " ".join(f"{x:.5f}" for x in g[0:5]))


def table(res):
    print("\n| run | MAP | seed alpha-1 mean | seed vM mean [Pa] | phi (seed mean, T* = 1) | phi psi (seed mean, T* = 1) |")
    print("|---|---|---|---|---|---|")
    for job, m in (("w8_CH_g1_s015", "old 2d8d99e"), ("w8_CH_g1_s015_raw", "new 6e18e0b"),
                   ("w8_CH_cal_k0882", "old 2d8d99e"), ("w8_CH_cal_k0882_raw", "new 6e18e0b")):
        s = seed_summary(res / f"{job}.json", res / f"{job}_pm.json")
        print(f"| {job} | {m} | {s['alpha_mean']:.4e} | {s['seqv_mean']:.4e} | "
              + " / ".join(f"{x:.3f}" for x in s["phi"]) + " | "
              + " / ".join(f"{x:.3f}" for x in s["phipsi"]) + " |")
    print("(species order " + ", ".join(SP) + ")")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=str(HERE / "results" / "2026-10-cal5"))
    res = Path(ap.parse_args().results)
    check_0d(res)
    table(res)


if __name__ == "__main__":
    main()

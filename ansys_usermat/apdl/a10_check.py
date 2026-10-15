#!/usr/bin/env python3
"""a10_check.py -- A-10 (CLOUD_TO_IKMHIWI03_2026-10-14_C.ja.md): growth from the point
model's own increase (prop(28) = 10, no kappa) in the element.

    python ansys_usermat/apdl/a10_check.py JOB [--ref REF] [--results DIR]

1. check: for the first Gauss point of the traced seed element, the element's
   ln J = 3 ln(1 + alpha) after every step against
   JAXFEM/growth_from_point_model.run_with_phi given the same phi_used series and
   the element's first state (gold); the largest difference should be ~1e-10 or
   below. Exact only without the local nutrient (prop(33) = 0): the Python model
   has no c_rel.
2. seed alpha - 1 mean and von Mises mean at T* = 1, next to REF (default the
   (a) run w8_<cond>_g1_s015);
3. the volume ratio alpha^3 of that Gauss point at Day 3, 6, 10, 15, 21
   (T* = 1 is Day 21, Day 1 at T* = 0), relative to Day 3, next to the measured
   CH volumes (0.222, 0.224, 0.344, 0.382, 0.400 x 10^6 um^3).
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
sys.path.insert(0, str(ROOT / "JAXFEM"))
from compare_cal5 import seed_summary  # noqa: E402

DAYS = (3, 6, 10, 15, 21)
CH_VOL = (0.222, 0.224, 0.344, 0.382, 0.400)


def seed_point_series(pm):
    """phi_used, alpha_new and the first gold state per step, GP 1 of the first seed element."""
    C = {k: i for i, k in enumerate(pm["columns"])}
    rows = pm["seed_trace"]
    e0 = rows[0][C["elem"]]
    by_step = {}
    for r in rows:
        if r[C["elem"]] == e0 and r[C["ip"]] == "1":
            by_step.setdefault(int(r[C["isubst"]]), r)      # first call of the step
    steps = sorted(by_step)
    phi = np.array([float(by_step[s][C["phi_used"]]) for s in steps])
    alpha = np.array([float(by_step[s][C["alpha_new"]]) for s in steps])
    g0 = np.array([float(by_step[steps[0]][C[f"gold_{i}"]]) for i in range(1, 13)])
    return e0, steps, phi, alpha, g0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("job")
    ap.add_argument("--ref", default="")
    ap.add_argument("--results", default=str(HERE / "results" / "2026-10-cal5"))
    ap.add_argument("--phi-cap", type=float, default=0.9)
    a = ap.parse_args()
    res = Path(a.results)
    pm = json.load(open(res / f"{a.job}_pm.json"))
    e0, steps, phi, alpha, g0 = seed_point_series(pm)
    lnj_el = 3 * np.log1p(alpha)
    print(f"{a.job}: element {e0} GP 1, {len(steps)} steps (substeps {steps[0]}..{steps[-1]})")

    import growth_from_point_model as gp  # noqa: E402  (imports JAX)
    lnj_py = gp.run_with_phi(phi, g0, np.array([float(x) for x in pm["theta"]]), phi_cap=a.phi_cap)[1:]
    d = np.abs(lnj_py - lnj_el)
    print(f"check: max |ln J (element) - ln J (python)| = {d.max():.3e} "
          f"(at step {int(d.argmax()) + 1}); ln J at T* = 1: element {lnj_el[-1]:.10f}, python {lnj_py[-1]:.10f}")

    cond = a.job.split("_")[1]
    ref = a.ref or f"w8_{cond}_g1_s015"
    print("\n| run | seed alpha-1 mean | seed vM mean [Pa] | vM/(alpha-1) [Pa] |")
    print("|---|---|---|---|")
    for job in (ref, a.job):
        s = seed_summary(res / f"{job}.json", res / f"{job}_pm.json")
        print(f"| {job} | {s['alpha_mean']:.4e} | {s['seqv_mean']:.4e} | {s['seqv_mean'] / s['alpha_mean']:.4f} |")

    # Day d at T* = (d - 1)/20; ln J at the step ends, linear between them
    t = np.concatenate([[0.0], np.arange(1, len(lnj_el) + 1) / len(lnj_el)])
    lj = np.concatenate([[0.0], lnj_el])
    jd = {d: float(np.exp(np.interp((d - 1) / 20.0, t, lj))) for d in DAYS}
    print("\n| Day | alpha^3 = J (GP 1) | J / J(Day 3) | measured CH [10^6 um^3] | measured / Day 3 |")
    print("|---|---|---|---|---|")
    for d, v in zip(DAYS, CH_VOL):
        print(f"| {d} | {jd[d]:.4f} | {jd[d] / jd[3]:.4f} | {v:.3f} | {v / CH_VOL[0]:.3f} |")


if __name__ == "__main__":
    main()

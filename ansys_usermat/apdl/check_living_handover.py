#!/usr/bin/env python3
"""check_living_handover.py -- does USSFin use the living share L that the material routine wrote? (addendum 9,
after the A-22 bug that check_living_trace.py could not see.)

Inputs: ussfin_L_trace.csv from the debug build F:\\biofilm_upf_living4dbg (USSFin writes, at every call, its time,
element, Gauss point, the L it uses and the point's coordinates, for elements 220 and 25), the point-model trace
<job>_pm.json (L written by the material routine at integration point 1 after each step) and fields_<job>.json
(element centroids).
For each USSFin call at the end of an increment, its L at integration point 1 is compared with the routine's L of the
same increment and of the increment before; the 8 points' mean coordinates are compared with the element centroid.

    python ansys_usermat/apdl/check_living_handover.py JOB CSV [--start-step 0.0025] [--dt 0.025]
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parent / "results" / "2026-10-living"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("job")
    ap.add_argument("csv")
    ap.add_argument("--start-step", type=float, default=0.0025)
    ap.add_argument("--dt", type=float, default=0.025)
    a = ap.parse_args()
    pm = json.load(open(RES / f"{a.job}_pm.json"))
    C = {k: i for i, k in enumerate(pm["columns"])}
    fl = json.load(open(RES / f"fields_{a.job}.json"))
    cen = {e: (x, y, z) for e, x, y, z in zip(fl["elem"], fl["cx"], fl["cy"], fl["cz"])}

    def t_end(ld, ist):
        return a.start_step if ld == 1 else a.start_step + a.dt * ist

    um = {}                                   # (elem, t) -> L after the step (last row of the increment)
    for r in pm["seed_trace"]:
        e, ld, ist = int(r[0]), int(r[C["ldstep"]]), int(r[C["isubst"]])
        um[(e, round(t_end(ld, ist), 9))] = float(r[C["L"]])
    us = {}                                   # (elem, gp, t) -> L used by USSFin; coordinates
    xyz = {}
    for r in csv.reader(open(a.csv)):
        t, e, g, L, x, y, z = float(r[0]), int(r[1]), int(r[2]), float(r[3]), float(r[4]), float(r[5]), float(r[6])
        us[(e, g, round(t, 9))] = L
        xyz.setdefault((e, g), (x, y, z))
    print(f"{a.job}: {len(us)} USSFin values, {len(um)} routine values")
    for e in sorted({k[0] for k in us}):
        m = np.mean([xyz[(e, g)] for g in range(1, 9) if (e, g) in xyz], axis=0)
        print(f"element {e}: mean of the 8 points' coordinates {np.round(m, 6)}, centroid {np.round(cen[e], 6)}, "
              f"difference {np.abs(m - np.array(cen[e])).max():.2e}")
        # USSFin's time may be the start or the end of its increment: try both readings and report the one that fits
        ts = sorted(t for (ee, g, t) in us if ee == e and g == 1)
        tm = sorted(t for (ee, t) in um if ee == e)
        for name, shift in (("USSFin time = end of its increment", 0), ("USSFin time = start of its increment", 1)):
            same, lag, n = 0.0, 0.0, 0
            for t in ts:
                # the routine's increments in time order; with shift 1 USSFin's time t is the end of the one before
                k = min(range(len(tm)), key=lambda j: abs(tm[j] - t)) if tm else None
                if k is None or abs(tm[k] - t) > 1e-9 and shift == 0:
                    continue
                if shift == 1:
                    nxt = [j for j in range(len(tm)) if tm[j] > t + 1e-9]
                    if not nxt:
                        continue
                    k = nxt[0]
                Lu, Lm = us[(e, 1, t)], um[(e, tm[k])]
                n += 1
                same = max(same, abs(Lu - Lm))
                if k > 0:
                    lag = max(lag, abs(Lu - um[(e, tm[k - 1])]))
            print(f"  {name}: {n} increments, max |L(USSFin) - L(routine, same increment)| = {same:.2e}, "
                  f"one increment behind = {lag:.2e}")


if __name__ == "__main__":
    main()

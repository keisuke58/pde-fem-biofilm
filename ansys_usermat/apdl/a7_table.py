#!/usr/bin/env python3
"""a7_table.py -- the A-7 report (archive/CLOUD_TO_IKMHIWI03_2026-10-14_C.ja.md):
compression slows the (c) growth, g = 1/(1 + max(p,0)/p_ref).

    python ansys_usermat/apdl/a7_table.py [--results DIR]

1. the off run (w8_CH_cal_k0882_gp0, exe F:\\biofilm_upf_cal5gp) against
   w8_CH_cal_k0882 (exe F:\\biofilm_upf_cal5sc): every number of all_stress
   equal, or the largest difference;
2. a table, seed elements at T* = 1: alpha - 1 mean (UMAT variable), von Mises
   mean [Pa], change of both against w8_CH_cal_k0882, the mean compression
   p = -tr(sigma)/3 [Pa], and alpha - 1 at the seed's centre and edge (the
   seed element nearest to / farthest from the seed's centroid, mean of the
   four nearest / farthest), with the lysis test w8_CH_cal_k0882_sc (14 Oct)
   as the last row when it is there.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

REF = "w8_CH_cal_k0882"
RUNS = ["w8_CH_cal_k0882_gp1", "w8_CH_cal_k0882_gp10", "w8_CH_cal_k0882_sc"]
MPA_TO_PA = 1e6


def load(res, job):
    f = res / f"{job}.json"
    return json.load(open(f)) if f.exists() else None


def seed_stats(d):
    a = d["all_stress"]
    el = np.array(a["elem"], dtype=int)
    seed = np.isin(el, np.array(d["seed_BIOFILM1"], dtype=int))
    al = np.array(a["alpha"])[seed]
    vm = np.array(a["seqv"])[seed] * MPA_TO_PA
    p = -(np.array(a["sx"]) + np.array(a["sy"]) + np.array(a["sz"]))[seed] / 3 * MPA_TO_PA
    xyz = np.stack([np.array(a[k])[seed] for k in ("cx", "cy", "cz")], axis=1)
    r = np.linalg.norm(xyz - xyz.mean(axis=0), axis=1)
    o = np.argsort(r)
    return dict(alpha=al.mean(), vm=vm.mean(), p=p.mean(),
                centre=al[o[:4]].mean(), edge=al[o[-4:]].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=str(Path(__file__).parent / "results" / "2026-10-cal5"))
    res = Path(ap.parse_args().results)
    ref, off = load(res, REF), load(res, REF + "_gp0")
    if off is None:
        print(f"{REF}_gp0: not there yet")
    else:
        worst = 0.0
        for k, v in ref["all_stress"].items():
            x, y = np.array(v, dtype=float), np.array(off["all_stress"][k], dtype=float)
            worst = max(worst, float(np.max(np.abs(x - y))))
        print(f"off vs {REF}: " + ("all_stress identical to the last digit" if worst == 0.0
                                   else f"largest difference {worst:.3e}"))
    s0 = seed_stats(ref)
    print("\n| run | seed alpha-1 mean | vs ref | seed vM mean [Pa] | vs ref | seed p mean [Pa] "
          "| alpha-1 centre | alpha-1 edge | edge/centre |")
    print("|---|---|---|---|---|---|---|---|---|")
    for job in [REF] + RUNS:
        d = load(res, job)
        if d is None:
            print(f"| {job} | not there yet |||||||")
            continue
        s = seed_stats(d)
        print(f"| {job} | {s['alpha']:.4e} | {s['alpha'] / s0['alpha'] - 1:+.1%} | {s['vm']:.4e} | "
              f"{s['vm'] / s0['vm'] - 1:+.1%} | {s['p']:.3e} | {s['centre']:.4e} | {s['edge']:.4e} | "
              f"{s['edge'] / s['centre']:.3f} |")


if __name__ == "__main__":
    main()

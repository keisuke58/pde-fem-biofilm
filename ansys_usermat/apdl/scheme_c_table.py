#!/usr/bin/env python3
"""scheme_c_table.py -- the report table of CLOUD_TO_IKMHIWI03_2026-10-14_C.ja.md
for scheme (c) (growth from the point model's living amount) and its reference.

    python ansys_usermat/apdl/scheme_c_table.py REF RUN [RUN ...]
        [--results ansys_usermat/apdl/results/2026-10-cal5]

REF and RUN are job names (w8_CH_g1_s015, w8_CH_cal_k0882_ab, ...). For the seed
elements (seed_BIOFILM1) at T* = 1: alpha - 1 mean (the UMAT variable), von Mises
mean [Pa], their ratio, the change of the von Mises mean against REF, and the sum
of psi_i and of phi_i psi_i at the first Gauss point of the traced seed elements
(compare_cal5.seed_summary). Prints a Markdown table.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_cal5 import seed_summary  # noqa: E402


def row(res, job):
    s = seed_summary(res / f"{job}.json", res / f"{job}_pm.json")
    p = json.load(open(res / f"{job}_pm.json"))
    seed = set(int(e) for e in json.load(open(res / f"{job}.json"))["seed_BIOFILM1"])
    c = p["columns"]
    i0 = c.index("gnew_1")
    g = np.array([[float(x) for x in r[i0:i0 + 12]] for r in p["last_rows"]
                  if int(r[0]) in seed and r[1] == "1"])
    s["sum_psi"] = g[:, 6:11].sum(axis=1).mean()
    s["sum_phipsi"] = (g[:, 0:5] * g[:, 6:11]).sum(axis=1).mean()
    s["alpha_star"] = p.get("constants", {}).get("alpha_star", "")
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--results", default=str(Path(__file__).parent / "results" / "2026-10-cal5"))
    a = ap.parse_args()
    res = Path(a.results)
    ref = row(res, a.ref)
    print("| run | alpha* | seed alpha-1 mean | seed vM mean [Pa] | vM/(alpha-1) [Pa] | vM vs ref | "
          "sum psi (seed GP) | sum phi psi (seed GP) |")
    print("|---|---|---|---|---|---|---|---|")
    for job in [a.ref] + a.runs:
        s = ref if job == a.ref else row(res, job)
        d = s["seqv_mean"] / ref["seqv_mean"] - 1
        print(f"| {job} | {s['alpha_star']} | {s['alpha_mean']:.4e} | {s['seqv_mean']:.4e} | "
              f"{s['seqv_mean'] / s['alpha_mean']:.4e} | {d:+.1%} | {s['sum_psi']:.4f} | {s['sum_phipsi']:.4f} |")


if __name__ == "__main__":
    main()

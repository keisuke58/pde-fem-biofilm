#!/usr/bin/env python3
"""Compare an ANSYS run on the 8^3 example model with the same run on 16^3.

    python ansys_usermat/apdl/compare_mesh.py all_stress_8.csv all_stress_16.csv
        [--grid8 8] [--track 220]

Both CSVs come from post_all_stress.mac (elem, seqv, sx, sy, sz, alpha,
centroid). A CSV without centroid columns needs --grid8 8 (the 8^3 one from
before 2 Oct; same assumed numbering as plot_3d.py). Stresses MPa -> Pa.

Printed per mesh: seed size (alpha above half its maximum), the seed averages
of alpha, von Mises and mean stress, the largest von Mises outside the seed,
and the 16^3 / 8^3 ratios next to the ratios of the Python study
(MESH_STUDY.md: 0.63 for the seed von Mises, 0.67 for the seed mean stress).
--track: one 8^3 element, compared with the average of the 16^3 elements
inside it (found by centroid).

What to look for:
  - the seed should have 32 elements on 8^3 and 256 on 16^3, with the same
    mean alpha (growth does not depend on the mesh);
  - the ratios near the Python ones: then the 8^3 stresses of the 5 Oct deck
    are too large by that factor, and the thesis quotes the 16^3 averages.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_3d import read  # noqa: E402

PY_RATIO = {"vm": 3.59e-5 / 5.70e-5, "p": 2.52e-5 / 3.77e-5}


def summary(col):
    ins = col["alpha"] > 0.5 * col["alpha"].max()
    vm, p = col["seqv"] * 1e6, col["p"] * 1e6
    return {"n": len(vm), "seed": int(ins.sum()), "alpha": col["alpha"][ins].mean(),
            "vm": vm[ins].mean(), "p": p[ins].mean(), "vm_out": vm[~ins].max()}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("csv8")
    ap.add_argument("csv16")
    ap.add_argument("--grid8", type=int, default=None)
    ap.add_argument("--track", type=int, default=220)
    a = ap.parse_args(argv)
    c8 = read(a.csv8, a.grid8, 2.0)
    c16 = read(a.csv16)
    s8, s16 = summary(c8), summary(c16)
    print(f"{'':22s} {'8^3':>11s} {'16^3':>11s} {'16/8':>7s} {'Python 16/8':>12s}")
    for k, name in (("n", "elements"), ("seed", "seed elements")):
        print(f"{name:22s} {s8[k]:11d} {s16[k]:11d}")
    for k, name in (("alpha", "seed alpha (column)"), ("vm", "seed von Mises [Pa]"),
                    ("p", "seed mean stress [Pa]"), ("vm_out", "max vM outside [Pa]")):
        py = f"{PY_RATIO[k]:12.2f}" if k in PY_RATIO else ""
        print(f"{name:22s} {s8[k]:11.3e} {s16[k]:11.3e} {s16[k] / s8[k]:7.2f} {py}")
    if s8["seed"] != 32 or s16["seed"] != 256:
        print("WARNING: seed size is not 32 / 256: check the alpha column and the deck")

    e = np.flatnonzero(c8["elem"].astype(int) == a.track)
    if e.size:
        c = np.array([c8["cx"][e[0]], c8["cy"][e[0]], c8["cz"][e[0]]])
        cen16 = np.stack([c16["cx"], c16["cy"], c16["cz"]], 1)
        sub = np.all(np.abs(cen16 - c) < 0.125 + 1e-9, axis=1)
        ids = c16["elem"][sub].astype(int).tolist()
        print(f"\nelement {a.track} (8^3) against the {sub.sum()} 16^3 elements inside it {ids}:")
        for k, f in (("alpha", 1), ("seqv", 1e6), ("p", 1e6)):
            v8, v16 = c8[k][e[0]] * f, c16[k][sub].mean() * f
            print(f"  {k:6s} {v8:11.3e} {v16:11.3e}  ratio {v16 / v8 if v8 else float('nan'):.2f}")


if __name__ == "__main__":
    main()

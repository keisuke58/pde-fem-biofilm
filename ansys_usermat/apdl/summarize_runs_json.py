#!/usr/bin/env python3
"""summarize_runs_json.py -- tables from the run JSON files of export_runs_json.py,
without ANSYS or F: (works in a cloud session).

    python ansys_usermat/apdl/summarize_runs_json.py [DIR] [--runs a,b,c]

DIR defaults to ansys_usermat/apdl/results/2026-10-05_ansys. For every run:
  stress  seed = BIOFILM1 of the deck. Seed mean von Mises and mean stress p,
          and p by layer around the seed (layer k = elements whose centroid
          lies k*0.25 mm outside the seed's elements, max norm; the 8^3 element
          size, so the layers mean the same on every mesh), with the share of
          elements in tension. Pa.
  alpha   the material routine's variable (alpha - 1 in the thesis): seed
          interior / surface means (surface = a face neighbour outside the
          seed) and the step between them in %.
  share   point-model runs only: phi_1/(phi_1+phi_2) where the point model ran
          (phi_1+phi_2 < 0.95; never-run points keep 0.5/0.5), in the seed and
          by layer, and the nutrient Nut1 in the seed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

DEFAULT = Path(__file__).resolve().parent / "results" / "2026-10-05_ansys"


def layers(cen, seed_mask, h_seed=0.25):
    """distance class of every element from the seed region (0 = seed)."""
    sc = cen[seed_mask]
    h = np.min(np.diff(np.unique(np.round(cen[:, 0], 9))))
    d = np.full(len(cen), np.inf)
    for b in sc:
        d = np.minimum(d, np.max(np.maximum(np.abs(cen - b) - h / 2, 0), axis=1))
    lay = np.ceil(d / h_seed - 1e-9).astype(int)
    lay[seed_mask] = 0
    return lay, h


def surface(cen, seed_mask, h):
    key = {tuple(np.round(c / h, 3)) for c in cen[seed_mask]}
    out = np.zeros(len(cen), bool)
    for i in np.flatnonzero(seed_mask):
        c = np.round(cen[i] / h, 3)
        for k in range(3):
            for s in (-1, 1):
                n = c.copy(); n[k] = round(n[k] + s, 3)
                if tuple(n) not in key:
                    out[i] = True
    return out


def summarize(rec):
    a = rec.get("all_stress")
    if not a or "cx" not in a:
        print(f"== {rec['run']}: no all_stress with centroids"); return
    elem = np.array(a["elem"], int)
    cen = np.stack([a["cx"], a["cy"], a["cz"]], 1)
    seed = np.isin(elem, rec.get("seed_BIOFILM1", []))
    lay, h = layers(cen, seed)
    vm = np.array(a["seqv"]) * 1e6
    p = (np.array(a["sx"]) + np.array(a["sy"]) + np.array(a["sz"])) / 3 * 1e6
    al = np.array(a["alpha"])
    surf = surface(cen, seed, h)
    inner = seed & ~surf
    print(f"== {rec['run']}: {len(elem)} elements, seed {seed.sum()}")
    print(f"   seed vM {vm[seed].mean():.3e}  p {p[seed].mean():.3e} Pa;  layer p (tension share): "
          + "  ".join(f"{k}: {p[lay == k].mean():+.2e} ({np.mean(p[lay == k] > 0):.2f})" for k in (1, 2, 3) if (lay == k).any()))
    ai = al[inner].mean() if inner.any() else np.nan
    print(f"   alpha seed interior {ai:.4e}  surface {al[surf].mean():.4e}  step {(ai - al[surf].mean()) / ai * 100:.2f} %")
    nf = rec.get("nut_field")
    if nf:
        ne = np.array(nf["elem"], int)
        idx = {e: i for i, e in enumerate(elem)}
        order = np.array([idx[e] for e in ne])
        p1, p2, nut = np.array(nf["phi1"]), np.array(nf["phi2"]), np.array(nf["nut1"])
        s = p1 + p2
        ran = (s > 1e-6) & (s < 0.95)
        share = np.where(ran, p1 / np.where(ran, s, 1), np.nan)
        L = lay[order]; S = seed[order]
        print(f"   point model ran at {ran.sum()} elements; Nut1 in seed {nut[S].min():.3f}..{nut[S].max():.3f}")
        print("   share phi1/(phi1+phi2): seed "
              f"{np.nanmin(share[S]):.4f}..{np.nanmax(share[S]):.4f} (mean {np.nanmean(share[S]):.4f})"
              + "".join(f";  layer {k} mean {np.nanmean(share[(L == k) & ran]):.4f}" for k in (1, 2) if ((L == k) & ran).any()))


def main(argv):
    runs = None
    if "--runs" in argv:
        i = argv.index("--runs"); runs = argv[i + 1].split(","); del argv[i:i + 2]
    d = Path(argv[0]) if argv else DEFAULT
    files = sorted(d.glob("*.json")) if runs is None else [d / f"{r}.json" for r in runs]
    for f in files:
        if f.exists():
            summarize(json.loads(f.read_text()))
        else:
            print(f"== {f.stem}: no file")


if __name__ == "__main__":
    main(sys.argv[1:])

#!/usr/bin/env python3
"""Klempt 2024 test case 4.2 "biofilm on nutrients" (Fig. 5-7, Table 4).

Set-up as the paper states it (sec. 4.2, Fig. 5): 20 um cube, 1 um mesh; the
nodes of the bottom face held at c = 1; in the plane directly above, a
centred disk of diameter 5 um set to phi = 1; Table 2 otherwise, with
g = 1e8 ("high", fills the cube by 13 % of the time) or g = 1e10 ("low").
Targets: Fig. 7 read off the PDF by colour (klempt2024_digitize.py), not by
eye; the earlier eye reading of the "high" phi curve was far too low early.

The same choices as test case 1 (klempt2024_case1_bc.py) are varied:
consumption g phi (Eq. 24/35 as printed) or g phi c (first order); the front
term as printed (an advection) or growth on both faces |grad phi . n_c|; the
seed an initial value or held. Everything else is klempt2024_quantitative.

A consistency check that needs no simulation: once the "high" cube is full
(phi = 1), the quasi-static nutrient with c = 1 at the bottom and zero flux
at the top has a closed form. First-order consumption gives
c = cosh((L - z)/lam)/cosh(L/lam), lam = sqrt(d/g) = 10 um, mean
tanh(2)/2 = 0.482; the paper's plateau is 0.489. Zero-order consumption gives
a negative c (mean -0.33 unclipped, 0.139 clipped). So Fig. 7 was computed
with first-order consumption and Table 2's d and g.

    python JAXFEM/klempt2024_case2.py
    python JAXFEM/klempt2024_case2.py --consumption first_order --growth abs --k-monod 0.01

Result (2026-10-02), largest |difference| to Fig. 7 over T* in [0, 1]:

  set-up                                        high phi  high c   low phi  low c
  Eq. 34/35 as printed                            0.98     0.49     0.29    0.37
  first order, Eq. 34 as printed (advection)      0.98     0.49     0.28    0.49
  first order, both-face growth                   0.90     0.46     0.23    0.43
  DIAGNOSTIC first order, both faces, k = 0.01    0.76     0.37     0.24    0.12
  DIAGNOSTIC the same, seed 9 layers (= 0.024)    0.53     0.18     0.27    0.34

  - Not reproduced. As printed, Eq. 34 moves the biofilm towards the nutrient
    (down, into the source), while Table 4 shows it growing upwards, away
    from it; only both-face growth grows upwards at all.
  - The "high" case fills 20 um by T* ~ 0.13, a front speed near 150 um/T*.
    The front term's speed is at most r c/(k+c) = 50 um/T* with Table 2
    (100 with k -> 0), so Table 2's r and k cannot produce it.
  - The paper's mean phi at T* = 0.01 is 0.021 (high) and 0.024 (low), ten
    times one node layer of the 5 um disk (0.0026): the paper's initial
    biofilm is larger than the text's one plane suggests.
  - What does hold is the nutrient plateau: 0.482 from first-order
    consumption with Table 2's d and g against the paper's 0.489, while the
    printed zero-order form gives 0.139. Together with test case 1 this
    points to first-order consumption in the paper's computations.

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402

RES = Path(__file__).resolve().parent / "klempt2024_results"
OUT = RES / "case2.json"
T_CMP = [0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.13, 0.15, 0.20, 0.30, 0.50, 1.00]
G = {"high": 1e8, "low": 1e10}


def targets(sub):
    d = json.loads((RES / "paper_curves_digitized.json").read_text())
    t = np.array(d["t"])
    out = {}
    for k in ("phi", "c"):
        y = np.array([np.nan if v is None else v for v in d["fig7"][f"{k}_{sub}"]])
        ok = ~np.isnan(y)
        out[k] = np.interp(T_CMP, t[ok], y[ok])
    return out


def setup(layers=1):
    phi0 = np.zeros_like(K.X)
    for k in range(1, layers + 1):
        phi0 += (np.isclose(K.Z, k * K.H) & ((K.X - 10) ** 2 + (K.Y - 10) ** 2 <= 2.5 ** 2))
    return np.clip(phi0, 0.0, 1.0), np.isclose(K.Z, 0.0)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", nargs="+", default=["high", "low"])
    ap.add_argument("--consumption", nargs="+", default=["printed", "first_order"])
    ap.add_argument("--growth", nargs="+", default=["printed", "abs"])
    ap.add_argument("--seed", nargs="+", default=["ic"])
    ap.add_argument("--t-end", type=float, default=1.0)
    ap.add_argument("--k-monod", type=float, default=1.0,
                    help="DIAGNOSTIC ONLY: Monod constant k (Table 2: 1)")
    ap.add_argument("--layers", type=int, default=1,
                    help="seed disk thickness in node layers (the paper: the plane above)")
    a = ap.parse_args(argv)
    try:
        results = json.loads(OUT.read_text())
    except (OSError, ValueError):
        results = {}
    K.K_M = a.k_monod
    phi0, mask = setup(a.layers)
    for sub in a.sub:
        tg = targets(sub)
        for cons in a.consumption:
            for gr in a.growth:
                for sd in a.seed:
                    t0 = time.time()
                    rec, phi, c = B.run_setup(phi0, mask, G[sub], sd, cons, gr,
                                              t_end=a.t_end)
                    tt = np.array(rec["t"])
                    sim = {k: np.interp(T_CMP, tt, rec[k]) for k in ("phi", "c")}
                    m = [i for i, t in enumerate(T_CMP) if t <= a.t_end + 1e-9]
                    diff = {k: float(np.max(np.abs(sim[k][m] - tg[k][m]))) for k in sim}
                    key = (f"{sub}/consumption={cons}/growth={gr}/seed={sd}"
                           + (f"/k={a.k_monod:g}" if a.k_monod != 1 else "")
                           + (f"/layers={a.layers}" if a.layers != 1 else ""))
                    results[key] = {"t": T_CMP, **{f"{k}_sim": [round(float(v), 4) for v in sim[k]]
                                                   for k in sim},
                                    **{f"{k}_paper": [round(float(v), 4) for v in tg[k]] for k in tg},
                                    **{f"max_abs_diff_{k}": round(diff[k], 4) for k in diff}}
                    print(f"\n== {key}  ({time.time() - t0:.0f}s)")
                    print("   t       " + " ".join(f"{t:5.2f}" for t in T_CMP))
                    for k in ("phi", "c"):
                        print(f"   {k:3s} sim  " + " ".join(f"{v:5.3f}" for v in sim[k]))
                        print(f"   {k:3s} pap  " + " ".join(f"{v:5.3f}" for v in tg[k]))
                    print(f"   max |diff| phi {diff['phi']:.3f}  c {diff['c']:.3f}", flush=True)
                    OUT.write_text(json.dumps(results, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()

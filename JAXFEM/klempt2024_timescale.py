#!/usr/bin/env python3
"""Klempt 2024 with Table 2 unchanged and one time scale per simulation.

The joint search (klempt2024_variant_search.py) found no single (k, r, g) for
test cases 4.1 and 4.2. Two readings of the paper, checked on 2026-10-03,
change the question:

  1. The nutrient source of 4.1. Fig. 2 draws a 5 um strip, but Table 3's
     diagonal cut (aspect sqrt(2), the corner with the nutrient at the top
     right) shows c ~ 1 only within 2-3 um of the corner, and the text says
     "In one of the corners". A strip w um wide would show as w*sqrt(2) along
     the top of the cut and w down its side. So the source is the edge line
     (one row of nodes), not a strip.
  2. The time axis. Sec. 4: "a normalized time T* = t/t_ref ... with a
     reference time t_ref which is chosen based on the conditions for growth
     at hand such that T* in [0, 1]". t_ref is chosen per simulation, so the
     rates of Table 2 can hold in solver time while each figure's T* is a
     different multiple of it. With c quasi-static, this is one factor s on
     every rate of Eq. 34/36 (r, beta, k_a); d and g enter only as d/g.

So: Table 2 as printed (d, g, k, r, beta, k_a), growth on both faces
|grad phi . n_c|, first-order consumption g phi c (Fig. 7's plateau), seed as
initial value, and one s per simulation, scanned.

    python JAXFEM/klempt2024_timescale.py        (a few minutes on 4 cores)
        -> JAXFEM/klempt2024_results/timescale.json

Result (2026-10-03), RMS difference to the digitised curves:

  simulation   best s   phi    c
  4.1 (edge)    1.5     0.03   0.05
  4.2 high      10      0.04   0.04
  4.2 low       4-5     0.09   0.05

  - All three curves come within 0.03-0.09 with Table 2 unchanged. The joint
    search, which varied k, r and g with one time axis, did not get below
    0.17 on its worst curve.
  - 4.2 needs a front speed about ten times Table 2's r c/(k+c) on the
    figure's axis; that is what s = 10 says. The early Table 4 panels agree:
    at T* = 0.01 the phi = 0.5 contour already stands about 3.5 um above the
    seed.
  - Weak point: "high" and "low" are drawn on one axis in Fig. 7, yet they
    want s = 10 and s = 4-5. With one s = 7 both are at about 0.09. Whether
    t_ref differed between the two runs is a question for the authors.
  - A 4.1 source 1 um wide (band1) is worse at every s (0.10 at best, c too
    high), which agrees with reading 1.
"""
from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_case2 as C  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402
import klempt2024_variant_search as VS  # noqa: E402

OUT = Path(__file__).resolve().parent / "klempt2024_results" / "timescale.json"
SCAN = {"4.1 edge": [1.0, 1.5, 2.0, 2.5],
        "4.1 band1": [1.0, 1.5, 2.0, 2.5],
        "4.2 high": [5.0, 6.0, 7.0, 8.0, 10.0, 12.0],
        "4.2 low": [4.0, 5.0, 6.0, 7.0, 8.0, 10.0]}


def one(job):
    case, s = job
    K.R, K.BETA, K.K_A = 100.0 * s, 2.0 * s, 1e-3 * s      # every rate x s
    tg = VS.digitized()
    if case.startswith("4.1"):
        phi0, _, _ = K.setup("fig4_edge")
        rec, _, _ = B.run_setup(phi0, B.nutrient_mask(case.split()[1]), 1e8,
                                "ic", "first_order", "abs")
        sc, rms = VS.score(rec, VS.T41, tg["4.1"])
    else:
        phi0, mask = C.setup(1)
        rec, _, _ = B.run_setup(phi0, mask, C.G[case.split()[1]], "ic", "first_order", "abs")
        sc, rms = VS.score(rec, C.T_CMP, tg[case])
    return case, s, {"score": round(sc, 4), "rms": {k: round(v, 4) for k, v in rms.items()},
                     "t": rec["t"], "phi": [round(v, 4) for v in rec["phi"]],
                     "c": [round(v, 4) for v in rec["c"]]}


def main():
    jobs = [(c, s) for c, ss in SCAN.items() for s in ss]
    res = {}
    with Pool(4) as pool:
        for case, s, e in pool.imap_unordered(one, jobs):
            res.setdefault(case, {})[f"{s:g}"] = e
            print(f"{case:10s} s = {s:4g}  score {e['score']:.3f}  phi {e['rms']['phi']:.3f}"
                  f"  c {e['rms']['c']:.3f}", flush=True)
    res = {c: dict(sorted(v.items(), key=lambda kv: float(kv[0]))) for c, v in res.items()}
    OUT.write_text(json.dumps(res, indent=1))
    print("\nbest s per simulation:")
    for c, v in res.items():
        s, e = min(v.items(), key=lambda kv: kv[1]["score"])
        print(f"  {c:10s} s = {s:4s}  phi {e['rms']['phi']:.3f}  c {e['rms']['c']:.3f}")
    print("wrote", OUT)


if __name__ == "__main__":
    main()

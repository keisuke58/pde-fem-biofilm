#!/usr/bin/env python3
"""Klempt 2024 test case 4.1 with the nutrient in a corner block.

Fig. 3 (45 %) shows the biofilm reaching a corner of the cube in a thin
tube, and the text says "In one of the corners". So the source is tried as
the nodes of an a-um cube at the corner (L, L, L), a = 1-4 um, with
Table 2 unchanged, growth on both faces, consumption g phi c and a time
scale s (klempt2024_timescale.py), scored against Fig. 4 (RMS, phi and c
averaged).

    python JAXFEM/klempt2024_corner_scan.py

Result (2026-10-03), best five:
  a = 2 um, s = 2     score 0.022  (phi 0.029, c 0.015)
  a = 2 um, s = 1.5   score 0.054
  a = 3 um, s = 1.5   score 0.056
  a = 1 um, s = 4     score 0.061
  a = 2 um, s = 3     score 0.072
The edge line (klempt2024_timescale.py) gave 0.040 (phi 0.029, c 0.051):
the corner block fits c three times better. Fig. 3's shape is not
reproduced, see klempt2024_fig3.py and KLEMPT2024_REPRODUCTION.md sec. 14.
"""
from multiprocessing import Pool
import numpy as np
import klempt2024_case1_bc as B, klempt2024_quantitative as K, klempt2024_variant_search as VS
def mask(a):
    return (K.X >= K.L - a - 1e-9) & (K.Y >= K.L - a - 1e-9) & (K.Z >= K.L - a - 1e-9)
def one(job):
    a, s = job
    K.R, K.BETA, K.K_A = 100*s, 2*s, 1e-3*s
    p,_,_ = K.setup("fig4_edge")
    rec,_,_ = B.run_setup(p, mask(a), 1e8, "ic", "first_order", "abs")
    sc, rms = VS.score(rec, VS.T41, VS.digitized()["4.1"])
    return a, s, sc, rms


def main():
    jobs = [(a, s) for a in (1, 2, 3, 4) for s in (1, 1.5, 2, 3, 4, 6)]
    with Pool(4) as p:
        for a, s, sc, rms in p.imap_unordered(one, jobs):
            print(f"a={a} s={s} score {sc:.3f} phi {rms['phi']:.3f} c {rms['c']:.3f}", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Reference for the four real Day-1 CLSM seeds (CS/CH/DS/DH) advanced by a
coarse mechanical increment, the way usermat_biofilm.f now handles one:
split into ceil(dTime/1e-4) ecology sub-steps, alpha += k_alpha * phi_int.

    python ecology_4region_reference_coarse.py [dTime] [n_increments]

Defaults to dTime=0.1, one increment -- the stepping of the partner deck
discussed in V222_PORT_INSTRUCTIONS.md. Seeds and k_alpha come from
ecology_4region_reference_all_real_clsm.py, so the two cannot drift.
"""
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "coupling"))

import ecology_jax  # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO  # noqa: E402
from ecology_4region_reference_all_real_clsm import (  # noqa: E402
    CONDITIONS, K_ALPHA, day1_seed,
)

DT_ECO_MAX = 1.0e-4   # must match usermat_biofilm.f


if __name__ == "__main__":
    dt = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
    n_inc = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    n_sub = max(1, int(np.ceil(dt / DT_ECO_MAX - 1e-9)))
    print(f"dTime={dt} x {n_inc} increment(s), n_sub={n_sub}, k_alpha={K_ALPHA}\n")
    for label, sheet, key in CONDITIONS:
        g, _ = day1_seed(sheet, key)
        alpha = 0.0
        for _ in range(n_inc):
            g, phi_int = ecology_jax.ecology_substeps(g, THETA_DEMO, dt, n_sub)
            alpha += K_ALPHA * phi_int
        g = np.asarray(g)
        print(f"--- {label} ---")
        print(f"alpha = {alpha:.10e}   stretch 1+alpha = {1 + alpha:.4f}")
        print("g(1:12): " + ", ".join(f"{v:.6e}" for v in g.tolist()))

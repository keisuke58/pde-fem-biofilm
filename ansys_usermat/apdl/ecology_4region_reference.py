#!/usr/bin/env python3
"""Independent reference for t_growth_cylinder_ecology_4region.dat's four
growth-layer materials (mat 2/3/4/5) -- extends
ecology_cylinder_reference.py's two-region check to a 2x2 (theta x Z) grid.
Seeds A/B are reused unchanged from the base ecology deck; C/D are new,
deliberately different-again compositions (not measured data -- same status
as A/B). Same chained ecology_step logic, dt=1e-5 x 10 substeps, THETA_DEMO.
"""
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "coupling"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import ecology_jax
from jax_hamilton_0d_5species_demo import THETA_DEMO

DT = 1.0e-5
N_SUBSTEPS = 10
K_ALPHA = 50.0

SEEDS = {
    "A (mat2, default zero-seed)": ecology_jax.default_initial_state(),
    "B (mat3, existing synthetic)": np.array([0.05, 0.05, 0.15, 0.25, 0.20, 0.30,
                                               0.999, 0.999, 0.999, 0.999, 0.999, 0.0]),
    "C (mat4, new: sparse uneven)": np.array([0.10, 0.02, 0.02, 0.02, 0.04, 0.10,
                                               0.99, 0.99, 0.99, 0.99, 0.99, 0.0]),
    "D (mat5, new: low uniform)":   np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.05,
                                               0.99, 0.99, 0.99, 0.99, 0.99, 0.0]),
}

if __name__ == "__main__":
    print(f"dt={DT}, N_SUBSTEPS={N_SUBSTEPS}, k_alpha={K_ALPHA}\n")
    for label, g0 in SEEDS.items():
        g = np.asarray(g0, dtype=np.float64)
        alpha = 0.0
        for _ in range(N_SUBSTEPS):
            g = np.asarray(ecology_jax.ecology_step(g, THETA_DEMO, DT))
            phi_tot = ecology_jax.living_fraction_total(g)
            alpha += DT * K_ALPHA * phi_tot
        ok = "OK" if alpha < 0.01 else "CHECK -- may not converge"
        print(f"{label}: alpha = {alpha:.10e}   (<0.01: {ok})")

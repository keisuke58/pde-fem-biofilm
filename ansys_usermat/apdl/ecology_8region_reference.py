import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, str(Path(r"c:\Users\nishioka\git\pde-fem-biofilm\ansys_usermat\coupling")))
sys.path.insert(0, str(Path(r"c:\Users\nishioka\git\pde-fem-biofilm")))

import ecology_jax
from jax_hamilton_0d_5species_demo import THETA_DEMO

DT = 1.0e-5
N_SUBSTEPS = 10
K_ALPHA = 50.0

SEEDS = {
    "A (mat2, outer, default zero-seed)": ecology_jax.default_initial_state(),
    "B (mat3, outer, existing synthetic)": np.array([0.05,0.05,0.15,0.25,0.20, 0.30,
                                                      0.999,0.999,0.999,0.999,0.999, 0.0]),
    "C (mat4, outer, sparse uneven)": np.array([0.10,0.02,0.02,0.02,0.04, 0.10,
                                                 0.99,0.99,0.99,0.99,0.99, 0.0]),
    "D (mat5, outer, low uniform)":   np.array([0.02,0.02,0.02,0.02,0.02, 0.05,
                                                 0.99,0.99,0.99,0.99,0.99, 0.0]),
    "E (mat6, inner, very sparse/early)": np.array([0.01,0.01,0.01,0.01,0.01, 0.02,
                                                     0.99,0.99,0.99,0.99,0.99, 0.0]),
    "F (mat7, inner, heavy base biofilm)": np.array([0.15,0.15,0.05,0.05,0.05, 0.40,
                                                      0.99,0.99,0.99,0.99,0.99, 0.0]),
    "G (mat8, inner, alternating)": np.array([0.03,0.08,0.03,0.08,0.03, 0.15,
                                               0.99,0.99,0.99,0.99,0.99, 0.0]),
    "H (mat9, inner, single dominant)": np.array([0.25,0.01,0.01,0.01,0.01, 0.05,
                                                   0.99,0.99,0.99,0.99,0.99, 0.0]),
}

for label, g0 in SEEDS.items():
    g = np.asarray(g0, dtype=np.float64)
    alpha = 0.0
    for _ in range(N_SUBSTEPS):
        g = np.asarray(ecology_jax.ecology_step(g, THETA_DEMO, DT))
        phi_tot = ecology_jax.living_fraction_total(g)
        alpha += DT * K_ALPHA * phi_tot
    print(f"{label}: alpha = {alpha:.10e}   (<0.01: {'OK' if alpha < 0.01 else 'CHECK'})")

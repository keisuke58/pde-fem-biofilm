#!/usr/bin/env python3
"""Independent reference for t_growth_cylinder_ecology_clsm.dat's two
growth-layer materials (mat 2 / mat 3) -- same two-region curved-shell
pattern as ecology_cylinder_reference.py, but region B's seed is REAL
Day-1 (Tag=1) Commensal/Static CLSM composition from
data/heine_species_distribution_biofilm.xlsx ("all cells" sheet), not a
synthetic stand-in. Same values clsm_phi_seed_reference.py already verified
on the single-element deck (t_growth_ecology_clsm_phi.dat) -- this script
chains the same real seed through 10 substeps (dt=1e-5 each) instead of 1,
matching t_growth_cylinder_ecology.dat's convention for comparing against
a real, unconstrained, curved geometry.

psi(1:5) is still the 0.999 placeholder here too -- see
t_growth_ecology_clsm_phi.dat's header comment for why the workbook's
"only living cells" ratio is not a valid substitute (values >1 for 2 of 5
species on this exact condition/day).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "coupling"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import ecology_jax  # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO  # noqa: E402
from clsm_phi_seed_reference import day1_phi  # noqa: E402

DT = 1.0e-5
N_SUBSTEPS = 10
K_ALPHA = 50.0


def real_clsm_seed():
    phi_raw = day1_phi()  # Commensal/Static, Day 1
    phi = phi_raw * (0.999999 / phi_raw.sum())
    phi0 = 1.0 - phi.sum()
    return np.concatenate([phi, [phi0], np.full(5, 0.999), [0.0]])


def run(label, g0):
    g = np.asarray(g0, dtype=np.float64)
    alpha = 0.0
    for _ in range(N_SUBSTEPS):
        g = np.asarray(ecology_jax.ecology_step(g, THETA_DEMO, DT))
        phi_tot = ecology_jax.living_fraction_total(g)
        alpha += DT * K_ALPHA * phi_tot
    print(f"--- {label} ---")
    print(f"alpha after {N_SUBSTEPS} steps: {alpha:.10e}")
    print("g(1:12):  " + ", ".join(f"{v:.10e}" for v in g.tolist()))
    return alpha, g


if __name__ == "__main__":
    print(f"dt={DT}, N_SUBSTEPS={N_SUBSTEPS}, k_alpha={K_ALPHA}\n")
    a_a, g_a = run(
        "region A (mat 2, all-zero seed -> INIT_ECO_IF_ZERO default)",
        ecology_jax.default_initial_state(),
    )
    a_b, g_b = run(
        "region B (mat 3, REAL Day-1 Commensal/Static CLSM seed)",
        real_clsm_seed(),
    )
    print(f"\nalpha ratio B/A: {a_b / a_a:.4f}")
    print(
        "both alphas well under the deck's known 0.01 convergence "
        f"threshold: {'OK' if max(a_a, a_b) < 0.01 else 'CHECK -- may not converge'}"
    )

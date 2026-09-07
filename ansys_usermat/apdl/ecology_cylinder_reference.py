#!/usr/bin/env python3
"""Independent reference for t_growth_cylinder_ecology.dat's ecology-driven
region A / region B growth-layer materials (mat 2 / mat 3).

Unlike t_growth_constrained.dat's closed form (F=I everywhere, so stress is
predictable), a two-layer curved-shell solve has no closed form -- but the
ecology *state trajectory* itself does, the same way the single/multi-element
ecology smoke tests verified ustatev(15:26) against ecology_jax.ecology_step.
This script chains that same call 10 times (one per committed ANSYS
substep, TIME=1e-4 / NSUBST,10 => dTime = 1e-5 each, matching the verified
single-step magnitude from t_growth_ecology.dat) for each material's initial
seed, and prints the resulting alpha and g(12) -- to be diffed against
SVAR(10) and SVAR(15:26) in the real ANSYS run's growth-layer elements,
split by material (mat 2 vs mat 3), the same role closed_form_reference.py
plays for the unit-cube decks.

Region B's seed is a synthetic stand-in for a spatially-varying CLSM
composition (coupling/README.md's still-open next step), NOT calibrated
data -- picked only to be visibly different from region A's all-zero seed
(which triggers INIT_ECO_IF_ZERO's default_initial_state()), so the two
regions provably grow at different rates.
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

REGION_B_SEED = np.array([
    0.05, 0.05, 0.15, 0.25, 0.20,   # phi(1:5)
    0.30,                            # phi0
    0.999, 0.999, 0.999, 0.999, 0.999,  # psi(1:5)
    0.0,                             # gamma
])


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
    a_a, g_a = run("region A (mat 2, all-zero seed -> INIT_ECO_IF_ZERO default)",
                    ecology_jax.default_initial_state())
    a_b, g_b = run("region B (mat 3, explicit synthetic seed)", REGION_B_SEED)
    print(f"\nalpha ratio B/A: {a_b / a_a:.4f}")
    print(f"both alphas well under the deck's known 0.01 convergence "
          f"threshold: {'OK' if max(a_a, a_b) < 0.01 else 'CHECK -- may not converge'}")

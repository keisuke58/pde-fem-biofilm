"""
clsm_phi_seed_reference.py -- independent Python reference for
t_growth_ecology_clsm_phi.dat's real-CLSM-seeded ustatev(15:26) TB,STATE
block (coupling/README.md next-steps #4's first open item, partial: real
phi, still-placeholder psi/theta -- see that deck's own header comment for
why psi is NOT derived from the workbook's "only living cells" sheet here).

Regenerates, from data/heine_species_distribution_biofilm.xlsx directly:
  1. the Day-1 (Tag=1) mean composition phi for Commensal/Static,
  2. the normalized TBDATA values the deck hardcodes,
  3. the one-step ecology_jax.ecology_step(..., THETA_DEMO, dt=1e-5) result
     and the resulting alpha/stress, matching what the deck's own header
     comment claims and what the real 2026-09-08 ANSYS run
     (growth_result_ecology_clsm_phi.txt) actually printed.

Run: python ansys_usermat/apdl/clsm_phi_seed_reference.py
Needs: openpyxl, jax (both already used elsewhere in this repo).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / "ansys_usermat" / "coupling"))

import openpyxl  # noqa: E402
import jax.numpy as jnp  # noqa: E402

from plot_heine_phi_psi import load  # noqa: E402
import ecology_jax as eco  # noqa: E402
import material_server as ms  # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO  # noqa: E402

XLSX = _REPO / "data" / "heine_species_distribution_biofilm.xlsx"
SPECIES = ["S. oralis", "A. naeslundii", "Veillonella", "F. nucleatum", "P. gingivalis"]
PSI_PLACEHOLDER = 0.999  # NOT from the workbook -- see module docstring


def day1_phi(condition="Static all cells", sheet="Commensal"):
    """Mean of Day-1 (Tag=1) per-species composition [%], normalized to a
    fraction, from the workbook's "all cells" sheet block."""
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    d = load(wb[sheet])
    rows = d[condition][1]  # Tag == 1
    phi_pct = np.nanmean(np.array(rows, dtype=float), axis=1)
    phi = phi_pct / 100.0
    return phi


def build_initial_state(phi_raw):
    """Same normalization ecology_jax.default_initial_state() uses for its
    own placeholder: scale phi to sum 0.999999, fill the remainder as phi0."""
    phi = phi_raw * (0.999999 / phi_raw.sum())
    phi0 = 1.0 - phi.sum()
    psi = np.full(5, PSI_PLACEHOLDER)
    g = jnp.zeros(12, dtype=jnp.float64)
    g = g.at[0:5].set(jnp.array(phi))
    g = g.at[5].set(phi0)
    g = g.at[6:11].set(jnp.array(psi))
    g = g.at[11].set(0.0)
    return g, phi, phi0, psi


def main():
    phi_raw = day1_phi()
    print("Day-1 Commensal/Static CLSM composition (raw fractions):")
    for name, v in zip(SPECIES, phi_raw):
        print(f"  {name:16s} {v:.6f}")
    print(f"  sum = {phi_raw.sum():.9f}\n")

    g0, phi, phi0, psi = build_initial_state(phi_raw)
    print("TBDATA,15 line (phi(1:5), phi0):")
    print("  " + ",".join(f"{v:.9G}" for v in list(phi) + [phi0]))
    print("TBDATA,21 line (psi(1:5), placeholder):")
    print("  " + ",".join(f"{v:.9G}" for v in psi))

    dt = 1e-5
    k_alpha = 50.0
    g1 = eco.ecology_step(g0, THETA_DEMO, dt)
    # alpha uses phi_tot of the *post-step* state -- matches
    # t_growth_ecology.dat's established convention, not the initial state.
    phi_tot_new = eco.living_fraction_total(g1)
    alpha_new = dt * k_alpha * phi_tot_new

    F = np.eye(3)
    sv, _, _ = ms.stress_core(F, F, alpha_new, 1.0, 0.0, 0.01, 0.0, 0.0, 1.0)

    print("\nOne ecology_step (THETA_DEMO, dt=1e-5) result:")
    print("  g_new(phi)  :", np.array(g1[0:5]))
    print("  g_new(phi0) :", float(g1[5]))
    print("  g_new(psi)  :", np.array(g1[6:11]))
    print("  g_new(gamma):", float(g1[11]))
    print(f"  phi_tot(g_new) = {phi_tot_new:.9f}")
    print(f"  alpha_new      = {alpha_new:.9e}")
    print(f"  SX=SY=SZ       = {sv[0]:.9e}")


if __name__ == "__main__":
    main()

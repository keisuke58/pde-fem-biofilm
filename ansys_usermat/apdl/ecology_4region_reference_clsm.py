#!/usr/bin/env python3
"""Independent reference for t_growth_cylinder_ecology_4region_clsm.dat's
four growth-layer materials (mat 2/3/4/5) -- extends
ecology_cylinder_reference_clsm.py's single-real-condition check to THREE
clinical conditions in data/heine_species_distribution_biofilm.xlsx
(region A / mat 2 stays the pre-existing all-zero -> INIT_ECO_IF_ZERO
default, kept unchanged as a regression baseline).

CORRECTION, 2026-09-08 (same day, later): an earlier version of this file
also computed DH_SEED (Dysbiotic/HOBIC) below and concluded it had to be
left out because F. nucleatum/P. gingivalis had "no measurement at all" at
Day-1. That conclusion was WRONG -- it was an artifact of a real bug in
plot_heine_phi_psi.load(), now fixed: that function used to derive the
species-block column width ONCE from the sheet's first header row and
reuse it everywhere, but the Dysbiotic sheet's "HOBIC ..." blocks are 9
columns/species wide while its "Static ..." blocks are 18 -- reusing
width=18 for HOBIC didn't just drop 2 species, it read species 2's real
data into the "species 1" slot, species 3's into "species 2", species 5's
into "species 3", and landed the last two slots on blank padding past the
real data. See t_growth_cylinder_ecology_4region_all_real_clsm.dat for the
corrected DH seed and all four real conditions on one deck.
"""
import sys
from pathlib import Path

import numpy as np
import openpyxl

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / "ansys_usermat" / "coupling"))

from plot_heine_phi_psi import load  # noqa: E402
import ecology_jax  # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO  # noqa: E402

XLSX = _REPO / "data" / "heine_species_distribution_biofilm.xlsx"
DT = 1.0e-5
N_SUBSTEPS = 10
K_ALPHA = 50.0

CONDITIONS = [
    ("B (mat3) CS - Commensal/Static", "Commensal", "Static all cells"),
    ("C (mat4) CH - Commensal/HOBIC", "Commensal", "HOBIC all cells"),
    ("D (mat5) DS - Dysbiotic/Static", "Dysbiotic", "Static all cells"),
]


def day1_seed(sheet_name, condition_key):
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    d = load(wb[sheet_name])
    rows = d[condition_key][1]  # Tag == 1
    phi_pct = np.nanmean(np.array(rows, dtype=float), axis=1)
    phi = phi_pct / 100.0
    phi = phi * (0.999999 / phi.sum())
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
    run("A (mat2, all-zero seed -> INIT_ECO_IF_ZERO default)",
        ecology_jax.default_initial_state())
    for label, sheet, key in CONDITIONS:
        run(label, day1_seed(sheet, key))

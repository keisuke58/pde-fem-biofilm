#!/usr/bin/env python3
"""Independent reference for
t_growth_cylinder_ecology_4region_all_real_clsm.dat's four growth-layer
materials (mat 2/3/4/5) -- ALL FOUR of this repo's clinical conditions
(RESEARCH_MODEL.md sec.0: CS/CH/DS/DH), each seeded with real Day-1
(Tag=1) CLSM composition from data/heine_species_distribution_biofilm.xlsx.
Retires the earlier "region A = default placeholder" convention
(ecology_cylinder_reference_clsm.py, ecology_4region_reference_clsm.py)
now that DH's Day-1 data is available too -- see this directory's
ecology_4region_reference_clsm.py for the bug (in
plot_heine_phi_psi.load(), now fixed) that made DH look unmeasured when it
first was not.
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
    ("A (mat2) CS - Commensal/Static", "Commensal", "Static all cells"),
    ("B (mat3) CH - Commensal/HOBIC", "Commensal", "HOBIC all cells"),
    ("C (mat4) DS - Dysbiotic/Static", "Dysbiotic", "Static all cells"),
    ("D (mat5) DH - Dysbiotic/HOBIC", "Dysbiotic", "HOBIC all cells"),
]


def day1_seed(sheet_name, condition_key):
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    d = load(wb[sheet_name])
    rows = d[condition_key][1]  # Tag == 1
    phi_pct = np.array([np.nanmean(np.array(r, dtype=float)) for r in rows])
    phi = phi_pct / 100.0
    phi = phi * (0.999999 / phi.sum())
    phi0 = 1.0 - phi.sum()
    return np.concatenate([phi, [phi0], np.full(5, 0.999), [0.0]]), phi_pct


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
    for label, sheet, key in CONDITIONS:
        g0, phi_pct = day1_seed(sheet, key)
        print(f"Day-1 composition (%): {np.round(phi_pct, 3).tolist()}  "
              f"sum={phi_pct.sum():.3f}")
        run(label, g0)

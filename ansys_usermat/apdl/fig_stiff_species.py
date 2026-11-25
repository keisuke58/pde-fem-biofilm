"""fig_stiff_species.py -- stiffness by species (--stiff-species mA,mB; Keio continuation, not the thesis, kept back):
E = E(phi_3D) (w_A m_A + w_B m_B), w_i = phi_i / sum_j phi_j of the point model (dead biomass included). Sensitivity only,
there are no measured moduli per species. 64^3, Fritsch et al. 2025 case I, antibiotic alpha*_0 = 100 from T* = 1/3,
outer step 0.0025, with dc/dt, minmod. Seed von Mises at T* = 1 for the separate and the mixed start; alpha - 1 does not
change (the stiffness enters only the stress solve). Table: papers/coupling_paper/EXPERIMENT_COMPARISON.ja.md section 4.

    python ansys_usermat/apdl/fig_stiff_species.py -> assets/fig_stiff_species.pdf/.png
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402

figstyle.apply()
R = HERE / "results" / "2026-10-3d"
P0 = "coupled_living_consumption_nf64_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv"
SETS = [("1-1", ""), ("2-1", "_stiff2-1"), ("1-2", "_stiff1-2"), ("0.5-1", "_stiff0.5-1"), ("1-0.5", "_stiff1-0.5")]
START = {"separate": ("", "separate_ab"), "mixed": ("_mixeq", "mixed_ab")}

vals = {}
for start, (mx, run) in START.items():
    for lab, tag in SETS:
        f = R / f"{P0}{mx}_a100_gab1_gc1_tab0.333333{tag}_transc1e+08_tvd_wgmirror_fast_{run}.json"
        r = json.load(open(f))["runs"][run]
        vals[(start, lab)] = r["seed_vm_Pa"]
        print(f"{start:8s} m={lab:6s} vM {r['seed_vm_Pa']*1e4:.2f}e-4 Pa  alpha-1 {r['seed_am1']*1e4:.2f}e-4  "
              f"ratio {r['seed_vm_Pa']/vals[(start, '1-1')]:.2f}")

fig, ax = plt.subplots(figsize=(4.2, 2.3), constrained_layout=True)
x = np.arange(len(SETS)); w = 0.38
for k, (start, col) in enumerate((("separate", "#1c7ed6"), ("mixed", "#e8590c"))):
    v = np.array([vals[(start, lab)] for lab, _ in SETS]) * 1e4
    b = ax.bar(x + (k - 0.5) * w, v, w, color=col, label=f"{start} start")
    for xi, vi, ri in zip(x + (k - 0.5) * w, v, v / v[0]):
        ax.text(xi, vi + 0.15, f"{ri:.2f}", ha="center", va="bottom", fontsize=6)
ax.set_xticks(x, [rf"$m_A, m_B$ = {lab.replace('-', ', ')}" for lab, _ in SETS], fontsize=6.5)
ax.set_ylabel(r"seed von Mises at $T^*=1$ [$10^{-4}$ Pa]", fontsize=7.5)
ax.tick_params(labelsize=7); ax.set_ylim(0, 13.5)
ax.legend(frameon=False, fontsize=7, loc="upper right")
A = ROOT / "assets"
fig.savefig(A / "fig_stiff_species.pdf"); fig.savefig(A / "fig_stiff_species.png", dpi=400)
print("wrote", A / "fig_stiff_species.png")

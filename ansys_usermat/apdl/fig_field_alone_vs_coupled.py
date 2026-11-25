"""fig_field_alone_vs_coupled.py -- the coupled model against the field model of Klempt 2024 alone (Oliver, 25 Nov 2026),
both without antibiotic, separate and mixed start, 64^3, outer step 0.025 (the runs of Table 1 / Fig. 7). Field model
alone = --living-one: the uptake g phi_3D c with sum_i phibar_i = 1 (Klempt 2024 Eq. 35, misprint corrected); the
coupled model uses g phi_3D (sum_i phibar_i) c. (a) Seed nutrient c, (b) volume reached by the biofilm (cells with phi_3D >= phi_min). Prints the seed
values at T* = 1 (c, alpha - 1, von Mises) and the volume.

    python ansys_usermat/apdl/fig_field_alone_vs_coupled.py -> assets/fig_field_alone_vs_coupled.pdf/.png
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
P = "coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_{mx}a100_gab1_gc1_tab0.333333{l1}_tvd_wgmirror_fast_{k}.json"
EX = [("separate", "", "separate", "#3b5bdb"), ("mixed", "mixeq_", "mixed", "#7048a8")]
fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.3), constrained_layout=True)
for ex, mx, k, col in EX:
    cp = json.load(open(R / P.format(mx=mx, l1="", k=k)))["runs"][k]
    fa = json.load(open(R / P.format(mx=mx, l1="_l1", k=k)))["runs"][k]
    for r, ls, lab in ((cp, "-", "coupled"), (fa, "--", "field model alone")):
        axs[0].plot(r["t"], r["c_seed"], color=col, ls=ls, lw=1.2, label=f"{ex}, {lab}")
        axs[1].plot(r["t"], np.array(r["colonised"]) * (2.0 / 64) ** 3, color=col, ls=ls, lw=1.2)
    axs[2].plot(cp["t"], cp["living_seed"], color=col, lw=1.2)
    d = lambda a, b: 100 * (b - a) / a   # noqa: E731
    print(f"{ex:8s} T*=1 seed, coupled / field alone: c {cp['c_seed'][-1]:.3f} / {fa['c_seed'][-1]:.3f} ({d(cp['c_seed'][-1], fa['c_seed'][-1]):+.1f} %), "
          f"alpha-1 {cp['seed_am1']*1e4:.3f} / {fa['seed_am1']*1e4:.3f} ({d(cp['seed_am1'], fa['seed_am1']):+.1f} %), "
          f"vM {cp['seed_vm_Pa']*1e4:.3f} / {fa['seed_vm_Pa']*1e4:.3f} ({d(cp['seed_vm_Pa'], fa['seed_vm_Pa']):+.1f} %), "
          f"reached volume {cp['colonised'][-1]*(2/64)**3:.3f} / {fa['colonised'][-1]*(2/64)**3:.3f} mm^3, at T*=0.2 {cp['colonised'][8]*(2/64)**3:.2f} / {fa['colonised'][8]*(2/64)**3:.2f}; c at T*=0.1: {cp['c_seed'][4]:.3f} / {fa['c_seed'][4]:.3f}")
for ax, (yl, ylim, t) in zip(axs, (("seed $c$", (0, 1.02), "(a) seed nutrient"), (r"volume with $\phi_{3D}\geq 0.001$ [mm$^3$]", (0, 8.2), "(b) volume reached"),
                                    (r"seed $\Sigma_i\bar\phi_i$", (0, 1.02), "(c) living fraction (coupled)"))):
    ax.set_xlim(0, 1); ax.set_ylim(*ylim); ax.set_xlabel(r"$T^*$", fontsize=8); ax.set_ylabel(yl, fontsize=8)
    ax.tick_params(labelsize=7); ax.set_title(t, fontsize=8, loc="left")
axs[0].legend(frameon=False, fontsize=5.8, loc="lower left")
A = ROOT / "assets"
fig.savefig(A / "fig_field_alone_vs_coupled.pdf"); fig.savefig(A / "fig_field_alone_vs_coupled.png", dpi=400)
print("wrote", A / "fig_field_alone_vs_coupled.png")

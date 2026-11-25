"""fig_species_volumes.py -- volume of each species in the whole cube over time (Oliver, 25 Nov 2026), from the 64^3
snapshots of the two coupled examples with antibiotic (the runs of the videos). V_i = int phi_3D phi_i/(phi_1+phi_2) dV
(all biomass of species i, dead included) and the living volume int phi_3D sum_i phibar_i dV, in mm^3 of the 8 mm^3 cube.
The runs solve one quarter (symmetry in x and z); the snapshots are stored for the full cube.

    python ansys_usermat/apdl/fig_species_volumes.py -> assets/fig_species_volumes.pdf/.png
"""
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
F = {"separate": "coupled_living_separate_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast_separate_ab.npz",
     "mixed": "coupled_living_mixed_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast_mixed_ab.npz"}
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.4), constrained_layout=True, sharey=True)
for ax, (ex, fn) in zip(axs, F.items()):
    z = np.load(ROOT / "figdata" / fn)
    ts = sorted({float(k[4:]) for k in z.files if k.startswith("phi_")})
    key = lambda q, t: next(k for k in z.files if k.startswith(q + "_") and abs(float(k[len(q) + 1:]) - t) < 1e-9)   # noqa: E731
    n = z[key("phi", ts[0])].shape[0]
    dv = (2.0 / n) ** 3                       # the snapshots are stored unfolded (full cube)
    VA, VB, VL = [], [], []
    for t in ts:
        phi = np.clip(z[key("phi", t)], 0, 1); s0 = z[key("shareA", t)]; on = np.isfinite(s0); sh = np.where(on, s0, 0.0)
        L = np.nan_to_num(z[key("L", t)])
        VA.append(float((phi * sh).sum() * dv)); VB.append(float((phi * (1 - sh) * on).sum() * dv)); VL.append(float((phi * L).sum() * dv))
    ax.plot(ts, VA, color="#3b5bdb", lw=1.3, label="species A (all)")
    ax.plot(ts, VB, color="#c92a2a", lw=1.3, label="species B (all)")
    ax.plot(ts, VL, color="#2b8a3e", lw=1.3, ls="--", label=r"living, $\int\phi_{3D}\Sigma_i\bar\phi_i\,dV$")
    ax.axvline(1 / 3, color="0.5", ls=(0, (1, 2)), lw=0.8)
    ax.set_title(f"({'ab'[list(F).index(ex)]}) {ex} start, with antibiotic", fontsize=8, loc="left")
    ax.set_xlim(0, 1); ax.set_xlabel(r"$T^*$", fontsize=8); ax.tick_params(labelsize=7)
    i3 = int(np.argmin(np.abs(np.array(ts) - 1 / 3)))
    print(f"{ex}: grid {n}^3 in file; V_A, V_B, V_living at T*={ts[0]:.3f}: {VA[0]:.3f} {VB[0]:.3f} {VL[0]:.3f}; "
          f"T*={ts[i3]:.3f}: {VA[i3]:.3f} {VB[i3]:.3f} {VL[i3]:.3f}; T*=1: {VA[-1]:.3f} {VB[-1]:.3f} {VL[-1]:.3f} mm^3")
axs[0].set_ylabel(r"volume [mm$^3$] (cube 8 mm$^3$)", fontsize=8)
axs[0].legend(frameon=False, fontsize=6.5, loc="upper left")
A = ROOT / "assets"
fig.savefig(A / "fig_species_volumes.pdf"); fig.savefig(A / "fig_species_volumes.png", dpi=400)
print("wrote", A / "fig_species_volumes.png")

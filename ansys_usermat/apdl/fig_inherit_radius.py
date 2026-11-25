"""fig_inherit_radius.py -- start of newly reached cells from all started cells within a radius R (Oliver's filter, 25 Nov
2026, --inherit-radius R) against the six face neighbours used so far. 64^3, outer step 0.025, with antibiotic, T* = 1.
Separate start: slice z = 0 (x-y plane through the seed and two of the B blocks), (a-c) share phi_1/(phi_1+phi_2) and
(d-f) living volume fraction where phi_3D > 0.5, for six neighbours, R = 0.1 mm and R = 0.2 mm. Prints the seed values,
the fraction of biofilm cells with a share between 0.3 and 0.7 (the mixed border) and, for the mixed start, the number of
living cells (sum phibar > 0.3) at T* = 1.

    python ansys_usermat/apdl/fig_inherit_radius.py -> assets/fig_inherit_radius.pdf/.png
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
TAG = "nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_{mx}a100_gab1_gc1_tab0.333333{r}_tvd_wgmirror_fast_{k}"
CASES = [("six neighbours", ""), ("$R = 0.1$ mm", "_inhr0.1"), ("$R = 0.2$ mm", "_inhr0.2")]


def load(k, mx, r):
    t = TAG.format(mx=mx, r=r, k=k)
    z = np.load(ROOT / "figdata" / f"coupled_living_{k}_{t}.npz")
    js = json.load(open(HERE / "results" / "2026-10-3d" / f"coupled_living_consumption_{t}.json"))["runs"][k]
    T = max((q[4:] for q in z.files if q.startswith("phi_")), key=float)
    return z, js, T


fig, axs = plt.subplots(2, 3, figsize=(6.3, 4.3), constrained_layout=True)
for j, (lab, r) in enumerate(CASES):
    z, js, T = load("separate_ab", "", r)
    phi, sh, L = z[f"phi_{T}"], z[f"shareA_{T}"], np.nan_to_num(z[f"L_{T}"])
    n = phi.shape[0]
    m = phi > 0.5
    mid = m & np.isfinite(sh) & (sh > 0.3) & (sh < 0.7)
    print(f"separate_ab {lab:16s} seed living {js['living_seed'][-1]:.3f} share {js['share_a_seed'][-1]:.3f} alpha-1 {js['seed_am1']*1e4:.3f} "
          f"vM {js['seed_vm_Pa']*1e4:.3f}; mixed border {100*mid.sum()/m.sum():.1f} % of biofilm cells")
    kz = n // 2                                 # cell just above z = 0
    for i, (fld, cmap, ttl) in enumerate(((sh, "coolwarm_r", r"$\phi_1/(\phi_1+\phi_2)$"), (L, "Greens", r"$\Sigma_i\bar\phi_i$"))):
        img = np.where(m[:, :, kz], fld[:, :, kz], np.nan).T
        ax = axs[i, j]
        im = ax.imshow(img, origin="lower", extent=(-1, 1, -1, 1), cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
        ax.set_title(f"({'abcdef'[3 * i + j]}) {lab}, {ttl}", fontsize=7.5, loc="left")
        ax.set_xlabel("$x$ [mm]", fontsize=7); ax.set_ylabel("$y$ [mm]", fontsize=7); ax.tick_params(labelsize=6.5)
        if j == 2:
            fig.colorbar(im, ax=axs[i, :], shrink=0.8).ax.tick_params(labelsize=6.5)
for lab, r in (CASES[0], CASES[2]):
    z, js, T = load("mixed_ab", "mixeq_", r)
    L, phi = np.nan_to_num(z[f"L_{T}"]), z[f"phi_{T}"]
    print(f"mixed_ab    {lab:16s} seed living {js['living_seed'][-1]:.3f} share {js['share_a_seed'][-1]:.3f} alpha-1 {js['seed_am1']*1e4:.3f} "
          f"vM {js['seed_vm_Pa']*1e4:.3f}; living cells (>0.3) {int(((L > 0.3) & (phi > 0.5)).sum())}")
A = ROOT / "assets"
fig.savefig(A / "fig_inherit_radius.pdf"); fig.savefig(A / "fig_inherit_radius.png", dpi=400)
print("wrote", A / "fig_inherit_radius.png")

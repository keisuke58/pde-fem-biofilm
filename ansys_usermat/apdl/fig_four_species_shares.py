"""fig_four_species_shares.py -- the four species of Klempt et al. 2026, case 4, in the coupled cube example with antibiotic
(64^3, the --snap-every --snap-stress run with the share of every species saved). Top: the share phi_i / sum_j phi_j of
each species where phi_3D > 0.5, at T* = 0.325 (just before the antibiotic) and T* = 1 (one octant removed, the camera
of all cube figures). Bottom: the share of each species averaged over the biofilm against T*. The antibiotic
sensitivities are b = (10, 2, 1, 0.01), so species 4 survives the antibiotic and takes over.

    python ansys_usermat/apdl/fig_four_species_shares.py [separate_ab|mixed_ab]   -> assets/fig_four_species_shares_<case>.png
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pyvista as pv

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply()
CASE = sys.argv[1] if len(sys.argv) > 1 else "separate_ab"
MIX = "_mixeq" if CASE.startswith("mixed") else ""
Z = np.load(ROOT / "figdata" / (f"coupled_living_{CASE}_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv{MIX}_a100_gab1_gc1"
                                 f"_tab0.333333_snapsvm_transc1e+08_k26c4_tvd_wgmirror_fast_{CASE}.npz"))
times = sorted({float(k.split("_", 1)[1]) for k in Z.files if k.startswith("sp0_")})
key = lambda n, t: f"{n}_{t}" if f"{n}_{t}" in Z.files else f"{n}_{t:g}"
T_SHOW = [min(times, key=lambda t: abs(t - 0.325)), times[-1]]
CAM = [(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)]
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)
CMAPS = ["Blues", "PuBu", "Reds", "Oranges"]
COL = ["#3b5bdb", "#74c0fc", "#c92a2a", "#f59f00"]
B = (10, 2, 1, 0.01)


def mesh(T, i):
    n = Z[key("phi", T)].shape[0]
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1.0,) * 3)
    cc = g.cell_centers().points
    octant = (cc[:, 0] > 0) & (cc[:, 1] > 0) & (cc[:, 2] > 0)
    g.cell_data["s"] = np.nan_to_num(Z[key(f"sp{i}", T)]).ravel(order="F")
    g.cell_data["bio"] = ((Z[key("phi", T)].ravel(order="F") > 0.5) & ~octant).astype(float)
    return g.threshold(0.5, scalars="bio")


def render(m, cmap):
    pl = pv.Plotter(off_screen=True, window_size=(900, 860))
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
    pl.add_mesh(m, scalars="s", cmap=cmap, clim=(0, 1), show_scalar_bar=False, **LIT)
    pl.camera_position = CAM
    pl.camera.zoom(0.95)
    img = pl.screenshot(return_img=True)
    pl.close()
    r = np.where((img < 250).any(axis=(1, 2)))[0]; c = np.where((img < 250).any(axis=(0, 2)))[0]
    return img[max(r[0] - 8, 0):r[-1] + 8, max(c[0] - 8, 0):c[-1] + 8]


fig = plt.figure(figsize=(6.3, 5.6))
W = 0.245
for j, T in enumerate(T_SHOW):
    for i in range(4):
        ax = fig.add_axes([0.02 + W * i, 0.70 - 0.29 * j, W - 0.01, 0.27])
        ax.imshow(render(mesh(T, i), CMAPS[i])); ax.set_axis_off()
        if j == 0:
            ax.set_title(rf"species {i + 1} ($b_{i + 1} = {B[i]:g}$)", fontsize=8.5, color=COL[i])
        if i == 0:
            ax.text(-0.03, 0.5, rf"$T^* = {T:.3f}$" if j == 0 else rf"$T^* = {T:.0f}$", transform=ax.transAxes,
                    rotation=90, va="center", ha="right", fontsize=8.5)
cb = fig.add_axes([0.3, 0.385, 0.4, 0.015])
c = fig.colorbar(plt.cm.ScalarMappable(cmap="Greys", norm=plt.Normalize(0, 1)), cax=cb, orientation="horizontal")
c.set_label(r"share $\phi_i/\Sigma_j\phi_j$ (each panel in the colour of its species)", fontsize=8); c.ax.tick_params(labelsize=7.5)
ax = fig.add_axes([0.1, 0.07, 0.62, 0.24])
for i in range(4):
    m = [float(np.nanmean(np.where(Z[key("phi", t)] > 0.5, Z[key(f"sp{i}", t)], np.nan))) for t in times]
    ax.plot(times, m, color=COL[i], lw=1.6, label=f"species {i + 1}")
ax.axvline(1 / 3, color="0.5", ls=":", lw=0.8); ax.text(0.345, 0.88, "antibiotic", fontsize=7.5, color="0.4")
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("$T^*$", fontsize=9); ax.set_ylabel("mean share in the biofilm", fontsize=8.5)
ax.tick_params(labelsize=8); ax.legend(frameon=False, fontsize=7.5, loc="center left", bbox_to_anchor=(1.01, 0.5))
d = -np.asarray(CAM[0], float); d /= np.linalg.norm(d); r_ = np.cross(d, CAM[2]); r_ /= np.linalg.norm(r_); u = np.cross(r_, d)
tx = fig.add_axes([0.86, 0.36, 0.12, 0.1])
for v, lab, col in (((1, 0, 0), "x", "#d9480f"), ((0, 1, 0), "y", "#2b8a3e"), ((0, 0, 1), "z", "#1c3faa")):
    e = np.array(v, float); px, py = e @ r_, e @ u
    tx.annotate("", xy=(px, py), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0, mutation_scale=7))
    tx.text(1.3 * px, 1.3 * py, f"${lab}$", color=col, fontsize=8, ha="center", va="center")
tx.set_xlim(-1.05, 1.05); tx.set_ylim(-0.6, 1.25); tx.set_aspect("equal"); tx.set_axis_off()
out = ROOT / "assets" / f"fig_four_species_shares_{CASE}.png"
fig.savefig(out, dpi=400)
print("wrote", out)

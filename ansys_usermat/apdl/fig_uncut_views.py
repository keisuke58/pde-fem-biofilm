"""fig_uncut_views.py -- the separate example with antibiotic at T* = 1 (64^3, the run of fig_coupled_3d_panels.py), drawn
without the removed octant from two opposite cameras, for the question of 2 Nov 2026 (Assoc. Prof. Muramatsu: the cut
3D views looked unsymmetric). Camera (3.6, 2.6, 4.4) as all 3D figures, and the camera mirrored in x and z,
(-3.6, 2.6, -4.4). The fields are mirror symmetric about x = 0 and z = 0 (fig_symmetry_sections.py), so the two views
are mirror images of each other. Rows: share of species A (where phi_3D > 0.5), von Mises stress (16^3).

    python ansys_usermat/apdl/fig_uncut_views.py   -> assets/fig_uncut_views.png
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
Z = np.load(ROOT / "figdata" / ("coupled_living_separate_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100"
                                 "_gab1_gc1_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast_separate_ab.npz"))
T = "1.0"
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)
CAMS = [[(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)], [(-3.6, 2.6, -4.4), (0, 0, 0), (0, 1, 0)]]
VMAX = float(2 * np.ceil(1e4 * np.percentile(Z[f"vm16_{T}"], 99.5) / 2))


def grid(m, name, values):
    g = pv.ImageData(dimensions=(m + 1,) * 3, spacing=(2.0 / m,) * 3, origin=(-1.0,) * 3)
    g.cell_data[name] = values.ravel(order="F")
    return g


phi = Z[f"phi_{T}"]
gA = grid(64, "shareA", np.where(phi > 0.5, np.nan_to_num(Z[f"shareA_{T}"]), -1.0)).threshold(-0.5, scalars="shareA")
gS = grid(Z[f"vm16_{T}"].shape[0], "vm", 1e4 * np.nan_to_num(Z[f"vm16_{T}"]))
ROWS = [(gA, "shareA", "coolwarm_r", (0, 1), r"$\phi_1/(\phi_1+\phi_2)$"),
        (gS, "vm", "magma_r", (0, VMAX), r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")]


def render(mesh, field, cmap, clim, cam):
    pl = pv.Plotter(off_screen=True, window_size=(1400, 1340))
    pl.enable_anti_aliasing("ssaa")
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey", line_width=3)
    pl.add_mesh(mesh, scalars=field, cmap=cmap, clim=clim, show_scalar_bar=False, **LIT)
    pl.camera_position = cam
    pl.camera.zoom(0.9)
    img = pl.screenshot(return_img=True)
    pl.close()
    rows = np.where((img < 250).any(axis=(1, 2)))[0]
    cols = np.where((img < 250).any(axis=(0, 2)))[0]
    return img[max(rows[0] - 12, 0):rows[-1] + 12, max(cols[0] - 12, 0):cols[-1] + 12]


def triad(ax, cam):
    d = -np.asarray(cam[0], float); d /= np.linalg.norm(d)
    r = np.cross(d, cam[2]); r /= np.linalg.norm(r); u = np.cross(r, d)
    for v, lab, col in (((1, 0, 0), "x", "#d9480f"), ((0, 1, 0), "y", "#2b8a3e"), ((0, 0, 1), "z", "#1c3faa")):
        e = np.array(v, float); px, py = e @ r, e @ u
        ax.annotate("", xy=(px, py), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0, mutation_scale=7))
        ax.text(1.3 * px, 1.3 * py, f"${lab}$", color=col, fontsize=8, ha="center", va="center")
    ax.set_xlim(-1.05, 1.05); ax.set_ylim(-0.7, 1.25); ax.set_aspect("equal"); ax.set_axis_off()


fig = plt.figure(figsize=(5.0, 5.0))
for j, cam in enumerate(CAMS):
    fig.text(0.27 + 0.42 * j, 0.975, ["(a) camera (3.6, 2.6, 4.4)", "(b) camera (-3.6, 2.6, -4.4)"][j],
             ha="center", va="top", fontsize=8.5)
    for i, (mesh, field, cmap, clim, lab) in enumerate(ROWS):
        ax = fig.add_axes([0.06 + 0.42 * j, 0.55 - 0.42 * i, 0.4, 0.4])
        ax.imshow(render(mesh, field, cmap, clim, cam), interpolation="none"); ax.set_axis_off()
        if j == 1:
            cb = fig.add_axes([0.9, 0.6 - 0.42 * i, 0.02, 0.3])
            c = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*clim)), cax=cb)
            c.set_label(lab, fontsize=8); c.ax.tick_params(labelsize=7)
    tx = fig.add_axes([0.06 + 0.42 * j + 0.15, 0.0, 0.1, 0.1]); triad(tx, cam)
o = ROOT / "assets" / "fig_uncut_views.png"
fig.savefig(o, dpi=400)
print("wrote", o)

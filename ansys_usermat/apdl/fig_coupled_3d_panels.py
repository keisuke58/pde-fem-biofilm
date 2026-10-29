"""fig_coupled_3d_panels.py -- figure for Chapter 4: the coupled model in 3D (Python, 64^3, TVD, the snapshots of the
--snap-stress run of the separate example with antibiotic) at T* = 1/3 (antibiotic onset) and T* = 1. Columns:
(a) share of group A, (b) living volume fraction, both where phi_3D > 0.5; (c) nutrient c; (d) antibiotic a;
(e) von Mises stress (16^3). One octant removed. Printed width 16 cm, Times (figstyle), renderings at 400 dpi; no title
(the caption carries it).

    python ansys_usermat/apdl/fig_coupled_3d_panels.py [separate_ab|mixed_ab] [inhliv|startpos]
        -> assets/fig_coupled_3d_panels_<case>.pdf / .png
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
# rule for newly reached points: startpos (by position, the A-75 runs) or inhliv (the neighbours' living composition,
# the base case from 25 Oct, A-82)
RULE = sys.argv[2] if len(sys.argv) > 2 else "inhliv"
MIX = "_mixeq" if CASE.startswith("mixed") else ""
Z = np.load(ROOT / "figdata" / (f"coupled_living_{CASE}_nf64_front1_phimin0.001_v2_fortran_ansysbc_{RULE}{MIX}_a100"
                                 f"_gab1_gc1_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast_{CASE}.npz"))
times = sorted({float(k.split("_", 1)[1]) for k in Z.files if k.startswith("phi_")})
T_SHOW = [min(times, key=lambda t: abs(t - 1 / 3)), times[-1]]
n = Z[f"phi_{times[0]}"].shape[0]
nm = Z[f"vm16_{times[0]}"].shape[0]
VMAX = float(2 * np.ceil(1e4 * max(np.percentile(Z[f"vm16_{T}"], 99.5) for T in T_SHOW) / 2))
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)
CAM = [(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)]
PANELS = [("shareA", "bio", "coolwarm_r", (0, 1), r"$\phi_1/(\phi_1+\phi_2)$"),
          ("L", "living", "Greens", (0, 1), r"$\Sigma_i\,\bar\phi_i$"),
          ("c", "cube", "Blues", (0.4, 1), r"$c$"),
          ("a", "ab", "Reds", (0, 1), r"$a$"),
          ("vm", "stress", "magma_r", (0, VMAX), r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")]


def cube(m):
    g = pv.ImageData(dimensions=(m + 1,) * 3, spacing=(2.0 / m,) * 3, origin=(-1.0,) * 3)
    cc = g.cell_centers().points
    return g, (cc[:, 0] > 0) & (cc[:, 1] > 0) & (cc[:, 2] > 0)


def meshes(T):
    g, octant = cube(n)
    for k in ("phi", "shareA", "L", "c", "a"):
        g.cell_data[k] = np.nan_to_num(Z[f"{k}_{T}"], nan=0.0).ravel(order="F")
    g.cell_data["bio"] = ((g.cell_data["phi"] > 0.5) & ~octant).astype(float)
    g.cell_data["cube"] = (~octant).astype(float)
    g.cell_data["living"] = g.cell_data["bio"] * (g.cell_data["L"] > 0.01)
    g.cell_data["ab"] = g.cell_data["cube"] * (g.cell_data["a"] > 0.01)
    s, octs = cube(nm)
    s.cell_data["vm"] = 1e4 * Z[f"vm16_{T}"].ravel(order="F")
    s.cell_data["stress"] = (~octs).astype(float)
    out = {m: g.threshold(0.5, scalars=m) for m in ("bio", "living", "cube", "ab")}
    out["stress"] = s.threshold(0.5, scalars="stress")
    return out


def render(mesh, field, cmap, clim):
    pl = pv.Plotter(off_screen=True, window_size=(1600, 1520))
    pl.enable_anti_aliasing("ssaa")
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey", line_width=3)
    if mesh.n_cells:
        pl.add_mesh(mesh, scalars=field, cmap=cmap, clim=clim, show_scalar_bar=False, **LIT)
    pl.camera_position = CAM
    pl.camera.zoom(0.9)   # the whole cube in view, the front corner included
    img = pl.screenshot(return_img=True)
    pl.close()
    rows = np.where((img < 250).any(axis=(1, 2)))[0]
    cols = np.where((img < 250).any(axis=(0, 2)))[0]
    p = 12
    return img[max(rows[0] - p, 0):rows[-1] + p, max(cols[0] - p, 0):cols[-1] + p]


fig = plt.figure(figsize=(6.3, 3.05))
W, H0 = 0.198, 0.375
for i, T in enumerate(T_SHOW):
    M = meshes(T)
    for k, (field, mesh, cmap, clim, label) in enumerate(PANELS):
        ax = fig.add_axes([0.01 + W * k, 0.555 - 0.415 * i, W - 0.004, H0])
        ax.imshow(render(M[mesh], field, cmap, clim), interpolation="none")   # PDF keeps the full render
        ax.set_axis_off()
        ax.text(0.0, 1.0, f"({'abcde'[k]}{i + 1})", transform=ax.transAxes, va="bottom", fontsize=8)
        if k == 0:
            ax.text(-0.02, 0.5, rf"$T^* = {T:.3f}$" if i == 0 else rf"$T^* = {T:.0f}$", transform=ax.transAxes,
                    rotation=90, va="center", ha="right", fontsize=8)
        if i == 1:
            cb = fig.add_axes([0.035 + W * k, 0.105, W - 0.05, 0.022])
            c = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*clim)), cax=cb,
                             orientation="horizontal")
            c.set_ticks(np.linspace(clim[0], clim[1], 3))
            c.ax.set_xticklabels([f"{v:g}" for v in np.linspace(clim[0], clim[1], 3)])
            c.ax.tick_params(labelsize=7.5, length=2, pad=1)
            cb.set_xlabel(label, fontsize=8, labelpad=1)
out = ROOT / "assets" / f"fig_coupled_3d_panels_{CASE}{'_startpos' if RULE == 'startpos' else ''}.pdf"   # the base case keeps the plain name
fig.savefig(out, dpi=400)
fig.savefig(out.with_suffix(".png"), dpi=400)
print("wrote", out, "times", T_SHOW, "stress range 0 ..", VMAX)

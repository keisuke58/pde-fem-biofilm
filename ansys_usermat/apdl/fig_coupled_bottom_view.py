"""fig_coupled_bottom_view.py -- the two coupled examples with antibiotic seen from below (Oliver, 25 Nov 2026: where does
species A survive?). T* = 1, 64^3, the snapshots of the videos. Rows: separate start, mixed start. Columns: (a) share
phi_1/(phi_1+phi_2) and (b) living volume fraction sum_i phibar_i where phi_3D > 0.5 (dead biofilm, sum_i phibar_i < 0.1,
as a translucent grey shell), as in the videos, but with the
camera below the cube (the nutrient enters through the bottom layer y < -0.75 mm). One octant (x, y, z > 0, at the top)
is removed as in the other views; it is hidden from below.

    python ansys_usermat/apdl/fig_coupled_bottom_view.py -> assets/fig_coupled_bottom_view.png
"""
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
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
CAM = [(3.6, -3.4, 4.4), (0, 0, 0), (0, 1, 0)]       # from below, y still up
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)
PAN = [("shareA", "coolwarm_r", r"(a) share $\phi_1/(\phi_1+\phi_2)$"), ("L", "Greens", r"(b) living $\Sigma_i\bar\phi_i$")]


def render(z, field, cmap, T):
    n = z[f"phi_{T}"].shape[0]
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1.0,) * 3)
    cc = g.cell_centers().points
    octant = (cc[:, 0] > 0) & (cc[:, 1] > 0) & (cc[:, 2] > 0)
    for k in ("phi", field):
        g.cell_data[k] = np.nan_to_num(z[f"{k}_{T}"], nan=0.0).ravel(order="F")
    bio = (g.cell_data["phi"] > 0.5) & ~octant
    if field == "L":                          # 25 Nov: dead biofilm (sum phibar < 0.1) as a translucent grey shell
        g.cell_data["bio"] = (bio & (g.cell_data["L"] >= 0.1)).astype(float)
        g.cell_data["dead"] = (bio & (g.cell_data["L"] < 0.1)).astype(float)
    else:
        g.cell_data["bio"] = bio.astype(float)
    m = g.threshold(0.5, scalars="bio")
    pl = pv.Plotter(off_screen=True, window_size=(900, 820))
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
    if m.n_cells:
        pl.add_mesh(m, scalars=field, cmap=cmap, clim=(0, 1), show_scalar_bar=False, **LIT)
    if field == "L":
        d = g.threshold(0.5, scalars="dead")
        if d.n_cells:
            pl.add_mesh(d.extract_surface(), color=(0.55, 0.55, 0.55), opacity=0.12, **LIT)
    pl.add_axes(line_width=3, labels_off=False)
    pl.camera_position = CAM
    pl.camera.zoom(0.98)
    img = pl.screenshot(return_img=True); pl.close()
    return img


fig = plt.figure(figsize=(6.3, 6.0))
for r, (ex, fn) in enumerate(F.items()):
    z = np.load(ROOT / "figdata" / fn)
    T = max((k[4:] for k in z.files if k.startswith("phi_")), key=float)
    for c, (field, cmap, title) in enumerate(PAN):
        ax = fig.add_axes([0.01 + 0.5 * c, 0.52 - 0.47 * r, 0.47, 0.42])
        ax.imshow(render(z, field, cmap, T)); ax.set_axis_off()
        ax.set_title(f"{ex} start, {title}", fontsize=8, loc="left")
for c, (field, cmap, _) in enumerate(PAN):
    cb = fig.add_axes([0.08 + 0.5 * c, 0.035, 0.34, 0.015])
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1)), cax=cb, orientation="horizontal").ax.tick_params(labelsize=7)
out = ROOT / "assets" / "fig_coupled_bottom_view.png"
fig.savefig(out, dpi=300)
print("wrote", out)

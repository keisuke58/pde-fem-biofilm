"""fig_thesis_bottom_view.py -- thesis figure (Chapter 4): where the living biomass survives, seen from below. Living
volume fraction sum_i phibar_i at T* = 1 of the separate (a) and the mixed (b) example with antibiotic, 64^3, outer step
0.025: the runs of Table 4.3 (their fields at T* = 1, figdata/, written by coupled_living_consumption.py). Biofilm phi_3D > 0.5; living part (sum phibar >= 0.1) in colour, dead biofilm as a
translucent grey shell; blue frame: the layer y < -0.75 mm where c is held at 1. Camera below the cube (the exception to
the common view of the other 3D figures, stated in the caption). Printed width 16 cm, Times, renders at 400 dpi.

    python ansys_usermat/apdl/fig_thesis_bottom_view.py -> assets/fig_thesis_bottom_view.pdf/.png
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
F = [("separate", "coupled_living_separate_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_separate_ab.npz"),
     ("mixed", "coupled_living_mixed_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_mixed_ab.npz")]
CAM = [(3.6, -3.4, 4.4), (0, 0, 0), (0, 1, 0)]       # from below, y up
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)
LIV_MIN = 0.1


def render(z):
    T = max((k[4:] for k in z.files if k.startswith("phi_")), key=float)
    n = z[f"phi_{T}"].shape[0]
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1.0,) * 3)
    cc = g.cell_centers().points
    octant = (cc[:, 0] > 0) & (cc[:, 1] > 0) & (cc[:, 2] > 0)
    phi = np.nan_to_num(z[f"phi_{T}"]).ravel(order="F"); L = np.nan_to_num(z[f"L_{T}"]).ravel(order="F")
    g.cell_data["L"] = L
    bio = (phi > 0.5) & ~octant
    g.cell_data["liv"] = (bio & (L >= LIV_MIN)).astype(float)
    g.cell_data["dead"] = (bio & (L < LIV_MIN)).astype(float)
    pl = pv.Plotter(off_screen=True, window_size=(1600, 1520))
    pl.enable_anti_aliasing("ssaa")
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey", line_width=3)
    pl.add_mesh(pv.Box((-1, 1, -1, -0.75, -1, 1)).outline(), color="#1c7ed6", line_width=5)
    m = g.threshold(0.5, scalars="liv")
    if m.n_cells:
        pl.add_mesh(m, scalars="L", cmap="Greens", clim=(0, 1), show_scalar_bar=False, **LIT)
    d = g.threshold(0.5, scalars="dead")
    if d.n_cells:
        pl.add_mesh(d.extract_surface(), color=(0.55, 0.55, 0.55), opacity=0.12, **LIT)
    pl.add_axes(line_width=4)
    pl.camera_position = CAM
    pl.camera.zoom(0.95)
    img = pl.screenshot(return_img=True); pl.close()
    rows = np.where((img < 250).any(axis=(1, 2)))[0]; cols = np.where((img < 250).any(axis=(0, 2)))[0]
    return img[max(rows[0] - 10, 0):rows[-1] + 10, max(cols[0] - 10, 0):cols[-1] + 10]


fig = plt.figure(figsize=(6.3, 3.2))
for k, (ex, fn) in enumerate(F):
    ax = fig.add_axes([0.02 + 0.49 * k, 0.17, 0.46, 0.8])
    ax.imshow(render(np.load(ROOT / "figdata" / fn)), interpolation="none"); ax.set_axis_off()
    ax.text(0.0, 1.0, f"({'ab'[k]})", transform=ax.transAxes, va="top", fontsize=9)
cb = fig.add_axes([0.3, 0.08, 0.4, 0.03])
c = fig.colorbar(plt.cm.ScalarMappable(cmap="Greens", norm=plt.Normalize(0, 1)), cax=cb, orientation="horizontal")
c.ax.tick_params(labelsize=7.5, length=2, pad=1); cb.set_xlabel(r"$\Sigma_i\,\bar\phi_i$", fontsize=8, labelpad=1)
out = ROOT / "assets" / "fig_thesis_bottom_view.pdf"
fig.savefig(out, dpi=400); fig.savefig(out.with_suffix(".png"), dpi=400)
print("wrote", out)

"""anim_coupled_tvd_3d.py -- 3D video of the coupled model (Python reference, 64^3, TVD front term, mirrored walls,
the main settings of results/2026-10-3d) from the snapshots of a --snap-every run: (a) the biofilm, surface
phi_3D = 0.5 coloured by the share of species A, phi_A/(phi_A + phi_B); (b) the living volume fraction sum_i phibar_i
on the same surface; (c) nutrient c and antibiotic a on the mid-planes z = 0 and x = 0 (c supplied from the bottom
layer y < -0.75 mm, a from the outer shell from T* = 1/3). A demonstration, not a thesis figure.

    python ansys_usermat/apdl/anim_coupled_tvd_3d.py separate_ab|mixed_ab
        -> assets/anim_coupled_tvd_<case>.mp4, assets/anim_coupled_tvd_<case>_end.png
"""
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASE = sys.argv[1]
MIX = "_mixeq" if CASE.startswith("mixed") else ""
TAG = (f"_nf64_front1_phimin0.001_v2_fortran_ansysbc_startpos{MIX}_a100_gab1_gc1_tab0.333333_snaps"
       "_tvd_wgmirror_fast")
src = ROOT / "figdata" / f"coupled_living_{CASE}{TAG}_{CASE}.npz"
Z = np.load(src)
times = sorted({float(k.split("_", 1)[1]) for k in Z.files if k.startswith("phi_")})
n = Z[f"phi_{times[0]}"].shape[0]
h = 2.0 / n
pv.global_theme.font.family = "times"


LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)   # light shading, colours close to the colour bar


def BAR(title, x, w):
    return dict(title=title, vertical=False, position_x=x, position_y=0.03, width=w, height=0.05, n_labels=3,
                fmt="%.1f", title_font_size=16, label_font_size=13)


def grid(T):
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(h,) * 3, origin=(-1.0,) * 3)
    for k in ("phi", "shareA", "L", "c", "a"):
        v = np.nan_to_num(Z[f"{k}_{T}"], nan=0.0)
        g.cell_data[k] = v.ravel(order="F")
    cc = g.cell_centers().points
    # shown: the biofilm phi_3D > 0.5 outside the front octant x > 0, y > 0, z > 0 (cut away so the inside shows)
    g.cell_data["show"] = ((g.cell_data["phi"] > 0.5) & ~((cc[:, 0] > 0) & (cc[:, 1] > 0) & (cc[:, 2] > 0))).astype(float)
    return g


out = ROOT / "assets" / f"anim_coupled_tvd_{CASE}.mp4"
import imageio.v2 as imageio  # noqa: E402

movie = imageio.get_writer(str(out), fps=8, quality=9)
cam = [(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)]
title = {"separate_ab": "separate start: A in the centre seed, B in the four edges along z; c supplied from y = -1",
         "mixed_ab": "mixed start: A and B together in the centre seed; c supplied from y = -1"}[CASE]
for T in times:
    g = grid(T)
    # a new plotter per frame: Plotter.clear() also removes the lights (flat, unshaded surfaces)
    pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(1920, 672), border=False)
    # the biofilm with the front octant cut away (grid())
    bio = g.threshold(0.5, scalars="show")
    for col, (field, label, cmap) in enumerate((("shareA", "share of A", "coolwarm_r"),
                                                ("L", "living volume fraction", "Greens"))):
        pl.subplot(0, col)
        pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
        if bio.n_cells:
            pl.add_mesh(bio, scalars=field, cmap=cmap, clim=(0, 1), **LIT, scalar_bar_args=BAR(label, 0.2, 0.6))
        pl.camera_position = cam
        pl.add_text(("(a) biofilm (phi > 0.5), front octant cut away" if col == 0
                     else "(b) living volume fraction"), font_size=11, position="upper_left")
    pl.subplot(0, 2)
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
    pl.add_mesh(g.cell_data_to_point_data().slice(normal="z", origin=(0, 0, 0)), scalars="c", cmap="Blues", clim=(0.4, 1), **LIT,
                scalar_bar_args=BAR("nutrient c (z = 0)", 0.05, 0.4))
    pl.add_mesh(g.cell_data_to_point_data().slice(normal="x", origin=(0, 0, 0)), scalars="a", cmap="Reds", clim=(0, 1), **LIT,
                scalar_bar_args=BAR("antibiotic a (x = 0)", 0.55, 0.4))
    pl.camera_position = cam
    pl.add_text(f"(c) T* = {T:.3f}", font_size=11, position="upper_left")
    pl.subplot(0, 0)
    pl.add_text(title, font_size=9, position="lower_left")
    img = pl.screenshot(return_img=True)
    for _ in range(9 if T == times[-1] else 1):
        movie.append_data(img)
    if T == times[-1]:
        imageio.imwrite(out.with_name(out.stem + "_end.png"), img)
    pl.close()
movie.close()
print("wrote", out)

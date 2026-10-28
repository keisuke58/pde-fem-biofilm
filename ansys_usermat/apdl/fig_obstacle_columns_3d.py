#!/usr/bin/env python3
"""fig_obstacle_columns_3d.py -- 3D view of growth against four rigid bars, in the style of Klempt et al. 2024,
Figs. 9 and 10 (Sect. 4.3), with the agreed coupling (16^3, main constants).

Column (a): set-up -- bars (dark grey), start biomass (red), nutrient layer y = -1 mm (green).
Columns (b)-(d): surface phi_3D = 0.5 at three times (orange), bars and the cube edges.
Rows: separate example (A in the seed, B in four corner blocks) and mixed example (A and B in the seed).
Lengths on the deck's 2 mm cube (1 mm = 10 um on the paper's 20 um cube); y is up, the nutrient comes from below.

Data: coupled_living_consumption.py --obstacle columns --nf 16 --snap-every 0.025 --only separate|mixed (main constants)
    python ansys_usermat/apdl/fig_obstacle_columns_3d.py -> assets/fig_obstacle_columns_3d.pdf/.png
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
import pyvista as pv  # noqa: E402

pv.OFF_SCREEN = True
# main settings (25 Oct): 32^3, outer step 0.0025, TVD, mirror walls, new points inherit the living composition
TAG = "_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv{mix}_a100_gab1_gc1_tab0.333333_snaps_tvd_wgmirror_fast_obstcolumns_{ex}"
TIMES = (0.05, 0.1, 0.2)
LEVEL = {"separate": 0.5, "mixed": 0.2}   # the mixed seed starts at phi_3D = 0.55 and stays below 0.5
BARS = []
for lo, hi in ((-0.75, -0.5), (0.5, 0.75)):
    BARS.append((-1.0, 1.0, -0.125, 0.125, lo, hi))        # along x
    BARS.append((-0.125, 0.125, lo, hi, -1.0, 1.0))        # along z
OUT = ROOT / "figdata" / "obst3d"
OUT.mkdir(parents=True, exist_ok=True)


def grid(u):
    n = u.shape[0]
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1, -1, -1))
    g.cell_data["u"] = u.ravel(order="F")
    return g


def scene(pl, add):
    pl.set_background("white")
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).extract_feature_edges(), color="#595959", line_width=2)
    for b in BARS:
        pl.add_mesh(pv.Box(b), color="#303030", smooth_shading=False)
    add(pl)
    pl.camera_position = [(3.6, 2.6, 4.6), (0, 0, 0), (0, 1, 0)]
    pl.camera.zoom(0.92)


def render(name, add):
    pl = pv.Plotter(off_screen=True, window_size=(900, 900))
    scene(pl, add)
    p = OUT / f"{name}.png"
    pl.screenshot(str(p), transparent_background=False)
    pl.close()
    return p


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figstyle.apply(size=8)
    rows = []
    for ex in ("separate", "mixed"):
        d = np.load(ROOT / "figdata" / f"coupled_living_{ex}{TAG.format(ex=ex, mix='_mixeq' if ex == 'mixed' else '')}.npz")
        n = d[f"phi_{TIMES[0]}"].shape[0]
        xc = -1 + (np.arange(n) + 0.5) * 2 / n
        X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
        import wp3_reference as W
        seed = W.seed_mask(n, W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json"))
        start = seed.copy()
        if ex == "separate":
            zmax = np.abs(Z[seed]).max() + 1e-9
            start |= (np.abs(X) > 0.75) & (np.abs(Y) > 0.75) & (np.abs(Z) <= zmax)

        def setup(pl, start=start):
            g = grid(start.astype(float))
            pl.add_mesh(g.threshold(0.5), color="#d62728", show_edges=False)
            pl.add_mesh(pv.Box((-1, 1, -1, -0.95, -1, 1)), color="#2ca02c", opacity=0.8)
        imgs = [render(f"{ex}_setup", setup)]
        for T in TIMES:
            u = d[f"phi_{T}"]

            def surf(pl, u=u, lev=LEVEL[ex]):
                # one layer of zeros outside the cube closes the surface on the walls, as in Klempt et al. Fig. 10
                up = np.pad(u, 1)
                n2 = up.shape[0]
                g = pv.ImageData(dimensions=(n2 + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1 - 2.0 / n,) * 3)
                g.cell_data["u"] = up.ravel(order="F")
                g = g.cell_data_to_point_data()
                s = g.contour([lev], scalars="u").smooth(n_iter=30)
                pl.add_mesh(s, color="#e8a33d", smooth_shading=True, specular=0.2)
            imgs.append(render(f"{ex}_{T}", surf))
        rows.append((ex, imgs))
    fig, axs = plt.subplots(2, 4, figsize=(6.3, 3.4), gridspec_kw=dict(wspace=0.02, hspace=0.08))
    titles = ["(a) set-up"] + [f"({'bcd'[i]}) $T^*={T:g}$" for i, T in enumerate(TIMES)]
    for r, (ex, imgs) in enumerate(rows):
        for c, p in enumerate(imgs):
            a = axs[r, c]
            im = plt.imread(p)
            a.imshow(im[60:-40, 60:-60])
            a.axis("off")
            if r == 0:
                a.set_title(titles[c], fontsize=8)
        axs[r, 0].text(-0.08, 0.5, f"{ex}, $\\phi_{{3D}}={LEVEL[ex]:g}$", rotation=90, va="center", ha="center", transform=axs[r, 0].transAxes,
                       fontsize=8)
    out = ROOT / "assets" / "fig_obstacle_columns_3d"
    fig.savefig(out.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), dpi=400, bbox_inches="tight")
    print("wrote", out.with_suffix(".png"))


if __name__ == "__main__":
    main()

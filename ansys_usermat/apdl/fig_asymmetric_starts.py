"""fig_asymmetric_starts.py -- thesis/paper figure of the asymmetric starts (separate example with antibiotic, T* = 1).

Rows: the symmetric example, species B only in the block along one edge (seed moved by one element towards -x), and
species B only in one corner cube. Columns: (a) the start and the boundary conditions (species A in the seed, species B,
the nutrient layer y < -0.75 mm held at c = 1, the antibiotic shell |x_i| > 0.75 mm held at a = 1 from T* = 1/3),
(b) growth alpha - 1 where phi_3D > 0.5 at T*, (c) von Mises stress at T* in the section z = 0 (16^3 stress mesh),
seed outline in white. One colour range per column for all rows. Camera and removed octant of (a), (b) as in every 3D
view of the cube examples (CLAUDE.md).

Needs the snapshots of the three runs (figdata/*.npz), e.g. for 64^3:
    B="--fortran --front 1 --gc 1 --gab 1 --alpha0 100 --tab 0.333333 --phi-min 0.001 --tvd --wall-grad mirror --v2 --ansys-bc --inherit-living"
    python3 coupled_living_consumption.py --nf 64 --amg  $B --transient-c 1e8 --snap-every 1.0 --snap-stress --b-one-edge --seed-shift -1 --only separate_ab
    python3 coupled_living_consumption.py --nf 64 --amg  $B --transient-c 1e8 --snap-every 1.0 --snap-stress --b-corner --only separate_ab
    python3 coupled_living_consumption.py --nf 64 --fast $B --transient-c 1e8 --snap-every 1.0 --snap-stress --only separate_ab
    python3 fig_asymmetric_starts.py [--nf 64]   -> thesis/figures/fig_asymmetric_starts.pdf (and .png)
"""
import argparse
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
import wp3_reference as WR  # noqa: E402

figstyle.apply()
ap = argparse.ArgumentParser()
ap.add_argument("--nf", type=int, default=64)
ap.add_argument("--t", type=float, default=1.0)
ap.add_argument("--out", default=str(ROOT / "thesis" / "figures" / "fig_asymmetric_starts"))
A = ap.parse_args()

BASE = (f"coupled_living_separate_ab_nf{A.nf}_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333"
        "_snapsvm_transc1e+08")
ROWS = [("symmetric", "_tvd_wgmirror_fast", "four", 0),
        ("B along one edge", "_b1edge_seedx-1_tvd_wgmirror_amg", "edge", -1),
        ("B in one corner", "_bcorner_tvd_wgmirror_amg", "corner", 0)]
CELLS8 = WR.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")


def start(m, where, shift):
    """Seed and species-B masks on m^3 cells, as coupled_living_consumption.run (--v2)."""
    seed = WR.seed_mask(m, CELLS8 + np.array([shift, 0, 0]))
    ax_ = -1 + (np.arange(m) + 0.5) * 2.0 / m
    X, Y, Z = np.meshgrid(ax_, ax_, ax_, indexing="ij")
    b = (np.abs(X) > 0.75) & (np.abs(Y) > 0.75) & (np.abs(Z) <= np.abs(Z[seed]).max() + 1e-9)
    if where == "edge":
        b &= (X > 0) & (Y < 0)
    if where == "corner":
        b = (X > 0.75) & (Y < -0.75) & (Z > 0.75)
    return seed, b, Y
CAM = [(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)]
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)


def load(tag):
    Z = np.load(ROOT / "figdata" / f"{BASE}{tag}_separate_ab.npz")
    t = min((float(k.split("_", 1)[1]) for k in Z.files if k.startswith("phi_")), key=lambda s: abs(s - A.t))
    return {k: Z[f"{k}_{t}"] for k in ("phi", "shareA", "alpha", "vm16")}


def cube(m):
    g = pv.ImageData(dimensions=(m + 1,) * 3, spacing=(2.0 / m,) * 3, origin=(-1.0,) * 3)
    cc = g.cell_centers().points
    return g, (cc[:, 0] > 0) & (cc[:, 1] > 0) & (cc[:, 2] > 0)


def render(mesh, field, cmap, clim):
    pl = pv.Plotter(off_screen=True, window_size=(900, 820))
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
    if mesh.n_cells:
        pl.add_mesh(mesh, scalars=field, cmap=cmap, clim=clim, show_scalar_bar=False, **LIT)
    pl.camera_position = CAM
    pl.camera.zoom(0.98)
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


D = [load(r[1]) for r in ROWS]
AMAX = float(np.ceil(max(np.percentile(1e4 * (d["alpha"] - 1.0), 99.5) for d in D)))
VMAX = float(2 * np.ceil(max(1e4 * np.nanpercentile(d["vm16"], 99.5) for d in D) / 2))
COLS = [("start", None, None, None),
        ("am1", "viridis", (0.0, AMAX), r"growth $\alpha-1$ [$10^{-4}$]"),
        ("vm", "magma_r", (0.0, VMAX), r"von Mises stress, $z=0$ [$10^{-4}$ Pa]")]
SEEDC, BC, NUTC, ABC = "#3b5bdb", "#c92a2a", "#2b8a3e", "#e8590c"


def render_start(seed, b, Y, oct_):
    m = seed.shape[0]
    g, _ = cube(m)
    pl = pv.Plotter(off_screen=True, window_size=(900, 820))
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
    shell = pv.Box((-1, 1, -1, 1, -1, 1)).extract_all_edges()
    pl.add_mesh(pv.Box((-0.75, 0.75, -0.75, 0.75, -0.75, 0.75)).outline(), color=ABC, line_width=1.5)
    for mask, col, op in (((Y < -0.75) & ~oct_.reshape(Y.shape, order="F") & ~seed & ~b, NUTC, 0.18), (seed, SEEDC, 1.0), (b, BC, 1.0)):
        g.cell_data["m"] = mask.ravel(order="F").astype(float)
        sub = g.threshold(0.5, scalars="m")
        if sub.n_cells:
            pl.add_mesh(sub, color=col, opacity=op, **LIT)
    pl.camera_position = CAM
    pl.camera.zoom(0.98)
    img = pl.screenshot(return_img=True)
    pl.close()
    return img

W, H = 6.3, 6.0                                     # printed width of the thesis text
fig = plt.figure(figsize=(W, H))
x0, cw, top, rh = 0.075, 0.305, 0.955, 0.275
for r, ((label, _, where, shift), d) in enumerate(zip(ROWS, D)):
    n = d["phi"].shape[0]
    g, oct_ = cube(n)
    g.cell_data["am1"] = (1e4 * (d["alpha"] - 1.0)).ravel(order="F")
    g.cell_data["bio"] = ((np.nan_to_num(d["phi"]).ravel(order="F") > 0.5) & ~oct_).astype(float)
    bio = g.threshold(0.5, scalars="bio")
    seed, b, Y = start(n, where, shift)
    vm = 1e4 * np.nan_to_num(d["vm16"])
    nm = vm.shape[0]
    sec = 0.5 * (vm[:, :, nm // 2 - 1] + vm[:, :, nm // 2])       # z = 0 lies between the two middle element layers
    for c, (field, cmap, clim, _) in enumerate(COLS):
        ax = fig.add_axes([x0 + c * cw, top - (r + 1) * rh, cw - 0.01, rh - 0.01])
        if field == "start":
            ax.imshow(render_start(seed, b, Y, oct_)); ax.set_axis_off()
        elif field == "am1":
            ax.imshow(render(bio, field, cmap, clim)); ax.set_axis_off()
        else:
            ax.imshow(sec.T, origin="lower", extent=(-1, 1, -1, 1), cmap=cmap, vmin=clim[0], vmax=clim[1], interpolation="nearest")
            sz = seed[:, :, n // 2]
            ax.contour(np.linspace(-1, 1, n), np.linspace(-1, 1, n), sz.T.astype(float), [0.5], colors="white", linewidths=0.8)
            ax.set_aspect("equal"); ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
            ax.tick_params(labelsize=7, length=2, pad=1)
            ax.set_xlabel("$x$ [mm]", fontsize=7.5, labelpad=0); ax.set_ylabel("$y$ [mm]", fontsize=7.5, labelpad=0)
            ax.set_position([x0 + c * cw + 0.045, top - (r + 1) * rh + 0.03, cw - 0.07, rh - 0.045])
        ax.text(0.02 if field != "vm" else -0.28, 0.98, f"({'abc'[c]}{r + 1})", transform=ax.transAxes, ha="left", va="top", fontsize=9)
    fig.text(0.035, top - (r + 0.5) * rh, label, rotation=90, ha="center", va="center", fontsize=9)
for c, (field, cmap, clim, lab) in enumerate(COLS):
    cb = fig.add_axes([x0 + c * cw + 0.03, top - 3 * rh - 0.045, cw - 0.07, 0.016])
    if field == "start":                      # legend instead of a colour bar
        cb.set_axis_off()
        items = ((SEEDC, "species A (seed)", True), (BC, "species B", True), (NUTC, "nutrient layer, $c=1$", True),
                 (ABC, r"antibiotic shell, $a=1$ from $T^*=1/3$", False))
        for i, (col, txt, fill) in enumerate(items):
            xx, yy = (0.0, 1.6 - 1.5 * i) if i else (0.0, 1.6)
            if i == 1:
                xx, yy = 0.66, 1.6
            if i >= 2:
                xx, yy = 0.0, 1.6 - 1.5 * (i - 1)
            cb.add_patch(plt.Rectangle((xx, yy - 0.45), 0.06, 0.9, facecolor=col if fill else "none", edgecolor=col,
                                       transform=cb.transAxes, clip_on=False, lw=1.2))
            cb.text(xx + 0.09, yy, txt, transform=cb.transAxes, fontsize=7.5, va="center")
        continue
    cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*clim)), cax=cb, orientation="horizontal")
    cbar.set_ticks(np.linspace(clim[0], clim[1], 3))
    cbar.ax.set_xticklabels([f"{v:g}" for v in np.linspace(clim[0], clim[1], 3)])
    cbar.ax.tick_params(labelsize=7.5)
    cb.set_xlabel(lab, fontsize=8.5, labelpad=2)
# coordinate axes as seen from the camera (y up), lower left
tx = fig.add_axes([0.0, 0.10, 0.08, 0.09])
dv = -np.asarray(CAM[0], float); dv /= np.linalg.norm(dv)
rv = np.cross(dv, CAM[2]); rv /= np.linalg.norm(rv); uv = np.cross(rv, dv)
for v, lab, col in (((1, 0, 0), "x", "#d9480f"), ((0, 1, 0), "y", "#2b8a3e"), ((0, 0, 1), "z", "#1c3faa")):
    e = np.array(v, float)
    tx.annotate("", xy=(e @ rv, e @ uv), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=1, mutation_scale=8))
    tx.text(1.35 * (e @ rv), 1.35 * (e @ uv), f"${lab}$", color=col, fontsize=8, ha="center", va="center")
tx.set_xlim(-1.1, 1.1); tx.set_ylim(-0.7, 1.3); tx.set_aspect("equal"); tx.set_axis_off()
out = Path(A.out)
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out.with_suffix(".pdf"), dpi=400)
fig.savefig(out.with_suffix(".png"), dpi=200)
print("wrote", out.with_suffix(".pdf"), f"alpha-1 up to {AMAX:g}e-4, vM up to {VMAX:g}e-4 Pa")

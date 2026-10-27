"""fig_real_tooth_coupled.py -- figure 9 of the coupling paper: the agreed coupling on the clinical molar crown
(real_tooth_coupled.py, run tooth_ab at NF^3 on the 10 mm box). (a) set-up: the crown with the start ring, species A
(x < 0) and B (x > 0); at T* = 1 with the half y > 0 cut away so the section through the film shows: (b) share
phi_1/(phi_1 + phi_2), (c) living volume fraction sum_i phibar_i, (d) nutrient c, (e) antibiotic a, all on the cells
with phi_3D > 0.5, (f) von Mises stress of the film (rigid crown, h = h_f). Printed width 16 cm, Times (figstyle),
renderings at full resolution inside the PDF; no title (the caption carries it).

    python papers/coupling_paper/fig_real_tooth_coupled.py [--nf 80]   -> papers/coupling_paper/fig_real_tooth_coupled.pdf/.png
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
sys.path.insert(0, str(ROOT / "ansys_usermat" / "apdl"))
import figstyle  # noqa: E402

figstyle.apply()
NF = int(sys.argv[sys.argv.index("--nf") + 1]) if "--nf" in sys.argv else 80
BOX = 10.0
H = BOX / NF
CASE = "tooth_ab"
OPTS = "front1_phimin0.001_v2_fortran_inhliv_a100_gab1_gc1_tab0.333333_tvd_wgmirror_obststl_box10_amg"
Z = np.load(ROOT / "figdata" / f"coupled_living_{CASE}_nf{NF}_{OPTS}_{CASE}.npz")
T = sorted({float(f.split("_", 1)[1]) for f in Z.files if f.startswith("phi_")})[-1]
ST = np.load(ROOT / "figdata" / f"real_tooth_coupled_stress_{CASE}_nf{NF}_T{T}_h1.npz")
OBST = Z[f"obst_{T}"].astype(bool)
FILM = (Z[f"phi_{T}"] > 0.5) & ~OBST
vm_mPa = np.nan_to_num(ST["vm"].astype(float)) * 1e3
VMAX = float(np.ceil(np.percentile(vm_mPa[ST["active"] & ~np.isnan(ST["vm"])], 99)))
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)
CAM = [(-7.0, -17.0, 12.0), (0.0, -0.5, -1.5), (0.0, 0.0, 1.0)]
PANELS = [("shareA", "coolwarm_r", (0, 1), r"$\phi_1/(\phi_1+\phi_2)$"),
          ("L", "Greens", (0, 1), r"$\Sigma_i\,\bar\phi_i$"),
          ("c", "Blues", (0, 1), r"$c$"),
          ("a", "Reds", (0, 1), r"$a$"),
          ("vm", "magma_r", (0, VMAX), r"$\sigma_\mathrm{vM}$ [mPa]")]


def crown(cut):
    """the crown surface placed as in the runs (coupled_living_consumption.tooth_sdf), optionally cut at y = 0"""
    m = pv.read(ROOT / "tier2b_real" / "_crown_real_watertight.stl").triangulate().clean()
    P = np.asarray(m.points)
    lo, hi = P.min(0), P.max(0)
    m.points = P - np.array([0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), lo[2] + BOX / 2])
    if cut:
        m = m.clip(normal=(0, 1, 0), origin=(0, 0, 0))
    return m


def voxels(values, keep):
    g = pv.ImageData(dimensions=(NF + 1,) * 3, spacing=(H,) * 3, origin=(-BOX / 2,) * 3)
    g.cell_data["v"] = np.asarray(values, float).ravel(order="F")
    g.cell_data["keep"] = keep.ravel(order="F").astype(float)
    return g.threshold(0.5, scalars="keep")


def stress_voxels(cut_mask):
    """the stress on the crop grid of the submodel (h = h_f: the field cells), masked like the fields"""
    o, h, box = ST["origin"], float(ST["h"]), ST["box"]
    vm = vm_mPa
    sub = cut_mask[box[0, 0]:box[0, 1], box[1, 0]:box[1, 1], box[2, 0]:box[2, 1]]
    g = pv.ImageData(dimensions=tuple(np.array(vm.shape) + 1), spacing=(h,) * 3, origin=tuple(o))
    g.cell_data["v"] = vm.ravel(order="F")
    g.cell_data["keep"] = (sub & ST["active"] & ~np.isnan(ST["vm"])).ravel(order="F").astype(float)
    return g.threshold(0.5, scalars="keep")


def render(add, zoom=1.0):
    pl = pv.Plotter(off_screen=True, window_size=(1500, 1300))
    pl.enable_anti_aliasing("ssaa")
    pl.set_background("white")
    pl.add_mesh(pv.Box((-BOX / 2, BOX / 2) * 3).outline(), color="grey", line_width=2)
    add(pl)
    pl.camera_position = CAM
    pl.camera.zoom(zoom)
    img = pl.screenshot(return_img=True)
    pl.close()
    rows = np.where((img < 250).any(axis=(1, 2)))[0]
    cols = np.where((img < 250).any(axis=(0, 2)))[0]
    p = 10
    return img[max(rows[0] - p, 0):rows[-1] + p, max(cols[0] - p, 0):cols[-1] + p]


ax_ = -BOX / 2 + (np.arange(NF) + 0.5) * H
Yc = ax_[None, :, None]
front = np.broadcast_to(Yc < 0, (NF,) * 3)             # the half y < 0 is kept, y > 0 cut away
tooth_full, tooth_cut = crown(False), crown(True)
seedA, seedB = Z[f"seed_{T}"].astype(bool), Z[f"edges_{T}"].astype(bool)
images = []


def setup(pl):
    pl.add_mesh(tooth_full, color="#d8d3c8", smooth_shading=True, ambient=0.35, diffuse=0.7, specular=0.1)
    pl.add_mesh(voxels(np.ones((NF,) * 3), seedA), color="#3b6fb6", **LIT)
    pl.add_mesh(voxels(np.ones((NF,) * 3), seedB), color="#c0392b", **LIT)


images.append(render(setup, zoom=1.35))
for field, cmap, clim, _ in PANELS:
    def add(pl, field=field, cmap=cmap, clim=clim):
        pl.add_mesh(tooth_cut, color="#d8d3c8", smooth_shading=True, ambient=0.35, diffuse=0.7, specular=0.1)
        if field == "vm":
            m = stress_voxels(FILM & front)
        else:
            m = voxels(np.nan_to_num(Z[f"{field}_{T}"], nan=0.0), FILM & front)
        if m.n_cells:
            pl.add_mesh(m, scalars="v", cmap=cmap, clim=clim, show_scalar_bar=False, **LIT)
    images.append(render(add))

fig = plt.figure(figsize=(6.3, 4.6))
W, H0 = 0.325, 0.40
for k, img in enumerate(images):
    i, j = divmod(k, 3)
    ax = fig.add_axes([0.012 + W * j, 0.575 - 0.47 * i, W - 0.012, H0])
    ax.imshow(img, interpolation="none")
    ax.set_axis_off()
    ax.text(0.0, 1.0, f"({'abcdef'[k]})", transform=ax.transAxes, va="bottom", fontsize=8)
    if k >= 1:
        field, cmap, clim, label = PANELS[k - 1]
        cb = fig.add_axes([0.05 + W * j, 0.575 - 0.47 * i - 0.04, W - 0.09, 0.018])
        c = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*clim)), cax=cb, orientation="horizontal")
        c.set_ticks(np.linspace(clim[0], clim[1], 3))
        c.ax.set_xticklabels([f"{v:g}" for v in np.linspace(clim[0], clim[1], 3)])
        c.ax.tick_params(labelsize=7.5, length=2, pad=1)
        cb.set_xlabel(label, fontsize=8, labelpad=1)
    else:
        ax.text(0.5, -0.06, r"start: A ($x<0$), B ($x>0$)", transform=ax.transAxes, ha="center", va="top", fontsize=8)
out = HERE / "fig_real_tooth_coupled.pdf"
fig.savefig(out, dpi=400)
fig.savefig(out.with_suffix(".png"), dpi=300)
print("wrote", out, "T* =", T, "vM range 0 ..", VMAX, "mPa")

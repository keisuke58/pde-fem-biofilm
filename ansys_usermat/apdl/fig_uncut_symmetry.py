"""fig_uncut_symmetry.py -- the stress of the separate start with antibiotic at T* = 1 without the octant cut (Assoc. Prof.
Muramatsu, 2 Nov 2026: the cut views looked unsymmetric). (a) the whole cube from the usual camera, (b) seen from above
(+y), (c) the slice y = 0 of the 16^3 stress mesh with the mirror lines x = 0 and z = 0. The set-up is mirror symmetric
in x and z (seed and B blocks), and so is the solution: prints the largest difference between mirrored elements.

    python ansys_usermat/apdl/fig_uncut_symmetry.py -> assets/fig_uncut_symmetry.png
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
Z = np.load(ROOT / "figdata" / "coupled_living_separate_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast_separate_ab.npz")
T = max((k[5:] for k in Z.files if k.startswith("vm16_")), key=float)
vm = 1e4 * Z[f"vm16_{T}"]                       # [10^-4 Pa], 16^3, axes (x, y, z)
nm = vm.shape[0]
VMAX = float(2 * np.ceil(np.percentile(vm, 99.5) / 2))
dx = np.abs(vm - vm[::-1, :, :]).max(); dz = np.abs(vm - vm[:, :, ::-1]).max()
print(f"T* = {T}: max von Mises {vm.max():.2f}e-4 Pa; largest mirror difference x: {dx:.2e}, z: {dz:.2e} (x 1e-4 Pa)")
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)


def render(cam):
    g = pv.ImageData(dimensions=(nm + 1,) * 3, spacing=(2.0 / nm,) * 3, origin=(-1.0,) * 3)
    g.cell_data["vm"] = vm.ravel(order="F")
    pl = pv.Plotter(off_screen=True, window_size=(1200, 1140))
    pl.add_mesh(g, scalars="vm", cmap="magma_r", clim=(0, VMAX), show_scalar_bar=False, **LIT)
    pl.add_axes(line_width=3)
    pl.camera_position = cam
    pl.camera.zoom(0.95)
    img = pl.screenshot(return_img=True); pl.close()
    return img


fig = plt.figure(figsize=(6.3, 2.5))
ax = fig.add_axes([0.0, 0.12, 0.32, 0.8]); ax.imshow(render([(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)])); ax.set_axis_off()
ax.set_title("(a) whole cube, usual view", fontsize=8, loc="left")
ax = fig.add_axes([0.33, 0.12, 0.32, 0.8]); ax.imshow(render([(0, 6.0, 0.001), (0, 0, 0), (0, 0, -1)])); ax.set_axis_off()
ax.set_title("(b) seen from above (+$y$)", fontsize=8, loc="left")
ax = fig.add_axes([0.71, 0.17, 0.25, 0.7])
j = nm // 2
im = ax.imshow(0.5 * (vm[:, j - 1, :] + vm[:, j, :]).T, origin="lower", extent=(-1, 1, -1, 1), cmap="magma_r", vmin=0, vmax=VMAX)
ax.axvline(0, color="#1c7ed6", lw=0.9, ls="--"); ax.axhline(0, color="#1c7ed6", lw=0.9, ls="--")
ax.set_xlabel("$x$ [mm]", fontsize=7); ax.set_ylabel("$z$ [mm]", fontsize=7); ax.tick_params(labelsize=6.5)
ax.set_title("(c) slice $y = 0$, mirror lines", fontsize=8, loc="left")
cb = fig.add_axes([0.08, 0.08, 0.5, 0.03])
fig.colorbar(plt.cm.ScalarMappable(cmap="magma_r", norm=plt.Normalize(0, VMAX)), cax=cb, orientation="horizontal").ax.tick_params(labelsize=6.5)
cb.set_xlabel(r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]", fontsize=7, labelpad=1)
out = ROOT / "assets" / "fig_uncut_symmetry.png"
fig.savefig(out, dpi=400)
print("wrote", out)

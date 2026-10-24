"""anim_coupled_tvd_3d.py -- 3D video of the coupled model (Python reference, 64^3, TVD front term, mirrored walls,
the main settings of results/2026-10-3d) from the snapshots of a --snap-every run. Four views of the 2 mm cube with one
octant removed: (a) share of species A, (b) living volume fraction, both where phi_3D > 0.5; (c) nutrient c,
(d) antibiotic a in the whole cube. The 3D views are rendered with pyvista; titles, colour bars and the caption are
laid out with matplotlib (figstyle: Times, mathtext). A demonstration, not a thesis figure.

    python ansys_usermat/apdl/anim_coupled_tvd_3d.py separate_ab|mixed_ab
        -> assets/anim_coupled_tvd_<case>.mp4, assets/anim_coupled_tvd_<case>_end.png
"""
import sys
from pathlib import Path

import imageio.v2 as imageio
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
CASE = sys.argv[1]
# inherit: newly colonised points take the phi-weighted mean state of their colonised neighbours (A spreads from the
# seed, B from the edge blocks); default: the start state by position, as in the ANSYS runs (A-19)
INHERIT = "inherit" in sys.argv[2:]
MIX = "_mixeq" if CASE.startswith("mixed") else ""
TAG = (f"_nf64_front1_phimin0.001_v2_fortran_ansysbc{'' if INHERIT else '_startpos'}{MIX}_a100_gab1_gc1_tab0.333333"
       "_snaps_tvd_wgmirror_fast")
Z = np.load(ROOT / "figdata" / f"coupled_living_{CASE}{TAG}_{CASE}.npz")
times = sorted({float(k.split("_", 1)[1]) for k in Z.files if k.startswith("phi_")})
n = Z[f"phi_{times[0]}"].shape[0]
h = 2.0 / n
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)   # light shading, colours close to the colour bar
CAM = [(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)]
# field, mesh ("bio": phi_3D > 0.5, "cube": all), colour map, range, title, colour-bar label
# mesh: "bio" = phi_3D > 0.5, "living" = bio and sum phibar > 0.01, "cube" = all, "ab" = a > 0.01 (values near 0
# are left out, i.e. transparent)
PANELS = [("shareA", "bio", "coolwarm_r", (0, 1), r"(a) share of species A", r"$\phi_A/(\phi_A+\phi_B)$"),
          ("L", "living", "Greens", (0, 1), r"(b) living volume fraction", r"$\Sigma_i\,\bar\phi_i$"),
          ("c", "cube", "Blues", (0.4, 1), r"(c) nutrient $c$", r"$c$"),
          ("a", "ab", "Reds", (0, 1), r"(d) antibiotic $a$", r"$a$")]
START = {"separate_ab": r"Separate start: species A in the centre seed, species B in four small blocks at the edges along $z$ "
                        r"($|x|, |y| > 0.75$ mm). " + (r"New biofilm takes the composition of its neighbours." if INHERIT
                        else r"New biofilm starts as B where $|x|, |y| > 0.56$ mm, else as A (rule of the ANSYS runs)."),
         "mixed_ab": r"Mixed start: species A and B in equal amounts together in the centre seed. "
                     + (r"New biofilm takes the composition of its neighbours." if INHERIT
                        else r"New biofilm starts from the mixed start state (rule of the ANSYS runs).")}[CASE]
CAPTION = "\n".join([
    r"Coupled growth model: Klempt et al. 2024, Eqs. 34$-$36, front term on ($r = 10$ mm/$T^*$), with the point model "
    r"of Klempt et al. 2026 (Fritsch et al. 2025, Case I, true parameters) at every point; 2 mm cube, $64^3$ finite "
    r"volumes, TVD front term.",
    START,
    r"(a) Share of species A (point-model species 1; B is species 2; blue: A, red: B) and (b) living volume fraction $\Sigma_i\,\bar\phi_i$, "
    r"$\bar\phi_i=\phi_i\psi_i$, both shown where the biofilm $\phi_{3D} > 0.5$.",
    r"(c) Nutrient $c$, held at 1 in the bottom layer $y < -0.75$ mm and consumed by the living biomass. "
    r"(d) Antibiotic $a$, held at 1 on the outer shell from $T^* = 1/3$ ($\alpha^*_0 = 100$). "
    r"One octant ($x, y, z > 0$) is removed to show the inside.",
    r"Constants: $\beta = 0.02$ mm$^2$/$T^*$, $r = 10$ mm/$T^*$, $K = 1$, $k_\alpha = 10^{-3}$/$T^*$, "
    r"$d = 10^8$ mm$^2$/$T^*$ for $c$ and $a$, $g\,l^2/d = 1$ for $c$ and $a$, "
    r"$\phi_{\min} = 10^{-3}$, $\alpha^*_0 = 100$ from $T^* = 1/3$.",
    r"Sources: $\beta, r, k_\alpha, d, g$ from Klempt et al. 2024, Table 2, converted to the 2 mm cube; $\alpha^*_0$ and the onset "
    r"from Klempt et al. 2026, case 4. $T^* = 1$ is the duration of the point model's data (Fritsch et al. 2025)."])


def grid(T):
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(h,) * 3, origin=(-1.0,) * 3)
    for k in ("phi", "shareA", "L", "c", "a"):
        g.cell_data[k] = np.nan_to_num(Z[f"{k}_{T}"], nan=0.0).ravel(order="F")
    cc = g.cell_centers().points
    octant = (cc[:, 0] > 0) & (cc[:, 1] > 0) & (cc[:, 2] > 0)
    g.cell_data["bio"] = ((g.cell_data["phi"] > 0.5) & ~octant).astype(float)
    g.cell_data["cube"] = (~octant).astype(float)
    g.cell_data["living"] = (g.cell_data["bio"] * (g.cell_data["L"] > 0.01)).astype(float)
    g.cell_data["ab"] = (g.cell_data["cube"] * (g.cell_data["a"] > 0.01)).astype(float)
    return g


def view(mesh, field, cmap, clim):
    pl = pv.Plotter(off_screen=True, window_size=(900, 820))
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
    if mesh.n_cells:
        pl.add_mesh(mesh, scalars=field, cmap=cmap, clim=clim, show_scalar_bar=False, **LIT)
    pl.camera_position = CAM
    pl.camera.zoom(0.98)
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


fig = plt.figure(figsize=(16, 6.2), dpi=160)                  # 2560 x 992 px (multiples of 16)
ax_img = [fig.add_axes([0.005 + 0.25 * k, 0.385, 0.24, 0.56]) for k in range(4)]
ax_bar = [fig.add_axes([0.04 + 0.25 * k, 0.355, 0.17, 0.019]) for k in range(4)]
out = ROOT / "assets" / f"anim_coupled_tvd_{CASE}{'_inherit' if INHERIT else ''}.mp4"
movie = imageio.get_writer(str(out), fps=8, quality=9)
for T in times:
    g = grid(T)
    meshes = {m: g.threshold(0.5, scalars=m) for m in ("bio", "living", "cube", "ab")}
    for k, (field, mesh, cmap, clim, title, label) in enumerate(PANELS):
        ax, cb = ax_img[k], ax_bar[k]
        ax.clear(); cb.clear()
        ax.imshow(view(meshes[mesh], field, cmap, clim))
        ax.set_axis_off()
        ax.set_title(title, fontsize=13, loc="left")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*clim))
        c = fig.colorbar(sm, cax=cb, orientation="horizontal")
        c.set_ticks(np.linspace(clim[0], clim[1], 3))
        c.ax.tick_params(labelsize=10)
        cb.set_xlabel(label, fontsize=12, labelpad=2)
    for t in list(fig.texts):
        t.remove()
    fig.text(0.995, 0.965, rf"$T^* = {T:.3f}$", ha="right", va="top", fontsize=15)
    fig.text(0.01, 0.255, CAPTION, ha="left", va="top", fontsize=11.5, linespacing=1.5)
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    for _ in range(9 if T == times[-1] else 1):
        movie.append_data(img)
    if T == times[-1]:
        imageio.imwrite(out.with_name(out.stem + "_end.png"), img)
movie.close()
print("wrote", out)

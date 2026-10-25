"""anim_coupled_tvd_3d.py -- 3D video of the coupled model (Python reference, 64^3, TVD front term, mirrored walls,
the main settings of results/2026-10-3d) from the snapshots of a --snap-every --snap-stress run. Five views of the
2 mm cube with one octant removed: (a) share of species A, (b) living volume fraction, both where phi_3D > 0.5;
(c) nutrient c, (d) antibiotic a (left out where a < 0.01), (e) von Mises stress on the 16^3 stress mesh. The 3D views
are rendered with pyvista; titles, colour bars and the caption are laid out with matplotlib (figstyle: Times,
mathtext). A demonstration, not a thesis figure.

    python ansys_usermat/apdl/anim_coupled_tvd_3d.py separate_ab|mixed_ab [inherit] [obst] [k26c4] [pulse]
        -> assets/anim_coupled_tvd_<case>[_inherit][_obst][_k26c4][_pulse]_5panel.mp4 and _end.png
  inherit: newly colonised points take the phi-weighted mean state of their colonised neighbours (default: the start
           state by position, the rule of the ANSYS runs); obst: the four rigid columns of Klempt 2024 Sect. 4.3;
  k26c4:   the four species of Klempt 2026 Table 2, case 4; pulse: the antibiotic in pulses of the same dose (A-71).
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
OPT = set(sys.argv[2:])
INHERIT, OBST, K26C4, PULSE = ("inherit" in OPT), ("obst" in OPT), ("k26c4" in OPT), ("pulse" in OPT)
INHLIV = "inhliv" in OPT          # the base rule from 25 Oct: the neighbours' living composition, fresh psi (A-82)
BOXES = [((-1.0, -0.125, lo), (1.0, 0.125, hi)) for lo, hi in ((-0.75, -0.5), (0.5, 0.75))] + \
        [((-0.125, lo, -1.0), (0.125, hi, 1.0)) for lo, hi in ((-0.75, -0.5), (0.5, 0.75))]
MIX = "_mixeq" if CASE.startswith("mixed") else ""
TAG = (f"_nf64_front1_phimin0.001_v2_fortran_ansysbc{'_inhliv' if INHLIV else ('' if INHERIT else '_startpos')}{MIX}_a100_gab1_gc1_tab0.333333"
       f"_snapsvm{'_k26c4' if K26C4 else ''}{'_abpulse0.2222' if PULSE else ''}_tvd_wgmirror_fast"
       f"{'_obstcolumns' if OBST else ''}")
Z = np.load(ROOT / ("figdata") / f"coupled_living_{CASE}{TAG}_{CASE}.npz")
times = sorted({float(k.split("_", 1)[1]) for k in Z.files if k.startswith("phi_")})
n = Z[f"phi_{times[0]}"].shape[0]
nm = Z[f"vm16_{times[0]}"].shape[0]
VMAX = 1e4 * max(np.percentile(Z[f"vm16_{T}"], 99.5) for T in times)   # one colour range over the video, 1e-4 Pa
VMAX = float(2 * np.ceil(VMAX / 2))                      # even, so the middle tick is a whole number
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)   # light shading, colours close to the colour bar
CAM = [(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)]
# field, mesh, colour map, range, title, colour-bar label; mesh: "bio" phi_3D > 0.5, "living" bio and sum phibar > 0.01,
# "cube" all, "ab" a > 0.01, "stress" the stress mesh (values near 0 are left out, i.e. transparent)
GROUP = "species 1, 2" if K26C4 else "species 1"
PANELS = [("shareA", "bio", "coolwarm_r", (0, 1), r"(a) share of group A", r"$\phi_A/(\phi_A+\phi_B)$"),
          ("L", "living", "Greens", (0, 1), r"(b) living volume fraction", r"$\Sigma_i\,\bar\phi_i$"),
          ("c", "cube", "Blues", (0.4, 1), r"(c) nutrient $c$", r"$c$"),
          ("a", "ab", "Reds", (0, 1), r"(d) antibiotic $a$", r"$a$"),
          ("vm", "stress", "magma_r", (0, VMAX), r"(e) von Mises stress", r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")]
B_WHERE = r"four small blocks at the edges along $z$ ($|x|, |y| > 0.75$ mm)"
if CASE.startswith("separate"):
    START = rf"Separate start: group A in the centre seed, group B in {B_WHERE}. "
else:
    START = r"Mixed start: all species in equal amounts together in the centre seed. "
START += (r"New biofilm takes the living composition of its neighbours and starts alive." if INHLIV else
          r"New biofilm takes the composition of its neighbours." if INHERIT else
          (r"New biofilm starts as B where $|x|, |y| > 0.56$ mm, else as A (rule of the ANSYS runs)."
           if CASE.startswith("separate") else r"New biofilm starts from the mixed start state (rule of the ANSYS runs)."))
MODEL = (r"Klempt et al. 2026, Table 2, case 4 (four species; group A = species 1, 2, group B = species 3, 4; antibiotic "
         r"sensitivities $b = 10, 2, 1, 0.01$)" if K26C4 else
         r"Fritsch et al. 2025, Case I, true parameters (two species; A = species 1, B = species 2)")
AB = (r"in pulses: $a = 2$ for the first half of each period $2/9$, 0 for the second half, the same dose as $a = 1$"
      if PULSE else r"held at 1")
CAPTION = "\n".join([
    r"Coupled growth model: Klempt et al. 2024, Eqs. 34$-$36, front term on ($r = 10$ mm/$T^*$), with the point model of "
    r"Klempt et al. 2026 at every point; 2 mm cube, $64^3$ finite volumes, TVD front term.",
    r"Point model: " + MODEL + ". " + START,
    r"(a) Share of group A (blue: A, red: B) and (b) living volume fraction $\Sigma_i\,\bar\phi_i$, $\bar\phi_i=\phi_i\psi_i$, "
    r"both where the biofilm $\phi_{3D} > 0.5$. (c) Nutrient $c$, held at 1 in the bottom layer $y < -0.75$ mm, "
    r"consumed by the living biomass.",
    r"(d) Antibiotic $a$ on the outer shell from $T^* = 1/3$ ($\alpha^*_0 = 100$), " + AB + r". "
    r"One octant ($x, y, z > 0$) is removed to show the inside.",
    r"(e) Stress. Growth: $\mathbf{F}_g = \alpha\mathbf{I}$, $\dot\alpha = k_\alpha\phi_{3D}$ (Eq. 36), small strain: eigenstrain "
    r"$(\alpha - 1)\mathbf{I}$. Constitutive law: $\sigma = E(\phi_{3D})\,\mathbb{C}(\nu) : (\varepsilon - (\alpha - 1)\mathbf{I})$, "
    r"$E(\phi_{3D}) = 10$ Pa $(\phi_{3D}^2 + 10^{-3})$, $\nu = 0.49$. Equilibrium: $\nabla\cdot\sigma = \mathbf{0}$, free cube "
    r"(rigid-body motion constrained), $16^3$ hexahedra. Shown: $\sigma_\mathrm{vM} = (\frac{3}{2}\,s : s)^{1/2}$."
    + (r" Grey: four rigid columns (Klempt et al. 2024, Sect. 4.3), no biofilm and no flux in them." if OBST else ""),
    r"Constants: $\beta = 0.02$ mm$^2$/$T^*$, $r = 10$ mm/$T^*$, $K = 1$, $k_\alpha = 10^{-3}$/$T^*$, "
    r"$d = 10^8$ mm$^2$/$T^*$ and $g\,l^2/d = 1$ for $c$ and $a$, $\phi_{\min} = 10^{-3}$ (Klempt et al. 2024, Table 2, "
    r"converted to the 2 mm cube); $T^* = 1$ is the duration of the point model's data."])


def cube_grid(m):
    g = pv.ImageData(dimensions=(m + 1,) * 3, spacing=(2.0 / m,) * 3, origin=(-1.0,) * 3)
    cc = g.cell_centers().points
    octant = (cc[:, 0] > 0) & (cc[:, 1] > 0) & (cc[:, 2] > 0)
    solid = np.zeros(len(cc), bool)
    if OBST:
        for lo, hi in BOXES:
            solid |= np.all((cc > np.array(lo)) & (cc < np.array(hi)), axis=1)
    return g, octant, solid


def grid(T):
    g, octant, solid = cube_grid(n)
    for k in ("phi", "shareA", "L", "c", "a"):
        g.cell_data[k] = np.nan_to_num(Z[f"{k}_{T}"], nan=0.0).ravel(order="F")
    g.cell_data["bio"] = ((g.cell_data["phi"] > 0.5) & ~octant).astype(float)
    g.cell_data["cube"] = (~octant & ~solid).astype(float)
    g.cell_data["living"] = (g.cell_data["bio"] * (g.cell_data["L"] > 0.01)).astype(float)
    g.cell_data["ab"] = (g.cell_data["cube"] * (g.cell_data["a"] > 0.01)).astype(float)
    s, octant_s, solid_s = cube_grid(nm)
    s.cell_data["vm"] = 1e4 * Z[f"vm16_{T}"].ravel(order="F")
    s.cell_data["stress"] = (~octant_s & ~solid_s).astype(float)
    return g, s


def view(mesh, field, cmap, clim):
    pl = pv.Plotter(off_screen=True, window_size=(900, 820))
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
    if mesh.n_cells:
        pl.add_mesh(mesh, scalars=field, cmap=cmap, clim=clim, show_scalar_bar=False, **LIT)
    if OBST:
        for lo, hi in BOXES:
            pl.add_mesh(pv.Box((lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])), color=(0.55, 0.55, 0.55), **LIT)
    pl.camera_position = CAM
    pl.camera.zoom(0.98)
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


fig = plt.figure(figsize=(20, 6.8), dpi=160)                  # 3200 x 1088 px (multiples of 16)
W = 0.2
ax_img = [fig.add_axes([0.004 + W * k, 0.40, W - 0.008, 0.545]) for k in range(5)]
ax_bar = [fig.add_axes([0.03 + W * k, 0.37, W - 0.06, 0.018]) for k in range(5)]
name = (f"anim_coupled_tvd_{CASE}{'_inhliv' if INHLIV else ''}{'_inherit' if INHERIT else ''}{'_obst' if OBST else ''}{'_k26c4' if K26C4 else ''}"
        f"{'_pulse' if PULSE else ''}_5panel")     # the earlier four-panel videos keep their names
out = ROOT / "assets" / f"{name}.mp4"
movie = imageio.get_writer(str(out), fps=8, quality=9)
for T in times:
    g, s = grid(T)
    meshes = {m: g.threshold(0.5, scalars=m) for m in ("bio", "living", "cube", "ab")}
    meshes["stress"] = s.threshold(0.5, scalars="stress")
    for k, (field, mesh, cmap, clim, title, label) in enumerate(PANELS):
        ax, cb = ax_img[k], ax_bar[k]
        ax.clear(); cb.clear()
        ax.imshow(view(meshes[mesh], field, cmap, clim))
        ax.set_axis_off()
        ax.set_title(title, fontsize=13, loc="left")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*clim))
        c = fig.colorbar(sm, cax=cb, orientation="horizontal")
        c.set_ticks(np.linspace(clim[0], clim[1], 3))
        if field == "vm":
            c.ax.set_xticklabels([f"{v:g}" for v in np.linspace(clim[0], clim[1], 3)])
        c.ax.tick_params(labelsize=10)
        cb.set_xlabel(label, fontsize=12, labelpad=2)
    for t in list(fig.texts):
        t.remove()
    fig.text(0.996, 0.975, rf"$T^* = {T:.3f}$", ha="right", va="top", fontsize=15)
    fig.text(0.008, 0.29, CAPTION, ha="left", va="top", fontsize=11.5, linespacing=1.5)
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    for _ in range(9 if T == times[-1] else 1):
        movie.append_data(img)
    if T == times[-1]:
        imageio.imwrite(out.with_name(out.stem + "_end.png"), img)
movie.close()
print("wrote", out)

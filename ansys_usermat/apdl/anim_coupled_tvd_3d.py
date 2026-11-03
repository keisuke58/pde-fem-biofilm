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
  k26c4:   the four species of Klempt 2026 Table 2, case 4; pulse: the antibiotic in pulses of the same dose (A-71);
  transc:  the runs with dc/dt in Eq. 35 (--transient-c 1e8), as the thesis figures; use with inhliv.
  fine:    coupling step 0.0025 and stress on 32^3 (A-90; run with --dtc 0.0025 --nm 32).
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
ALPHA = None                      # set after loading: alpha - 1 replaces the antibiotic panel whenever alpha is saved at every snapshot (the standard from 2 Nov; "noalpha" keeps a)
INHLIV = "inhliv" in OPT          # the base rule from 25 Oct: the neighbours' living composition, fresh psi (A-82)
BOXES = [((-1.0, -0.125, lo), (1.0, 0.125, hi)) for lo, hi in ((-0.75, -0.5), (0.5, 0.75))] + \
        [((-0.125, lo, -1.0), (0.125, hi, 1.0)) for lo, hi in ((-0.75, -0.5), (0.5, 0.75))]
MIX = "_mixeq" if CASE.startswith("mixed") else ""
# asymmetric runs (2-3 Nov, full domain with --amg): b1edge = B only in the edge block x > 0.75, y < -0.75 mm; bcorner = B only
# in the corner cube x > 0.75, y < -0.75, z > 0.75 mm; seedx-1 = the seed moved by one element (0.25 mm) towards -x
B1EDGE, BCORNER = ("b1edge" in OPT), ("bcorner" in OPT)
SEEDX = next((o[5:] for o in OPT if o.startswith("seedx")), "")
ASYM = ("_b1edge" if B1EDGE else "") + ("_bcorner" if BCORNER else "") + (f"_seedx{SEEDX}" if SEEDX else "")
FULL = bool(ASYM)
TRANSC = "transc" in OPT          # Eq. 35 with dc/dt (--transient-c 1e8), the runs of the thesis figures (29 Oct)
FINE = "fine" in OPT              # A-90: coupling step 0.0025 (A-89) and the stress on a 32^3 mesh
TAG = (f"_nf64_front1{'_dtc0.0025' if FINE else ''}_phimin0.001_v2_fortran_ansysbc{'_inhliv' if INHLIV else ('' if INHERIT else '_startpos')}{MIX}_a100_gab1_gc1_tab0.333333"
       f"_snapsvm{'_nm32' if FINE else ''}{'_transc1e+08' if TRANSC else ''}{'_k26c4' if K26C4 else ''}{'_abpulse0.2222' if PULSE else ''}{ASYM}_tvd_wgmirror"
       f"{'' if FULL else '_fast'}{'_obstcolumns' if OBST else ''}{'_amg' if FULL else ''}")
Z = np.load(ROOT / ("figdata") / f"coupled_living_{CASE}{TAG}_{CASE}.npz")
times = sorted({float(k.split("_", 1)[1]) for k in Z.files if k.startswith("phi_")})
ALPHA = ("noalpha" not in OPT) and all(f"alpha_{t}" in Z.files for t in times)
DOM = K26C4 and all(f"sp0_{t}" in Z.files for t in times)   # 2 Nov: panel (a) = dominant species of four, coloured
SPC = np.array([[59, 91, 219], [116, 192, 252], [201, 42, 42], [245, 159, 0]], float)   # species 1..4
n = Z[f"phi_{times[0]}"].shape[0]
nm = Z[f"vm16_{times[0]}"].shape[0]
VMAX = 1e4 * max(np.nanpercentile(Z[f"vm16_{T}"], 99.5) for T in times)   # nan: the elements of the rigid bars   # one colour range over the video, 1e-4 Pa
VMAX = float(2 * np.ceil(VMAX / 2))                      # even, so the middle tick is a whole number
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)   # light shading, colours close to the colour bar
CAM = [(3.6, 2.6, 4.4), (0, 0, 0), (0, 1, 0)]
# field, mesh, colour map, range, title, colour-bar label; mesh: "bio" phi_3D > 0.5, "living" bio and sum phibar > 0.01,
# "cube" all, "ab" a > 0.01, "stress" the stress mesh (values near 0 are left out, i.e. transparent)
GROUP = "species 1, 2" if K26C4 else "species 1"
GRP = "group" if K26C4 else "species"     # two species: the notation of the thesis and the slides, phi_1/(phi_1+phi_2)
PANELS = [("shareA", "bio", "coolwarm_r", (0, 1), rf"(a) share of {GRP} A",
           r"$(\phi_1+\phi_2)/\Sigma_i\phi_i$" if K26C4 else r"$\phi_1/(\phi_1+\phi_2)$"),
          ("L", "living", "Greens", (0, 1), r"(b) living volume fraction", r"$\Sigma_i\,\bar\phi_i$"),
          ("c", "cube", "Blues", (0.4, 1), r"(c) nutrient $c$", r"$c$"),
          ("a", "ab", "Reds", (0, 1), r"(d) antibiotic $a$", r"$a$"),
          ("vm", "stress", "magma_r", (0, VMAX), r"(e) von Mises stress", r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")]
if DOM:
    PANELS[0] = ("dom", "bio", None, None, r"(a) dominant species", r"darker: larger share $\phi_i/\Sigma_j\phi_j$")
if ALPHA:
    AMAX = float(np.ceil(1e4 * np.percentile(Z[f"alpha_{times[-1]}"] - 1.0, 99.5)))
    PANELS[3] = ("am1", "bio", "viridis", (0, AMAX), r"(d) growth $\alpha - 1$", r"$\alpha - 1$ [$10^{-4}$]")   # replaces a
B_WHERE = (r"one block along the edge $x > 0.75$, $y < -0.75$ mm" if B1EDGE else
           r"one cube of 0.25 mm at the corner $x > 0.75$, $y < -0.75$, $z > 0.75$ mm" if BCORNER else
           r"four small blocks at the edges along $z$ ($|x|, |y| > 0.75$ mm)")
if SEEDX:
    B_WHERE += rf" (seed moved by {abs(int(SEEDX)) * 0.25:g} mm towards {'-' if int(SEEDX) < 0 else '+'}$x$)"
if CASE.startswith("separate"):
    START = rf"Separate start: {GRP} A in the centre seed, {GRP} B in {B_WHERE}. "
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
    r"Klempt et al. 2026 at every point; 2 mm cube, $64^3$ finite volumes, TVD front term"
    + (r", coupling step $\Delta T^* = 0.0025$." if FINE else "."),
    *((r"Point model: " + MODEL + ".", START) if (K26C4 or FULL) else (r"Point model: " + MODEL + ". " + START,)),
    (r"(a) Dominant species (colour), darker where its share is larger," if DOM else rf"(a) Share of {GRP} A (blue: A, red: B)") + rf" and (b) living volume fraction $\Sigma_i\,\bar\phi_i$, $\bar\phi_i=\phi_i\psi_i$, "
    r"both where the biofilm $\phi_{3D} > 0.5$. (c) Nutrient $c$, held at 1 in the layer on $y = -1$ mm ($y < -0.75$ mm), "
    r"consumed by the living biomass: $\dot c = d\,\nabla^2 c - g\,\phi_{3D}(\Sigma_i\,\bar\phi_i)\,c$ (Eq. 35, backward Euler).",
    r"Inside that layer $\nabla c = 0$, so the front term does not act and biofilm enters it by diffusion only (the gap at the bottom).",
    (r"(d) Growth $\alpha - 1$ where $\phi_{3D} > 0.5$. Antibiotic $a$ (not shown) on the outer shell from $T^* = 1/3$ ($\alpha^*_0 = 100$), "
     if ALPHA else r"(d) Antibiotic $a$ on the outer shell from $T^* = 1/3$ ($\alpha^*_0 = 100$), ") + AB + r". "
    r"One octant ($x, y, z > 0$) is removed to show the inside.",
    r"(e) Stress. Growth: $\mathbf{F}_g = \alpha\mathbf{I}$, $\dot\alpha = k_\alpha\phi_{3D}$ (Eq. 36), small strain: eigenstrain "
    r"$(\alpha - 1)\mathbf{I}$. Constitutive law: $\sigma = E(\phi_{3D})\,\mathbb{C}(\nu) : (\varepsilon - (\alpha - 1)\mathbf{I})$, "
    r"$E(\phi_{3D}) = 10$ Pa $(\phi_{3D}^2 + 10^{-3})$, $\nu = 0.49$. Equilibrium: $\nabla\cdot\sigma = \mathbf{0}$, "
    + (r"frictionless contact with the rigid columns (no other supports), $" if OBST else r"free cube (rigid-body motion constrained), $")
    + str(nm) + r"^3$ hexahedra. Shown: $\sigma_\mathrm{vM} = (\frac{3}{2}\,s : s)^{1/2}$."
    + ("\n" + r"Growth follows all biomass $\phi_{3D}$, dead included, so the stress keeps rising after the antibiotic."
       if FINE else "")
    + ("\n" + r"Grey: four rigid columns (Klempt et al. 2024, Sect. 4.3), no biofilm and no flux in them." if OBST else ""),
    r"Constants: $\beta = 0.02$ mm$^2$/$T^*$, $r = 10$ mm/$T^*$, $k = 1$, $k_\alpha = 10^{-3}$/$T^*$, "
    r"$d = 10^8$ mm$^2$/$T^*$ and $g\,\ell^2/d = 1$ for $c$ and $a$, $\phi_{\min} = 10^{-3}$ (Klempt et al. 2024, Table 2, "
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
    if DOM:
        sh = np.stack([np.nan_to_num(Z[f"sp{i}_{T}"]).ravel(order="F") for i in range(4)], 1)
        top = sh.argmax(1); w = np.clip((sh.max(1) - 0.25) / 0.75, 0, 1)[:, None]   # 0.25 (all equal) .. 1 (one species)
        g.cell_data["dom"] = (255 - w * (255 - SPC[top])).astype(np.uint8)
    if ALPHA:
        g.cell_data["am1"] = (1e4 * (Z[f"alpha_{T}"] - 1.0)).ravel(order="F")
    g.cell_data["bio"] = ((g.cell_data["phi"] > 0.5) & ~octant).astype(float)
    g.cell_data["cube"] = (~octant & ~solid).astype(float)
    g.cell_data["living"] = (g.cell_data["bio"] * (g.cell_data["L"] > 0.01)).astype(float)
    g.cell_data["ab"] = (g.cell_data["cube"] * (g.cell_data["a"] > 0.01)).astype(float)
    s, octant_s, solid_s = cube_grid(nm)
    s.cell_data["vm"] = 1e4 * np.nan_to_num(Z[f"vm16_{T}"]).ravel(order="F")
    s.cell_data["stress"] = (~octant_s & ~solid_s).astype(float)
    return g, s


def view(mesh, field, cmap, clim):
    pl = pv.Plotter(off_screen=True, window_size=(900, 820))
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).outline(), color="grey")
    if mesh.n_cells:
        if field == "dom":
            pl.add_mesh(mesh, scalars=field, rgb=True, show_scalar_bar=False, **LIT)
        else:
            pl.add_mesh(mesh, scalars=field, cmap=cmap, clim=clim, show_scalar_bar=False, **LIT)
    if OBST:
        for lo, hi in BOXES:
            pl.add_mesh(pv.Box((lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])), color=(0.55, 0.55, 0.55), **LIT)
    pl.camera_position = CAM
    pl.camera.zoom(0.98)
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


NP = len(PANELS)
HIST = "hist" in OPT               # 2 Nov (option): a small time-history plot of the seed means at the right, drawn up to the frame
fig = plt.figure(figsize=((24 if NP == 6 else 20) + (5 if HIST else 0), 6.8), dpi=160)    # 3200 x 1088 px (4000 with hist)
W = (0.8 if HIST else 1.0) / NP
ax_img = [fig.add_axes([0.004 + W * k, 0.40, W - 0.008, 0.545]) for k in range(NP)]
# coordinate axes as seen from the panels' camera, once, in the free lower right corner (y up)
_tx = fig.add_axes([0.905, 0.01, 0.09, 0.2])
_d = -np.asarray(CAM[0], float); _d /= np.linalg.norm(_d)
_r = np.cross(_d, CAM[2]); _r /= np.linalg.norm(_r); _u = np.cross(_r, _d)
for _v, _lab, _col in (((1, 0, 0), "x", "#d9480f"), ((0, 1, 0), "y", "#2b8a3e"), ((0, 0, 1), "z", "#1c3faa")):
    _e = np.array(_v, float)
    _tx.annotate("", xy=(_e @ _r, _e @ _u), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=_col, lw=2, mutation_scale=16))
    _tx.text(1.3 * (_e @ _r), 1.3 * (_e @ _u), f"${_lab}$", color=_col, fontsize=16, ha="center", va="center")
_tx.set_xlim(-1.05, 1.05); _tx.set_ylim(-0.6, 1.25); _tx.set_aspect("equal"); _tx.set_axis_off()
ax_bar = [fig.add_axes([0.03 + W * k, 0.37, W - 0.06 * (0.8 if HIST else 1.0), 0.018]) for k in range(NP)]
if HIST:                                  # seed means over all frames (from the snapshots), normalised where they are not fractions
    import wp3_reference as _W
    _c8 = _W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
    _s64, _s16 = _W.seed_mask(n, _c8), _W.seed_mask(nm, _c8)
    _mean = lambda k, T, m: float(np.nanmean(np.where(m, np.nan_to_num(Z[f"{k}_{T}"]), np.nan)))
    HS = {"living": [_mean("L", T, _s64) for T in times], "c": [_mean("c", T, _s64) for T in times],
          "share": [_mean("shareA", T, _s64) for T in times], "vm": [_mean("vm16", T, _s16) for T in times]}
    if ALPHA:
        HS["am1"] = [_mean("alpha", T, _s64) - 1.0 for T in times]
    ax_h = fig.add_axes([0.835, 0.47, 0.15, 0.38])
name = (f"anim_coupled_tvd_{CASE}{'_inhliv' if INHLIV else ''}{'_inherit' if INHERIT else ''}{'_obst' if OBST else ''}{'_k26c4' if K26C4 else ''}"
        f"{'_pulse' if PULSE else ''}{ASYM}{'_fine' if FINE else ''}{'_transc' if TRANSC else ''}{'_alpha' if ALPHA else ''}{'_hist' if HIST else ''}_{NP}panel")     # the earlier four-panel videos keep their names
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
        if field == "dom":                 # four colour boxes instead of a colour bar
            for i in range(4):
                cb.add_patch(plt.Rectangle((i, 0), 0.92, 1, color=SPC[i] / 255))
                cb.text(i + 0.46, 0.5, f"{i + 1}", ha="center", va="center", fontsize=10, color="white")
            cb.set_xlim(0, 4); cb.set_ylim(0, 1); cb.set_xticks([]); cb.set_yticks([])
            for sp_ in cb.spines.values():
                sp_.set_visible(False)
            cb.set_xlabel(f"({'abcde'[k]})  species  (" + label + ")", fontsize=12, labelpad=2)
            continue
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*clim))
        c = fig.colorbar(sm, cax=cb, orientation="horizontal")
        c.set_ticks(np.linspace(clim[0], clim[1], 3))
        if field == "vm":
            c.ax.set_xticklabels([f"{v:g}" for v in np.linspace(clim[0], clim[1], 3)])
        c.ax.tick_params(labelsize=10)
        cb.set_xlabel(f"({'abcde'[k]})  " + label, fontsize=12, labelpad=2)   # panel letter with the colour bar (2 Nov)
    for t in list(fig.texts):
        t.remove()
    fig.text(0.996, 0.975, rf"$T^* = {T:.3f}$", ha="right", va="top", fontsize=15)
    if CASE.endswith("_ab") and T >= 1 / 3 - 1e-9:   # 2 Nov: mark the antibiotic onset (its panel is replaced by alpha - 1)
        fig.text(0.996, 0.935, (r"antibiotic in pulses from $T^* = 1/3$" if PULSE else r"antibiotic $a = 1$ on the outer shell from $T^* = 1/3$"),
                 ha="right", va="top", fontsize=12, color="#a61e1e")
    fig.text(0.008, 0.29, CAPTION, ha="left", va="top", fontsize=10.5 if K26C4 else 11.5, linespacing=1.4 if K26C4 else 1.5)   # four species: one line more
    if HIST:
        ax_h.clear(); kk = times.index(T) + 1
        for key, lab, col, scale in (("living", r"$\Sigma_i\bar\phi_i$", "#2b8a3e", 1.0), ("share", "share of A", "#495057", 1.0), ("c", "$c$", "#1c7ed6", 1.0),
                                     ("am1", r"$(\alpha-1)/\max$", "#7950f2", None), ("vm", r"$\sigma_\mathrm{vM}/\max$", "#c2255c", None)):
            if key not in HS:
                continue
            v = np.array(HS[key]); v = v / (np.nanmax(v) if scale is None else scale)
            ax_h.plot(times[:kk], v[:kk], color=col, lw=1.8, label=lab)
        ax_h.axvline(1 / 3, color="0.5", ls=":", lw=0.8)
        ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1.05); ax_h.set_xlabel("$T^*$", fontsize=12); ax_h.tick_params(labelsize=10)
        ax_h.set_title("seed means", fontsize=13, loc="left"); ax_h.legend(frameon=False, fontsize=9.5, loc="upper left", bbox_to_anchor=(0, -0.22), ncol=2)
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    for _ in range(9 if T == times[-1] else 1):
        movie.append_data(img)
    if T == times[-1]:
        imageio.imwrite(out.with_name(out.stem + "_end.png"), img)
movie.close()
print("wrote", out)

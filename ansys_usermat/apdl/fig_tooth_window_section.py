"""fig_tooth_window_section.py -- quick check figure of the cervical window (papers/coupling_paper/REAL_TOOTH_WINDOW.ja.md):
the section y = 0 of a window run (coupled_living_consumption.py --obstacle stl --window mixed|layered --stl-scale 100
--box 30 --only window_ab) and of its film stress (real_tooth_coupled.py --window). (a) the start film on the crown
(species A, B of the layered start), (b) phi_3D at T* = 1 (and the film at the antibiotic onset as a contour),
(c) von Mises stress of the film at T* = 1 (h = h_f). Lengths in um (1 code mm = 10 um). Times (figstyle), no title.

    python ansys_usermat/apdl/fig_tooth_window_section.py --nf 64 --window mixed [--case window_ab] [--opts <tag>]
        -> assets/fig_tooth_window_section_<window>_nf<NF>[_noab].png
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(size=8)


def _opt(name, default, cast=float):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


NF = _opt("--nf", 64, int)
WINDOW = _opt("--window", "mixed", str)
CASE = _opt("--case", "window_ab", str)
BOX = _opt("--box", 30.0)
OPTS = _opt("--opts", "front1_phimin0.001_v2_fortran_inhliv_a100_gab1_gc1_tab0.333333_transc1e+08_tvd_superbee_wgmirror"
            f"_obststl_win{WINDOW}_stlx100_box30_amg", str)
UM = 10.0                                            # um per code mm
Z = np.load(ROOT / "figdata" / f"coupled_living_{CASE}_nf{NF}_{OPTS}_{CASE}.npz")
times = sorted({float(f.split("_", 1)[1]) for f in Z.files if f.startswith("phi_")})
T_END = times[-1]
obst = Z[f"obst_{T_END}"].astype(bool)
seed, edges = Z[f"seed_{T_END}"].astype(bool), Z[f"edges_{T_END}"].astype(bool)
H = BOX / NF
j = NF // 2                                          # the cell layer just above y = 0
e = (-BOX / 2 + np.arange(NF + 1) * H) * UM          # cell edges in um
st = np.load(ROOT / "figdata" / f"real_tooth_coupled_stress_{CASE}_nf{NF}_T{T_END}_h1.npz")
vm = np.where(st["active"], np.nan_to_num(st["vm"].astype(float)), np.nan) * 1e3   # mPa
o, hs, bx = st["origin"], float(st["h"]), st["box"]
jv = int(round((0.0 - o[1]) / hs))                   # stress elements above y = 0
ev = [(o[k] + np.arange(vm.shape[k] + 1) * hs) * UM for k in range(3)]

fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.45), constrained_layout=True)
grey = ListedColormap(["white", "#c9c4b8"])
for ax in axs:
    ax.pcolormesh(e, e, obst[:, j, :].T, cmap=grey, vmin=0, vmax=1, shading="flat", rasterized=True)
    ax.set_aspect("equal"); ax.grid(False)
    ax.set_xlim(e[0], e[-1]); ax.set_ylim(e[0], e[-1])
    ax.set_xlabel(r"$x$ [$\mu$m]")
    ax.tick_params(labelsize=7.5)
axs[0].set_ylabel(r"$z$ [$\mu$m]")
# (a) start film
ax = axs[0]
film0 = np.ma.masked_where(~(seed | edges)[:, j, :], np.where(edges[:, j, :], 2.0, 1.0))
ax.pcolormesh(e, e, film0.T, cmap=ListedColormap(["#3b6fb6", "#c0392b"]), vmin=1, vmax=2, shading="flat", rasterized=True)
ax.set_title("(a) start, 20 $\\mu$m: " + ("A inside, B outside" if edges.any() else "mixed"), fontsize=8, loc="left")
ax.text(-140, 100, "saliva", fontsize=7.5, color="0.3")
ax.text(60, 100, "crown", fontsize=7.5, color="0.3")
# (b) phi at T_END with the onset film as a contour
ax = axs[1]
phi = Z[f"phi_{T_END}"][:, j, :]
im = ax.pcolormesh(e, e, np.ma.masked_where(obst[:, j, :], phi).T, cmap="Greens", vmin=0, vmax=1, shading="flat", rasterized=True)
if len(times) > 1:
    phi_on = Z[f"phi_{times[0]}"][:, j, :]
    cc = (e[:-1] + e[1:]) / 2
    ax.contour(cc, cc, phi_on.T, levels=[0.5], colors=["k"], linewidths=0.6, linestyles="--")
    ax.plot([], [], "k--", lw=0.6, label=f"$T^*={times[0]:g}$")
    ax.legend(loc="upper left", fontsize=6.5, frameon=False, handlelength=1.4, borderaxespad=0.2)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label=r"$\phi_{3D}$")
ax.set_title(f"(b) $\\phi_{{3D}}$ at $T^* = {T_END:g}$", fontsize=8, loc="left")
# (c) von Mises at T_END
ax = axs[2]
vmax = float(np.ceil(np.nanpercentile(vm, 99)))
im = ax.pcolormesh(ev[0], ev[2], vm[:, jv, :].T, cmap="magma_r", vmin=0, vmax=vmax, shading="flat", rasterized=True)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label=r"$\sigma_\mathrm{vM}$ [mPa]")
ax.set_title(f"(c) von Mises at $T^* = {T_END:g}$", fontsize=8, loc="left")
out = ROOT / "assets" / f"fig_tooth_window_section_{WINDOW}_nf{NF}{'' if CASE.endswith('_ab') else '_noab'}.png"
fig.savefig(out, dpi=300)
print("wrote", out)

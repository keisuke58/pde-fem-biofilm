"""anim_cal5_ch.py -- colloquium animation of the Chapter 4 run with the calibrated five-species point model
(condition CH, commensal under flow, final MAP of the measured-data calibration of Chapter 3, Day-1 start), in the
element with the one-way coupling of the chapter (Eq. 36), 16^3 elements of the 2 mm cube, ANSYS.

Data (IKMHIWI03, A-18 of 20 Oct 2026):
  results/2026-10-cal5/fields_w16_CH_g1_s015_final_d1_dense.json   phi, alpha - 1, von Mises (MPa) per element at
                                                                  T* = 0.025, 0.05, ..., 1.0
  results/2026-10-cal5/w16_CH_g1_s015_final_d1_pm.json            point-model state at a seed Gauss point
Panels: the mid-plane z = 0 through the seed (mean of the two element layers next to it): phi, alpha - 1, von Mises;
below, the point model at a seed Gauss point: shares phi_i / sum phi_j and living volume fractions phi_i psi_i,
drawn up to the current time. Fixed colour scales (the maximum over all frames).

    python ansys_usermat/apdl/anim_cal5_ch.py   -> assets/anim_cal5_ch.gif, .mp4, .png (last frame), .pdf
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
from fig_cal5_appendix import NAMES, trace  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter  # noqa: E402

RES = HERE / "results" / "2026-10-cal5"
F = json.load(open(RES / "fields_w16_CH_g1_s015_final_d1_dense.json"))
N = 16
t = np.array(F["t"])
cx, cy, cz = (np.array(F[k]) for k in ("cx", "cy", "cz"))
ix = np.rint((cx + 1) * N / 2 - 0.5).astype(int)
iy = np.rint((cy + 1) * N / 2 - 0.5).astype(int)
iz = np.rint((cz + 1) * N / 2 - 0.5).astype(int)


def grid(v):
    g = np.zeros((N, N, N))
    g[ix, iy, iz] = v
    return g


def mid(g):
    return 0.5 * (g[:, :, N // 2 - 1] + g[:, :, N // 2])


PHI = [mid(grid(np.array(v))) for v in F["phi"]]
AM1 = [mid(grid(np.array(v))) * 1e3 for v in F["alpha"]]          # 1e-3
VM = [mid(grid(np.array(v))) * 1e9 for v in F["seqv"]]            # MPa -> mPa
tp, phi, psi = trace("w16_CH_g1_s015_final_d1_pm.json")
assert len(tp) == len(t) and np.allclose(tp, t), (len(tp), len(t))
sh = phi / phi.sum(1, keepdims=True)
liv = phi * psi

fig = plt.figure(figsize=(7.2, 4.9))
gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.82], hspace=0.42, wspace=0.9)
axm = [fig.add_subplot(gs[0, 2 * j:2 * j + 2]) for j in range(3)]
axl = [fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]
e = np.linspace(-1, 1, N + 1)
specs = [(PHI, "(a) $\\phi$", "Greens", 1.0), (AM1, "(b) $\\alpha-1$  ($10^{-3}$)", "Oranges", max(a.max() for a in AM1)),
         (VM, "(c) von Mises (mPa)", "viridis", max(v.max() for v in VM))]
ims = []
for a, (data, lab, cm, top) in zip(axm, specs):
    im = a.pcolormesh(e, e, data[0].T, cmap=cm, vmin=0, vmax=top, shading="flat", rasterized=True)
    cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=7.5)
    a.set_title(lab, fontsize=10, pad=3)
    a.set_aspect("equal"); a.set_xticks([-1, 0, 1]); a.set_yticks([-1, 0, 1])
    a.tick_params(labelsize=8); a.grid(False)
    a.set_xlabel("$x$ (mm)", fontsize=9)
    ims.append((im, data))
axm[0].set_ylabel("$y$ (mm)", fontsize=9)

cols = plt.rcParams["axes.prop_cycle"].by_key()["color"]
lines = []
for i in range(5):
    l0, = axl[0].plot([], [], lw=1.2, color=cols[i], label=NAMES[i])
    l1, = axl[1].plot([], [], lw=1.2, color=cols[i])
    lines += [(l0, sh[:, i]), (l1, liv[:, i])]
ls, = axl[1].plot([], [], lw=1.5, color="#222", label="$\\Sigma_i\\phi_i\\psi_i$")
lines.append((ls, liv.sum(1)))
axl[0].set(xlabel="$T^*$", ylabel="$\\phi_i/\\Sigma_j\\phi_j$", xlim=(0, 1), ylim=(0, 0.85))
axl[1].set(xlabel="$T^*$", ylabel="$\\phi_i\\psi_i$", xlim=(0, 1), ylim=(0, 0.9))
axl[0].set_title("(d) composition at a seed point", fontsize=10, pad=3)
axl[1].set_title("(e) living volume fractions", fontsize=10, pad=3)
axl[0].legend(ncol=2, fontsize=7.5, frameon=False, loc="center right", handlelength=1.2, columnspacing=0.8)
axl[1].legend(fontsize=7.5, frameon=False, loc="upper right", handlelength=1.2)
for a in axl:
    a.tick_params(labelsize=8)
txt = fig.text(0.5, 0.012, "", ha="center", fontsize=8.5)
fig.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.135)


def draw(i):
    i = min(i, len(t) - 1)
    for im, data in ims:
        im.set_array(data[i].T.ravel())
    for ln, y in lines:
        ln.set_data(t[:i + 1], y[:i + 1])
    txt.set_text(f"$T^*$ = {t[i]:.3f}   calibrated five-species point model (commensal, flow), one-way coupling; "
                 f"mid-plane $z$ = 0 of the 2 mm cube, seed in the centre")
    return [im for im, _ in ims] + [ln for ln, _ in lines] + [txt]


nfr = len(t) + 12                                                   # hold the last frame for 1.5 s
anim = FuncAnimation(fig, draw, frames=nfr)
out = ROOT / "assets" / "anim_cal5_ch"
anim.save(out.with_suffix(".gif"), writer=PillowWriter(fps=8), dpi=200)
try:
    import imageio_ffmpeg
    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    anim.save(out.with_suffix(".mp4"), writer=FFMpegWriter(fps=8, bitrate=4000), dpi=250)
except Exception as exc:                                            # no ffmpeg: the gif is enough
    print("mp4 skipped:", exc)
draw(len(t) - 1)
fig.savefig(out.with_suffix(".png"), dpi=400); fig.savefig(out.with_suffix(".pdf"))
s = np.isin(np.array(F["elem"]), json.load(open(RES / "w16_CH_g1_s015_final_d1.json"))["seed_BIOFILM1"])
print("seed at T* = 1: alpha-1 %.4e  vM %.4e Pa" % (np.mean(np.array(F["alpha"][-1])[s]), np.mean(np.array(F["seqv"][-1])[s]) * 1e6))
print("shares end", dict(zip(NAMES, sh[-1].round(3))), "living", liv.sum(1)[[0, -1]].round(3))
print("wrote", out.with_suffix(".gif"))

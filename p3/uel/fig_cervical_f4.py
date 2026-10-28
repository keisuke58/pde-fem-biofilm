"""fig_cervical_f4.py -- still figure for the synthetic cervical runs (anim_cervical_3d.py; the four-species case of
Fritsch et al. 2025, no experimental data): (a) composition of the seed over time (phibar_i = phi_i psi_i of species
1-4, as shares of their sum, element means of the point model over the seed layers), (b), (c) von Mises / E on the
section theta = 45 degrees (r-z plane, deformed, true scale) at the end, without and with the homeostatic-pressure law
(p_h/E = 0.05), with the largest value marked.

    python p3/uel/fig_cervical_f4.py   -> assets/fig_cervical_f4.pdf, .png
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(9)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

A = np.load(ROOT / "figdata/anim_cervical_f4_0.npz")
B = np.load(ROOT / "figdata/anim_cervical_f4_0.05.npz")
NR = 4
nt = 3 * NR
nodes, elems, layer = A["nodes"], A["elems"], A["layer"]
jk = np.array([((e // NR) % nt, e // (NR * nt)) for e in range(len(elems))])
SEC = jk[:, 0] == nt // 2                      # elements with theta from 45 to 52.5 degrees
seed = layer < NR // 2

fig = plt.figure(figsize=(6.3, 2.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.45)

# (a) composition
ax = fig.add_subplot(gs[0])
g = A["g"][:, seed, :].mean(1)                 # (time, 12)
pb = g[:, 0:4] * g[:, 6:10]
sh = pb / pb.sum(1, keepdims=True)
cols = ["#1e50a0", "#2a9d8f", "#e09f3e", "#c0392b"]
for i in range(4):
    ax.plot(A["t"], sh[:, i], color=cols[i], lw=1.4, label=f"species {i + 1}")
ax.set(xlabel="time $t$", ylabel="share of $\\bar\\phi_i$ in the seed", ylim=(0, 0.6), xlim=(0, A["t"][-1]))
ax.legend(frameon=False, fontsize=7.5, loc="upper left", ncol=2, handlelength=1.4, columnspacing=0.8)
ax.text(-0.3, 1.02, "(a)", transform=ax.transAxes, fontsize=9)

# (b), (c) sections
vmax = max(A["vm"][-1][SEC].max(), B["vm"][-1][SEC].max())
face = [0, 1, 5, 4]                            # i,j,k  i+1,j,k  i+1,j,k+1  i,j,k+1
for k, (D, lab) in enumerate(((A, "(b) without the law"), (B, "(c) with the law"))):
    ax = fig.add_subplot(gs[k + 1])
    x = nodes + D["u"][-1][:, :3]
    rz = np.stack([np.hypot(x[:, 0], x[:, 1]), x[:, 2]], 1)
    polys = [rz[elems[e][face]] for e in np.where(SEC)[0]]
    v = D["vm"][-1][SEC]
    pc = PolyCollection(polys, array=v, cmap="viridis", clim=(0, vmax), edgecolors="0.85", linewidths=0.2)
    ax.add_collection(pc)
    r0 = nodes[:, 0][::NR + 1] * 0 + np.hypot(nodes[::NR + 1, 0], nodes[::NR + 1, 1])
    zc = nodes[::NR + 1, 2]
    o = np.argsort(zc)
    ax.fill_betweenx(zc[o], 1.8, r0[o], color="#c8c3b8", lw=0)          # the collar
    im = int(np.argmax(v))
    c = polys[im].mean(0)
    ax.annotate(f"max {v[im]:.2f}", c, xytext=(2.0, 1.75) if c[1] < 1 else (1.88, 0.15), fontsize=7.5,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.grid(False)
    ax.set(xlim=(1.8, 3.8), ylim=(-0.05, 2.15), xlabel="$r$ (mm)", aspect="equal")
    if k == 0:
        ax.set_ylabel("$z$ (mm)")
    ax.text(0.0, 1.02, lab, transform=ax.transAxes, fontsize=8.5)
cb = fig.colorbar(pc, ax=fig.axes[1:], fraction=0.05, pad=0.08)
cb.set_label("von Mises / $E$")
fig.text(0.99, -0.02, "$p_h/E = 0.05$; deformed, true scale; section at 45$^\\circ$", ha="right", fontsize=7)
out = ROOT / "assets" / "fig_cervical_f4"
fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
print("seed shares end", sh[-1].round(3), "max vM/E", A["vm"][-1][SEC].max(), B["vm"][-1][SEC].max())

"""
plot_at2_thesis.py
==================
Thesis-quality figure for AT2 phase-field biofilm fracture (fully coupled Q4 FEM).

Run on fifa:
    cd ~/IKM_Hiwi/FEM
    PATH=~/texlive/2025/bin/x86_64-linux:$PATH \\
        python3 plot_at2_thesis.py [--save out.pdf]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

sys.path.insert(0, str(Path(__file__).parent))           # FEM root (phase_field_at2_2d_sparse.py)
sys.path.insert(1, str(Path(__file__).parent / "JAXFEM"))  # thesis_style.py
import phase_field_at2_2d_sparse as pf
from thesis_style import use, PALETTE, clean_ax

# ── CLI ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--save", default="phase_field_at2_2d_thesis.pdf")
parser.add_argument("--nx",   type=int, default=120)
parser.add_argument("--nz",   type=int, default=80)
parser.add_argument("--nsteps", type=int, default=120)
args = parser.parse_args()

NX, NZ, NSTEPS = args.nx, args.nz, args.nsteps

# ── Simulation ────────────────────────────────────────────────────────────────
cases = [
    (pf.K_EFF_CH, "z_gradient",   "CH"),
    (pf.K_EFF_DH, "z_gradient",   "DH-base"),
    (pf.K_EFF_DH, "pg_substrate", "DH-Pg"),
]

all_res = {}
for k_base, profile, label in cases:
    t0 = time.time()
    out = pf.run_coupled(k_base, profile, NX, NZ, n_steps=NSTEPS, label=label)
    all_res[label] = out
    print(f"  {label}: {time.time()-t0:.1f} s")

tc = {lbl: pf.t_crit(out["results"]) for lbl, out in all_res.items()}
r1 = tc["DH-base"] / tc["CH"] if tc["DH-base"] and tc["CH"] else float("nan")
r2 = tc["DH-Pg"]   / tc["CH"] if tc["DH-Pg"]   and tc["CH"] else float("nan")

hx   = pf.W / (NX - 1)
hz   = pf.L / (NZ - 1)
x_um = np.linspace(0, pf.W * 1e6, NX)   # (Nx,)  µm
z_um = np.linspace(0, pf.L * 1e6, NZ)   # (Nz,)  µm

# ── Style ─────────────────────────────────────────────────────────────────────
figsize = use(width_frac=1.0, aspect=0.70)

CMAP  = "inferno"                # dark=intact, yellow=fractured (standard for AT2)
NORM  = Normalize(0, 1)
C_CH  = PALETTE["ch"]            # #2196A6  teal
C_DH  = PALETTE["dh"]            # #D95F02  orange
C_PG  = "#8B1A1A"                # dark red  for DH-Pg

fig = plt.figure(figsize=figsize)

# GridSpec ─────────────────────────────────────────────────────────────────────
#   cols: [map0, map1, map2, gap, cbar, gap2, curves(2), gap3, bar]
gs = gridspec.GridSpec(
    3, 7,
    figure=fig,
    height_ratios=[1, 1, 1.0],
    width_ratios=[1, 1, 1, 0.03, 0.04, 1.6, 0.7],
    hspace=0.55, wspace=0.22,
)

# ── Damage maps (rows 0 & 1) ──────────────────────────────────────────────────
map_rows = [
    ("CH",    r"\textbf{CH}",                 C_CH),
    ("DH-Pg", r"\textbf{DH}$+$\textit{Pg}",  C_PG),
]
stages  = ["Early", "Mid", "Final"]
fracs   = [0.40, 0.72, 1.0]   # later snapshots show more developed damage

for row, (lbl, rlabel, rcol) in enumerate(map_rows):
    res = all_res[lbl]["results"]
    n   = len(res)

    for col, (frac, stage) in enumerate(zip(fracs, stages)):
        idx = max(0, min(n - 1, int(n * frac) - 1))
        ax  = fig.add_subplot(gs[row, col])
        d   = res[idx]["d_2d"]
        t_s = res[idx]["t"]
        dm  = res[idx]["d_max"]

        ax.pcolormesh(x_um, z_um, d.T, cmap=CMAP, norm=NORM,
                      rasterized=True, shading="nearest")
        ax.contour(x_um, z_um, d.T, levels=[0.5],
                   colors=["white"], linewidths=0.8, linestyles="--")

        ax.set_aspect("equal")
        ax.set_xlim(0, pf.W * 1e6)
        ax.set_ylim(0, pf.L * 1e6)

        ax.set_title(
            rf"\textit{{{stage}}},\ $t={t_s:.0f}$\,s" + "\n"
            + rf"$d_{{\max}}={dm:.2f}$",
            fontsize=7, pad=2,
        )

        if col == 0:
            ax.set_ylabel(rf"$z$\,[$\mu$m]", fontsize=8, labelpad=2)
            # Row label as text to the left
            ax.text(
                -0.38, 0.5, rlabel,
                transform=ax.transAxes,
                fontsize=8, rotation=90, va="center", ha="center", color=rcol,
            )
        else:
            ax.set_yticklabels([])

        if row == 1:
            ax.set_xlabel(rf"$x$\,[$\mu$m]", fontsize=8)
        else:
            ax.set_xticklabels([])

        ax.tick_params(labelsize=7)

# Shared colorbar spanning both map rows
cax = fig.add_subplot(gs[0:2, 4])
cb  = fig.colorbar(ScalarMappable(norm=NORM, cmap=CMAP), cax=cax)
cb.set_label(r"Damage $d$", fontsize=8, labelpad=3)
cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cb.ax.tick_params(labelsize=7)
cb.ax.axhline(0.5, color="white", lw=1.2, ls="--")     # mark d=0.5 on bar

# ── Time curves (row 2, cols 0-4 merged via span) ────────────────────────────
ax_t   = fig.add_subplot(gs[2, 0:5])
ax_bar = fig.add_subplot(gs[2, 6])

curve_info = [
    ("CH",      C_CH, "-",  "o", "CH (z-gradient)"),
    ("DH-base", C_DH, "--", "s", "DH-base (z-gradient)"),
    ("DH-Pg",   C_PG, "-",  "^", "DH+Pg (substrate)"),
]

for key, col, ls, mk, lab in curve_info:
    res  = all_res[key]["results"]
    ts   = [r["t"]     for r in res]
    dmx  = [r["d_max"] for r in res]
    step = max(1, len(ts) // 18)
    ax_t.plot(ts, dmx, color=col, ls=ls, lw=1.0,
              marker=mk, markevery=step, ms=2.5, label=lab)
    if tc[key]:
        ax_t.axvline(tc[key], color=col, ls=":", lw=0.7, alpha=0.7)
        ax_t.text(tc[key] + 0.6, 0.04, rf"${tc[key]:.0f}$\,s",
                  color=col, fontsize=6.5, rotation=90, va="bottom", ha="right")

ax_t.axhline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
ax_t.text(1.0, 0.52, r"$d=0.5$", color="gray", fontsize=6.5)
ax_t.set_xlabel(r"Physical time $t$\,[s]", fontsize=8)
ax_t.set_ylabel(r"$d_{\max}(t)$", fontsize=8)
ax_t.set_ylim(-0.02, 1.08)
ax_t.legend(fontsize=7, loc="lower right", framealpha=0.85,
            handlelength=1.5, labelspacing=0.3)
ax_t.tick_params(labelsize=7)
clean_ax(ax_t)

# Bar chart: t_crit
lbls_bar  = ["CH", "DH-base", "DH-Pg"]
disp_lbls = ["CH", "DH", "DH+Pg"]
vals_bar  = [tc[l] or 0 for l in lbls_bar]
bcols     = [C_CH, C_DH, C_PG]

xb   = np.arange(len(lbls_bar))
bars = ax_bar.bar(xb, vals_bar, color=bcols, alpha=0.85,
                  edgecolor="k", lw=0.5, width=0.60)
ax_bar.set_xticks(xb)
ax_bar.set_xticklabels(disp_lbls, fontsize=6.5, rotation=35, ha="right")
ax_bar.set_ylabel(r"$t_{\mathrm{crit}}$\,[s]", fontsize=8)
for bar, v in zip(bars, vals_bar):
    ax_bar.text(bar.get_x() + bar.get_width() / 2, v + 1.5,
                rf"${v:.0f}$", ha="center", va="bottom", fontsize=7)
ax_bar.tick_params(labelsize=7)
clean_ax(ax_bar)

# ── Suptitle ──────────────────────────────────────────────────────────────────
fig.suptitle(
    rf"AT2 phase-field fracture -- fully coupled Q4 FEM"
    rf"\ ($E={pf.E_BIO:.0f}$\,Pa,\ $\nu={pf.NU}$,\ "
    rf"$G_c={pf.G_C:.0e}$\,J/m$^2$,\ $\ell={pf.ELL*1e6:.0f}$\,$\mu$m,\ "
    rf"$h/\ell={max(hx, hz)/pf.ELL:.2f}$,\ "
    rf"$k_{{\mathrm{{ratio}}}}={pf.K_RATIO:.3f}$,\ "
    rf"$t_{{\mathrm{{crit}}}}^{{\mathrm{{DH}}}}/t_{{\mathrm{{crit}}}}^{{\mathrm{{CH}}}}={r1:.3f}$)",
    fontsize=8, y=1.01,
)

out_path = Path(__file__).parent / args.save
plt.savefig(str(out_path), bbox_inches="tight")
print(f"Saved: {out_path}")

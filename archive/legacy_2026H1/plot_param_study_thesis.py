"""
plot_param_study_thesis.py
==========================
Thesis-quality parameter study figure for AT2 biofilm fracture.
Sweeps (ell, G_c) at fixed growth rate; shows alpha_c, t_crit, ell sensitivity.

Run on fifa:
    cd ~/IKM_Hiwi/FEM
    PATH=~/texlive/2025/bin/x86_64-linux:$PATH \\
        python3 plot_param_study_thesis.py [--save param_study_thesis.pdf]
"""
import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(1, str(Path(__file__).parent / "JAXFEM"))
import phase_field_at2_2d_sparse as pf
from thesis_style import use, PALETTE, clean_ax

# ── CLI ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--save", default="param_study_thesis.pdf")
args = parser.parse_args()

# ── Baseline ──────────────────────────────────────────────────────────────────
ELL_BASE = 5.0e-6
GC_BASE  = 1.0e-5
E_BIO    = pf.E_BIO
NU       = pf.NU
NX, NZ   = 120, 80
N_STEPS  = 200

ALPHA_C_BASE = math.sqrt(GC_BASE * (1.0 - NU) / (2.0 * E_BIO * ELL_BASE))
K_EFF_BASE   = ALPHA_C_BASE * 1.6 / 120

print(f"Baseline alpha_c={ALPHA_C_BASE:.4f}  K_eff={K_EFF_BASE:.4e} 1/s")

# ── Parameter grid ────────────────────────────────────────────────────────────
ELL_VALS = np.array([3.0, 5.0, 8.0, 12.0]) * 1e-6
GC_VALS  = np.array([0.3, 1.0, 3.0, 10.0]) * 1e-5
n_ell    = len(ELL_VALS)
n_gc     = len(GC_VALS)

t_crit_grid  = np.full((n_ell, n_gc), np.nan)
alpha_c_grid = np.zeros((n_ell, n_gc))

print("Running (ell, Gc) grid ...")
for i, ell in enumerate(ELL_VALS):
    for j, gc in enumerate(GC_VALS):
        pf.ELL     = ell
        pf.G_C     = gc
        pf.ALPHA_C = math.sqrt(gc * (1.0 - NU) / (2.0 * E_BIO * ell))
        alpha_c_grid[i, j] = pf.ALPHA_C
        out = pf.run(K_EFF_BASE, "z_gradient", NX, NZ, n_steps=N_STEPS, label="")
        tc  = pf.t_crit(out["results"])
        t_crit_grid[i, j] = tc if tc else np.nan
        print(f"  ell={ell*1e6:.0f}um Gc={gc:.1e} -> alpha_c={pf.ALPHA_C:.4f} t={tc:.1f}s")

print("Running ell sweep (Gc fixed) ...")
ell_curves = {}
for ell in ELL_VALS:
    pf.ELL = ell; pf.G_C = GC_BASE
    pf.ALPHA_C = math.sqrt(GC_BASE * (1.0 - NU) / (2.0 * E_BIO * ell))
    out = pf.run(K_EFF_BASE, "z_gradient", NX, NZ, n_steps=N_STEPS, label="")
    ell_curves[ell] = out["results"]

# Restore globals
pf.ELL = ELL_BASE; pf.G_C = GC_BASE; pf.ALPHA_C = ALPHA_C_BASE

# ── Figure ────────────────────────────────────────────────────────────────────
figsize = use(width_frac=1.0, aspect=0.38)
fig = plt.figure(figsize=figsize)

gs = gridspec.GridSpec(
    1, 5,
    figure=fig,
    width_ratios=[1, 0.07, 1, 0.07, 1.15],
    wspace=0.10,
)

# Discrete grid positions for pcolormesh (edges halfway between points)
# x-axis: G_c index 0..n_gc, y-axis: ell index 0..n_ell
xi = np.arange(n_gc + 1) - 0.5   # edges: -0.5, 0.5, 1.5, 2.5, 3.5
yi = np.arange(n_ell + 1) - 0.5

ell_labels = [rf"${v*1e6:.0f}$" for v in ELL_VALS]
gc_labels  = [rf"${v*1e5:.1f}$" for v in GC_VALS]

# ── Panel A: alpha_c ──────────────────────────────────────────────────────────
ax_a  = fig.add_subplot(gs[0, 0])
cax_a = fig.add_subplot(gs[0, 1])

NORM_A = Normalize(alpha_c_grid.min(), alpha_c_grid.max())
im_a   = ax_a.pcolormesh(xi, yi, alpha_c_grid, cmap="viridis", norm=NORM_A,
                          shading="flat")
cb_a   = fig.colorbar(im_a, cax=cax_a)
cb_a.set_label(r"$\alpha_c$", fontsize=8, labelpad=2)
cb_a.ax.tick_params(labelsize=7)

# Annotate cells
for i in range(n_ell):
    for j in range(n_gc):
        v = alpha_c_grid[i, j]
        ax_a.text(j, i, f"{v:.3f}", ha="center", va="center",
                  fontsize=6.5,
                  color="white" if NORM_A(v) < 0.6 else "black")

# Baseline marker
ell_base_idx = list(ELL_VALS).index(ELL_BASE) if ELL_BASE in ELL_VALS else None
gc_base_idx  = list(GC_VALS).index(GC_BASE)  if GC_BASE  in GC_VALS  else None
if ell_base_idx is not None and gc_base_idx is not None:
    ax_a.add_patch(plt.Rectangle(
        (gc_base_idx - 0.5, ell_base_idx - 0.5), 1, 1,
        fill=False, edgecolor="white", lw=1.8, ls="--"))

ax_a.set_xticks(range(n_gc));  ax_a.set_xticklabels(gc_labels, fontsize=7)
ax_a.set_yticks(range(n_ell)); ax_a.set_yticklabels(ell_labels, fontsize=7)
ax_a.set_xlabel(r"$G_c\ [10^{-5}\ \mathrm{J/m}^2]$", fontsize=8)
ax_a.set_ylabel(r"$\ell\ [\mu\mathrm{m}]$", fontsize=8)
ax_a.set_title(r"Critical eigenstrain $\alpha_c$", fontsize=8, pad=4)
ax_a.set_xlim(-0.5, n_gc - 0.5); ax_a.set_ylim(-0.5, n_ell - 0.5)

# ── Panel B: t_crit ───────────────────────────────────────────────────────────
ax_b  = fig.add_subplot(gs[0, 2])
cax_b = fig.add_subplot(gs[0, 3])

NORM_B = Normalize(np.nanmin(t_crit_grid), np.nanmax(t_crit_grid))
im_b   = ax_b.pcolormesh(xi, yi, t_crit_grid, cmap="plasma", norm=NORM_B,
                          shading="flat")
cb_b   = fig.colorbar(im_b, cax=cax_b)
cb_b.set_label(r"$t_\mathrm{crit}\ [\mathrm{s}]$", fontsize=8, labelpad=2)
cb_b.ax.tick_params(labelsize=7)

for i in range(n_ell):
    for j in range(n_gc):
        v = t_crit_grid[i, j]
        if not np.isnan(v):
            ax_b.text(j, i, f"{v:.0f}", ha="center", va="center",
                      fontsize=6.5,
                      color="white" if NORM_B(v) < 0.55 else "black")

if ell_base_idx is not None and gc_base_idx is not None:
    ax_b.add_patch(plt.Rectangle(
        (gc_base_idx - 0.5, ell_base_idx - 0.5), 1, 1,
        fill=False, edgecolor="white", lw=1.8, ls="--"))

ax_b.set_xticks(range(n_gc));  ax_b.set_xticklabels(gc_labels, fontsize=7)
ax_b.set_yticks(range(n_ell)); ax_b.set_yticklabels(ell_labels, fontsize=7)
ax_b.set_xlabel(r"$G_c\ [10^{-5}\ \mathrm{J/m}^2]$", fontsize=8)
ax_b.set_ylabel(r"$\ell\ [\mu\mathrm{m}]$", fontsize=8)
ax_b.set_title(
    r"Fracture time $t_\mathrm{crit}$ (CH, fixed $\dot{\alpha}$)",
    fontsize=8, pad=4)
ax_b.set_xlim(-0.5, n_gc - 0.5); ax_b.set_ylim(-0.5, n_ell - 0.5)

# ── Panel C: ell sensitivity curves ──────────────────────────────────────────
ax_c = fig.add_subplot(gs[0, 4])

cmap_ell = plt.cm.viridis(np.linspace(0.15, 0.85, n_ell))
markers   = ["o", "s", "^", "D"]

for k, ell in enumerate(ELL_VALS):
    res  = ell_curves[ell]
    ts   = [r["t"]     for r in res]
    dmx  = [r["d_max"] for r in res]
    ac   = math.sqrt(GC_BASE * (1.0 - NU) / (2.0 * E_BIO * ell))
    tc   = pf.t_crit(res)
    step = max(1, len(ts) // 14)
    ax_c.plot(ts, dmx, color=cmap_ell[k], lw=1.0,
              marker=markers[k], markevery=step, ms=2.5,
              label=rf"$\ell={ell*1e6:.0f}\,\mu$m ($\alpha_c={ac:.3f}$)")
    if tc:
        ax_c.axvline(tc, color=cmap_ell[k], ls=":", lw=0.7, alpha=0.7)

ax_c.axhline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
ax_c.text(2, 0.52, r"$d=0.5$", color="gray", fontsize=6.5)
ax_c.set_xlabel(r"Physical time $t\ [\mathrm{s}]$", fontsize=8)
ax_c.set_ylabel(r"$d_\mathrm{max}(t)$", fontsize=8)
ax_c.set_ylim(-0.02, 1.08)
ax_c.set_title(
    rf"$\ell$ sensitivity ($G_c={GC_BASE:.0e}$\,J/m$^2$ fixed)",
    fontsize=8, pad=4)
ax_c.legend(fontsize=6.5, loc="lower right", framealpha=0.88,
            handlelength=1.3, labelspacing=0.3)
ax_c.tick_params(labelsize=7)
clean_ax(ax_c)

# ── Suptitle ──────────────────────────────────────────────────────────────────
fig.suptitle(
    rf"AT2 phase-field sensitivity -- prescribed eigenstrain (CH), "
    rf"$N_x\times N_z={NX}\times{NZ}$, "
    rf"fixed $\dot{{\alpha}}={K_EFF_BASE:.2e}$\,s$^{{-1}}$, "
    rf"$\alpha_c=\sqrt{{G_c(1-\nu)/(2E\ell)}}$; "
    rf"dashed box: baseline $\ell={ELL_BASE*1e6:.0f}\,\mu$m, "
    rf"$G_c={GC_BASE:.0e}$\,J/m$^2$",
    fontsize=7.5, y=1.02,
)

out_path = Path(__file__).parent / args.save
plt.savefig(str(out_path), bbox_inches="tight")
print(f"Saved: {out_path}")

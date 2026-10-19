"""
param_study_at2.py
==================
Parameter sensitivity of AT2 phase-field biofilm fracture.

Sweeps (ell, G_c) at fixed physical growth rate K_EFF_CH_BASE.
Uses prescribed-eigenstrain mode (fast: ~0.1 s / run).

Run on fifa:
    cd ~/IKM_Hiwi/FEM
    PATH=~/texlive/2025/bin/x86_64-linux:$PATH \\
        python3 param_study_at2.py [--save param_study_at2.pdf]
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
from matplotlib.colors import LogNorm

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(1, str(Path(__file__).parent / "JAXFEM"))
import phase_field_at2_2d_sparse as pf
from thesis_style import use, PALETTE, clean_ax

# ── CLI ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--save", default="param_study_at2.pdf")
args = parser.parse_args()

# ── Baseline (fixed physical growth rate) ────────────────────────────────────
ELL_BASE  = 5.0e-6     # [m]
GC_BASE   = 1.0e-5     # [J/m²]
E_BIO     = pf.E_BIO   # 1000 Pa
NU        = pf.NU       # 0.45
NX, NZ    = 120, 80
N_STEPS   = 200        # generous; early-stop when d_base > 0.95

ALPHA_C_BASE = math.sqrt(GC_BASE * (1.0 - NU) / (2.0 * E_BIO * ELL_BASE))
K_EFF_BASE   = ALPHA_C_BASE * 1.6 / 120   # fixed physical growth rate [1/s]

print(f"Baseline: ell={ELL_BASE*1e6:.0f}um  G_c={GC_BASE:.0e}  "
      f"alpha_c={ALPHA_C_BASE:.4f}  K_eff={K_EFF_BASE:.4e}")

# ── Parameter grid ────────────────────────────────────────────────────────────
ELL_VALS = np.array([3.0, 5.0, 8.0, 12.0]) * 1e-6   # [m]
GC_VALS  = np.array([0.3, 1.0, 3.0, 10.0]) * 1e-5   # [J/m²]

n_ell = len(ELL_VALS)
n_gc  = len(GC_VALS)

# ── Run grid ──────────────────────────────────────────────────────────────────
t_crit_grid = np.full((n_ell, n_gc), np.nan)
alpha_c_grid = np.zeros((n_ell, n_gc))

t0_all = time.time()
for i, ell in enumerate(ELL_VALS):
    for j, gc in enumerate(GC_VALS):
        # Patch module globals
        pf.ELL    = ell
        pf.G_C    = gc
        pf.ALPHA_C = math.sqrt(gc * (1.0 - NU) / (2.0 * E_BIO * ell))
        alpha_c_grid[i, j] = pf.ALPHA_C

        # Rebuild phase-field matrix constant
        # (LAP is rebuilt inside run() per call, so just need globals correct)
        t0 = time.time()
        out = pf.run(K_EFF_BASE, "z_gradient", NX, NZ,
                     n_steps=N_STEPS, label="")
        tc = pf.t_crit(out["results"])
        t_crit_grid[i, j] = tc if tc else np.nan

        dt = time.time() - t0
        print(f"  ell={ell*1e6:.0f}um  G_c={gc:.1e}  "
              f"alpha_c={pf.ALPHA_C:.4f}  t_crit={tc:.1f}s  ({dt:.2f}s)")

print(f"\nTotal: {time.time()-t0_all:.1f} s")

# Restore baseline globals
pf.ELL     = ELL_BASE
pf.G_C     = GC_BASE
pf.ALPHA_C = ALPHA_C_BASE

# ── Also: d_max(t) curves for ell sweep (fixed G_c = baseline) ───────────────
print("\nRunning ell sweep (fixed G_c)...")
ell_curves = {}
gc_fixed   = GC_BASE
for ell in ELL_VALS:
    pf.ELL     = ell
    pf.G_C     = gc_fixed
    pf.ALPHA_C = math.sqrt(gc_fixed * (1.0 - NU) / (2.0 * E_BIO * ell))
    out = pf.run(K_EFF_BASE, "z_gradient", NX, NZ, n_steps=N_STEPS, label="")
    ell_curves[ell] = out["results"]
    print(f"  ell={ell*1e6:.0f}um  alpha_c={pf.ALPHA_C:.4f}  "
          f"t_crit={pf.t_crit(out['results']):.1f}s")

# Restore
pf.ELL = ELL_BASE; pf.G_C = GC_BASE; pf.ALPHA_C = ALPHA_C_BASE

# ── Plot ──────────────────────────────────────────────────────────────────────
figsize = use(width_frac=1.0, aspect=0.42)
fig, axes = plt.subplots(1, 3, figsize=figsize,
                         gridspec_kw={"wspace": 0.40})

# Panel A: alpha_c heatmap
ax = axes[0]
im = ax.imshow(alpha_c_grid, aspect="auto", origin="lower",
               cmap="viridis",
               extent=[GC_VALS[0]*1e5, GC_VALS[-1]*1e5,
                       ELL_VALS[0]*1e6, ELL_VALS[-1]*1e6])
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label(r"$\alpha_c$", fontsize=8)
cb.ax.tick_params(labelsize=7)
# Mark baseline
ax.axhline(ELL_BASE*1e6, color="white", ls="--", lw=0.8)
ax.axvline(GC_BASE*1e5,  color="white", ls="--", lw=0.8)
ax.set_xlabel(r"$G_c$ [$10^{-5}$ J/m$^2$]", fontsize=8)
ax.set_ylabel(r"$\ell$ [$\mu$m]", fontsize=8)
ax.set_title(r"Critical eigenstrain $\alpha_c$", fontsize=8)
ax.set_xticks(GC_VALS * 1e5)
ax.set_yticks(ELL_VALS * 1e6)
ax.tick_params(labelsize=7)
clean_ax(ax)

# Panel B: t_crit heatmap
ax = axes[1]
im2 = ax.imshow(t_crit_grid, aspect="auto", origin="lower",
                cmap="plasma",
                extent=[GC_VALS[0]*1e5, GC_VALS[-1]*1e5,
                        ELL_VALS[0]*1e6, ELL_VALS[-1]*1e6])
cb2 = fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)
cb2.set_label(r"$t_\mathrm{crit}$ [s]", fontsize=8)
cb2.ax.tick_params(labelsize=7)
ax.axhline(ELL_BASE*1e6, color="white", ls="--", lw=0.8)
ax.axvline(GC_BASE*1e5,  color="white", ls="--", lw=0.8)
ax.set_xlabel(r"$G_c$ [$10^{-5}$ J/m$^2$]", fontsize=8)
ax.set_ylabel(r"$\ell$ [$\mu$m]", fontsize=8)
ax.set_title(r"Fracture time $t_\mathrm{crit}$ (CH, fixed $\dot{\alpha}$)", fontsize=8)
ax.set_xticks(GC_VALS * 1e5)
ax.set_yticks(ELL_VALS * 1e6)
ax.tick_params(labelsize=7)
# Annotate each cell
for ii in range(n_ell):
    for jj in range(n_gc):
        v = t_crit_grid[ii, jj]
        if not np.isnan(v):
            ax.text(GC_VALS[jj]*1e5, ELL_VALS[ii]*1e6,
                    f"{v:.0f}", ha="center", va="center",
                    fontsize=6, color="white")
clean_ax(ax)

# Panel C: d_max(t) curves for ell sweep
ax = axes[2]
cmap_ell = plt.cm.viridis(np.linspace(0.1, 0.9, n_ell))
for k, ell in enumerate(ELL_VALS):
    res  = ell_curves[ell]
    ts   = [r["t"]     for r in res]
    dmx  = [r["d_max"] for r in res]
    step = max(1, len(ts) // 15)
    tc   = pf.t_crit(res)
    ax.plot(ts, dmx, color=cmap_ell[k], lw=1.0,
            marker="o", markevery=step, ms=2.5,
            label=rf"$\ell={ell*1e6:.0f}\,\mu$m ($\alpha_c={math.sqrt(gc_fixed*(1-NU)/(2*E_BIO*ell)):.3f}$)")
    if tc:
        ax.axvline(tc, color=cmap_ell[k], ls=":", lw=0.7, alpha=0.6)

ax.axhline(0.5, color="gray", ls=":", lw=0.8)
ax.set_xlabel(r"Physical time $t$ [s]", fontsize=8)
ax.set_ylabel(r"$d_\mathrm{max}(t)$", fontsize=8)
ax.set_ylim(-0.02, 1.08)
ax.set_title(r"$\ell$ sensitivity ($G_c$ fixed)", fontsize=8)
ax.legend(fontsize=6.5, loc="lower right", framealpha=0.85,
          handlelength=1.3, labelspacing=0.3)
ax.tick_params(labelsize=7)
clean_ax(ax)

fig.suptitle(
    rf"AT2 sensitivity: $\ell$, $G_c$ (prescribed eigenstrain, CH, Nx={NX}, Nz={NZ},"
    rf" fixed $\dot{{\alpha}}={K_EFF_BASE:.2e}$\,s$^{{-1}}$,"
    rf" baseline $\ell={ELL_BASE*1e6:.0f}\,\mu$m, $G_c={GC_BASE:.0e}$\,J/m$^2$"
    rf" shown as dashed white)",
    fontsize=7.5, y=1.02,
)

out_path = Path(__file__).parent / args.save
plt.savefig(str(out_path), bbox_inches="tight")
print(f"Saved: {out_path}")

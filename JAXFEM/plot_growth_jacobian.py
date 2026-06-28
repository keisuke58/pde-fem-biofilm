"""
plot_growth_jacobian.py
=======================
Thesis-quality figure: composition sensitivity dσ/dφ_i (growth Jacobian).

Two-panel layout:
  Left  — absolute sensitivity ∂σ/∂φ_i  [relative stress / fraction]
  Right — normalised elasticity ε_i = (φ_i/σ) ∂σ/∂φ_i  [dimensionless]

Both CH (commensal) and DH (dysbiotic) conditions shown side by side.
Pg dominance is the key message.

Usage:
    PATH=~/texlive/2025/bin/x86_64-linux:$PATH \
        python JAXFEM/plot_growth_jacobian.py [--save fig_jacobian.pdf]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, str(Path(__file__).parent))
from thesis_style import use, PALETTE, clean_ax

# ── JAX / model (same as growth_jacobian.py) ─────────────────────────────────
import jax
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

K_B     = jnp.array([2.0, 0.8, 0.2])
SPECIES = [r"$P_g$", r"$F_n$", "Other"]
B_EXP   = 2.68

PHI_CH = jnp.array([0.42, 0.33, 0.25])
PHI_DH = jnp.array([0.21, 0.48, 0.31])


def sigma_hat(phi):
    return jnp.dot(phi, K_B) ** B_EXP

grad_sigma = jax.grad(sigma_hat)


def jacobian_data(phi):
    sig  = float(sigma_hat(phi))
    dsig = np.array(grad_sigma(phi))
    elas = np.array(phi) * dsig / sig
    return sig, dsig, elas


# ── Compute ───────────────────────────────────────────────────────────────────
sig_ch, dsig_ch, elas_ch = jacobian_data(PHI_CH)
sig_dh, dsig_dh, elas_dh = jacobian_data(PHI_DH)
ratio = sig_ch / sig_dh

# ── Figure ────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--save", default="fig_growth_jacobian.pdf")
args = parser.parse_args()

figsize = use(width_frac=1.0, aspect=0.38)
fig, axes = plt.subplots(1, 2, figsize=figsize)

x     = np.arange(len(SPECIES))
width = 0.35
c_ch  = PALETTE["ch"]
c_dh  = PALETTE["dh"]

# ── Left: absolute sensitivity ────────────────────────────────────────────────
ax = axes[0]
ax.bar(x - width/2, dsig_ch, width, color=c_ch, label="CH")
ax.bar(x + width/2, dsig_dh, width, color=c_dh, label="DH")
ax.set_xticks(x)
ax.set_xticklabels(SPECIES)
ax.set_ylabel(r"$\partial\hat{\sigma}/\partial\phi_i$")
ax.set_title(r"Absolute sensitivity")
clean_ax(ax)
ax.legend(frameon=False)

# ── Right: elasticity ─────────────────────────────────────────────────────────
ax = axes[1]
ax.bar(x - width/2, elas_ch, width, color=c_ch, label="CH")
ax.bar(x + width/2, elas_dh, width, color=c_dh, label="DH")
ax.set_xticks(x)
ax.set_xticklabels(SPECIES)
ax.set_ylabel(r"Elasticity $\varepsilon_i = (\phi_i/\hat{\sigma})\,\partial\hat{\sigma}/\partial\phi_i$")
ax.set_title(r"Normalised sensitivity")
clean_ax(ax)
ax.legend(frameon=False)

# ratio annotation on right panel
ax.text(0.97, 0.97,
        r"$\hat{\sigma}_\mathrm{CH}/\hat{\sigma}_\mathrm{DH} = %.2f$" % ratio,
        transform=ax.transAxes, ha="right", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.5))

fig.tight_layout()
out = Path(args.save)
fig.savefig(out)
print(f"Saved: {out}")

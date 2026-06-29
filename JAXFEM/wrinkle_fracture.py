"""
wrinkle_fracture.py
===================
Phase diagram: buckle-delamination vs. AT2 fracture for a thin biofilm on rigid substrate.

A biofilm of thickness h grows under equibiaxial strain α(t) = k_eff · t on a rigid
titanium substrate (E_sub >> E_film).  Two competing failure modes:

  (A) AT2 fracture (delamination):
      nucleates when α ≥ α_c = sqrt(G_c / (E · ℓ))

  (B) Buckle-delamination:
      a pre-existing interface defect of half-length a_0 buckles when
      the compressive stress exceeds the Euler buckling load:
          ε_buckle(a_0, h) = π² h² / (12 (1 - ν²) a_0²)   [plane-strain Euler column]

Phase boundary: α_c = ε_buckle  →  a_0*(h) = π h / sqrt(12 (1-ν²) α_c)

Below boundary (small defects): AT2 fracture first  → delamination without buckling
Above boundary (large defects): buckle-delamination first  → wrinkling/blistering
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase_field_at2_2d_sparse import G_C, ELL, E_BIO, K_EFF_CH, K_EFF_DH
from thesis_style import use, PALETTE, TEXTWIDTH_IN, clean_ax

NU = 0.45  # biofilm Poisson ratio (nearly incompressible soft hydrogel)


def alpha_c(gc: float, E: float = E_BIO, ell: float = ELL) -> float:
    """AT2 critical growth strain for damage nucleation."""
    return np.sqrt(gc / (E * ell))


def eps_buckle(a0: np.ndarray, h: float) -> np.ndarray:
    """Critical compressive strain for Euler buckling of a delaminated flap."""
    return np.pi**2 * h**2 / (12.0 * (1.0 - NU**2) * a0**2)


def a0_boundary(gc: float, h: float) -> float:
    """Phase boundary defect size a_0* for given G_c and film thickness h."""
    return np.pi * h / np.sqrt(12.0 * (1.0 - NU**2) * alpha_c(gc))


def make_figure(out: str = "fig_wrinkle_fracture.pdf"):
    use()
    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH_IN, TEXTWIDTH_IN * 0.46))

    h_arr  = np.logspace(-6, -3, 400)   # 1 μm – 1 mm
    a0_arr = np.logspace(-6, -2, 400)   # 1 μm – 10 mm

    G_CASES = [
        (G_C,  r"model $G_c\!=\!10^{-5}\,\mathrm{J\,m^{-2}}$", "0.55"),
        (5e-3, r"lit lo $\Gamma\!=\!5\,\mathrm{mJ\,m^{-2}}$",   PALETTE["ch"]),
        (1e-2, r"lit hi $\Gamma\!=\!10\,\mathrm{mJ\,m^{-2}}$",  PALETTE["dh"]),
    ]

    # ── left panel: phase diagram (h, a_0) ───────────────────────────────────
    ax = axes[0]
    for gc, label, col in G_CASES:
        a0_b = np.array([a0_boundary(gc, h) for h in h_arr])
        ax.loglog(h_arr * 1e6, a0_b * 1e6, lw=1.6, color=col, label=label)

    # Shade AT2-first region (model G_c boundary as reference)
    a0_b_model = np.array([a0_boundary(G_C, h) for h in h_arr])
    ax.fill_between(h_arr * 1e6, a0_b_model * 1e6, 1e4,
                    alpha=0.07, color=PALETTE["dh"])
    ax.fill_between(h_arr * 1e6, 1e0, a0_b_model * 1e6,
                    alpha=0.07, color=PALETTE["ch"])
    ax.text(2.5, 3000, "buckle-delam.\nfirst", fontsize=6, color=PALETTE["dh"])
    ax.text(300,    2, "AT2 fracture\nfirst",  fontsize=6, color=PALETTE["ch"])

    # Realistic biofilm parameter box
    ax.axvspan(10, 500, alpha=0.10, color="0.5", zorder=0)
    ax.axhspan(1,  100, alpha=0.10, color="0.5", zorder=0)
    ax.text(12, 120, "biofilm\nrange", fontsize=5.5, color="0.35")

    ax.set_xlabel(r"Film thickness $h$ [μm]")
    ax.set_ylabel(r"Defect half-length $a_0$ [μm]")
    ax.set_xlim(1, 1000)
    ax.set_ylim(1, 1e4)
    ax.legend(fontsize=5.5, loc="upper left")
    ax.set_title(r"(a) Phase diagram: buckle-delam. vs AT2", fontsize=7, pad=4)
    clean_ax(ax)

    # ── right panel: t_crit vs G_c for CH and DH ─────────────────────────────
    ax2 = axes[1]
    gc_sweep = np.logspace(-7, -1, 400)
    ac_sweep = np.array([alpha_c(g) for g in gc_sweep])
    t_ch = ac_sweep / K_EFF_CH
    t_dh = ac_sweep / K_EFF_DH

    ax2.loglog(gc_sweep, t_ch, color=PALETTE["ch"], lw=1.6, label="CH")
    ax2.loglog(gc_sweep, t_dh, color=PALETTE["dh"], lw=1.6, ls="--", label="DH-unif")

    # Literature Γ band
    ax2.axvspan(5e-3, 1e-2, alpha=0.13, color="0.7", zorder=0)
    ax2.text(5.5e-3, 8e3, r"$\Gamma_\mathrm{lit}$", fontsize=5.5, color="0.4")

    # Literature t_crit band
    ax2.axhline(48 * 3600, color="0.55", lw=0.8, ls=":")
    ax2.axhline(72 * 3600, color="0.55", lw=0.8, ls=":")
    ax2.text(1.2e-7, 74 * 3600, "48–72 h", fontsize=5.5, color="0.4")

    ax2.set_xlabel(r"$G_c$ [J m$^{-2}$]")
    ax2.set_ylabel(r"$t_\mathrm{crit}$ [s]")
    ax2.set_xlim(1e-7, 1e-1)
    ax2.legend(fontsize=6)
    ax2.set_title(r"(b) $t_\mathrm{crit}(G_c)$ for CH and DH", fontsize=7, pad=4)
    clean_ax(ax2)

    fig.tight_layout(pad=0.8)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")

    print("\nPhase boundary summary (a_0* at h = 100 μm):")
    for gc, label, _ in G_CASES:
        ac = alpha_c(gc)
        a0s = a0_boundary(gc, h=100e-6)
        print(f"  {label:50s}  α_c={ac:.4f}  a_0*={a0s*1e6:.1f} μm")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="fig_wrinkle_fracture.pdf")
    args = ap.parse_args()
    make_figure(args.out)

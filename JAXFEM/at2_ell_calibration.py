"""
at2_ell_calibration.py
======================
Self-consistency check: is AT2 regularisation length ℓ = ELL = 5 μm physically justified?

AT2 (second-order) 1-D localisation solution (Tanne et al. 2018, JMPS):
    σ_c = sqrt( 3·E·Gc / (8·ℓ) )

Inverse: given E, Gc, σ_c  →  ℓ_phys = 3·E·Gc / (8·σ_c²)

If ℓ_phys ≈ ELL for (E_BIO, Gc_lit, σ_c_lit), the choice ELL = 5 μm is not
a free numerical parameter but is consistent with biofilm material data.

Panel (a): ℓ_phys(σ_c) for several G_c values — where does the ELL=5μm line land?
Panel (b): (G_c, σ_c) contour map of ℓ_phys; ELL=5μm contour highlighted.

Usage:
    python at2_ell_calibration.py
    python at2_ell_calibration.py --out fig_at2_ell_calib.pdf
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase_field_at2_2d_sparse import ELL, E_BIO, G_C
from thesis_style import use, PALETTE, TEXTWIDTH_IN, clean_ax


def sigma_c(E: float, gc: float, ell: float) -> float:
    """AT2 cohesive strength from Tanne et al. 2018 (1-D bar localisation)."""
    return np.sqrt(3.0 * E * gc / (8.0 * ell))


def ell_phys(E: float, gc: float, sc: float) -> float:
    """Process-zone length implied by (E, Gc, σ_c) via AT2 formula."""
    return 3.0 * E * gc / (8.0 * sc**2)


def make_figure(out: str = "fig_at2_ell_calib.pdf"):
    use()
    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH_IN, TEXTWIDTH_IN * 0.46))

    sc_arr = np.logspace(1, 3, 400)   # σ_c : 10 – 1000 Pa

    # G_c cases: model value + literature band endpoints
    GC_CASES = [
        (G_C,  r"model $G_c\!=\!10^{-5}$",      "0.55",        ":"),
        (5e-3, r"lit lo $\Gamma\!=\!5\,$mJ/m²",  PALETTE["ch"], "-"),
        (1e-2, r"lit hi $\Gamma\!=\!10\,$mJ/m²", PALETTE["dh"], "--"),
    ]

    # ── panel (a): ℓ_phys vs σ_c ──────────────────────────────────────────────
    ax = axes[0]
    for gc, label, col, ls in GC_CASES:
        ell_arr = np.array([ell_phys(E_BIO, gc, sc) for sc in sc_arr])
        ax.loglog(sc_arr, ell_arr * 1e6, lw=1.6, color=col, ls=ls, label=label)

    # ELL reference line
    ax.axhline(ELL * 1e6, color="k", lw=1.0, ls="-.",
               label=rf"ELL $= {ELL*1e6:.0f}\,\mu$m")

    # Literature σ_c band for oral biofilm (Hall-Stoodley et al.; Stoodley et al.)
    ax.axvspan(200, 700, alpha=0.10, color="0.5", zorder=0)
    ax.text(210, 0.05, r"$\sigma_c^\mathrm{lit}$" "\n(oral BF)", fontsize=5.5, color="0.35")

    ax.set_xlabel(r"Cohesive strength $\sigma_c$ [Pa]")
    ax.set_ylabel(r"Process-zone length $\ell_\mathrm{phys}$ [μm]")
    ax.set_xlim(10, 1000)
    ax.set_ylim(0.01, 1e4)
    ax.legend(fontsize=5.5, loc="upper right")
    ax.set_title(r"(a) $\ell_\mathrm{phys}(\sigma_c)$ for several $G_c$", fontsize=7, pad=4)
    clean_ax(ax)

    # ── panel (b): contour map ℓ_phys(G_c, σ_c) ──────────────────────────────
    ax2 = axes[1]
    gc_2d = np.logspace(-6, -1, 300)
    sc_2d = np.logspace(1, 3, 300)
    GC2, SC2 = np.meshgrid(gc_2d, sc_2d)
    ELL2 = ell_phys(E_BIO, GC2, SC2) * 1e6  # μm

    # Filled contours
    levels = [0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 500]
    cs = ax2.contourf(gc_2d, sc_2d, ELL2,
                      levels=levels, norm=plt.matplotlib.colors.LogNorm(),
                      cmap="YlOrBr_r", alpha=0.85)
    fig.colorbar(cs, ax=ax2, label=r"$\ell_\mathrm{phys}$ [μm]", pad=0.02,
                 ticks=[0.1, 1, 5, 10, 50, 100])

    # ELL = 5 μm contour
    ct = ax2.contour(gc_2d, sc_2d, ELL2,
                     levels=[ELL * 1e6],
                     colors=["k"], linewidths=1.4, linestyles=["-"])
    ax2.clabel(ct, fmt=rf"$\ell={ELL*1e6:.0f}\,\mu$m", fontsize=6)

    # Mark current model point (G_C, σ_c implied by ELL)
    sc_model = sigma_c(E_BIO, G_C, ELL)
    ax2.scatter([G_C], [sc_model], color="k", s=30, zorder=5,
                label=rf"model: $\sigma_c={sc_model:.0f}$ Pa")

    # Mark literature G_c range × σ_c range box
    ax2.axvspan(5e-3, 1e-2, alpha=0.12, color=PALETTE["ch"], zorder=0)
    ax2.axhspan(200, 700,   alpha=0.12, color=PALETTE["dh"], zorder=0)
    ax2.text(5.5e-3, 720, r"$\Gamma_\mathrm{lit}$" r" × $\sigma_c^\mathrm{lit}$",
             fontsize=5.5, color="0.3")

    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel(r"$G_c$ [J m$^{-2}$]")
    ax2.set_ylabel(r"Cohesive strength $\sigma_c$ [Pa]")
    ax2.set_xlim(1e-6, 1e-1)
    ax2.set_ylim(10, 1000)
    ax2.legend(fontsize=5.5, loc="upper left")
    ax2.set_title(r"(b) Contour: $\ell_\mathrm{phys}(G_c,\,\sigma_c)$", fontsize=7, pad=4)
    clean_ax(ax2)

    fig.tight_layout(pad=0.8)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")

    # ── console summary ────────────────────────────────────────────────────────
    print(f"\nAT2 self-consistency check (E = {E_BIO:.0f} Pa, ELL = {ELL*1e6:.0f} μm):")
    print(f"  σ_c implied by (G_C_model, ELL) = {sigma_c(E_BIO, G_C, ELL):.1f} Pa")
    for gc, label, _, _ in GC_CASES:
        sc_self = sigma_c(E_BIO, gc, ELL)
        print(f"  σ_c for {label:40s}  Gc={gc:.1e} → σ_c = {sc_self:.1f} Pa")
    print()
    print("  ℓ_phys for Gc_lit = 5 mJ/m²:")
    for sc in [100, 300, 500, 700]:
        lp = ell_phys(E_BIO, 5e-3, sc) * 1e6
        note = " ← ELL ✓" if 3 < lp < 15 else ""
        print(f"    σ_c = {sc:4d} Pa  →  ℓ_phys = {lp:.1f} μm{note}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="fig_at2_ell_calib.pdf")
    args = ap.parse_args()
    make_figure(args.out)

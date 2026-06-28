"""
gc_identification.py
====================
Identify biofilm fracture energy G_c from observed delamination time.

AT2 analytical inverse (uniform-field 0D):
    t_crit = sqrt(G_c / (E * ell)) / k_eff^b
    G_c    = (k_eff^b * t_crit)^2 * E * ell

The formula is exact for the uniform case because damage nucleates when
    H = 1/2 * E * alpha^2 = G_c / (2*ell)
with alpha = k_eff^b * t, giving the closed-form above.

Usage
-----
    python gc_identification.py                              # sweep + figure
    python gc_identification.py --t-exp 75 --condition CH   # validate vs sim
    python gc_identification.py --t-exp 86400 --condition CH --k-eff-phys 1e-5
    python gc_identification.py --t-exp 3600 --t-exp-std 360 --condition DH

Note on time units
------------------
K_EFF_CH/DH in phase_field_at2_2d_sparse.py are numerically scaled so the
simulation covers ~200 steps.  To identify G_c from real Debener data, pass
a physically calibrated k_eff via --k-eff-phys (units: 1/s or same as t_exp).
"""

import math
import argparse
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase_field_at2_2d_sparse import (
    G_C as GC_MODEL, ELL, E_BIO, K_EFF_CH, K_EFF_DH, K_RATIO, ALPHA_C,
    run, t_crit as compute_tcrit,
)
from thesis_style import use, PALETTE, TEXTWIDTH_IN, clean_ax


# ── analytical formulas ──────────────────────────────────────────────────────

def gc_from_tcrit(t_exp: float, k_eff: float,
                  E: float = E_BIO, ell: float = ELL) -> float:
    """Identify G_c [J/m²] from observed delamination time t_exp [s]."""
    return (k_eff * t_exp) ** 2 * E * ell


def tcrit_from_gc(Gc: float, k_eff: float,
                  E: float = E_BIO, ell: float = ELL) -> float:
    """Predict t_crit [s] for given G_c [J/m²] and k_eff [1/s]."""
    return math.sqrt(Gc / (E * ell)) / k_eff


def gc_uncertainty(t_exp: float, dt_exp: float, k_eff: float) -> float:
    """
    Propagate ±dt_exp [s] into ΔG_c [J/m²].
    G_c = (k*t)²*E*ℓ  →  ΔG_c/G_c = 2*Δt/t  (first-order)
    """
    return 2.0 * gc_from_tcrit(t_exp, k_eff) * dt_exp / t_exp


# ── validation against AT2 simulation ────────────────────────────────────────

def validate(nx: int = 80, nz: int = 50, n_steps: int = 200) -> None:
    """
    Run AT2 simulations at K_EFF_CH and K_EFF_DH (module defaults),
    recover G_c from t_crit, compare to GC_MODEL.
    """
    for label, k in [("CH", K_EFF_CH), ("DH-unif", K_EFF_DH)]:
        res = run(k, "uniform", nx, nz, n_steps, label=f"validate-{label}")
        tc  = compute_tcrit(res["results"])
        if tc is None:
            print(f"  {label}: t_crit not reached in simulation")
            continue
        gc_id = gc_from_tcrit(tc, k)
        err   = abs(gc_id - GC_MODEL) / GC_MODEL * 100
        print(f"  {label:8s}: t_crit={tc:.1f}s  "
              f"G_c_id={gc_id:.3e}  G_c_model={GC_MODEL:.3e}  err={err:.1f}%")


# ── G_c sweep ────────────────────────────────────────────────────────────────

def sweep(n: int = 300,
          gc_lo: float = 1e-8,
          gc_hi: float = 1e-1) -> tuple:
    """
    Return (gc_arr, tc_ch, tc_dh) over log-spaced G_c values.
    k_eff is held fixed (module K_EFF_CH / K_EFF_DH).
    """
    gc_arr = np.logspace(math.log10(gc_lo), math.log10(gc_hi), n)
    tc_ch  = np.array([tcrit_from_gc(g, K_EFF_CH) for g in gc_arr])
    tc_dh  = np.array([tcrit_from_gc(g, K_EFF_DH) for g in gc_arr])
    return gc_arr, tc_ch, tc_dh


# ── figure ────────────────────────────────────────────────────────────────────

def make_figure(t_exp=None, t_exp_std=None, condition="CH", k_eff_phys=None,
                out="fig_gc_identification.pdf"):
    use()
    import matplotlib.pyplot as plt

    gc_arr, tc_ch, tc_dh = sweep()
    PURPLE = "#8B6BB1"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(TEXTWIDTH_IN, TEXTWIDTH_IN * 0.42))

    # ── left panel: G_c vs t_crit identification curve ──────────────────────
    ax1.loglog(tc_ch, gc_arr, color=PALETTE["ch"], lw=1.8, label="CH")
    ax1.loglog(tc_dh, gc_arr, color=PALETTE["dh"], lw=1.8, ls="--", label="DH-unif")

    # current model reference
    ax1.axhline(GC_MODEL, color="0.65", lw=0.8, ls=":")
    ax1.scatter([tcrit_from_gc(GC_MODEL, K_EFF_CH)],
                [GC_MODEL], color=PALETTE["ch"], s=30, zorder=5)
    ax1.scatter([tcrit_from_gc(GC_MODEL, K_EFF_DH)],
                [GC_MODEL], color=PALETTE["dh"], s=30, zorder=5)
    ax1.text(tcrit_from_gc(GC_MODEL, K_EFF_CH) * 1.12, GC_MODEL * 1.5,
             r"model ($G_c={10^{-5}}$)", fontsize=5.5, color="0.45")

    # physical k_eff identification curve (if provided)
    if k_eff_phys is not None:
        tc_phys = np.array([tcrit_from_gc(g, k_eff_phys) for g in gc_arr])
        ax1.loglog(tc_phys, gc_arr, color=PURPLE, lw=1.2, ls=":", label="phys $k$")

    # experimental point
    if t_exp is not None:
        k = k_eff_phys if k_eff_phys is not None else (
            K_EFF_CH if condition == "CH" else K_EFF_DH)
        col = PALETTE["ch"] if condition == "CH" else PALETTE["dh"]
        gc_id = gc_from_tcrit(t_exp, k)
        ax1.scatter([t_exp], [gc_id], marker="*", s=80, color=col,
                    zorder=6, label=rf"Debener $t^\mathrm{{exp}}$={t_exp:.0f} s")
        if t_exp_std is not None:
            dGc = gc_uncertainty(t_exp, t_exp_std, k)
            ax1.errorbar([t_exp], [gc_id],
                         yerr=[[dGc], [dGc]],
                         xerr=[[t_exp_std], [t_exp_std]],
                         fmt="none", color=col, capsize=3, lw=1.0)
        print(f"\nIdentification ({condition}): "
              f"t_exp={t_exp:.1f}s → G_c = {gc_id:.3e} J/m²")
        if t_exp_std is not None:
            dGc = gc_uncertainty(t_exp, t_exp_std, k)
            print(f"  ±{t_exp_std:.1f}s → ±{dGc:.3e} J/m²  "
                  f"({dGc/gc_id*100:.0f}%)")

    ax1.set_xlabel(r"$t_\mathrm{crit}$ [s]")
    ax1.set_ylabel(r"$G_c$ [J\,m$^{-2}$]")
    ax1.legend(fontsize=6, loc="upper left")
    ax1.set_title(r"(a) $G_c$–$t_\mathrm{crit}$ identification", fontsize=7, pad=4)

    # ── right panel: uncertainty propagation ΔGc/Gc = 2 Δt/t ─────────────
    dt_rel  = np.linspace(0.01, 0.50, 200)
    dgc_rel = 2.0 * dt_rel

    ax2.plot(dt_rel * 100, dgc_rel * 100, "k-", lw=1.8)
    # reference: ±20 % in time → ±40 % in G_c
    ax2.axvline(20, color="0.65", lw=0.8, ls=":")
    ax2.axhline(40, color="0.65", lw=0.8, ls=":")
    ax2.fill_between([0, 20], [0, 0], [40, 40], alpha=0.12, color="steelblue")
    ax2.text(1, 43, r"$\pm$20\% $t$ $\to$ $\pm$40\% $G_c$",
             fontsize=6, color="0.35")

    ax2.set_xlabel(r"$\Delta t_\mathrm{exp}/t_\mathrm{exp}$ [\%]")
    ax2.set_ylabel(r"$\Delta G_c/G_c$ [\%]")
    ax2.set_xlim(0, 50)
    ax2.set_ylim(0, 100)
    ax2.set_title(r"(b) Uncertainty propagation", fontsize=7, pad=4)

    fig.tight_layout(pad=0.8)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t-exp",      type=float, default=None,
                    help="Experimental delamination time [s]")
    ap.add_argument("--t-exp-std",  type=float, default=None,
                    help="1-sigma uncertainty on t_exp [s]")
    ap.add_argument("--condition",  choices=["CH", "DH"], default="CH",
                    help="Biofilm condition of the experiment")
    ap.add_argument("--k-eff-phys", type=float, default=None,
                    help="Physically calibrated k_eff^b [1/s]; overrides model default")
    ap.add_argument("--validate",   action="store_true",
                    help="Run AT2 sim to validate analytical inverse")
    ap.add_argument("--out",        default="fig_gc_identification.pdf")
    args = ap.parse_args()

    print(f"E={E_BIO:.0f} Pa  G_c(model)={GC_MODEL:.1e} J/m²  "
          f"ell={ELL*1e6:.0f} um  K_ratio={K_RATIO:.2f}")
    print(f"k_eff_CH={K_EFF_CH:.3e} /s  k_eff_DH={K_EFF_DH:.3e} /s")

    if args.validate:
        print("\nRunning validation simulations …")
        validate()

    make_figure(
        t_exp=args.t_exp,
        t_exp_std=args.t_exp_std,
        condition=args.condition,
        k_eff_phys=args.k_eff_phys,
        out=args.out,
    )


if __name__ == "__main__":
    main()

"""
replicon_at2_bridge.py
======================
Bridge: ecological Jacobian stability ↔ mechanical fracture time.

For each TMCMC posterior sample θ = (A, b):
  1. Integrate gLV ODE to φ* (steady state)
  2. k_eff(φ*) = Σ_i φ_i* · k_alpha_i   (species-weighted effective growth rate)
  3. t_frac = α_c / k_eff  (time to reach AT2 fracture threshold; proxy)
  4. λ_max(Re J(φ*))       (community Jacobian stability margin)

Scatter: x = λ_max(J), y = t_frac, color = condition.

Interpretation:
  DH: λ_max(J) → 0 (marginal stability) AND low k_eff (Pg has k_alpha=0.3)
      → ecologically unstable, but mechanically slow to fracture
  CH: λ_max(J) spread (robust), higher k_eff
      → ecologically stable, but mechanically faster to fracture
  → ecological–mechanical trade-off

Usage:
    python replicon_at2_bridge.py
    python replicon_at2_bridge.py --n-ss 200 --out fig_bridge.pdf
"""

import argparse
import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent.parent   # ~/IKM_Hiwi/
sys.path.insert(0, str(Path(__file__).parent))

from phase_field_at2_2d_sparse import G_C, ELL, E_BIO, K_EFF_CH, K_EFF_DH
from thesis_style import use, PALETTE, TEXTWIDTH_IN, clean_ax
from replicon_analysis import (
    load_samples, theta_to_matrices, community_jacobian, steady_state,
    ULTIMATE_BASE, CH_10K, DH_10K, DS_10K, CS_10K,
    C_CONST, K_HILL, N_HILL, T_STEADY, SPECIES,
)

K_ALPHA = np.array([1.0, 0.8, 0.4, 0.6, 0.3])   # So, An, Vd, Fn, Pg
ALPHA_C = np.sqrt(G_C / (E_BIO * ELL))


def k_eff_from_phi(phi: np.ndarray) -> float:
    """Effective growth rate from steady-state composition."""
    phi = np.maximum(phi, 0.0)
    total = phi.sum()
    if total < 1e-6:
        return 0.0
    return float(np.dot(phi / total, K_ALPHA))


def analyze_bridge(samples: np.ndarray, n_ss: int = 100):
    """
    Returns arrays: lam_max_J, k_eff, t_frac.
    Only samples 0..n_ss-1 have Jacobian computed (ODE is slow).
    """
    N = min(n_ss, len(samples))
    lam_J  = np.full(N, np.nan)
    keff   = np.zeros(N)
    t_frac = np.full(N, np.nan)

    phi0 = np.array([0.35, 0.15, 0.20, 0.20, 0.10])

    for k in range(N):
        A, b = theta_to_matrices(samples[k])
        phi_eq = steady_state(A, b, phi0.copy())
        if phi_eq.sum() < 0.01:
            continue
        J = community_jacobian(phi_eq, A, b)
        lam_J[k] = np.max(np.real(np.linalg.eigvals(J)))
        keff[k]  = k_eff_from_phi(phi_eq)
        if keff[k] > 1e-9:
            t_frac[k] = ALPHA_C / keff[k]

    mask = ~np.isnan(lam_J) & ~np.isnan(t_frac) & (keff > 0)
    return {"lam_J": lam_J[mask], "k_eff": keff[mask], "t_frac": t_frac[mask]}


def make_figure(results: dict, out: str = "fig_replicon_at2_bridge.pdf"):
    """
    results: {label: {lam_J, k_eff, t_frac}}
    """
    use()
    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH_IN, TEXTWIDTH_IN * 0.46))

    COLOR_MAP = {lbl: (PALETTE["ch"] if lbl.upper().startswith("C") else PALETTE["dh"])
                 for lbl in results}

    # ── left: λ_max(J) vs t_frac ─────────────────────────────────────────────
    ax = axes[0]
    for lbl, res in results.items():
        col = COLOR_MAP[lbl]
        ax.scatter(res["lam_J"], res["t_frac"],
                   s=8, color=col, alpha=0.5, label=lbl, zorder=3)

    ax.axvline(0, color="0.6", lw=0.8, ls="--")
    ax.set_xlabel(r"$\lambda_\mathrm{max}(\mathrm{Re}\,J)$")
    ax.set_ylabel(r"$t_\mathrm{frac} = \alpha_c / k_\mathrm{eff}$ [–]")
    ax.set_title(r"(a) Ecological stability vs fracture time", fontsize=7, pad=4)
    ax.legend(fontsize=6, markerscale=2)
    clean_ax(ax)

    # ── right: k_eff distribution per condition ────────────────────────────
    ax2 = axes[1]
    data   = [res["k_eff"] for res in results.values()]
    labels = list(results.keys())
    cols   = [COLOR_MAP[lbl] for lbl in labels]

    vp = ax2.violinplot(data, positions=range(len(labels)),
                        showmedians=True, showextrema=False)
    for body, col in zip(vp["bodies"], cols):
        body.set_facecolor(col); body.set_alpha(0.35); body.set_edgecolor(col)
    vp["cmedians"].set_color("k"); vp["cmedians"].set_linewidth(1.2)

    rng = np.random.default_rng(0)
    for xi, (vals, col) in enumerate(zip(data, cols)):
        pts = vals if len(vals) <= 200 else rng.choice(vals, 200, replace=False)
        jitter = rng.uniform(-0.08, 0.08, len(pts))
        ax2.scatter(xi + jitter, pts, s=2, color=col, alpha=0.4, zorder=3)

    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, fontsize=7)
    ax2.set_ylabel(r"$k_\mathrm{eff}(\phi^*)$")
    ax2.set_title(r"(b) Effective growth rate at $\phi^*$", fontsize=7, pad=4)
    clean_ax(ax2)

    meds = {lbl: np.median(res["k_eff"]) for lbl, res in results.items()}
    for lbl, med in meds.items():
        print(f"  {lbl}: k_eff median = {med:.4f}  "
              f"t_frac median = {np.median(results[lbl]['t_frac']):.4f}")

    fig.tight_layout(pad=0.8)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ch-run", default=CH_10K)
    ap.add_argument("--dh-run", default=DH_10K)
    ap.add_argument("--cs-run", default=None)
    ap.add_argument("--ds-run", default=None)
    ap.add_argument("--n-ss",   type=int, default=100,
                    help="# samples for which to compute ODE + Jacobian")
    ap.add_argument("--out",    default="fig_replicon_at2_bridge.pdf")
    args = ap.parse_args()

    load_queue = [("CH", args.ch_run)]
    if args.cs_run:
        load_queue.append(("CS", args.cs_run))
    if args.ds_run:
        load_queue.append(("DS", args.ds_run))
    load_queue.append(("DH", args.dh_run))

    results = {}
    for label, run in load_queue:
        print(f"Loading {label}: {run}")
        s = load_samples(run)
        print(f"  {s.shape[0]} samples  →  analyzing {args.n_ss} for bridge ...")
        results[label] = analyze_bridge(s, n_ss=args.n_ss)

    print(f"\nα_c = {ALPHA_C:.4f}  (G_c={G_C:.1e}, E={E_BIO:.0f}Pa, ℓ={ELL*1e6:.0f}μm)")
    make_figure(results, out=args.out)


if __name__ == "__main__":
    main()

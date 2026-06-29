"""
replicon_species_contrib.py
===========================
Species-level decomposition of the effective fracture-growth rate k_eff.

k_eff(phi*) = sum_i  (phi_i* / sum_j phi_j*)  * K_ALPHA_i

For each TMCMC posterior sample, integrate the gLV ODE to phi*, then compute
how much each of the 5 species contributes to k_eff.

Panels:
  (a) Stacked mean contribution  sum_i c_i   per condition (bar chart).
      Colour-coded by species, width = k_eff median.
  (b) Violin of total k_eff per condition (same as bridge panel (b), repeated
      here for completeness) with species median contribution overlay as text.

Key interpretation:
  - So (K_ALPHA=1.0): high-growth commensal, accelerates fracture.
  - Pg (K_ALPHA=0.3): slow-growth pathobiont, retards fracture.
  DH has higher Pg fraction → k_eff lower → t_frac longer despite being dysbiotic.

Usage:
    python replicon_species_contrib.py
    python replicon_species_contrib.py --n-ss 200 --out fig_species_contrib.pdf
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).parent))

from phase_field_at2_2d_sparse import G_C, ELL, E_BIO
from thesis_style import use, PALETTE, TEXTWIDTH_IN, clean_ax
from replicon_analysis import (
    load_samples, theta_to_matrices, steady_state,
    CH_10K, DH_10K, CS_10K, DS_10K, SPECIES,
)
from replicon_at2_bridge import K_ALPHA, ALPHA_C, k_eff_from_phi

# Species colours: distinct within the teal/orange palette family
SP_COLORS = ["#2a9d8f", "#57cc99", "#a8dadc",   # So, An, Vd  (cool)
             "#e9c46a", "#e76f51"]               # Fn, Pg      (warm)


def analyze_contrib(samples: np.ndarray, n_ss: int = 100):
    """
    Returns:
        contribs : (N, 5) — species contribution to k_eff for each sample.
        phi_stars : (N, 5) — steady-state composition for each sample.
    """
    N = min(n_ss, len(samples))
    contribs  = np.full((N, 5), np.nan)
    phi_stars = np.full((N, 5), np.nan)
    phi0 = np.array([0.35, 0.15, 0.20, 0.20, 0.10])

    for k in range(N):
        A, b = theta_to_matrices(samples[k])
        phi_eq = steady_state(A, b, phi0.copy())
        total = phi_eq.sum()
        if total < 0.01:
            continue
        frac = np.maximum(phi_eq, 0.0) / total
        phi_stars[k] = frac
        contribs[k]  = frac * K_ALPHA   # species contribution to k_eff

    # Drop samples where ODE failed
    mask = ~np.isnan(contribs[:, 0])
    return contribs[mask], phi_stars[mask]


def make_figure(results: dict, out: str = "fig_species_contrib.pdf"):
    """
    results: {label: (contribs (N,5), phi_stars (N,5))}
    """
    use()
    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH_IN, TEXTWIDTH_IN * 0.50))

    labels = list(results.keys())
    n_cond = len(labels)
    COLOR_MAP = {lbl: (PALETTE["ch"] if lbl.upper().startswith("C") else PALETTE["dh"])
                 for lbl in labels}

    # ── panel (a): stacked mean contribution bar ──────────────────────────────
    ax = axes[0]
    x = np.arange(n_cond)
    bar_w = 0.55
    bottoms = np.zeros(n_cond)
    sp_patches = []

    for si, sp in enumerate(SPECIES):
        means = np.array([np.nanmean(results[lbl][0][:, si]) for lbl in labels])
        bars = ax.bar(x, means, bar_w, bottom=bottoms,
                      color=SP_COLORS[si], label=sp, zorder=3)
        bottoms += means
        sp_patches.append(mpatches.Patch(color=SP_COLORS[si],
                                          label=rf"{sp} ($\kappa$={K_ALPHA[si]:.1f})"))

    # Annotate total k_eff (= sum of contributions) above each bar
    for xi, lbl in enumerate(labels):
        total_k = np.nanmean(results[lbl][0].sum(axis=1))
        ax.text(xi, bottoms[xi] + 0.01, f"{total_k:.2f}",
                ha="center", va="bottom", fontsize=6, color="0.3")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel(r"Mean species contribution to $k_\mathrm{eff}$")
    ax.set_title(r"(a) Species decomposition of $k_\mathrm{eff}(\phi^*)$", fontsize=7, pad=4)
    ax.legend(handles=sp_patches, fontsize=5.5, loc="upper right",
              framealpha=0.85, ncol=1)
    clean_ax(ax)

    # ── panel (b): violin of k_eff with median labels ─────────────────────────
    ax2 = axes[1]
    keff_data = [results[lbl][0].sum(axis=1) for lbl in labels]
    cols = [COLOR_MAP[lbl] for lbl in labels]

    vp = ax2.violinplot(keff_data, positions=range(n_cond),
                        showmedians=True, showextrema=False)
    for body, col in zip(vp["bodies"], cols):
        body.set_facecolor(col); body.set_alpha(0.35); body.set_edgecolor(col)
    vp["cmedians"].set_color("k"); vp["cmedians"].set_linewidth(1.2)

    rng = np.random.default_rng(42)
    for xi, (vals, col) in enumerate(zip(keff_data, cols)):
        pts = vals if len(vals) <= 200 else rng.choice(vals, 200, replace=False)
        jitter = rng.uniform(-0.08, 0.08, len(pts))
        ax2.scatter(xi + jitter, pts, s=2, color=col, alpha=0.4, zorder=3)

    # Dominant species annotation per condition
    for xi, lbl in enumerate(labels):
        mean_c = results[lbl][0].mean(axis=0)
        dom_sp = SPECIES[np.argmax(mean_c)]
        dom_frac = mean_c.max() / mean_c.sum() * 100
        med = np.median(keff_data[xi])
        ax2.text(xi, med + 0.03, f"{dom_sp}\n{dom_frac:.0f}%",
                 ha="center", va="bottom", fontsize=5, color="0.3")

    ax2.set_xticks(range(n_cond))
    ax2.set_xticklabels(labels, fontsize=7)
    ax2.set_ylabel(r"$k_\mathrm{eff}(\phi^*)$")
    ax2.set_title(r"(b) $k_\mathrm{eff}$ distribution per condition", fontsize=7, pad=4)
    clean_ax(ax2)

    fig.tight_layout(pad=0.8)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")

    # ── console summary ────────────────────────────────────────────────────────
    print(f"\nSpecies K_ALPHA: " + "  ".join(f"{s}={k:.1f}" for s, k in zip(SPECIES, K_ALPHA)))
    print(f"{'Cond':6s}  {'k_eff med':>10s}  " +
          "  ".join(f"{s:>6s}" for s in SPECIES))
    for lbl in labels:
        c = results[lbl][0]
        med_k = np.median(c.sum(axis=1))
        mean_c = c.mean(axis=0)
        print(f"{lbl:6s}  {med_k:>10.4f}  " +
              "  ".join(f"{v:>6.3f}" for v in mean_c))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ch-run", default=CH_10K)
    ap.add_argument("--dh-run", default=DH_10K)
    ap.add_argument("--cs-run", default=None)
    ap.add_argument("--ds-run", default=None)
    ap.add_argument("--n-ss",   type=int, default=150,
                    help="# samples for ODE integration (default 150)")
    ap.add_argument("--out",    default="fig_species_contrib.pdf")
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
        print(f"  {s.shape[0]} samples → integrating {args.n_ss} ODEs ...")
        results[label] = analyze_contrib(s, n_ss=args.n_ss)
        n_ok = (~np.isnan(results[label][0][:, 0])).sum()
        print(f"  {n_ok} samples converged")

    make_figure(results, out=args.out)


if __name__ == "__main__":
    main()

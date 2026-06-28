"""
replicon_analysis.py
====================
Stability analysis of the 5-species gLV (Hamilton model) via the
disordered-gLV lens (eLife 105948, Pearce et al. 2024/2025).

For each TMCMC posterior sample:
  - Reconstruct A(5,5) and b(5) via theta_to_matrices
  - Compute interaction heterogeneity sigma = std(A_ij)
  - Compute eigenvalues of A (symmetric, so real)
  - Compute Jacobian at fixed point and its stability margin

Compare CH (commensal HOBIC) vs DH (dh_baseline) to test:
    H: DH has lower sigma and/or lambda_max closer to 0 (marginal stability)
    consistent with eLife 105948 (dysbiosis = reduced interaction heterogeneity).

Usage:
    python replicon_analysis.py
    python replicon_analysis.py --out fig_replicon.pdf
"""

import argparse
import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent.parent  # ~/IKM_Hiwi/
RUNS = ROOT / "data_5species" / "_runs"
sys.path.insert(0, str(Path(__file__).parent))
from thesis_style import use, PALETTE, TEXTWIDTH_IN, clean_ax

# ── species ─────────────────────────────────────────────────────────────────
SPECIES = ["So", "An", "Vd", "Fn", "Pg"]
C_CONST = 25.0
K_HILL  = 0.05
N_HILL  = 4.0
T_STEADY = 50.0


def theta_to_matrices(theta: np.ndarray):
    """20-param θ → A(5,5), b(5). Identical layout to attractor_analysis.py."""
    A = np.zeros((5, 5))
    b = np.zeros(5)
    A[0, 0] = theta[0];  A[0, 1] = A[1, 0] = theta[1];  A[1, 1] = theta[2]
    b[0] = theta[3];     b[1] = theta[4]
    A[2, 2] = theta[5];  A[2, 3] = A[3, 2] = theta[6];  A[3, 3] = theta[7]
    b[2] = theta[8];     b[3] = theta[9]
    A[0, 2] = A[2, 0] = theta[10];  A[0, 3] = A[3, 0] = theta[11]
    A[1, 2] = A[2, 1] = theta[12];  A[1, 3] = A[3, 1] = theta[13]
    A[4, 4] = theta[14]; b[4] = theta[15]
    A[0, 4] = A[4, 0] = theta[16];  A[1, 4] = A[4, 1] = theta[17]
    A[2, 4] = A[4, 2] = theta[18];  A[3, 4] = A[4, 3] = theta[19]
    return A, b


def load_samples(run_name: str) -> np.ndarray:
    """Load theta_full posterior samples; shape (N, 20).

    run_name can be:
      - a short name resolved under RUNS (data_5species/_runs/)
      - an absolute path string to a directory containing samples.npy
    """
    run_path = Path(run_name)
    if not run_path.is_absolute():
        run_path = RUNS / run_name
    raw = np.load(run_path / "samples.npy")
    map_path = run_path / "theta_MAP.json"
    with open(map_path) as f:
        d = json.load(f)

    # Two theta_MAP.json layouts:
    #   (A) {"theta_full": [...], "active_indices": [...], ...}
    #   (B) {"0": v0, "1": v1, ..., "19": v19}  — ultimate_10000p format
    if "theta_full" in d:
        n_full = len(d["theta_full"])
        if raw.shape[1] == n_full:
            return raw
        active = np.array(d["active_indices"])
        full = np.zeros((raw.shape[0], n_full))
        full[:, active] = raw[:, :len(active)]
        return full
    else:
        # Layout B: keys are stringified indices; samples already full
        return raw


def ode_rhs(t, phi, A, b):
    hill = K_HILL**N_HILL / (K_HILL**N_HILL + phi**N_HILL)
    return phi * (C_CONST * (A @ phi) - b) * hill


def steady_state(A, b, phi0=None):
    """Integrate ODE to T_STEADY; return final φ."""
    if phi0 is None:
        phi0 = np.ones(5) * 0.2
    sol = solve_ivp(ode_rhs, [0, T_STEADY], phi0, args=(A, b),
                    method="RK45", rtol=1e-6, atol=1e-8, dense_output=False)
    return np.maximum(sol.y[:, -1], 0.0)


def community_jacobian(phi_eq, A, b):
    """
    Linearize ODE around fixed point phi_eq.
    J_ij = d(dφ_i/dt)/dφ_j at φ = φ_eq.
    (hill factor approximated as constant at fixed point value)
    """
    hill = K_HILL**N_HILL / (K_HILL**N_HILL + phi_eq**N_HILL)
    # Full derivative including hill factor is complex; use simplified form:
    # dφ_i/dt ≈ φ_i * (C * Σ A_ij φ_j - b_i) * h_i
    # J_ii = (C * Σ A_ij φ_j* - b_i) * h_i + φ_i * C * A_ii * h_i  (at eq: first term ≈ 0)
    J = np.zeros((5, 5))
    net = C_CONST * (A @ phi_eq) - b  # ≈ 0 at true equilibrium
    for i in range(5):
        for j in range(5):
            J[i, j] = phi_eq[i] * C_CONST * A[i, j] * hill[i]
        J[i, i] += net[i] * hill[i]
    return J


def analyze_samples(samples: np.ndarray, n_ss: int = 50):
    """
    For each posterior sample compute stability metrics.
    Returns dict of arrays (shape N or (N, 5)).
    n_ss: number of samples for which to compute steady-state (slow).
    """
    N = len(samples)
    sigma     = np.zeros(N)   # std of A_ij elements
    lam_max_A = np.zeros(N)   # largest eigenvalue of A
    lam_max_J = np.full(N, np.nan)  # largest real eigenvalue of Jacobian

    phi0 = np.array([0.35, 0.15, 0.20, 0.20, 0.10])  # typical CH-like init

    for k, theta in enumerate(samples):
        A, b = theta_to_matrices(theta)
        sigma[k]     = np.std(A)
        lam_max_A[k] = np.max(np.linalg.eigvalsh(A))  # A symmetric → real eigs

        if k < n_ss:
            phi_eq = steady_state(A, b, phi0.copy())
            if phi_eq.sum() > 0.01:
                J = community_jacobian(phi_eq, A, b)
                lam_max_J[k] = np.max(np.real(np.linalg.eigvals(J)))

    return {"sigma": sigma, "lam_max_A": lam_max_A, "lam_max_J": lam_max_J}


def make_figure(results_dict: dict, out="fig_replicon.pdf"):
    """
    results_dict: {label: res_dict} where label is e.g. "CH-300p", "DH-10k"
    Colors follow PALETTE["ch"]/["dh"]; alpha lighter for 10k.
    """
    use()
    fig, axes = plt.subplots(1, 3, figsize=(TEXTWIDTH_IN, TEXTWIDTH_IN * 0.44))

    panels = [
        ("sigma",     r"$\sigma(A_{ij})$",          r"(a) Interaction heterogeneity"),
        ("lam_max_A", r"$\lambda_\mathrm{max}(A)$",  r"(b) Largest eigenvalue of $A$"),
        ("lam_max_J", r"$\lambda_\mathrm{max}(\mathrm{Re}\,J)$", r"(c) Community Jacobian"),
    ]

    labels = list(results_dict.keys())
    cols = []
    alphas = []
    for lbl in labels:
        base = "ch" if lbl.upper().startswith("CH") else "dh"
        cols.append(PALETTE[base])
        alphas.append(0.25 if "10k" in lbl else 0.45)

    for ax, (key, ylabel, title) in zip(axes, panels):
        data = []
        for lbl, res in results_dict.items():
            v = res[key]
            if key == "lam_max_J":
                v = v[~np.isnan(v)]
            data.append(v)

        positions = list(range(len(labels)))
        vp = ax.violinplot(data, positions=positions,
                           showmedians=True, showextrema=False)
        for body, col, alp in zip(vp["bodies"], cols, alphas):
            body.set_facecolor(col)
            body.set_alpha(alp + 0.15)
            body.set_edgecolor(col)
        vp["cmedians"].set_color("k")
        vp["cmedians"].set_linewidth(1.2)

        # jitter scatter (subsample 10k to 200 pts)
        for xi, (vals, col) in enumerate(zip(data, cols)):
            rng = np.random.default_rng(xi)
            pts = vals if len(vals) <= 300 else rng.choice(vals, 200, replace=False)
            jitter = rng.uniform(-0.08, 0.08, len(pts))
            ax.scatter(xi + jitter, pts, s=2, color=col, alpha=0.35, zorder=3)

        if key in ("lam_max_A", "lam_max_J"):
            ax.axhline(0, color="0.5", lw=0.8, ls="--")

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=6, rotation=15, ha="right")
        ax.set_ylabel(ylabel, fontsize=7)
        ax.set_title(title, fontsize=7, pad=4)
        clean_ax(ax)

        meds = [f"{lbl}={np.median(v):.3f}" for lbl, v in zip(labels, data)]
        print(f"{key}: " + "  ".join(meds))

    fig.tight_layout(pad=0.8)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")


ULTIMATE_BASE = Path.home() / "IKM_Hiwi" / "nife" / "results" / "ultimate_10000p"

# Canonical paths: ultimate_10000p is the converged standard (10 000 TMCMC samples).
# dh_baseline/commensal_hobic_posterior (300 samples) were prior-dominated and
# should not be used for inference — kept only for comparison via --compare-300p.
CH_10K = str(ULTIMATE_BASE / "commensal_hobic")
DH_10K = str(ULTIMATE_BASE / "dh_baseline")
CH_300P = "commensal_hobic_posterior"
DH_300P = "dh_baseline"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ch-run", default=CH_10K,
                    help="CH run dir (default: ultimate_10000p/commensal_hobic)")
    ap.add_argument("--dh-run", default=DH_10K,
                    help="DH run dir (default: ultimate_10000p/dh_baseline)")
    ap.add_argument("--compare-300p", action="store_true",
                    help="Also load 300-sample posteriors for comparison")
    ap.add_argument("--n-ss",   type=int, default=50,
                    help="# samples for which to compute Jacobian (slow)")
    ap.add_argument("--out",    default="fig_replicon.pdf")
    args = ap.parse_args()

    results = {}

    print(f"Loading CH: {args.ch_run}")
    ch = load_samples(args.ch_run)
    print(f"  {ch.shape[0]} samples")
    print("Analyzing CH ...")
    results["CH"] = analyze_samples(ch, n_ss=args.n_ss)

    print(f"Loading DH: {args.dh_run}")
    dh = load_samples(args.dh_run)
    print(f"  {dh.shape[0]} samples")
    print("Analyzing DH ...")
    results["DH"] = analyze_samples(dh, n_ss=args.n_ss)

    if args.compare_300p:
        print(f"Loading CH-300p: {CH_300P}")
        ch300 = load_samples(CH_300P)
        print(f"  {ch300.shape[0]} samples")
        print("Analyzing CH-300p ...")
        results["CH-300p"] = analyze_samples(ch300, n_ss=args.n_ss)

        print(f"Loading DH-300p: {DH_300P}")
        dh300 = load_samples(DH_300P)
        print(f"  {dh300.shape[0]} samples")
        print("Analyzing DH-300p ...")
        results["DH-300p"] = analyze_samples(dh300, n_ss=args.n_ss)

    make_figure(results, out=args.out)


if __name__ == "__main__":
    main()

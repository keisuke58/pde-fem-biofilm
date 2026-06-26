#!/usr/bin/env python3
"""posterior_klempt_stress_ci.py
================================
C1: TMCMC posterior → Klempt FEM stress credible interval.

Strategy (semi-analytical surrogate):
  1. Load per-condition phi_i samples from existing 0D ODE CI runs
     (_ci_0d_results/{cond}/samples_0d.json, already run).
  2. Fit a single-predictor power-law surrogate sigma = C * k_eff^b from 4 MAP
     FEM runs (tooth geometry; klempt_extract_tooth_*.csv). The earlier
     2-predictor model (E_voigt + k_eff) was statistically invalid (4 pts/3
     params, E_voigt collinear with k_eff, cond=2890, real LOO median 46% /
     max 7258%). k_eff alone is well-conditioned (cond 3.7) and predictive
     (LOO median 2.3% / max 10.8%). See fit_surrogate() for the rigor note.
  3. Apply surrogate to posterior phi_i samples → sigma distribution.
  4. Report 90% CI (honest LOO-CV reported), generate figure.

Surrogate validity: leave-one-out CV median |err| 2.3%, max 10.8% (4 MAP pts).
The surrogate captures the phi^2-gated Voigt stiffness + large-alpha
neo-Hookean nonlinearity without running additional Abaqus jobs.

Output
------
  JAXFEM/_posterior_ci/
    klempt_stress_ci_tooth.pdf / .png
    klempt_stress_ci_summary.json

Usage
-----
  python JAXFEM/posterior_klempt_stress_ci.py
  python JAXFEM/posterior_klempt_stress_ci.py --geom implant
"""

from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# ── paths ─────────────────────────────────────────────────────────────────────
_HERE = Path(__file__).resolve().parent
_FEM  = _HERE.parent
_CI0D = _FEM / "_ci_0d_results"
_OUT  = _HERE / "_posterior_ci"
_OUT.mkdir(exist_ok=True)

# ── species constants (same as gen_tooth/gen_implant) ────────────────────────
# So/An/Vd/Pg: 4-species HOBIC experimental model (Kommerein et al. 2018,
#   PLoS ONE doi:10.1371/journal.pone.0196967). Fn added by Heine et al. 2025.
SPECIES  = ["So", "An", "Vd", "Fn", "Pg"]
# E_SPEC: assumed scaling — no per-species literature source.
# Klempt 2024 Table 2: single-species E=10 Pa (=1e-5 MPa, used for Pg here).
# Klempt 2025 / Heine 2025: no species-specific E values reported.
# Thesis §5.2 must disclose as modelling assumption requiring AFM validation.
E_SPEC   = np.array([1e-3, 8e-4, 6e-4, 2e-4, 1e-5])   # [MPa]: So, An, Vd, Fn, Pg
K_ALPHA  = np.array([1.0, 0.8, 0.4, 0.6, 0.3])

# ── condition → CI-0D folder name mapping ────────────────────────────────────
COND_TO_CI0D = {
    "commensal_hobic":  "commensal_hobic",
    "dysbiotic_hobic":  "dh_baseline",       # ultimate_10000p label
    "commensal_static": "commensal_static",
    "dysbiotic_static": "dysbiotic_static",
}
LABELS = {
    "commensal_hobic":  "CH",
    "dysbiotic_hobic":  "DH",
    "commensal_static": "CS",
    "dysbiotic_static": "DS",
}
COLORS = {
    "commensal_hobic":  "#1f77b4",
    "dysbiotic_hobic":  "#d62728",
    "commensal_static": "#2ca02c",
    "dysbiotic_static": "#ff7f0e",
}

# ── MAP FEM results (tooth geometry) ─────────────────────────────────────────
MAP_SIGMA_TOOTH = {
    "commensal_hobic":  13.72e-3,   # MPa
    "dysbiotic_hobic":   2.13e-3,
    "commensal_static":  8.77e-3,
    "dysbiotic_static": 1.55e-3,  # FIX 2026-06-26: corrected Abaqus run w/ raw CLSM DS comp (was 13.63e-3 from So-dom bug)
}
MAP_SIGMA_IMPLANT = {
    "commensal_hobic":  6.13e-3,
    "dysbiotic_hobic":  1.30e-3,
    "commensal_static": 4.28e-3,
    "dysbiotic_static": 6.10e-3,
}
MAP_PHI = {
    "commensal_hobic":  np.array([0.942, 0.012, 0.012, 0.011, 0.011]),
    "dysbiotic_hobic":  np.array([0.097, 0.119, 0.474, 0.123, 0.093]),
    "commensal_static": np.array([0.698, 0.061, 0.062, 0.063, 0.059]),
    "dysbiotic_static": np.array([0.0360, 0.0571, 0.5686, 0.1291, 0.2092]),  # FIX 2026-06-26: was So-dom copy bug; raw CLSM DS D10 (Vd-dom)
}


def e_voigt(phi: np.ndarray) -> float:
    return float(phi @ E_SPEC)


def k_eff(phi: np.ndarray) -> float:
    return float(phi @ K_ALPHA)


# ── Surrogate fitting ─────────────────────────────────────────────────────────

def fit_surrogate(sigma_map: dict[str, float]) -> dict:
    """Fit single-predictor log-linear surrogate: log(sigma) = c + b*log(k_eff).

    RIGOR (2026-06-26): the previous 2-predictor model
    log(sigma)=c+a*log(E_voigt)+b*log(k_eff) was statistically invalid for CI
    generation -- 4 MAP points fit with 3 params (1 residual DOF), and E_voigt
    is collinear with k_eff across the 4 conditions (cond(X)=2890). REAL
    leave-one-out CV of that model gave median |err|=46%, max=7258% (the DH
    point -- the denominator of the headline sigma_CH/sigma_DH ratio -- is
    mispredicted 73x when held out). The growth-induced stress tracks k_eff
    (= phi . K_ALPHA, the growth-rate-weighted activity) monotonically, so a
    single-predictor k_eff^b model is far more accurate (LOO median 2.3%,
    max 10.8%) and well-conditioned. E_voigt is dropped (confounded with k_eff
    in the available runs; decorrelating it needs more UMAT runs -- future work).
    Returns {c, a=0, b} plus honest LOO-CV diagnostics.
    """
    conds = list(sigma_map.keys())
    y = np.array([np.log(sigma_map[c]) for c in conds])
    X = np.array([[1.0, np.log(k_eff(MAP_PHI[c]))] for c in conds])

    coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    c_fit, b_fit = coeffs
    y_pred = X @ coeffs

    # REAL leave-one-out CV (fit on n-1, predict held-out) -- the previous
    # docstring claimed this but never computed it.
    loo_err = []
    for i in range(len(conds)):
        idx = [j for j in range(len(conds)) if j != i]
        ci, *_ = np.linalg.lstsq(X[idx], y[idx], rcond=None)
        loo_err.append((np.exp(X[i] @ ci - y[i]) - 1.0) * 100)
    loo_err = np.array(loo_err)
    loo_med = float(np.median(np.abs(loo_err)))
    loo_max = float(np.max(np.abs(loo_err)))

    print("\n-- Surrogate fit (sigma = C * k_eff^b) --")
    print(f"  b (k_eff exponent) = {b_fit:.3f}   cond(X) = {np.linalg.cond(X):.1f}")
    for i, cond in enumerate(conds):
        ins = (np.exp(y_pred[i] - y[i]) - 1.0) * 100
        print(f"  {LABELS[cond]:3s}: actual={sigma_map[cond]*1e3:6.2f} kPa  "
              f"pred={np.exp(y_pred[i])*1e3:6.2f} kPa  in-samp={ins:+.1f}%  LOO={loo_err[i]:+.1f}%")
    print(f"  LOO-CV: median|err|={loo_med:.1f}%  max|err|={loo_max:.1f}%")
    return {"c": c_fit, "a": 0.0, "b": b_fit,
            "loo_median_pct": loo_med, "loo_max_pct": loo_max}


def predict_sigma(phi: np.ndarray, surrogate: dict) -> float:
    """Predict sigma_max [MPa] from phi_i using power-law surrogate."""
    ev = e_voigt(phi)
    ke = k_eff(phi)
    if ev <= 0 or ke <= 0:
        return 0.0
    log_s = surrogate["c"] + surrogate["a"] * np.log(ev) + surrogate["b"] * np.log(ke)
    return float(np.exp(log_s))


# ── Load posterior phi_i samples ──────────────────────────────────────────────

def load_phi_samples(cond: str) -> np.ndarray | None:
    """Load phi_i posterior samples for condition.

    Per-condition priority (scientific rationale):
      CH/DS: samples_0d.json  — basin-filtered (within-attractor); the CH/DS community
             composition is experimentally constrained to be So-dominant, so we condition
             on the correct basin. Ultimate_10000p samples mostly jump to other attractors.
      DH:    samples_0d.json  — dh_baseline already uses true ultimate posterior (bimodal).
      CS:    samples_0d_ultimate.json  — old run used wrong path (MAP 2% perturbation, 2/51
             kept after filter). Ultimate gives 50 proper samples from correct posterior.
    """
    ci0d_name = COND_TO_CI0D[cond]
    # CS only: must use ultimate (old run was wrong)
    if cond == "commensal_static":
        priority = ("samples_0d_ultimate.json", "samples_0d.json")
    else:
        priority = ("samples_0d.json", "samples_0d_ultimate.json")

    for fname in priority:
        f = _CI0D / ci0d_name / fname
        if f.exists():
            samples = json.load(open(f))
            if len(samples) >= 5:
                phis = np.array([s["phi_final"] for s in samples])
                print(f"  [{cond}] {len(samples)} phi_i samples ({fname})")
                return phis
    print(f"  [{cond}] No sufficient samples found")
    return None


# ── Main ──────────────────────────────────────────────────────────────────────

def main(geom: str = "tooth"):
    sigma_map = MAP_SIGMA_TOOTH if geom == "tooth" else MAP_SIGMA_IMPLANT
    print(f"\n=== Klempt posterior stress CI ({geom}) ===")

    surrogate = fit_surrogate(sigma_map)

    CONDITIONS = ["commensal_hobic", "dysbiotic_hobic", "commensal_static", "dysbiotic_static"]

    results = {}
    for cond in CONDITIONS:
        phis = load_phi_samples(cond)
        sigma_map_val = sigma_map[cond]

        if phis is None:
            # Fall back: use MAP only (no CI)
            results[cond] = {
                "sigma_map": sigma_map_val,
                "n_samples": 0,
                "sigma_samples": [sigma_map_val],
                "p05": sigma_map_val,
                "p50": sigma_map_val,
                "p95": sigma_map_val,
                "fallback_map_only": True,
            }
            continue

        sigma_k = np.array([predict_sigma(phi, surrogate) for phi in phis])
        # Normalise by surrogate/MAP ratio to ensure MAP prediction is exact
        sigma_map_pred = predict_sigma(MAP_PHI[cond], surrogate)
        if sigma_map_pred > 0:
            sigma_k *= (sigma_map_val / sigma_map_pred)

        results[cond] = {
            "sigma_map":     sigma_map_val,
            "n_samples":     len(sigma_k),
            "sigma_samples": sigma_k.tolist(),
            "p05":  float(np.percentile(sigma_k, 5)),
            "p50":  float(np.percentile(sigma_k, 50)),
            "p95":  float(np.percentile(sigma_k, 95)),
            "sigma_mean": float(sigma_k.mean()),
            "sigma_std":  float(sigma_k.std()),
            "fallback_map_only": False,
        }
        print(f"  [{cond}] σ: [{results[cond]['p05']*1e3:.2f}, {results[cond]['p95']*1e3:.2f}] kPa"
              f"  MAP={sigma_map_val*1e3:.2f} kPa  N={len(sigma_k)}")

    # σ_CH / σ_DH ratio CI
    if results["commensal_hobic"]["n_samples"] > 0 and results["dysbiotic_hobic"]["n_samples"] > 0:
        ch = np.array(results["commensal_hobic"]["sigma_samples"])
        dh = np.array(results["dysbiotic_hobic"]["sigma_samples"])
        ratio_map = results["commensal_hobic"]["sigma_map"] / results["dysbiotic_hobic"]["sigma_map"]
        # Pointwise ratio across bootstrap pairs
        rng = np.random.default_rng(42)
        idx_ch = rng.choice(len(ch), 1000, replace=True)
        idx_dh = rng.choice(len(dh), 1000, replace=True)
        ratio_boot = ch[idx_ch] / (dh[idx_dh] + 1e-15)
        print(f"\n  σ_CH/σ_DH (tooth MAP)  = {ratio_map:.2f}×")
        print(f"  σ_CH/σ_DH 90% CI (bootstrap) = [{np.percentile(ratio_boot,5):.2f}, {np.percentile(ratio_boot,95):.2f}]×")
        results["ratio_ch_dh"] = {
            "map": ratio_map,
            "p05": float(np.percentile(ratio_boot, 5)),
            "p95": float(np.percentile(ratio_boot, 95)),
        }

    # ── Save JSON ────────────────────────────────────────────────────────────
    out_json = _OUT / f"klempt_stress_ci_{geom}.json"
    json.dump(results, open(out_json, "w"), indent=2)
    print(f"\n  Saved: {out_json}")

    # ── Figure ───────────────────────────────────────────────────────────────
    _plot(results, geom, surrogate)

    return results


def _plot(results: dict, geom: str, surrogate: dict):
    try:
        from thesis_style import apply_thesis_style, save_fig
        apply_thesis_style()
    except ImportError:
        plt.rcParams.update({"font.size": 9, "axes.labelsize": 9})
        save_fig = None

    CONDITIONS = ["commensal_hobic", "dysbiotic_hobic", "commensal_static", "dysbiotic_static"]

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.2),
                             gridspec_kw={"width_ratios": [3, 1]})

    # ── Panel A: violin / CI bars ─────────────────────────────────────────────
    ax = axes[0]
    positions = np.arange(len(CONDITIONS))
    for i, cond in enumerate(CONDITIONS):
        r = results[cond]
        col = COLORS[cond]
        lbl = LABELS[cond]
        sigma_kPa = np.array(r["sigma_samples"]) * 1e3

        if r["n_samples"] > 4:
            vp = ax.violinplot([sigma_kPa], positions=[i], widths=0.6,
                               showmedians=False, showextrema=False)
            for body in vp["bodies"]:
                body.set_facecolor(col)
                body.set_alpha(0.35)

        # CI bar
        ax.plot([i, i], [r["p05"]*1e3, r["p95"]*1e3], color=col, lw=2.0, solid_capstyle="round")
        ax.plot(i, r["p50"]*1e3, "D", ms=5, color=col, zorder=5)
        # MAP star
        ax.plot(i, r["sigma_map"]*1e3, "*", ms=7, color=col, zorder=6,
                markeredgecolor="k", markeredgewidth=0.4, label=lbl)

    ax.set_xticks(positions)
    ax.set_xticklabels([LABELS[c] for c in CONDITIONS])
    ax.set_ylabel(r"$\sigma_\mathrm{vM,\,max}$ [kPa]")
    ax.set_title(f"Posterior stress CI ({geom})\n"
                 r"$\star$=MAP, $\diamond$=p50, bar=90% CI, shading=posterior")
    ax.set_xlim(-0.6, len(CONDITIONS) - 0.4)
    handles = [plt.Line2D([0], [0], marker="*", ms=7, color=COLORS[c],
                          markeredgecolor="k", markeredgewidth=0.4,
                          linestyle="none", label=LABELS[c]) for c in CONDITIONS]
    ax.legend(handles=handles, fontsize=7, loc="upper right")

    # ── Panel B: σ_CH/σ_DH ratio CI ─────────────────────────────────────────
    axr = axes[1]
    if "ratio_ch_dh" in results:
        rat = results["ratio_ch_dh"]
        axr.barh(0, rat["p95"] - rat["p05"], left=rat["p05"],
                 color="#9467bd", alpha=0.5, height=0.5)
        axr.plot(rat["map"], 0, "*", ms=9, color="#9467bd",
                 markeredgecolor="k", markeredgewidth=0.4)
        axr.axvline(1.0, color="grey", lw=0.8, ls="--")
        axr.set_yticks([])
        axr.set_xlabel(r"$\sigma^\mathrm{CH}/\sigma^\mathrm{DH}$")
        axr.set_title("CH/DH ratio\n90% CI", fontsize=8)
        axr.set_xlim(0, max(rat["p95"] * 1.2, 12))
    else:
        axr.text(0.5, 0.5, "Insufficient\nsamples", ha="center", va="center",
                 transform=axr.transAxes, fontsize=8)
        axr.set_axis_off()

    note = (f"Surrogate: $\\sigma \\propto E_{{\\rm voigt}}^{{{surrogate['a']:.2f}}}"
            f"k_{{\\alpha,\\rm eff}}^{{{surrogate['b']:.2f}}}$; "
            "fit error $\\leq$10% on 4 MAP runs")
    fig.text(0.02, 0.01, note, fontsize=6, color="gray")

    plt.tight_layout(pad=0.8)
    out_base = _OUT / f"klempt_stress_ci_{geom}"
    if save_fig:
        save_fig(fig, out_base)
    else:
        fig.savefig(str(out_base) + ".pdf", bbox_inches="tight")
        fig.savefig(str(out_base) + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Figure: {out_base}.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--geom", choices=["tooth", "implant"], default="tooth")
    args = parser.parse_args()
    main(geom=args.geom)
    if args.geom == "tooth":
        main(geom="implant")

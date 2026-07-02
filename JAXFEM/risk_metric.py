#!/usr/bin/env python3
"""risk_metric.py
================
T2.3 — Clinical risk metric  P[sigma > threshold]  from the TMCMC posterior.

This is the terminal stage of the mechano-ecological pipeline
(posterior -> PDE alpha -> FEM -> **risk**). It turns the per-condition
posterior stress distribution produced by ``posterior_klempt_stress_ci.py``
(``_posterior_ci/klempt_stress_ci_{geom}.json``, key ``sigma_samples``) into a
single clinically interpretable number: the probability that the
growth-induced peak von Mises stress exceeds a detachment-relevant threshold.

Why an exceedance probability
-----------------------------
The headline ``sigma_CH/sigma_DH`` ratio is a point/central statistic. For a
risk read-out we instead ask, per condition, *what fraction of the posterior
puts the biofilm above a mechanically meaningful stress level* -- i.e.
``P[sigma_max > tau]``. Because the choice of ``tau`` is not itself well
established (biofilm cohesive/detachment strength spans a wide band), we do NOT
report a single threshold in isolation: the survival function ``P[sigma > tau]``
is reported over a sweep, and point read-outs are given at several reference
thresholds. The exceedance probability carries a bootstrap CI, since the
per-condition posterior sample counts are small (n ~ 22-51).

Units
-----
Stresses in the CI json are Abaqus MPa (mm/N/MPa). Thresholds are accepted and
reported in kPa (clinically more legible); 1 kPa = 1e-3 MPa.

Outputs
-------
  JAXFEM/_risk/
    risk_summary_{geom}.json     (tracked; numbers)
    risk_survival_{geom}.png     (gitignored; survival curves)
    risk_bars_{geom}.png         (gitignored; per-condition risk at tau)

Usage
-----
  python JAXFEM/risk_metric.py                       # tooth, default thresholds
  python JAXFEM/risk_metric.py --geom implant
  python JAXFEM/risk_metric.py --threshold-kpa 5.0   # headline tau
  python JAXFEM/risk_metric.py --no-plot             # numbers only
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_CI_DIR = _HERE / "_posterior_ci"
_OUT = _HERE / "_risk"

MPA_PER_KPA = 1.0e-3

# Canonical condition order + display labels (matches the CI json / report).
CONDITIONS = ["commensal_hobic", "commensal_static", "dysbiotic_static", "dysbiotic_hobic"]
COND_LABEL = {
    "commensal_hobic": "CH",
    "commensal_static": "CS",
    "dysbiotic_static": "DS",
    "dysbiotic_hobic": "DH",
}

# Reference thresholds (kPa) at which point read-outs are always reported, so
# the summary is informative regardless of the single --threshold-kpa choice.
REFERENCE_THRESHOLDS_KPA = [2.5, 5.0, 10.0]


def exceedance_prob(samples_mpa: np.ndarray, tau_mpa: float) -> float:
    """Empirical P[sigma > tau] from posterior samples (strict exceedance)."""
    samples_mpa = np.asarray(samples_mpa, dtype=float)
    if samples_mpa.size == 0:
        return float("nan")
    return float(np.mean(samples_mpa > tau_mpa))


def bootstrap_ci(
    samples_mpa: np.ndarray,
    tau_mpa: float,
    n_boot: int = 4000,
    alpha: float = 0.10,
    seed: int = 0,
) -> tuple[float, float]:
    """Bootstrap (1-alpha) CI for the exceedance probability.

    Resamples the posterior stress samples with replacement. Returns
    (low, high) percentile bounds. With degenerate (near-constant) samples the
    interval collapses to a point, which is the honest answer.
    """
    samples_mpa = np.asarray(samples_mpa, dtype=float)
    n = samples_mpa.size
    if n == 0:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = np.mean(samples_mpa[idx] > tau_mpa, axis=1)
    lo = float(np.quantile(boot, alpha / 2.0))
    hi = float(np.quantile(boot, 1.0 - alpha / 2.0))
    return (lo, hi)


def survival_curve(samples_mpa: np.ndarray, taus_mpa: np.ndarray) -> np.ndarray:
    """Vectorised survival function P[sigma > tau] over a threshold grid."""
    samples_mpa = np.asarray(samples_mpa, dtype=float)
    taus_mpa = np.asarray(taus_mpa, dtype=float)
    if samples_mpa.size == 0:
        return np.full_like(taus_mpa, np.nan, dtype=float)
    # (n_tau, n_samp) comparison -> mean over samples
    return np.mean(samples_mpa[None, :] > taus_mpa[:, None], axis=1)


def load_ci(geom: str) -> dict:
    """Load the committed posterior CI json for a geometry."""
    path = _CI_DIR / f"klempt_stress_ci_{geom}.json"
    if not path.exists():
        raise FileNotFoundError(
            f"CI json not found: {path}\n"
            f"Run:  python JAXFEM/posterior_klempt_stress_ci.py --geom {geom}"
        )
    with open(path) as f:
        return json.load(f)


def _condition_samples(entry: dict) -> np.ndarray:
    """Extract stress samples (MPa) from a CI-json condition entry.

    Falls back to the MAP value as a single sample if no posterior samples are
    present (fallback_map_only conditions), so the metric still returns a value.
    """
    samples = entry.get("sigma_samples")
    if samples:
        return np.asarray(samples, dtype=float)
    if "sigma_map" in entry:
        return np.asarray([entry["sigma_map"]], dtype=float)
    return np.asarray([], dtype=float)


def compute_risk(
    ci: dict,
    threshold_kpa: float,
    reference_kpa: list[float] | None = None,
    seed: int = 0,
) -> dict:
    """Compute the risk summary for all conditions present in a CI dict.

    Returns a JSON-serialisable dict:
      { headline_threshold_kpa, conditions: {cond: {...}}, ... }
    """
    reference_kpa = list(reference_kpa or REFERENCE_THRESHOLDS_KPA)
    # Ensure the headline threshold is among the reported reference points.
    ref_set = sorted({round(t, 6) for t in reference_kpa} | {round(threshold_kpa, 6)})

    tau_mpa = threshold_kpa * MPA_PER_KPA
    out_conditions: dict[str, dict] = {}

    for cond in [c for c in CONDITIONS if c in ci]:
        entry = ci[cond]
        samples = _condition_samples(entry)
        n = int(samples.size)
        degenerate = bool(n <= 1 or np.std(samples) < 1e-9)

        p = exceedance_prob(samples, tau_mpa)
        lo, hi = bootstrap_ci(samples, tau_mpa, seed=seed)

        ref_points = {}
        for t_kpa in ref_set:
            ref_points[f"{t_kpa:g}"] = round(
                exceedance_prob(samples, t_kpa * MPA_PER_KPA), 4
            )

        out_conditions[cond] = {
            "label": COND_LABEL.get(cond, cond),
            "n_samples": n,
            "degenerate": degenerate,
            "sigma_map_kpa": round(entry.get("sigma_map", float("nan")) / MPA_PER_KPA, 4),
            "sigma_mean_kpa": round(float(np.mean(samples)) / MPA_PER_KPA, 4),
            "sigma_p50_kpa": round(float(np.median(samples)) / MPA_PER_KPA, 4),
            "sigma_p95_kpa": round(float(np.quantile(samples, 0.95)) / MPA_PER_KPA, 4),
            "risk": round(p, 4),
            "risk_ci90": [round(lo, 4), round(hi, 4)],
            "risk_by_threshold_kpa": ref_points,
        }

    return {
        "headline_threshold_kpa": threshold_kpa,
        "reference_thresholds_kpa": ref_set,
        "metric": "P[sigma_max > threshold]  (empirical posterior exceedance)",
        "note": (
            "Threshold is a modelling choice, not an established constant; the "
            "survival curve (see risk_survival_*.png) is the primary read-out. "
            "Degenerate conditions have near-constant posterior stress "
            "(fallback / collapsed CI) and give a step-function risk."
        ),
        "conditions": out_conditions,
    }


def plot_risk(ci: dict, summary: dict, geom: str, out_dir: Path) -> list[Path]:
    """Survival curves + per-condition risk bars. Returns written paths."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    conds = [c for c in CONDITIONS if c in ci]
    colors = {
        "commensal_hobic": "#2a9d8f",
        "commensal_static": "#8ab17d",
        "dysbiotic_static": "#e9c46a",
        "dysbiotic_hobic": "#e76f51",
    }

    # ── survival curves ──────────────────────────────────────────────
    all_samples = np.concatenate(
        [_condition_samples(ci[c]) for c in conds if _condition_samples(ci[c]).size]
    )
    tau_max_kpa = float(np.quantile(all_samples, 0.999)) / MPA_PER_KPA * 1.05
    taus_kpa = np.linspace(0.0, max(tau_max_kpa, 1.0), 300)
    taus_mpa = taus_kpa * MPA_PER_KPA

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for c in conds:
        s = _condition_samples(ci[c])
        surv = survival_curve(s, taus_mpa)
        ax.plot(taus_kpa, surv, lw=2.0, color=colors.get(c, None), label=COND_LABEL.get(c, c))
    thr = summary["headline_threshold_kpa"]
    ax.axvline(thr, color="0.4", ls="--", lw=1.0)
    ax.text(thr, 1.02, f"τ={thr:g} kPa", ha="center", va="bottom", fontsize=8, color="0.3")
    ax.set_xlabel("threshold  τ  (kPa)")
    ax.set_ylabel("P[σ$_{max}$ > τ]")
    ax.set_ylim(-0.02, 1.08)
    ax.set_title(f"Growth-stress exceedance risk — {geom}")
    ax.legend(frameon=False, title="condition")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    p1 = out_dir / f"risk_survival_{geom}.png"
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    written.append(p1)

    # ── per-condition risk bars at headline threshold ────────────────
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    labels = [COND_LABEL.get(c, c) for c in conds]
    risks = [summary["conditions"][c]["risk"] for c in conds]
    cis = [summary["conditions"][c]["risk_ci90"] for c in conds]
    yerr = np.array(
        [[r - lo for r, (lo, hi) in zip(risks, cis)], [hi - r for r, (lo, hi) in zip(risks, cis)]]
    )
    yerr = np.clip(yerr, 0.0, None)
    bars = ax.bar(labels, risks, color=[colors.get(c, None) for c in conds], alpha=0.9)
    ax.errorbar(range(len(conds)), risks, yerr=yerr, fmt="none", ecolor="0.25", capsize=4, lw=1.2)
    for b, r in zip(bars, risks):
        ax.text(b.get_x() + b.get_width() / 2, min(r + 0.03, 1.0), f"{r:.2f}", ha="center", fontsize=9)
    ax.set_ylabel(f"P[σ$_{{max}}$ > {thr:g} kPa]")
    ax.set_ylim(0, 1.12)
    ax.set_title(f"Detachment-stress risk by condition — {geom}\n(90% bootstrap CI)")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    p2 = out_dir / f"risk_bars_{geom}.png"
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    written.append(p2)

    return written


def _print_table(summary: dict, geom: str) -> None:
    thr = summary["headline_threshold_kpa"]
    print(f"\nRisk metric — geometry: {geom}   headline τ = {thr:g} kPa")
    print(f"  metric: {summary['metric']}")
    print(f"  {'cond':<5} {'n':>3} {'σ_map':>8} {'σ_p95':>8} {'P[σ>τ]':>8}  {'90% CI':>16}  flag")
    print(f"  {'-'*62}")
    for cond, d in summary["conditions"].items():
        lo, hi = d["risk_ci90"]
        flag = "degenerate" if d["degenerate"] else ""
        print(
            f"  {d['label']:<5} {d['n_samples']:>3} {d['sigma_map_kpa']:>8.3f} "
            f"{d['sigma_p95_kpa']:>8.3f} {d['risk']:>8.3f}  [{lo:>5.2f},{hi:>5.2f}]{'':>4}  {flag}"
        )
    print(f"  {'-'*62}")
    ref = summary["reference_thresholds_kpa"]
    print(f"  reference thresholds (kPa): {ref}")
    print("  P[σ>τ] by threshold:")
    for cond, d in summary["conditions"].items():
        row = "  ".join(f"{t}:{d['risk_by_threshold_kpa'][t]:.2f}" for t in map(lambda x: f"{x:g}", ref))
        print(f"    {d['label']:<5} {row}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="P[sigma > threshold] risk metric from posterior stress")
    ap.add_argument("--geom", default="tooth", choices=["tooth", "implant"])
    ap.add_argument("--threshold-kpa", type=float, default=5.0, help="Headline threshold τ (kPa)")
    ap.add_argument(
        "--reference-kpa",
        type=float,
        nargs="+",
        default=None,
        help=f"Reference thresholds to tabulate (default {REFERENCE_THRESHOLDS_KPA})",
    )
    ap.add_argument("--seed", type=int, default=0, help="Bootstrap RNG seed")
    ap.add_argument("--out-dir", default=None, help="Output dir (default JAXFEM/_risk)")
    ap.add_argument("--no-plot", action="store_true", help="Skip figure generation")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir) if args.out_dir else _OUT
    out_dir.mkdir(parents=True, exist_ok=True)

    ci = load_ci(args.geom)
    summary = compute_risk(
        ci,
        threshold_kpa=args.threshold_kpa,
        reference_kpa=args.reference_kpa,
        seed=args.seed,
    )
    summary["geom"] = args.geom

    summary_path = out_dir / f"risk_summary_{args.geom}.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    _print_table(summary, args.geom)
    print(f"\n  summary → {summary_path}")

    if not args.no_plot:
        try:
            figs = plot_risk(ci, summary, args.geom, out_dir)
            for p in figs:
                print(f"  figure  → {p}")
        except Exception as exc:  # figures are non-essential; never fail the metric
            print(f"  [plot] skipped: {exc}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

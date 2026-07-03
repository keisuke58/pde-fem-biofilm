#!/usr/bin/env python3
"""validate_composition.py — model composition vs Heine measured composition.

Closes the model <-> experiment loop for the one condition where both sides are
available: **dysbiotic, static**.

  model      : the pipeline's 0D posterior composition
               (_ci_0d_results/dysbiotic_static/samples_0d.json, phi_So..phi_Pg)
  experiment : Heine CLSM measured composition, Dysbiotic sheet, "Static all
               cells" (data/heine_species_distribution_biofilm.xlsx)

Both are the same five species in the same order (S. oralis, A. naeslundii,
Veillonella, F. nucleatum, P. gingivalis), so we can compare them directly.
The experiment is time-resolved (days 1..21); the model gives a single mature
composition, so we compare against the Heine time-mean (with the measured
day-to-day range shown) and report:

  * per-species absolute error (percentage points)
  * mean absolute error (MAE, pp)
  * total variation distance (TVD, 0..1)

Outputs a comparison figure and a metrics JSON, both consumed by
tests/test_validation_composition.py as a regression guard.

  python validate_composition.py     # -> assets/validation_composition_dysbiotic.png
                                      #    _validation/composition_metrics.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpyxl

_HERE = Path(__file__).resolve().parent
_MODEL = _HERE / "_ci_0d_results" / "dysbiotic_static" / "samples_0d.json"
_XLSX = _HERE / "data" / "heine_species_distribution_biofilm.xlsx"
_OUT_FIG = _HERE / "assets" / "validation_composition_dysbiotic.png"
_OUT_JSON = _HERE / "_validation" / "composition_metrics.json"

SPECIES = ["S. oralis", "A. naeslundii", "Veillonella", "F. nucleatum", "P. gingivalis"]
MODEL_KEYS = ["phi_So", "phi_An", "phi_Vd", "phi_Fn", "phi_Pg"]
TAGS = [1, 3, 6, 10, 15, 21]
C_MODEL, C_EXP = "#0072B2", "#D55E00"       # two series (Okabe-Ito), identity by series
SURFACE = "#fcfcfb"


def model_composition() -> np.ndarray:
    """(n_samples, 5) posterior composition, normalised to %."""
    s = json.loads(_MODEL.read_text())
    M = np.array([[x[k] for k in MODEL_KEYS] for x in s], float)
    return 100.0 * M / M.sum(1, keepdims=True)


def heine_composition() -> np.ndarray:
    """(len(TAGS), 5) measured dysbiotic-static all-cells composition, %."""
    wb = openpyxl.load_workbook(_XLSX, data_only=True, read_only=True)
    rows = list(wb["Dysbiotic"].iter_rows(values_only=True))
    hdr = next(r for r in rows if r[1] and str(r[1]).startswith("S. oralis"))
    sp = [j for j in range(1, len(hdr)) if hdr[j]]
    width = sp[1] - sp[0]
    cond, blocks = None, {}
    for r in rows:
        v = r[0]
        if isinstance(v, str) and "cells" in v:
            cond = v.strip()
            blocks[cond] = {}
            continue
        if cond and v in TAGS:
            row = []
            for c in sp:
                vals = [r[c + k] for k in range(width)
                        if isinstance(r[c + k], (int, float))]
                row.append(np.mean(vals) if vals else np.nan)
            blocks[cond][v] = row
    A = blocks["Static all cells"]
    mat = np.array([A[t] for t in TAGS], float)
    return 100.0 * mat / mat.sum(1, keepdims=True)


def metrics(model_mean, exp_mean):
    err = model_mean - exp_mean
    return {
        "per_species_abs_error_pp": {SPECIES[i]: float(abs(err[i])) for i in range(5)},
        "mae_pp": float(np.abs(err).mean()),
        "max_abs_error_pp": float(np.abs(err).max()),
        "tvd": float(0.5 * np.abs(err).sum() / 100.0),
        "model_mean_pct": {SPECIES[i]: float(model_mean[i]) for i in range(5)},
        "exp_mean_pct": {SPECIES[i]: float(exp_mean[i]) for i in range(5)},
        "condition": "dysbiotic_static",
        "note": "model = 0D posterior composition; experiment = Heine CLSM "
                "all-cells, static; comparison vs Heine time-mean (days 1-21).",
    }


def main():
    M = model_composition()
    H = heine_composition()
    model_mean = M.mean(0)
    exp_mean = H.mean(0)
    m = metrics(model_mean, exp_mean)

    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 10, "axes.edgecolor": "#888",
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    })
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    x = np.arange(5)
    w = 0.38

    # model bars (posterior spread as whisker; near-degenerate here)
    m_lo, m_hi = np.percentile(M, 5, 0), np.percentile(M, 95, 0)
    ax.bar(x - w / 2, model_mean, width=w, color=C_MODEL, edgecolor=SURFACE,
           linewidth=1.2, label="Model (0D posterior)", zorder=3)
    ax.errorbar(x - w / 2, model_mean, yerr=[model_mean - m_lo, m_hi - model_mean],
                fmt="none", ecolor="#04466e", elinewidth=1.2, capsize=2, zorder=4)
    # experiment bars (mean + measured day-to-day range as whisker)
    e_lo, e_hi = H.min(0), H.max(0)
    ax.bar(x + w / 2, exp_mean, width=w, color=C_EXP, edgecolor=SURFACE,
           linewidth=1.2, label="Heine measured (mean, days 1–21)", zorder=3)
    ax.errorbar(x + w / 2, exp_mean, yerr=[exp_mean - e_lo, e_hi - exp_mean],
                fmt="none", ecolor="#7a3208", elinewidth=1.2, capsize=2, zorder=4)

    for xi in x:
        ax.text(xi - w / 2, model_mean[xi] + 1.2, f"{model_mean[xi]:.0f}",
                ha="center", va="bottom", fontsize=7.5, color=C_MODEL)
        ax.text(xi + w / 2, exp_mean[xi] + 1.2, f"{exp_mean[xi]:.0f}",
                ha="center", va="bottom", fontsize=7.5, color=C_EXP)

    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES, fontsize=9)
    ax.set_ylabel("relative composition (%)", fontsize=10)
    ax.set_ylim(0, max(model_mean.max(), exp_mean.max()) + 12)
    ax.tick_params(length=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_title(
        "Model vs experiment — dysbiotic/static composition\n"
        f"MAE = {m['mae_pp']:.1f} pp   ·   TVD = {m['tvd']:.3f}   "
        f"(worst: {max(m['per_species_abs_error_pp'], key=m['per_species_abs_error_pp'].get)})",
        fontsize=11, fontweight="bold")

    fig.tight_layout()
    _OUT_FIG.parent.mkdir(exist_ok=True)
    fig.savefig(_OUT_FIG, dpi=200, bbox_inches="tight")
    _OUT_JSON.parent.mkdir(exist_ok=True)
    _OUT_JSON.write_text(json.dumps(m, indent=2))
    print(f"wrote {_OUT_FIG}")
    print(f"wrote {_OUT_JSON}")
    print(f"  MAE = {m['mae_pp']:.2f} pp   TVD = {m['tvd']:.3f}   "
          f"max = {m['max_abs_error_pp']:.2f} pp")
    return m


if __name__ == "__main__":
    main()

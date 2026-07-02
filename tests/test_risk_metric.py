#!/usr/bin/env python3
"""
Unit tests for the T2.3 risk metric  P[sigma > threshold].

Run from the repo root:
    python -m pytest tests/test_risk_metric.py -v
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

_FEM_DIR = Path(__file__).resolve().parent.parent
_JAXFEM = _FEM_DIR / "JAXFEM"
for p in (str(_FEM_DIR), str(_JAXFEM)):
    if p not in sys.path:
        sys.path.insert(0, p)

import risk_metric as rm  # noqa: E402

MPA_PER_KPA = rm.MPA_PER_KPA


def test_exceedance_prob_basic():
    samples = np.array([1.0, 2.0, 3.0, 4.0]) * MPA_PER_KPA  # kPa -> MPa
    # strictly greater than 2 kPa -> {3,4} -> 0.5
    assert rm.exceedance_prob(samples, 2.0 * MPA_PER_KPA) == pytest.approx(0.5)
    # above the max -> 0
    assert rm.exceedance_prob(samples, 10.0 * MPA_PER_KPA) == 0.0
    # below the min -> 1
    assert rm.exceedance_prob(samples, 0.0) == 1.0


def test_exceedance_empty():
    assert np.isnan(rm.exceedance_prob(np.array([]), 1.0))


def test_survival_monotone_decreasing():
    rng = np.random.default_rng(1)
    samples = rng.lognormal(mean=-4.0, sigma=0.5, size=500)  # MPa
    taus = np.linspace(0, samples.max() * 1.1, 50)
    surv = rm.survival_curve(samples, taus)
    # survival function is non-increasing
    assert np.all(np.diff(surv) <= 1e-12)
    # bounded in [0, 1]
    assert surv.min() >= 0.0 and surv.max() <= 1.0
    # at tau=0 essentially all samples exceed
    assert surv[0] == pytest.approx(1.0)


def test_bootstrap_ci_brackets_point_estimate():
    rng = np.random.default_rng(2)
    samples = rng.normal(loc=5.0, scale=1.0, size=200) * MPA_PER_KPA
    tau = 5.0 * MPA_PER_KPA
    p = rm.exceedance_prob(samples, tau)
    lo, hi = rm.bootstrap_ci(samples, tau, n_boot=2000, seed=0)
    assert 0.0 <= lo <= p <= hi <= 1.0
    assert hi - lo > 0.0  # non-degenerate sample -> finite width


def test_bootstrap_ci_degenerate_collapses():
    # all samples identical -> exceedance is 0 or 1 deterministically, CI width 0
    samples = np.full(30, 2.0 * MPA_PER_KPA)
    lo, hi = rm.bootstrap_ci(samples, 1.0 * MPA_PER_KPA, n_boot=1000, seed=0)
    assert lo == 1.0 and hi == 1.0
    lo, hi = rm.bootstrap_ci(samples, 3.0 * MPA_PER_KPA, n_boot=1000, seed=0)
    assert lo == 0.0 and hi == 0.0


def test_compute_risk_synthetic():
    # Two conditions: one clearly above tau, one clearly below.
    ci = {
        "commensal_hobic": {
            "sigma_map": 10.0 * MPA_PER_KPA,
            "sigma_samples": list(np.full(20, 10.0 * MPA_PER_KPA)),
        },
        "dysbiotic_static": {
            "sigma_map": 1.0 * MPA_PER_KPA,
            "sigma_samples": list(np.full(20, 1.0 * MPA_PER_KPA)),
        },
    }
    summ = rm.compute_risk(ci, threshold_kpa=5.0)
    assert summ["conditions"]["commensal_hobic"]["risk"] == 1.0
    assert summ["conditions"]["dysbiotic_static"]["risk"] == 0.0
    # headline threshold is folded into the reference set
    assert 5.0 in summ["reference_thresholds_kpa"]


def test_compute_risk_fallback_map_only():
    # No samples -> falls back to MAP as a single sample; still returns a value.
    ci = {"commensal_hobic": {"sigma_map": 8.0 * MPA_PER_KPA}}
    summ = rm.compute_risk(ci, threshold_kpa=5.0)
    d = summ["conditions"]["commensal_hobic"]
    assert d["n_samples"] == 1
    assert d["degenerate"] is True
    assert d["risk"] == 1.0


@pytest.mark.skipif(
    not (_JAXFEM / "_posterior_ci" / "klempt_stress_ci_tooth.json").exists(),
    reason="committed posterior CI json required",
)
def test_real_ci_json_tooth():
    ci = rm.load_ci("tooth")
    summ = rm.compute_risk(ci, threshold_kpa=5.0)
    conds = summ["conditions"]
    # all reported risks are valid probabilities
    for d in conds.values():
        assert 0.0 <= d["risk"] <= 1.0
        lo, hi = d["risk_ci90"]
        assert 0.0 <= lo <= hi <= 1.0
    # commensal-HOBIC (early, So-dominant) carries the highest growth-stress:
    # its exceedance at a mid threshold should dominate dysbiotic-static.
    assert conds["commensal_hobic"]["risk"] >= conds["dysbiotic_static"]["risk"]

"""Tests for JAXFEM/psi_spread_sensitivity.py.

The script's claim -- that composition reaches the stress only once the
species' viabilities differ, and by how much -- is only worth anything if it
is standing on the same ground as the runs the thesis quotes. So the tests
pin it against two figures that were established independently of it:

  * at zero psi spread the four conditions' alpha differ by ~0.6%, the
    1.006x reported for the 4-condition deck after the calibration constants
    were unified;
  * over the longer horizon that rises to ~1.2% (it was ~2.4% before the
    Hill-gate fix of 2026-10-01: with the gate off, species 5's interaction
    used to be multiplied by 0; coupling/ODE_TMCMC_CROSSCHECK.md).

If either drifts, the sweep is measuring something else and its answer about
what psi would have to be is not about this model.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("jax")
pytest.importorskip("openpyxl")

_ROOT = Path(__file__).resolve().parents[1]
_SRC = _ROOT / "JAXFEM" / "psi_spread_sensitivity.py"

if not (_ROOT / "data" / "heine_species_distribution_biofilm.xlsx").exists():
    pytest.skip("Heine CLSM workbook not present", allow_module_level=True)

_spec = importlib.util.spec_from_file_location("psi_spread_sensitivity", _SRC)
pss = importlib.util.module_from_spec(_spec)
sys.modules["psi_spread_sensitivity"] = pss
_spec.loader.exec_module(pss)


@pytest.fixture(scope="module")
def comps():
    return pss.day1_compositions()


def _spread(comps, psi, horizon):
    dt_h, n_sub = pss.HORIZONS[horizon]
    return pss.spread_of([pss.alpha(comps[c], psi, dt_h, n_sub)
                          for c, _, _ in pss.CONDITIONS])


def test_all_four_conditions_are_real_compositions(comps):
    assert set(comps) == {"CS", "CH", "DS", "DH"}
    for name, phi in comps.items():
        assert phi.shape == (5,), name
        assert np.isclose(phi.sum(), 0.999999, atol=1e-6), name
        assert np.all(phi >= 0.0), name


def test_uniform_psi_reproduces_the_reported_alpha_spread(comps):
    """The baseline this whole sweep is measured against: with psi identical
    across species the four conditions differ by the 1.006x (~0.6%) the
    4-condition deck reports."""
    s = _spread(comps, np.full(5, pss.PSI_MAX), "deck")
    assert 0.005 < s < 0.008, f"expected ~0.6%, got {s:.4%}"


def test_the_spread_is_larger_over_the_longer_horizon(comps):
    """~1.2% near t = 1e-2 (2.4% before the Hill-gate fix, 2026-10-01)."""
    s = _spread(comps, np.full(5, pss.PSI_MAX), "peak")
    assert 0.010 < s < 0.015, f"expected ~1.2%, got {s:.4%}"


def test_giving_the_species_different_viabilities_brings_the_spread_back(comps):
    """The point of the study, as a test: if a psi spread changed nothing there
    would be no measurement to make, and a regression that ignored psi would
    pass every other check here."""
    uniform = _spread(comps, np.full(5, pss.PSI_MAX), "deck")
    tilted = _spread(comps, pss.PSI_MAX * (1.0 - 0.5 * np.linspace(0, 1, 5)),
                     "deck")
    assert tilted > 5 * uniform, f"uniform {uniform:.4%}, tilted {tilted:.4%}"


def test_pattern_spans_the_unit_interval():
    """s is defined as (psi_max - psi_min)/psi_max, which only holds if the
    random shape really spans [0, 1]."""
    rng = np.random.default_rng(3)
    for _ in range(20):
        v = pss.pattern(rng)
        assert v.shape == (5,)
        assert v.min() == pytest.approx(0.0) and v.max() == pytest.approx(1.0)

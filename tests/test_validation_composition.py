"""Regression guard for the model <-> experiment composition validation.

Locks in the agreement between the pipeline's dysbiotic-static composition and
the Heine CLSM measurement (validate_composition.py). If a pipeline change
silently degrades that agreement, these bounds fail.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))

pytest.importorskip("openpyxl")
_DATA = _ROOT / "data" / "heine_species_distribution_biofilm.xlsx"
_MODEL = _ROOT / "_ci_0d_results" / "dysbiotic_static" / "samples_0d.json"
pytestmark = pytest.mark.skipif(
    not (_DATA.exists() and _MODEL.exists()),
    reason="Heine workbook or 0D posterior samples not present")

import validate_composition as vc  # noqa: E402


@pytest.fixture(scope="module")
def result():
    M = vc.model_composition()
    H = vc.heine_composition()
    return vc.metrics(M.mean(0), H.mean(0)), M, H


def test_compositions_are_valid_simplex(result):
    _, M, H = result
    assert np.allclose(M.sum(1), 100.0, atol=1e-6)
    assert np.allclose(H.sum(1), 100.0, atol=1e-6)
    assert (M >= 0).all() and (H >= 0).all()


def test_agreement_within_tolerance(result):
    m, _, _ = result
    # current: MAE 4.2 pp, TVD 0.105, max 10.5 pp — bounds give headroom but
    # would catch a real regression in the calibrated composition.
    assert m["mae_pp"] < 8.0, m
    assert m["tvd"] < 0.20, m
    assert m["max_abs_error_pp"] < 15.0, m


def test_dysbiotic_structure_reproduced(result):
    """Both model and experiment must be Veillonella-dominated (the dysbiotic
    hallmark), and P. gingivalis present at double digits — a structural check
    beyond the aggregate error."""
    m, _, _ = result
    model, exp = m["model_mean_pct"], m["exp_mean_pct"]
    assert model["Veillonella"] == max(model.values())
    assert exp["Veillonella"] == max(exp.values())
    assert model["P. gingivalis"] > 10.0

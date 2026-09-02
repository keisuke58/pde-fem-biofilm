"""Pin the saturation of the production phi->E bridge.

See E_SATURATION_FINDING.md. compute_E_composite clips to [E_MIN_PA, E_MAX_PA],
and at the current calibration that clip is active over much of the healthy
half of composition space -- so distinct conditions can report identical
stiffness, silently, in a study whose headline is a comparison between
conditions.

The clamp is deliberate and is not changed here. These tests exist so that if
the calibration moves, or the bound is raised, it is a visible decision rather
than a quiet change to every reported number.
"""
import numpy as np
import pytest

import material_models as mm

EVEN = np.full(5, 0.2)
COMMENSAL = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
TWO_PRODUCERS = np.array([0.5, 0.5, 0.0, 0.0, 0.0])
MOSTLY_PG = np.array([0.05, 0.05, 0.05, 0.05, 0.8])


def _E(phi):
    return float(mm.compute_E_composite(np.asarray(phi, dtype=float)))


@pytest.mark.parametrize("phi,label", [(COMMENSAL, "all commensal"),
                                       (EVEN, "even mix"),
                                       (TWO_PRODUCERS, "two EPS producers")])
def test_healthy_side_compositions_are_on_the_upper_bound(phi, label):
    """Documents the reach of the clamp. If one of these ever comes back below
    the cap, the calibration or the bound has moved -- update the finding and
    every number derived from this bridge."""
    assert _E(phi) == pytest.approx(mm.E_MAX_PA), (
        f"{label} is no longer saturating; E_SATURATION_FINDING.md is stale")


def test_distinct_compositions_collapse_to_the_same_stiffness():
    """The consequence that matters: the comparison loses a real difference."""
    assert _E(EVEN) == _E(TWO_PRODUCERS), "no longer collapsing"


def test_the_dysbiotic_side_is_not_clamped():
    """The clamp is not global -- the pathogen-dominated end reports the
    model's own value, so a comparison across the dysbiosis axis does carry
    signal. Only the healthy half is flattened."""
    e = _E(MOSTLY_PG)
    assert mm.E_MIN_PA < e < mm.E_MAX_PA, e


def test_the_bounds_themselves_have_not_moved():
    assert (mm.E_MIN_PA, mm.E_MAX_PA) == (10.0, 1000.0), (
        "E_MIN_PA/E_MAX_PA changed — every number derived from this bridge "
        "changes with them; see E_SATURATION_FINDING.md §'If it turns out to bind'")

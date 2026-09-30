"""The repo's 0D Hamilton integrator (the code behind ecology_jax.py and the
ANSYS ecology bridge) reproduces every numerical example of Klempt et al.,
Arch Appl Mech 96:164 (2026) / arXiv:2509.01274 -- see
JAXFEM/klempt2026_cases.py for the case data and the two places where the
printed equations disagree with the paper's own figures.

Reference values are read off the paper's figures by eye, so the
tolerances are figure-reading tolerances, not numerical ones.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("jax")

_JAXFEM = Path(__file__).resolve().parents[1] / "JAXFEM"
sys.path.insert(0, str(_JAXFEM))

import klempt2026_reproduction as rep  # noqa: E402
from klempt2026_cases import CASES  # noqa: E402

PHIBAR_TOL = 0.035


def _obs(name, **kw):
    return rep.observables(rep.simulate(CASES[name], **kw))


@pytest.mark.parametrize("name", sorted(CASES))
def test_case_matches_paper_figure(name):
    case, obs = CASES[name], _obs(name)
    assert obs["finite"]
    assert obs["sum_end"] == pytest.approx(1.0, abs=1e-8)
    ref = case["paper"]
    if "phibar_end" in ref:
        np.testing.assert_allclose(obs["phibar_end"], ref["phibar_end"], atol=PHIBAR_TOL)
    if "phi_end" in ref:
        np.testing.assert_allclose(obs["phi_end"], ref["phi_end"], atol=PHIBAR_TOL)
    if "phi0_zero_step" in ref:
        want = ref["phi0_zero_step"]
        assert obs["phi0_zero_step"] is not None
        assert abs(obs["phi0_zero_step"] - want) <= max(60, 0.12 * want)


def test_eq19_as_printed_does_not_reproduce_fig7():
    """Pins the finding that the figures were computed without Eq. 19's 1/2:
    with it, the 4-species case never takes off within the plotted window."""
    obs = _obs("4sp_case1", as_printed=True)
    assert max(obs["phibar_end"]) < 0.1
    assert obs["phi0_zero_step"] is None


def test_eq17_plus_gamma_does_not_reproduce_fig3():
    """Pins the finding that Eq. 17's printed +gamma was not what produced
    the figures (the variation w.r.t. psi has no gamma either)."""
    obs = _obs("2sp_case3", variant="psi_gamma")
    assert abs(obs["phibar_end"][0] - CASES["2sp_case3"]["paper"]["phibar_end"][0]) > 0.2

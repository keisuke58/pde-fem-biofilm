"""With nothing varying in space, the 1D and 2D Hamilton+nutrient PDEs must
reproduce the 0D integrator (the one checked against Klempt et al. 2026)
node by node, at the calibration c* from ecology_constants.py. This pins the
reaction coupling and the nutrient scaling of both PDE codes to the 0D model:
a PDE that silently ran at a different c* (as the 1D code did before
2026-09-29, at c* = 1) fails the first assertion and passes the second."""
import sys
from pathlib import Path

import pytest

pytest.importorskip("jax")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "JAXFEM"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pde_uniform_consistency as puc  # noqa: E402


@pytest.fixture(scope="module")
def rows():
    return {r[0]: r for r in puc.compare()[1]}


@pytest.mark.parametrize("label", ["1D (default)", "2D (default scale)"])
def test_pde_matches_0d_at_calibration_c_star(rows, label):
    _, spread, c_spread, d_cstar, d_one = rows[label]
    assert spread == 0.0          # a uniform state stays uniform
    assert c_spread == 0.0        # no consumption: nutrient stays at the boundary value
    assert d_cstar < 1e-9         # every node follows the 0D chain at C_STAR
    assert d_one > 1e-3           # ... and demonstrably not the old c* = 1 one

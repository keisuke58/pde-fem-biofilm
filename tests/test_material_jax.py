"""JAX material core + exact (AD) tangent — cross-validated against the
verified NumPy/finite-difference path.

material_jax.py mirrors material_server.py's core so the tangent can be taken
by forward-mode AD instead of a finite difference with a hand-picked step.
This file holds the two together and, just as importantly, pins what the
exact tangent is and is NOT worth — several claims that are tempting to
assume turn out to be false when measured:

  * the FD tangent at the USERMAT's PERT=1e-7 is already accurate to ~1e-7
    relative, so the exact tangent is not a convergence fix;
  * the FD step-size optimum is NOT material-dependent (sigma scales with
    C10, which cancels out of the relative error), so one global PERT is fine
    across the ~31x stiffness range the composition path introduces.

What AD does uniquely buy is d(sigma)/d(theta) for posterior propagation,
which a finite-difference path could only match with extra solves per Gauss
point.

Requires jax; skipped where it is not installed.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("jax")

_COUP = Path(__file__).resolve().parents[1] / "ansys_usermat" / "coupling"
sys.path.insert(0, str(_COUP))

import jax.numpy as jnp                     # noqa: E402
import material_jax as mj                   # noqa: E402
import material_server as ms                # noqa: E402

F_TEST = np.array([[1.05, 0.02, 0.0], [0.0, 0.98, 0.01], [0.0, 0.0, 0.97]])
I3 = np.eye(3)

# (label, alpha, C10, C01, D1, eta, mtype, dt)
CASES = [
    ("elastic neo-Hookean", (0.2, 1.0, 0.0, 0.01, 0.0, 0.0, 1.0)),
    ("viscous", (0.3, 1.0, 0.0, 0.01, 5.0, 0.0, 0.1)),
    ("Mooney-Rivlin", (0.1, 0.8, 0.3, 0.02, 0.0, 1.0, 1.0)),
    ("viscous + Mooney-Rivlin", (0.4, 0.6, 0.2, 0.015, 3.0, 1.0, 0.05)),
]

# Composition-derived constants (composition_to_material.py), the stiffest and
# softest of the four clinical conditions.
P_CH = (0.2, 166.388, 24.958, 0.00241, 72.5, 1.0, 0.01)
P_DS = (0.2, 5.351, 0.803, 0.075, 432.5, 1.0, 0.01)


def _rel(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b)))
                 / max(float(np.max(np.abs(np.asarray(b)))), 1e-30))


@pytest.mark.parametrize("label,params", CASES, ids=[c[0] for c in CASES])
def test_jax_core_matches_numpy_core(label, params):
    """The AD tangent is only trustworthy if the function it differentiates is
    the same one the verified path evaluates."""
    sv_np, fv_np, det_np = ms.stress_core(F_TEST, I3, *params)
    sv_jx, fv_jx, det_jx = mj.stress_core_jax(
        jnp.asarray(F_TEST), jnp.asarray(I3), *params)
    np.testing.assert_allclose(np.asarray(sv_jx), sv_np, rtol=1e-12, atol=1e-13)
    np.testing.assert_allclose(np.asarray(fv_jx), fv_np, rtol=1e-12, atol=1e-13)
    assert float(det_jx) == pytest.approx(det_np, rel=1e-12)


@pytest.mark.parametrize("label,params", CASES, ids=[c[0] for c in CASES])
def test_exact_tangent_agrees_with_fd_within_truncation(label, params):
    """Two independent routes to the same 6x6 must agree to the FD's own
    truncation level. This is what independently validates the finite
    -difference tangent -- and so the Fortran one, which is 0-ULP identical
    to the NumPy core (test_coupling_vs_fortran.py)."""
    D_fd = ms.dsde_perturbation(F_TEST, I3, params)
    D_ex = np.asarray(mj.dsde_exact(F_TEST, I3, params))
    assert _rel(D_fd, D_ex) < 1e-6


def test_fd_tangent_is_already_accurate_enough_for_newton():
    """Pins the negative result: at the USERMAT's PERT the FD tangent is
    within ~1e-7 of exact, far tighter than Newton's convergence rate
    responds to. So the tangent does not explain the cylinder-shell case
    that stops converging at alpha=0.015 -- look elsewhere."""
    for params in (P_CH, P_DS):
        D_fd = ms.dsde_perturbation(F_TEST, I3, params, h=1.0e-7)
        D_ex = np.asarray(mj.dsde_exact(F_TEST, I3, params))
        assert _rel(D_fd, D_ex) < 1e-6


def test_fd_step_optimum_is_not_material_dependent():
    """Also a negative result, recorded because the opposite is easy to
    assume: sigma scales with C10, so C10 cancels out of the relative error
    and the best step is the same for the stiffest and softest condition
    despite a ~31x spread in C10. One global PERT is adequate."""
    steps = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

    def best_step(params):
        D_ex = np.asarray(mj.dsde_exact(F_TEST, I3, params))
        errs = [(_rel(ms.dsde_perturbation(F_TEST, I3, params, h=h), D_ex), h)
                for h in steps]
        return min(errs)[1]

    assert P_CH[1] / P_DS[1] > 20.0            # the spread really is large
    assert best_step(P_CH) == best_step(P_DS)


@pytest.mark.parametrize("label,params", CASES, ids=[c[0] for c in CASES])
def test_tangent_has_no_step_to_tune(label, params):
    """dsde_exact takes no step argument at all -- the property that makes it
    worth having. Guards against someone reintroducing one."""
    import inspect
    assert "h" not in inspect.signature(mj.dsde_exact).parameters
    D = np.asarray(mj.dsde_exact(F_TEST, I3, params))
    assert D.shape == (6, 6)
    assert np.all(np.isfinite(D))


def test_elastic_limit_is_differentiable_not_nan():
    """eta = 0 divides by zero in the viscous arm; the guard must keep the
    derivative finite, not merely the value. A jnp.where over a NaN arm
    silently poisons the gradient, so value-only checks would miss this."""
    params = (0.2, 1.0, 0.0, 0.01, 0.0, 0.0, 1.0)          # eta = 0
    D = np.asarray(mj.dsde_exact(F_TEST, I3, params))
    assert np.all(np.isfinite(D)), "elastic-limit tangent contains NaN/inf"
    S = np.asarray(mj.dsigma_dparams(F_TEST, I3, params))
    assert np.all(np.isfinite(S)), "elastic-limit parameter sensitivity is NaN"


@pytest.mark.parametrize("j,name", list(enumerate(mj.SENS_PARAMS)))
def test_parameter_sensitivity_matches_central_difference(j, name):
    """d(sigma)/d(theta) is the capability FD-over-a-socket cannot practically
    match; validate every column against a central difference on the
    independent NumPy core."""
    params = P_CH
    S = np.asarray(mj.dsigma_dparams(F_TEST, I3, params))
    assert S.shape == (6, len(mj.SENS_PARAMS))

    h = 1e-6 * max(abs(params[j]), 1.0)
    up, dn = list(params), list(params)
    up[j] += h
    dn[j] -= h
    fd = (ms.stress_core(F_TEST, I3, *up)[0]
          - ms.stress_core(F_TEST, I3, *dn)[0]) / (2.0 * h)
    assert _rel(S[:, j], fd) < 1e-6, f"d(sigma)/d({name}) disagrees"


def test_server_tangent_backend_defaults_to_fd():
    """The default must stay the finite difference: that is what keeps
    kUsePy=1 vs kUsePy=0 an exact equivalence check rather than an
    approximate one."""
    assert ms.TANGENT_BACKEND == "fd"


def test_server_tangent_backend_switch():
    original = ms.TANGENT_BACKEND
    try:
        ms.set_tangent_backend("jax")
        D_jax = ms._tangent(F_TEST, I3, P_CH)
        ms.set_tangent_backend("fd")
        D_fd = ms._tangent(F_TEST, I3, P_CH)
        assert _rel(D_jax, np.asarray(mj.dsde_exact(F_TEST, I3, P_CH))) < 1e-12
        assert _rel(D_jax, D_fd) < 1e-6          # same quantity, both routes
    finally:
        ms.set_tangent_backend(original)
    with pytest.raises(ValueError, match="unknown tangent backend"):
        ms.set_tangent_backend("nope")

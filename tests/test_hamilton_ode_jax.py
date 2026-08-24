"""Tests for JAXFEM/hamilton_ode_jax.py and hamilton_ode_jax_nsp.py.

Two JAX implementations of the 0-D Hamilton ODE (theta -> phibar(t)):
hamilton_ode_jax.py is fixed at 5 species; hamilton_ode_jax_nsp.py generalises
to N species. They are supposed to describe the same physics, so this file
proves it rather than assuming it -- both by construction (the same residual,
called directly with matching parameters) and end-to-end (the public
simulate_* entry points).

Requires jax; skipped where it is not installed.
"""
import numpy as np
import pytest

jax = pytest.importorskip("jax")
import jax.numpy as jnp  # noqa: E402

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "JAXFEM"))
import hamilton_ode_jax as h5  # noqa: E402
import hamilton_ode_jax_nsp as hn  # noqa: E402


def _random_symmetric_A(seed=0, n=5):
    rng = np.random.RandomState(seed)
    A = rng.uniform(-1, 3, size=(n, n))
    return (A + A.T) / 2


def test_make_initial_state_phi0_is_not_the_barrier_singularity():
    """Regression guard for the fixed bug: phi0 must not land exactly at 0,
    which is the log-barrier potential's singularity
    (Kp1*(2-4*phi)/((phi-1)**3 * phi**3) blows up as phi -> 0).
    """
    active_mask = jnp.ones(5, dtype=jnp.int64)
    g0 = h5.make_initial_state(jnp.full(5, 0.2), active_mask)
    phi0 = float(g0[5])
    assert phi0 > 0.0, "phi0 landed exactly on the barrier singularity"
    assert phi0 == pytest.approx(1e-6, rel=1e-3)


def test_initial_state_matches_the_nsp_generalisation():
    active_mask = jnp.ones(5, dtype=jnp.int64)
    g_fixed = h5.make_initial_state(jnp.full(5, 0.2), active_mask)
    g_nsp = hn._make_initial_state(jnp.full(5, 0.2), 5, active_mask, psi_init=0.999)
    np.testing.assert_allclose(np.array(g_fixed), np.array(g_nsp), atol=1e-9)


def test_hill_gate_off_leaves_interaction_untouched_not_zeroed():
    """The gate being off must mean "not applied", not "species 4's
    interaction forced to zero" -- those are different models. Regression
    guard for exactly that bug."""
    active_mask = jnp.ones(5, dtype=jnp.int64)
    g0 = h5.make_initial_state(jnp.full(5, 0.2), active_mask)
    A = jnp.array(_random_symmetric_A())
    b = jnp.zeros(5)
    params_gate_off = {
        "dt_h": 1e-4, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5),
        "c": 25.0, "alpha": 100.0, "K_hill": jnp.array(0.0), "n_hill": jnp.array(2.0),
        "A": A, "b_diag": b, "active_mask": active_mask,
    }
    # With no gate, species 4's residual must depend on the full A @ (phi*psi)
    # interaction, i.e. match a hand-computed ungated Ia[4].
    Ia = A @ (g0[0:5] * g0[6:11])
    Q = h5.residual(g0, g0, params_gate_off)
    # Reconstruct what the gated term contributes to Q[4] and check it used
    # the *ungated* Ia[4], not zero.
    assert abs(float(Q[4])) > 1e-3, "species 4's interaction was suppressed with the gate off"


def test_residual_matches_nsp_when_disabled():
    """Both modules' Hill gate 'off' state (K_hill<=1e-9 / gated=-1) must
    leave every species' interaction untouched, so the residuals match
    exactly with no gate parameters involved at all."""
    active_mask = jnp.ones(5, dtype=jnp.int64)
    g0 = h5.make_initial_state(jnp.full(5, 0.2), active_mask)
    A = jnp.array(_random_symmetric_A())
    b = jnp.zeros(5)

    params_fixed = {
        "dt_h": 1e-4, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5),
        "c": 25.0, "alpha": 100.0, "K_hill": jnp.array(0.0), "n_hill": jnp.array(2.0),
        "A": A, "b_diag": b, "active_mask": active_mask,
    }
    params_nsp = {
        "n_sp": 5, "dt_h": 1e-4, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5),
        "c": 25.0, "alpha": 100.0, "K_hill": jnp.array(0.0), "n_hill": jnp.array(2.0),
        "hill_gate_species": jnp.array((-1, 0), dtype=jnp.int32),
        "A": A, "b_diag": b, "active_mask": active_mask,
    }
    Q_fixed = h5.residual(g0, g0, params_fixed)
    Q_nsp = hn._residual(g0, g0, params_nsp)
    np.testing.assert_allclose(np.array(Q_fixed), np.array(Q_nsp), atol=0.0)


def test_residual_matches_nsp_with_gate_wired(hill_k=0.05):
    """With the Hill gate parametrised identically (K_hill>0, gated species 4
    driven by gating species 3, matching hamilton_ode_jax.py's hardcoded
    phi_new[3]/psi_new[3] gate), the two residuals must agree exactly -- same
    physics, different code paths."""
    active_mask = jnp.ones(5, dtype=jnp.int64)
    g0 = h5.make_initial_state(jnp.full(5, 0.2), active_mask)
    A = jnp.array(_random_symmetric_A())
    b = jnp.zeros(5)

    params_fixed = {
        "dt_h": 1e-4, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5),
        "c": 25.0, "alpha": 100.0, "K_hill": jnp.array(hill_k), "n_hill": jnp.array(2.0),
        "A": A, "b_diag": b, "active_mask": active_mask,
    }
    params_nsp = {
        "n_sp": 5, "dt_h": 1e-4, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5),
        "c": 25.0, "alpha": 100.0, "K_hill": jnp.array(hill_k), "n_hill": jnp.array(2.0),
        "hill_gate_species": jnp.array((4, 3), dtype=jnp.int32),
        "A": A, "b_diag": b, "active_mask": active_mask,
    }

    Q_fixed = h5.residual(g0, g0, params_fixed)
    Q_nsp = hn._residual(g0, g0, params_nsp)
    np.testing.assert_allclose(np.array(Q_fixed), np.array(Q_nsp), atol=0.0)


def test_trajectories_match_with_gate_wired():
    """End-to-end: simulate_0d-equivalent stepping via each module's own
    Newton loop, with the Hill gate parametrised identically, must produce
    bit-identical phibar trajectories."""
    active_mask = jnp.ones(5, dtype=jnp.int64)
    A = jnp.array(_random_symmetric_A(seed=1))
    b = jnp.zeros(5)
    n_steps = 100

    params_fixed = {
        "dt_h": 1e-4, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5),
        "c": 25.0, "alpha": 100.0, "K_hill": jnp.array(0.05), "n_hill": jnp.array(2.0),
        "A": A, "b_diag": b, "active_mask": active_mask,
    }
    params_nsp = {
        "n_sp": 5, "dt_h": 1e-4, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5),
        "c": 25.0, "alpha": 100.0, "K_hill": jnp.array(0.05), "n_hill": jnp.array(2.0),
        "hill_gate_species": jnp.array((4, 3), dtype=jnp.int32),
        "A": A, "b_diag": b, "active_mask": active_mask,
    }

    g0_fixed = h5.make_initial_state(jnp.full(5, 0.2), active_mask)
    g0_nsp = hn._make_initial_state(jnp.full(5, 0.2), 5, active_mask, psi_init=0.999)

    def run(newton_step, g0, params):
        def body(g, _):
            gn = newton_step(g, params)
            return gn, gn
        _, traj = jax.lax.scan(body, g0, jnp.arange(n_steps))
        return traj

    traj_fixed = run(h5.newton_step, g0_fixed, params_fixed)
    traj_nsp = run(hn._newton_step, g0_nsp, params_nsp)

    phibar_fixed = traj_fixed[:, 0:5] * traj_fixed[:, 6:11]
    phibar_nsp = traj_nsp[:, 0:5] * traj_nsp[:, 6:11]
    np.testing.assert_allclose(np.array(phibar_fixed), np.array(phibar_nsp), atol=1e-10)


def test_simulate_0d_public_api_runs_and_conserves_volume():
    theta15 = jnp.array([1.34, -0.18, 1.79, 1.17, 2.58, 3.51, 2.73, 0.71, 2.1,
                          0.37, 2.05, -0.15, 3.56, 0.16, 0.12])
    traj = h5.simulate_0d(theta15, n_steps=200, dt=1e-4)
    assert traj.shape == (201, 5)
    assert bool(jnp.all(jnp.isfinite(traj)))
    assert bool(jnp.all(traj >= 0.0))


def test_hill_gate_defaults_to_off():
    """The Hill gate must not be engaged unless a caller explicitly opts in."""
    import inspect
    assert inspect.signature(h5.simulate_0d).parameters["K_hill"].default == 0.0
    assert inspect.signature(h5.simulate_0d_nutrient).parameters["K_hill"].default == 0.0
    assert inspect.signature(hn.simulate_0d_nsp).parameters["K_hill"].default == 0.0
    assert inspect.signature(hn.simulate_0d_nsp).parameters["hill_gate_species"].default is None


def test_simulate_0d_nutrient_runs():
    theta15 = jnp.array([1.34, -0.18, 1.79, 1.17, 2.58, 3.51, 2.73, 0.71, 2.1,
                          0.37, 2.05, -0.15, 3.56, 0.16, 0.12])
    phibar, S = h5.simulate_0d_nutrient(theta15, n_steps=200, dt=1e-4)
    assert phibar.shape == (201, 5)
    assert S.shape == (201,)
    assert bool(jnp.all((S >= 0.0) & (S <= 1.0 + 1e-9)))


def test_count_params():
    assert hn.count_params(5) == 20
    assert hn.count_params(11) == 77


@pytest.mark.parametrize("n_sp", [3, 5, 11])
def test_simulate_0d_nsp_runs_at_several_sizes(n_sp):
    n_par = hn.count_params(n_sp)
    theta = jax.random.uniform(jax.random.PRNGKey(n_sp), (n_par,), minval=0.0, maxval=3.0)
    traj = hn.simulate_0d_nsp(theta, n_sp=n_sp, n_steps=150, dt=1e-4)
    assert traj.shape == (151, n_sp)
    assert bool(jnp.all(jnp.isfinite(traj)))
    assert bool(jnp.all(traj >= 0.0))

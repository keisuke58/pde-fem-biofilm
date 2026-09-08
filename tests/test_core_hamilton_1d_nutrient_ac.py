"""Tests for JAXFEM/core_hamilton_1d_nutrient_ac.py.

core_hamilton_1d_nutrient_ac.py adds per-species Allen-Cahn interface
sharpening and chemotaxis on top of core_hamilton_1d_nutrient.py's existing
Hamilton-reaction + nutrient-diffusion model -- see that module's docstring
for exactly which term comes from which source (an original synthesis, not
a literal reproduction of either Klempt paper). The critical property this
file checks: with the new terms switched off (Gamma_ac=0, R_chemo=0), the
extended module must reproduce the existing, already-used
core_hamilton_1d_nutrient.simulate_hamilton_1d_nutrient bit-for-bit -- so
this file cannot have silently changed behaviour for every caller already
depending on it (multiscale_coupling_1d.py, extract_alpha_field_1d.py, ...).

Requires jax; skipped where it is not installed.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

jax = pytest.importorskip("jax")
import jax.numpy as jnp  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from JAXFEM.core_hamilton_1d import THETA_DEMO  # noqa: E402
from JAXFEM.core_hamilton_1d_nutrient import (  # noqa: E402
    simulate_hamilton_1d_nutrient,
)
from JAXFEM.core_hamilton_1d_nutrient_ac import (  # noqa: E402
    simulate_hamilton_1d_nutrient_ac,
)

_SMALL = dict(n_macro=4, n_react_sub=3, n_sub_c=3, N=8, L=1.0, dt_h=1e-5)


def test_zeroed_ac_terms_reproduce_the_nutrient_only_model_exactly():
    D_eff = jnp.full(5, 1e-3)
    G_ref, c_ref = simulate_hamilton_1d_nutrient(THETA_DEMO, D_eff, **_SMALL)
    G_ac, c_ac, alpha_ac = simulate_hamilton_1d_nutrient_ac(
        THETA_DEMO,
        D_eff,
        Gamma_ac=jnp.zeros(5),
        R_chemo=jnp.zeros(5),
        **_SMALL,
    )
    np.testing.assert_allclose(np.array(G_ac), np.array(G_ref), atol=0.0)
    np.testing.assert_allclose(np.array(c_ac), np.array(c_ref), atol=0.0)


def test_alpha_field_is_nonnegative_and_nondecreasing_with_zero_ac():
    """Eq.36 diagnostic field: alpha_dot = k_alpha*phi_total >= 0 always
    (phi_total >= 0 by construction), so alpha(x,t) must be monotone
    nondecreasing along the trajectory at every node."""
    D_eff = jnp.full(5, 1e-3)
    _, _, alpha_ac = simulate_hamilton_1d_nutrient_ac(
        THETA_DEMO,
        D_eff,
        Gamma_ac=jnp.zeros(5),
        R_chemo=jnp.zeros(5),
        **_SMALL,
    )
    alpha_np = np.array(alpha_ac)
    assert np.all(alpha_np >= 0.0)
    assert np.all(np.diff(alpha_np, axis=0) >= -1e-12)


def test_runs_finite_and_conserves_the_simplex_with_ac_terms_on():
    D_eff = jnp.full(5, 1e-3)
    Gamma_ac = jnp.full(5, 2.0)
    R_chemo = jnp.full(5, 1.0)
    G_ac, c_ac, alpha_ac = simulate_hamilton_1d_nutrient_ac(
        THETA_DEMO, D_eff, Gamma_ac, R_chemo, **_SMALL
    )
    assert bool(jnp.all(jnp.isfinite(G_ac)))
    assert bool(jnp.all(jnp.isfinite(c_ac)))
    assert bool(jnp.all(jnp.isfinite(alpha_ac)))

    phi = G_ac[:, :, 0:5]
    phi0 = G_ac[:, :, 5]
    assert bool(jnp.all(phi >= -1e-9))
    assert bool(jnp.all(phi0 >= -1e-9))
    simplex_sum = jnp.sum(phi, axis=2) + phi0
    np.testing.assert_allclose(np.array(simplex_sum), 1.0, atol=1e-6)


def test_ac_terms_actually_change_the_trajectory():
    """Sanity guard against a wiring bug that leaves Gamma_ac/R_chemo
    inert -- turning them on must change phi relative to the off case."""
    D_eff = jnp.full(5, 1e-3)
    G_off, _, _ = simulate_hamilton_1d_nutrient_ac(
        THETA_DEMO, D_eff, jnp.zeros(5), jnp.zeros(5), **_SMALL
    )
    G_on, _, _ = simulate_hamilton_1d_nutrient_ac(
        THETA_DEMO, D_eff, jnp.full(5, 5.0), jnp.full(5, 3.0), **_SMALL
    )
    assert not bool(jnp.allclose(G_on[:, :, 0:5], G_off[:, :, 0:5]))

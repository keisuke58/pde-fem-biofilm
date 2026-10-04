"""Reduced one- and two-species runs through the unchanged five-species bridge.

The ANSYS side speaks a fixed five-species interface (g(12), theta(20),
ustatev(72:83), prop(8:27)). A one- or two-species run is obtained by masking
species off in the ecology model, set per server with --active-species. These
tests pin that this IS the n-species model -- the masked five-species step
reproduces JAXFEM/hamilton_ode_jax_nsp at n = 1 and n = 2 -- and that the
naive alternative, starting the extra species at zero, is not.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

jax = pytest.importorskip("jax")
import jax.numpy as jnp                                   # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(_ROOT / "ansys_usermat" / "coupling"), str(_ROOT),
                str(_ROOT / "JAXFEM")]

import ecology_jax as eco                                 # noqa: E402
import hamilton_ode_jax_nsp as nsp                        # noqa: E402
import material_server as ms                              # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO      # noqa: E402
from ecology_constants import C_STAR, ALPHA_STAR, K_HILL, N_HILL  # noqa: E402

DT, N = 1.0e-4, 400
T = np.asarray(THETA_DEMO, dtype=float)


def _theta5(n):
    """n-species parameters placed where the five-species encoding keeps
    them: A's upper triangle first (column-major), b from index 15."""
    th = np.zeros(20)
    na = n * (n + 1) // 2
    th[:na] = T[:na]
    th[15:15 + n] = T[15:15 + n]
    return th


def _nsp_run(n, g0):
    th = jnp.array(list(T[:n * (n + 1) // 2]) + list(T[15:15 + n]))
    A, b = nsp.theta_to_matrices(th, n)
    p = {"n_sp": n, "Kp1": 1e-4, "A": A, "b_diag": b, "dt_h": DT,
         "Eta": jnp.ones(n), "EtaPhi": jnp.ones(n), "c": C_STAR,
         "alpha": ALPHA_STAR, "K_hill": jnp.array(K_HILL, float),
         "n_hill": jnp.array(N_HILL, float),
         "hill_gate_species": jnp.array((-1, 0), dtype=jnp.int32),
         "active_mask": jnp.ones(n, dtype=jnp.int64)}
    step = jax.jit(lambda g: nsp._newton_step(g, p))
    g = jnp.asarray(g0)
    for _ in range(N):
        g = step(g)
    return np.asarray(g)


def _map(n):
    """indices of (phi_1..n, phi0, psi_1..n, gamma) in the five-species g"""
    return list(range(n)) + [5] + list(range(6, 6 + n)) + [11]


@pytest.mark.parametrize("n", [1, 2])
def test_masked_five_species_is_the_n_species_model(n):
    seed = np.asarray(eco.default_initial_state())       # what Fortran seeds
    g5, _ = eco.ecology_substeps(seed, _theta5(n), DT * N, N, n_active=n)
    g5 = np.asarray(g5)
    s = eco.seed_inactive(seed, n)
    gn = _nsp_run(n, s[_map(n)])
    np.testing.assert_allclose(g5[_map(n)], gn, rtol=1e-10, atol=1e-12)
    off = list(range(n, 5)) + list(range(6 + n, 11))
    assert np.all(g5[off] == 0.0), "switched-off species must stay exactly 0"


def test_starting_extra_species_at_zero_is_not_the_two_species_model():
    """The clip revives a zero species: this is why the mask exists."""
    seed = eco.seed_inactive(np.asarray(eco.default_initial_state()), 2)
    g_naive, _ = eco.ecology_substeps(seed, _theta5(2), DT * N, N, n_active=5)
    g_mask, _ = eco.ecology_substeps(seed, _theta5(2), DT * N, N, n_active=2)
    assert np.max(np.abs(np.asarray(g_naive)[2:5])) > 1e-6
    assert np.all(np.asarray(g_mask)[2:5] == 0.0)


def test_five_species_path_is_bit_identical_to_before():
    seed = np.asarray(eco.default_initial_state())
    a, pa = eco.ecology_substeps(seed, THETA_DEMO, 1e-3, 10)
    b, pb = eco.ecology_substeps(seed, THETA_DEMO, 1e-3, 10, n_active=5)
    assert np.array_equal(np.asarray(a), np.asarray(b)) and pa == pb


def test_seeding_zeroes_switched_off_species_once():
    seed = np.asarray(eco.default_initial_state())
    s = eco.seed_inactive(seed, 2)
    assert np.all(s[[2, 3, 4, 8, 9, 10]] == 0.0)
    assert s[5] == pytest.approx(1.0 - s[0:2].sum(), abs=1e-15)
    assert np.array_equal(eco.seed_inactive(s, 2), s)     # idempotent


def test_server_uses_the_configured_species_count():
    seed = list(map(float, np.asarray(eco.default_initial_state())))
    req = {"g": seed, "theta": list(_theta5(1)), "dt_h": 1e-3, "n_sub": 10}
    try:
        ms.set_active_species(1)
        g = json.loads(ms.evaluate_ecology(req))["g_new"]
        assert all(x == 0.0 for x in g[1:5] + g[7:11])
        g2 = json.loads(ms.evaluate_ecology({**req, "n_active": 5}))["g_new"]
        assert any(x != 0.0 for x in g2[1:5])            # per-request override
    finally:
        ms.set_active_species(5)


def test_out_of_range_species_count_is_refused():
    with pytest.raises(ValueError):
        ms.set_active_species(0)
    with pytest.raises(ValueError):
        eco.active_mask(6)

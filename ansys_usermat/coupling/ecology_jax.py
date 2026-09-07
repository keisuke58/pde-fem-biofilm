#!/usr/bin/env python3
"""ecology_jax.py -- Gauss-point bridge to the verified 0D Hamilton ecology ODE.

Wraps jax_hamilton_0d_5species_demo.py's Newton-per-step integrator (the
extended-Hamilton-principle 0D ODE system in (phi, phi0, psi, gamma),
RESEARCH_MODEL.md sec.1) as a single-step "ecology_step" callable, matching
the shape of the material bridge next to it (material_server.py): one Gauss
point, one increment, one call. This is coupling/README.md's next-steps
item 4 -- unlike the growth PDE (JAXFEM/klempt_pde_jax.py), a per-Gauss-point
Fortran call cannot see neighbouring points, so only the *local* reaction
dynamics (no diffusion) are representable this way; that is what "0D" means
here, not an approximation of the 1D/2D reaction-diffusion PDE.

State vector g (12,): phi[0:5], phi0[5], psi[6:11], gamma[11] -- same layout
jax_hamilton_0d_5species_demo.integrate_0d uses.
theta (20,): the 15 independent A_ij plus 5 b_i entries, same encoding as
theta_to_matrices there (the TMCMC-calibrated interaction parameters,
RESEARCH_MODEL.md sec.1).

No physics is duplicated here -- newton_step_jit/theta_to_matrices are
imported from the demo module directly, so this cannot silently drift from
the verified 0D reference.
"""
from __future__ import annotations

import sys
from pathlib import Path

import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from jax_hamilton_0d_5species_demo import (  # noqa: E402
    newton_step_jit, theta_to_matrices,
)

G_DIM = 12
THETA_DIM = 20


def default_hparams(dt_h: float) -> dict:
    """Fixed ecology-model hyperparameters, exactly as integrate_0d builds
    them (Kp1, Eta, EtaPhi, c, the Hill-saturation 'alpha', K_hill, n_hill,
    active_mask) -- everything the per-call (g, theta) pair does not carry.

    Note: this dict's "alpha" key is the ODE's own Hill-function coefficient
    (theta_to_matrices/residual's naming), unrelated to the mechanical
    growth driver alpha used elsewhere in this repo (Fg=(1+alpha)I). Kept
    under the demo module's own name for a byte-identical import rather than
    renamed and risking drift from the verified reference.
    """
    return {
        "dt_h": dt_h,
        "Kp1": 1e-4,
        "Eta": jnp.ones(5),
        "EtaPhi": jnp.ones(5),
        "c": 100.0,
        "alpha": 100.0,
        "K_hill": 0.05,
        "n_hill": 4.0,
        "active_mask": jnp.ones(5, dtype=jnp.int64),
    }


def default_initial_state() -> jnp.ndarray:
    """The same initial composition integrate_0d starts from -- used as a
    fallback when a Gauss point has no measured composition of its own yet.
    Real per-point initial conditions should come from CLSM data, the same
    way composition_to_material.py seeds C10/C01/D1/eta from phi; this is a
    placeholder, not a substitute for that."""
    phi_init = jnp.array([0.12, 0.12, 0.08, 0.05, 0.0])
    phi_init = phi_init * (0.999999 / jnp.sum(phi_init))
    phi0_init = 1.0 - jnp.sum(phi_init)
    psi_init = jnp.ones(5) * 0.999
    g = jnp.zeros(12, dtype=jnp.float64)
    g = g.at[0:5].set(phi_init)
    g = g.at[5].set(phi0_init)
    g = g.at[6:11].set(psi_init)
    g = g.at[11].set(0.0)
    return g


def ecology_step(g_prev, theta, dt_h: float) -> jnp.ndarray:
    """Advance the 0D Hamilton ODE state by one increment dt_h.

    g_prev, theta: array-like, shapes (12,) and (20,). Returns g_new (12,).
    Uses the same fixed-6-Newton-iteration scheme as integrate_0d -- one
    call here is exactly one iteration of that function's inner scan body,
    so a Gauss point stepped once per solver increment reproduces
    integrate_0d exactly when dt_h is held constant across increments.
    """
    A, b_diag = theta_to_matrices(jnp.asarray(theta, dtype=jnp.float64))
    params = default_hparams(dt_h)
    params["A"] = A
    params["b_diag"] = b_diag
    return newton_step_jit(jnp.asarray(g_prev, dtype=jnp.float64), params)


def living_fraction_total(g) -> float:
    """phi_tot = sum_i phi_i * psi_i -- the living-volume-fraction driver for
    the growth reaction term k_alpha * phi_tot (RESEARCH_MODEL.md sec.2), the
    only piece of that PDE representable at a single Gauss point without the
    diffusion term."""
    g = jnp.asarray(g, dtype=jnp.float64)
    phi = g[0:5]
    psi = g[6:11]
    return float(jnp.sum(phi * psi))

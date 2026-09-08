"""
JAXFEM/core_hamilton_1d_nutrient_ac.py
=======================================
Hamilton 1D multi-species + nutrient diffusion, extended with per-species
Allen-Cahn interface sharpening and chemotaxis ("ac" = Allen-Cahn +
chemotaxis).

**This is an original synthesis, not a literal reproduction of either Klempt
paper.** Built 2026-09-08 after `CITATION_AUDIT.md` F1d established that this
repo has two, structurally different, independently-verified Klempt models:

  - Klempt, Geisler, Soleimani, Junker (2026), *A continuum multi-species
    biofilm model with a novel interaction scheme* (arXiv:2509.01274 / AAM
    96, 164) -- Eq. 10, 16-18: a 0-D, N-species, Hamilton-derived reaction
    model, `Ia = A @ (phi*psi)` bilinear/additive, symmetric A, no space, no
    nutrient transport. This is `core_hamilton_1d.py`'s reaction step,
    unmodified here.
  - Klempt, Soleimani, Wriggers, Junker (2024), *A Hamilton principle-based
    model for diffusion-driven biofilm growth* (Biomech Model Mechanobiol
    23, 2091-2113) -- Eq. 34-36: a single scalar field phi (biofilm vs
    void, no species), with spatial diffusion, an Allen-Cahn double well,
    logistic-Monod growth, chemotaxis toward a diffusing/consumed nutrient
    c(x,t), and growth kinematics alpha via alpha_dot = k_alpha*phi.
    `felix_complete_reproduction.py` reproduces this one exactly (PASS).

Neither paper has a model with BOTH multi-species interaction AND spatial
transport/chemotaxis at once. What follows below combines them by taking,
from Eq. 34, only the two spatial mechanisms the Hamilton reaction model has
no analogue for -- interface sharpening and chemotaxis -- and applying one
independent copy of each per species:

    d(phi_i)/dt =   D_eff_i * lap(phi_i)                              [transport, generalizes Eq.34's beta*lap term; already core_hamilton_1d.py's diffusion_step]
                   - Gamma_ac_i * phi_i*(1-phi_i)*(1-2*phi_i)         [per-species interface sharpening -- OUR generalization of Eq.34's double well]
                   - R_chemo_i * monod(c) * d(phi_i)/dx * sign(dc/dx) [per-species chemotaxis -- OUR generalization of Eq.34's chemotaxis term,
                                                                        same discretization as felix_complete_reproduction.py's verified step_phi]

**Deliberately NOT included**: Eq.34's own logistic-Monod growth term
(`K_LOG*phi*(1-phi)*monod`). The Hamilton reaction step already supplies
growth/decay via a different, already-verified mechanism (Eq.16-18); adding
Eq.34's growth term on top would double-count growth with no derivation
justifying the combination. Also not included: `felix_complete_reproduction
.py`'s own `K_ALPHA*alpha` feed-forward term (an extension beyond the
literal Eq.34-36 set) -- the growth-kinematics field `alpha_field` here is a
pure diagnostic output (Eq.36, `alpha_dot = k_alpha * phi_total`), matching
`klempt_pde_jax.py`'s existing role, with no feedback into the reaction
dynamics.

**Naming note, easy to get wrong**: `params["alpha"]` inherited from
`core_hamilton_1d_nutrient.py` is the paper's antibiotic concentration
alpha* (Eq. 17's `alpha*.psi_i.b_i` term) -- unrelated to the mechanical
growth-kinematics alpha (`Fg=(1+alpha)I`) tracked here as `alpha_field`.
Keep them distinct; this module never conflates the two.

**Regression guard**: with `Gamma_ac = 0` and `R_chemo = 0` for every
species, this module's `simulate_hamilton_1d_nutrient_ac` must reproduce
`core_hamilton_1d_nutrient.py`'s `simulate_hamilton_1d_nutrient` bit-for-bit
(same reaction step, same nutrient step, same diffusion, extra terms zeroed)
-- see `tests/test_core_hamilton_1d_nutrient_ac.py`.
"""
from __future__ import annotations

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

from .core_hamilton_1d import (
    clip_state,
    make_initial_state,
    theta_to_matrices,
)
from .core_hamilton_1d_nutrient import (
    nutrient_step,
    reaction_step_c,
)

# ---------------------------------------------------------------------------
# Spatial step: D_eff diffusion + per-species Allen-Cahn + per-species
# chemotaxis, all explicit, folded into one step (they act on the same
# phi field over the same macro-step dt).
# ---------------------------------------------------------------------------


def interface_chemotaxis_step(G, c_field, params):
    """Advance phi_i(x) by one macro-step under diffusion + Allen-Cahn +
    chemotaxis. Mirrors core_hamilton_1d.diffusion_step's clip/renormalize
    convention so the (phi, phi0) simplex constraint stays consistent with
    what the next reaction step's Newton solve expects.
    """
    D_eff = params["D_eff"]  # (5,)
    Gamma_ac = params["Gamma_ac"]  # (5,)
    R_chemo = params["R_chemo"]  # (5,)
    k_monod = params["k_monod"]
    dx = params["dx"]
    dt = params["dt_h"] * params["n_react_sub"]
    eps = 1e-12

    phi = G[:, 0:5]  # (N, 5)
    N = phi.shape[0]

    # Laplacian, interior nodes only (Neumann/no-flux at both ends, matching
    # core_hamilton_1d.diffusion_step's existing convention).
    lap = jnp.zeros_like(phi)
    lap_interior = (phi[0 : N - 2, :] + phi[2:N, :] - 2.0 * phi[1 : N - 1, :]) / (dx * dx)
    lap = lap.at[1 : N - 1, :].set(lap_interior)

    # Per-species gradient, interior nodes only (zero at the boundaries ->
    # no chemotactic flux through either end, consistent with the no-flux
    # diffusion convention above).
    grad_phi = jnp.zeros_like(phi)
    grad_phi_interior = (phi[2:N, :] - phi[0 : N - 2, :]) / (2.0 * dx)
    grad_phi = grad_phi.at[1 : N - 1, :].set(grad_phi_interior)

    grad_c = jnp.zeros_like(c_field)
    grad_c_interior = (c_field[2:N] - c_field[0 : N - 2]) / (2.0 * dx)
    grad_c = grad_c.at[1 : N - 1].set(grad_c_interior)
    sign_grad_c = jnp.sign(grad_c)  # (N,)

    monod_c = c_field / (k_monod + c_field + eps)  # (N,)

    double_well = -Gamma_ac[jnp.newaxis, :] * phi * (1.0 - phi) * (1.0 - 2.0 * phi)
    chemo = (
        R_chemo[jnp.newaxis, :]
        * monod_c[:, jnp.newaxis]
        * grad_phi
        * sign_grad_c[:, jnp.newaxis]
    )

    phi_dot = D_eff[jnp.newaxis, :] * lap + double_well - chemo
    phi_new = phi + dt * phi_dot
    phi_new = jnp.clip(phi_new, 0.0, 1.0)
    phi_sum = jnp.sum(phi_new, axis=1)
    phi_sum = jnp.minimum(phi_sum, 1.0)
    phi0_new = 1.0 - phi_sum

    G_new = G.at[:, 0:5].set(phi_new)
    G_new = G_new.at[:, 5].set(phi0_new)
    return G_new


def alpha_field_step(alpha_field, phi_total, params):
    """Eq.36 (Klempt 2024), diagnostic-only: alpha_dot = k_alpha * phi_total.
    No feedback into the Hamilton reaction dynamics -- see module docstring.
    """
    k_alpha = params["k_alpha"]
    dt = params["dt_h"] * params["n_react_sub"]
    return alpha_field + dt * k_alpha * phi_total


# ---------------------------------------------------------------------------
# Main simulation entry point
# ---------------------------------------------------------------------------


def simulate_hamilton_1d_nutrient_ac(
    theta,
    D_eff,
    Gamma_ac,
    R_chemo,
    D_c=1.0,
    g_eff=50.0,
    k_monod=1.0,
    k_alpha=0.05,
    n_macro=60,
    n_react_sub=20,
    n_sub_c=10,
    N=30,
    L=1.0,
    dt_h=1e-5,
):
    """Hamilton 1D (5 species) + nutrient c(x,t) + per-species Allen-Cahn
    interface sharpening + per-species chemotaxis + diagnostic growth field
    alpha(x,t). See module docstring for exactly which term comes from which
    source.

    Parameters
    ----------
    theta     : (20,)  TMCMC-calibrated Hamilton interaction parameters.
    D_eff     : (5,)   per-species diffusion coefficient.
    Gamma_ac  : (5,)   per-species Allen-Cahn double-well strength (0 = off).
    R_chemo   : (5,)   per-species chemotaxis strength (0 = off).
    D_c, g_eff, k_monod : nutrient diffusion/consumption, as in
                          core_hamilton_1d_nutrient.simulate_hamilton_1d_nutrient.
    k_alpha   : growth-kinematics accumulation rate (Eq.36).

    Returns
    -------
    G_all     : (n_macro+1, N, 12)  Hamilton state trajectory.
    c_all     : (n_macro+1, N)      nutrient field trajectory.
    alpha_all : (n_macro+1, N)      growth-kinematics field trajectory.
    """
    dx = L / (N - 1)
    A, b_diag = theta_to_matrices(theta)
    active_mask = jnp.ones(5, dtype=jnp.int64)

    params = {
        "dt_h": dt_h,
        "Kp1": 1e-4,
        "Eta": jnp.ones(5),
        "EtaPhi": jnp.ones(5),
        "alpha": 100.0,  # antibiotic alpha* (Eq.17) -- NOT alpha_field below
        "K_hill": 0.05,
        "n_hill": 4.0,
        "A": A,
        "b_diag": b_diag,
        "active_mask": active_mask,
        "n_react_sub": n_react_sub,
        "D_eff": jnp.asarray(D_eff, dtype=jnp.float64),
        "dx": dx,
        "D_c": D_c,
        "g_eff": g_eff,
        "k_monod": k_monod,
        "n_sub_c": n_sub_c,
        "Gamma_ac": jnp.asarray(Gamma_ac, dtype=jnp.float64),
        "R_chemo": jnp.asarray(R_chemo, dtype=jnp.float64),
        "k_alpha": k_alpha,
    }

    G0 = make_initial_state(N, active_mask)
    c0 = jnp.ones(N, dtype=jnp.float64)
    alpha0 = jnp.zeros(N, dtype=jnp.float64)

    def body(carry, _):
        G, c, alpha_field = carry
        phi_total = G[:, 0:5].sum(axis=1)
        G = reaction_step_c(G, c, params)
        G = interface_chemotaxis_step(G, c, params)
        c = nutrient_step(c, phi_total, params)
        alpha_field = alpha_field_step(alpha_field, phi_total, params)
        return (G, c, alpha_field), (G, c, alpha_field)

    _, (G_traj, c_traj, alpha_traj) = jax.lax.scan(
        body, (G0, c0, alpha0), jnp.arange(n_macro)
    )

    G_all = jnp.concatenate([G0[jnp.newaxis], G_traj], axis=0)
    c_all = jnp.concatenate([c0[jnp.newaxis], c_traj], axis=0)
    alpha_all = jnp.concatenate([alpha0[jnp.newaxis], alpha_traj], axis=0)
    return G_all, c_all, alpha_all

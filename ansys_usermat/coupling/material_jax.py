#!/usr/bin/env python3
"""material_jax.py — JAX mirror of the material core, with an EXACT tangent.

`material_server.py` builds the 6x6 material Jacobian the same way
`usermat_biofilm.f` does: a forward finite difference of the stress under
spatial perturbations of F, with a hand-picked step (`PERT = 1.0d-7`). That is
what the ANSYS run was verified against, so it stays the default.

This module computes the *same* quantity by forward-mode AD instead
(`jax.jacfwd`), which removes the step entirely.

What that is and is not worth, measured rather than assumed
(tests/test_material_jax.py pins all of these):

  * It does NOT rescue Newton convergence. At PERT=1e-7 the FD tangent
    already agrees with the exact one to ~3e-8 relative, which is far tighter
    than anything Newton's convergence *rate* responds to. In particular this
    rules out the tangent as the explanation for the cylinder-shell case that
    stops converging at alpha=0.015 (ansys_usermat/README.md) -- that has to
    be looked for elsewhere.
  * The step-size sensitivity is also NOT material-dependent, which is worth
    recording because the opposite is easy to assume: sigma scales with C10,
    so C10 cancels out of the relative error, and the optimum sits at h=1e-8
    for both the stiffest (CH, C10=166 Pa) and softest (DS, C10=5.4 Pa)
    condition. One global PERT is fine.
  * What it does buy: no tuned magic number and no truncation error at all;
    an independent check on the FD tangent (and so on the Fortran one, which
    is 0-ULP identical to the NumPy core); and -- the real prize --
    differentiability with respect to the material parameters, so d(sigma)/d(theta)
    for posterior/UQ propagation comes out of the same code instead of
    needing a finite-difference sweep over a socket.

`stress_core_jax` is a line-by-line mirror of `material_server._sigma_and_fv`,
so the tests can hold the two to machine precision and the AD tangent inherits
that provenance.

Voigt order: Abaqus 11,22,33,12,13,23 (same as material_server.py).
"""
from __future__ import annotations

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

I3 = jnp.eye(3)
VOIGT = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))   # Abaqus order


def _sigma_and_fv_jax(F, Fv_old, alpha, C10, C01, D1, eta, mtype, dt):
    """Cauchy stress (3x3), updated Fv (3x3), detFe.

    Mirrors material_server._sigma_and_fv. The two `if` branches there
    (Mooney-Rivlin, viscous vs elastic) become jnp.where so the function stays
    traceable and differentiable; both arms are cheap.
    """
    Fg_inv = I3 / jnp.maximum(1.0 + alpha, 1e-15)
    Ftrial = F @ Fg_inv

    def elastic(Fv):
        Fe = Ftrial @ jnp.linalg.inv(Fv)
        Be = Fe @ Fe.T
        detFe = jnp.maximum(jnp.linalg.det(Fe), 1e-15)
        tmp1 = detFe ** (-2.0 / 3.0)
        i1b = tmp1 * jnp.trace(Be)
        return Be, detFe, tmp1, i1b

    def mr_extra(Be, tmp1, i1b):                 # C01 deviatoric contribution
        T3 = i1b * tmp1 * Be - tmp1 ** 2 * (Be @ Be)
        return T3 - (jnp.trace(T3) / 3.0) * I3

    # trial state -> viscous flow driver (deviatoric Kirchhoff)
    Be, detFe, tmp1, i1b = elastic(Fv_old)
    tau = 2.0 * C10 * tmp1 * (Be - (i1b / 3.0) * I3)
    tau = tau + jnp.where(mtype > 0.5, 2.0 * C01 * mr_extra(Be, tmp1, i1b), 0.0)

    # eta may be exactly 0 (elastic); guard the division so the unused arm of
    # the where does not produce a NaN that would poison the derivative.
    eta_safe = jnp.where(eta > 1e-20, eta, 1.0)
    Fv_flow = (I3 + dt / (2.0 * eta_safe * detFe) * tau) @ Fv_old
    Fv_new = jnp.where(eta > 1e-20, Fv_flow, Fv_old)

    # recompute with updated Fv -> Cauchy stress
    Be, detFe, tmp1, i1b = elastic(Fv_new)
    press = (2.0 / D1) * (detFe - 1.0) * detFe
    sig = (2.0 * C10 * tmp1 * (Be - (i1b / 3.0) * I3) + press * I3) / detFe
    sig = sig + jnp.where(mtype > 0.5,
                          2.0 * C01 * mr_extra(Be, tmp1, i1b) / detFe, 0.0)
    return sig, Fv_new, detFe


def stress_core_jax(F, Fv_old, alpha, C10, C01, D1, eta, mtype, dt):
    """Return (stress[6] Abaqus Voigt, Fv_new[3,3], detFe)."""
    sig, Fv_new, detFe = _sigma_and_fv_jax(
        F, Fv_old, alpha, C10, C01, D1, eta, mtype, dt)
    sv = jnp.array([sig[i, j] for i, j in VOIGT])
    return sv, Fv_new, detFe


def _perturbation(v):
    """Symmetric spatial perturbation dE from a 6-vector, matching the
    convention material_server.dsde_perturbation and usermat_biofilm.f both
    use: component k contributes 1/2 to each of (a,b) and (b,a), so a normal
    component (a==b) lands at 1.0 and a shear at 1/2 on each off-diagonal.
    Differentiating w.r.t. v at v=0 therefore reproduces exactly what their
    (sigma(h) - sigma(0))/h estimates."""
    dE = jnp.zeros((3, 3))
    for k, (a, b) in enumerate(VOIGT):
        dE = dE.at[a, b].add(0.5 * v[k])
        dE = dE.at[b, a].add(0.5 * v[k])
    return dE


def dsde_exact(F, Fv_old, params):
    """6x6 material Jacobian by forward-mode AD -- no finite-difference step.

    Same definition as material_server.dsde_perturbation (spatial perturbation
    of F, Fv held at its entry value), evaluated exactly.
    """
    F = jnp.asarray(F, dtype=jnp.float64)
    Fv_old = jnp.asarray(Fv_old, dtype=jnp.float64)

    def sigma_of(v):
        sv, _, _ = stress_core_jax((I3 + _perturbation(v)) @ F, Fv_old, *params)
        return sv

    return jax.jacfwd(sigma_of)(jnp.zeros(6))


# jit once; params are traced values so a single compile serves every call.
dsde_exact_jit = jax.jit(dsde_exact, static_argnums=())


# Parameters d(sigma)/d(theta) is taken with respect to, in this order.
SENS_PARAMS = ("alpha", "C10", "C01", "D1", "eta")


def dsigma_dparams(F, Fv_old, params):
    """d(sigma_k)/d(theta_j) as a (6, 5) array, theta = SENS_PARAMS.

    This is the capability the finite-difference path cannot practically
    match: propagating the TMCMC posterior needs the stress response to the
    calibrated parameters, and getting that by differencing would mean a
    sweep of extra solves (or extra socket round-trips) per Gauss point. Here
    it falls out of the same evaluation.

    mtype and dt are held fixed -- mtype is a discrete model switch, and dt is
    the solver's increment rather than a material parameter.
    """
    F = jnp.asarray(F, dtype=jnp.float64)
    Fv_old = jnp.asarray(Fv_old, dtype=jnp.float64)
    alpha, C10, C01, D1, eta, mtype, dt = params

    def sigma_of(theta):
        sv, _, _ = stress_core_jax(F, Fv_old, theta[0], theta[1], theta[2],
                                   theta[3], theta[4], mtype, dt)
        return sv

    return jax.jacfwd(sigma_of)(jnp.array([alpha, C10, C01, D1, eta]))

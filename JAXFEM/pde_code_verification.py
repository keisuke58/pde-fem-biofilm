#!/usr/bin/env python3
"""Verify the PDE discretisations against exact solutions.

This is *code verification* in the V&V sense: it asks whether the discretisation
solves the equations it claims to solve, and at what order. It needs no
experimental data and no ANSYS, which is why it is worth doing -- every
validation claim in this repository rests on the solver being right, and until
now the only check on the nutrient PDE was
`tests/test_pde_uniform_consistency.py`, which compares the 1D and 2D codes to
the 0D integrator on a *uniform* field. A uniform field makes the Laplacian
vanish identically, so that test cannot see the spatial operator or either
boundary condition at all.

Nothing here imports the solver's discretisation to build its reference. The
exact solutions below are derived from the continuum statement of the problem,
so the comparison is a real one (the same reasoning as
`ansys_usermat/apdl/closed_form_reference.py`).

------------------------------------------------------------------------------
The problem
------------------------------------------------------------------------------

`core_hamilton_1d_nutrient.nutrient_step` advances

    dc/dt = D_c c'' - g_eff * phi * c/(k_monod + c)      on 0 <= x <= L
    c'(0) = 0        (Neumann, the tooth surface)
    c(L)  = 1        (Dirichlet, the gingival crevicular fluid)

The Monod factor has no closed-form steady state, but its two asymptotic limits
do, and both are reachable exactly as written by choosing `k_monod`.

**First-order limit** (`c << k_monod`, so `c/(k+c) -> c/k`). The steady problem
is `D c'' = (g phi/k) c`, i.e. `c'' = (theta/L)^2 c` with the Thiele modulus

    theta = L sqrt(g_eff phi / (k_monod D_c))

The solution satisfying `c'(0) = 0` and `c(L) = 1` is

    c(x) = cosh(theta x / L) / cosh(theta)

This is the classical reaction-diffusion slab profile that the Thiele modulus is
defined by -- and `core_hamilton_1d_nutrient`'s own docstring already quotes a
Thiele number for its default parameters, without ever having checked the
profile.

**Zeroth-order limit** (`c >> k_monod`, so `c/(k+c) -> 1`). The steady problem is
`D c'' = g phi`, a constant, giving the parabola

    c(x) = 1 - (g_eff phi / (2 D_c)) (L^2 - x^2)

This limit is not a convenience: zeroth order in `c` is **Klempt et al. (2024)
Eq. 35** (`c_dot - d grad^2 c + g phi = 0`, fixed by their Eq. 24 as "the
simplest possible functional dependency, a linear relation"). So the parabola is
the exact steady solution of the nutrient equation of the paper this work
reproduces. See `KLEMPT2024_REPRODUCTION.md`.

------------------------------------------------------------------------------
What this found, and what was done about it
------------------------------------------------------------------------------

Both solvers wrote their zero-flux row as `(u1 - u0)/h^2`, half the correct
second difference, which places the wall half a cell outside the node. It cost
one order in 1D and made the 2D operator **inconsistent** at the wall -- the
error did not shrink with `h` at all -- leaking about 5 % of the species mass
through a boundary declared zero-flux.

Fixed 2026-09-30 in all three places it appeared:

  * `core_hamilton_1d_nutrient.nutrient_step`          (1D nutrient, x=0)
  * `core_hamilton_2d_nutrient.laplacian_2d_neumann`   (2D species, 4 walls)
  * `core_hamilton_2d_nutrient._make_nutrient_step_mixed` (2D nutrient, 3 walls)

The `_half` variants below preserve the pre-fix arithmetic so the measurement
can still be reproduced, and `_mirror` variants are asserted to reproduce the
now-committed routines bit for bit. `PDE_VERIFICATION_FINDINGS.md` has the
numbers and `JAXFEM/boundary_fix_impact.py` the effect on the reported 2D
condition spread; `tests/test_pde_code_verification.py` pins all of it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import jax                                                    # noqa: E402
import jax.numpy as jnp                                       # noqa: E402

jax.config.update("jax_enable_x64", True)

from JAXFEM import core_hamilton_1d_nutrient as p1            # noqa: E402

L = 1.0


# ---------------------------------------------------------------------------
# exact solutions -- derived above, independent of the implementation
# ---------------------------------------------------------------------------

def exact_first_order(x, D_c, g_eff, k_monod, phi):
    """cosh profile; theta is the Thiele modulus."""
    theta = L * np.sqrt(g_eff * phi / (k_monod * D_c))
    return np.cosh(theta * x / L) / np.cosh(theta)


def exact_zeroth_order(x, D_c, g_eff, phi):
    """parabola; this is Klempt 2024 Eq. 35's steady state."""
    return 1.0 - (g_eff * phi / (2.0 * D_c)) * (L * L - x * x)


# ---------------------------------------------------------------------------
# driving the repository's own step to steady state
# ---------------------------------------------------------------------------

def _params(N, D_c, g_eff, k_monod, cfl=0.25):
    dx = L / (N - 1)
    dt = cfl * dx * dx / D_c
    return {"D_c": D_c, "g_eff": g_eff, "k_monod": k_monod, "dx": dx,
            "dt_h": dt, "n_react_sub": 1, "n_sub_c": 1}, dx, dt


def steady_state(N, D_c, g_eff, k_monod, phi, step=None, n_diff_times=20.0):
    """Iterate to steady state and return (x, c, residual).

    `step` defaults to the repository's `nutrient_step`; passing the variant
    below swaps only the Neumann coefficient.

    The loop is compiled once per grid with `lax.fori_loop`. Calling the step
    from Python instead re-traces it every iteration, which is what made an
    earlier version of this script take longer than the solve it was timing.

    `n_diff_times` is in units of the diffusion time `L^2/D_c`, so the step
    count scales with the grid the way the explicit stability limit demands.
    The returned residual `max|c_{n+1} - c_n|/dt` reports whether steady state
    was actually reached, rather than assuming it.
    """
    if step is None:
        step = p1.nutrient_step
    params, dx, dt = _params(N, D_c, g_eff, k_monod)
    n_steps = int(n_diff_times * L * L / D_c / dt)

    @jax.jit
    def run(c0):
        return jax.lax.fori_loop(0, n_steps, lambda _, c: step(c, phi, params),
                                 c0)

    c = run(jnp.ones(N))
    residual = float(jnp.max(jnp.abs(step(c, phi, params) - c))) / dt
    return np.linspace(0.0, L, N), np.asarray(c), residual


# ---------------------------------------------------------------------------
# the candidate fix, as one changed coefficient
# ---------------------------------------------------------------------------

def nutrient_step_mirror(c_field, phi_total, params):
    """`nutrient_step`'s Neumann row at its ghost-node value -- now the
    committed one, which a test below asserts bit for bit.

    A ghost node `c_{-1} = c_1` (the reflecting wall `c'(0) = 0`) gives

        c''(0) ~ (c_{-1} - 2 c_0 + c_1)/dx^2 = 2 (c_1 - c_0)/dx^2

    and the finite-volume reading of the same boundary -- a half-width control
    volume at node 0, flux `D (c_1 - c_0)/dx` in over a cell of width `dx/2` --
    gives the identical factor 2. `nutrient_step` writes
    `(c_1 - c_0)/dx^2`, half of both, while its comment says "ghost node
    approach".

    Everything else here is copied unchanged, and
    `tests/test_nutrient_verification.py` asserts that setting `NEUMANN = 1.0`
    reproduces `nutrient_step` bit for bit -- so this really is that routine
    plus one coefficient, not a second opinion about the whole scheme.
    """
    return _step_impl(c_field, phi_total, params, 2.0)


def nutrient_step_half(c_field, phi_total, params):
    """The same body with the pre-2026-09-30 coefficient, kept so the
    measurement that motivated the fix stays reproducible."""
    return _step_impl(c_field, phi_total, params, 1.0)


def _step_impl(c_field, phi_total, params, neumann):
    D_c = params["D_c"]
    g_eff = params["g_eff"]
    k_monod = params["k_monod"]
    dx = params["dx"]
    dt = params["dt_h"] * params["n_react_sub"]
    n_sub_c = params["n_sub_c"]
    dt_c = dt / n_sub_c

    def sub_step(c, _):
        lap = jnp.zeros_like(c)
        interior = (c[:-2] + c[2:] - 2.0 * c[1:-1]) / (dx * dx)
        lap = lap.at[1:-1].set(interior)
        lap = lap.at[0].set(neumann * (c[1] - c[0]) / (dx * dx))
        consumption = g_eff * phi_total * c / (k_monod + c + 1e-12)
        c_new = c + dt_c * (D_c * lap - consumption)
        c_new = jnp.clip(c_new, 0.0, None)
        c_new = c_new.at[-1].set(1.0)
        return c_new, None

    c_final, _ = jax.lax.scan(sub_step, c_field, jnp.arange(n_sub_c))
    return c_final


# ---------------------------------------------------------------------------
# the two studies
# ---------------------------------------------------------------------------

CASES = {
    # theta = 2: a profile with real curvature, nutrient reaching the wall
    "first-order (Thiele theta=2)": dict(
        D_c=1.0, k_monod=1.0e6, g_eff=4.0e6, phi=1.0, kind="first"),
    # Klempt 2024 Eq. 35's own law; c(0) = 0.5 keeps it positive
    "zeroth-order (Klempt Eq. 35)": dict(
        D_c=1.0, k_monod=1.0e-8, g_eff=1.0, phi=1.0, kind="zeroth"),
}

GRIDS = (21, 41, 81, 161)


def _exact(x, cfg):
    if cfg["kind"] == "first":
        return exact_first_order(x, cfg["D_c"], cfg["g_eff"], cfg["k_monod"],
                                 cfg["phi"])
    return exact_zeroth_order(x, cfg["D_c"], cfg["g_eff"], cfg["phi"])


def study(step=None, grids=None):
    """sup-norm error against the exact solution, per case and grid."""
    out = {}
    for label, cfg in CASES.items():
        rows = []
        for N in (grids or GRIDS):
            x, c, res = steady_state(N, cfg["D_c"], cfg["g_eff"],
                                     cfg["k_monod"], cfg["phi"], step=step)
            ce = _exact(x, cfg)
            rows.append((N, L / (N - 1),
                         float(np.max(np.abs(c - ce))),
                         float(abs(c[0] - ce[0])),
                         float(np.max(np.abs(c[1:] - ce[1:]))), res))
        out[label] = rows
    return out


def orders(rows):
    """observed order between successive grids, from the sup-norm column."""
    return [np.log(a[2] / b[2]) / np.log(a[1] / b[1])
            for a, b in zip(rows[:-1], rows[1:])]


def main():
    for name, step in (("as committed (wall on the node)", None),
                       ("pre-fix (wall half a cell out)", nutrient_step_half)):
        print(f"\n=== {name} " + "=" * (52 - len(name)))
        res = study(step)
        for label, rows in res.items():
            print(f"\n{label}")
            print(f"  {'N':>5} {'h':>9} {'sup err':>11} {'err@x=0':>11} "
                  f"{'sup err x>0':>12} {'residual':>10}")
            for N, h, e, e0, ei, res in rows:
                print(f"  {N:>5} {h:>9.5f} {e:>11.3e} {e0:>11.3e} "
                      f"{ei:>12.3e} {res:>10.1e}")
            print("  observed order: " +
                  ", ".join(f"{p:.2f}" for p in orders(rows)))


if __name__ == "__main__":
    main()
    main_2d()


# ---------------------------------------------------------------------------
# the 2D operator carries the same boundary treatment, on all four walls
# ---------------------------------------------------------------------------
#
# `core_hamilton_2d_nutrient.laplacian_2d_neumann` -- which drives the *species*
# diffusion in the 2D code (`diffusion_step`), and so the field
# `condition_spread_2d.py` measures -- writes the identical rows:
#
#     lap_x[0, :]  = (u[1, :] - u[0, :]) / dx^2        and three more
#
# Its comment says "ghost u[-1] = u[0]", which is exactly what the arithmetic
# does, so this is a deliberate choice rather than a slip. It *is* zero-flux --
# but for a wall sitting half a cell outside the node, which costs an order.
# Reflecting about the node instead (`ghost u[-1] = u[1]`) puts the wall on the
# node and restores second order.
#
# Verified on the operator directly, against an exact eigenfunction: for
# zero-flux walls on [0, Lx] x [0, Ly] the eigenfunctions are
# cos(n pi x/Lx) cos(m pi y/Ly) with eigenvalue -((n pi/Lx)^2 + (m pi/Ly)^2).

def laplacian_2d_neumann_mirror(u, dx, dy):
    """`laplacian_2d_neumann` with the walls reflected about the node."""
    lap_x = jnp.zeros_like(u)
    lap_x = lap_x.at[1:-1, :].set(
        (u[:-2, :] + u[2:, :] - 2.0 * u[1:-1, :]) / (dx * dx))
    lap_x = lap_x.at[0, :].set(2.0 * (u[1, :] - u[0, :]) / (dx * dx))
    lap_x = lap_x.at[-1, :].set(2.0 * (u[-2, :] - u[-1, :]) / (dx * dx))
    lap_y = jnp.zeros_like(u)
    lap_y = lap_y.at[:, 1:-1].set(
        (u[:, :-2] + u[:, 2:] - 2.0 * u[:, 1:-1]) / (dy * dy))
    lap_y = lap_y.at[:, 0].set(2.0 * (u[:, 1] - u[:, 0]) / (dy * dy))
    lap_y = lap_y.at[:, -1].set(2.0 * (u[:, -2] - u[:, -1]) / (dy * dy))
    return lap_x + lap_y


def laplacian_2d_neumann_half(u, dx, dy):
    """`laplacian_2d_neumann` as it stood before 2026-09-30: ghost = boundary
    node, so every wall row carries half the correct second difference."""
    lap_x = jnp.zeros_like(u)
    lap_x = lap_x.at[1:-1, :].set(
        (u[:-2, :] + u[2:, :] - 2.0 * u[1:-1, :]) / (dx * dx))
    lap_x = lap_x.at[0, :].set((u[1, :] - u[0, :]) / (dx * dx))
    lap_x = lap_x.at[-1, :].set((u[-2, :] - u[-1, :]) / (dx * dx))
    lap_y = jnp.zeros_like(u)
    lap_y = lap_y.at[:, 1:-1].set(
        (u[:, :-2] + u[:, 2:] - 2.0 * u[:, 1:-1]) / (dy * dy))
    lap_y = lap_y.at[:, 0].set((u[:, 1] - u[:, 0]) / (dy * dy))
    lap_y = lap_y.at[:, -1].set((u[:, -2] - u[:, -1]) / (dy * dy))
    return lap_x + lap_y


def study_2d_operator(op=None, grids=GRIDS):
    """sup-norm error of the 2D Laplacian on an exact Neumann eigenfunction."""
    from JAXFEM import core_hamilton_2d_nutrient as p2
    if op is None:
        op = p2.laplacian_2d_neumann
    rows = []
    for N in grids:
        h = L / (N - 1)
        g = np.linspace(0.0, L, N)
        X, Y = np.meshgrid(g, g, indexing="ij")
        u = np.cos(np.pi * X / L) * np.cos(np.pi * Y / L)
        exact = -2.0 * (np.pi / L) ** 2 * u
        got = np.asarray(op(jnp.asarray(u), h, h))
        rows.append((N, h, float(np.max(np.abs(got - exact))),
                     float(np.max(np.abs(got[1:-1, 1:-1] - exact[1:-1, 1:-1])))))
    return rows


def main_2d():
    from JAXFEM import core_hamilton_2d_nutrient as p2
    for name, op in (("as committed (mirrored about the node)",
                      p2.laplacian_2d_neumann),
                     ("pre-fix (ghost = boundary node)",
                      laplacian_2d_neumann_half)):
        rows = study_2d_operator(op)
        print(f"\n2D Laplacian, Neumann eigenfunction -- {name}")
        print(f"  {'N':>5} {'h':>9} {'sup err':>11} {'interior only':>14}")
        for N, h, e, ei in rows:
            print(f"  {N:>5} {h:>9.5f} {e:>11.3e} {ei:>14.3e}")
        ords = [np.log(a[2] / b[2]) / np.log(a[1] / b[1])
                for a, b in zip(rows[:-1], rows[1:])]
        print("  observed order: " + ", ".join(f"{p:.2f}" for p in ords))

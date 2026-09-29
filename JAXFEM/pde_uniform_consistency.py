#!/usr/bin/env python3
"""pde_uniform_consistency.py -- does the spatial Hamilton+nutrient PDE reduce
to the verified 0D integrator when nothing varies in space?

With a spatially uniform initial state and no nutrient consumption, the
diffusion terms vanish and every node must follow the 0D material-point
trajectory (jax_hamilton_0d_5species_demo.newton_step, the integrator
klempt2026_reproduction.py checks against Klempt et al. 2026) with the same
c*. This runs the PDE codes' own building blocks, in their own step order
(reaction sub-steps -> species diffusion -> nutrient), from a uniform state,
and compares every node with the 0D chain at the calibration c* and at c* = 1.

  1D  core_hamilton_1d_nutrient : reaction sees C_STAR * c (c = 1 at the boundary)
  2D  core_hamilton_2d_nutrient : reaction sees c_hamilton_scale * c (default C_STAR)

Before 2026-09-29 the 1D code passed c itself (c* = 1) and the 2D default
was 1 with one caller passing 100; both now take c* from ecology_constants.py.

    python JAXFEM/pde_uniform_consistency.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE.parent / "ansys_usermat" / "coupling"))

import ecology_jax  # noqa: E402
from ecology_constants import C_STAR, K_HILL, N_HILL  # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO, newton_step_jit, theta_to_matrices  # noqa: E402
from JAXFEM import core_hamilton_1d_nutrient as p1  # noqa: E402
from JAXFEM import core_hamilton_2d_nutrient as p2  # noqa: E402

DT_H = 1e-5
N_SUB = 20
N_MACRO = 10


def _hamilton_params(c):
    A, b = theta_to_matrices(jnp.asarray(THETA_DEMO, dtype=jnp.float64))
    return {"dt_h": DT_H, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5),
            "c": c, "alpha": 100.0, "K_hill": K_HILL, "n_hill": N_HILL,
            "A": A, "b_diag": b, "active_mask": jnp.ones(5, dtype=jnp.int64)}


def zero_d(c, n_steps):
    g = jnp.asarray(ecology_jax.default_initial_state())
    p = _hamilton_params(c)
    for _ in range(n_steps):
        g = newton_step_jit(g, p)
    return np.asarray(g)


def run_1d(n_nodes=11):
    """simulate_hamilton_1d_nutrient's loop body, uniform start, g_eff = 0.
    (params["c"] is ignored by the 1D reaction, which uses C_STAR * c_node.)"""
    p = _hamilton_params(C_STAR)
    p.update({"n_react_sub": N_SUB, "D_eff": jnp.full(5, 1e-3),
              "dx": 1.0 / (n_nodes - 1), "D_c": 1.0, "g_eff": 0.0,
              "k_monod": 1.0, "n_sub_c": 10})
    G = jnp.tile(jnp.asarray(ecology_jax.default_initial_state()), (n_nodes, 1))
    c = jnp.ones(n_nodes)
    for _ in range(N_MACRO):
        phi_total = G[:, 0:5].sum(axis=1)
        G = p1.reaction_step_c(G, c, p)
        G = p1.diffusion_step(G, p)
        c = p1.nutrient_step(c, phi_total, p)
    return np.asarray(G), np.asarray(c)


def run_2d(c_scale=None, nx=5, ny=4):
    """run_simulation_coupled's loop body, uniform start, zero consumption.
    c_scale=None uses run_simulation_coupled's own default."""
    if c_scale is None:
        import inspect
        c_scale = inspect.signature(p2.run_simulation_coupled).parameters["c_hamilton_scale"].default
    cfg = p2.Config2D(Nx=nx, Ny=ny, dt_h=DT_H, n_react_sub=N_SUB, n_macro=N_MACRO,
                      g_consumption=np.zeros(5))
    p = _hamilton_params(cfg.c_hamilton)
    reaction = p2._make_reaction_step_c(cfg.n_react_sub, cfg.newton_iters)
    nutrient = p2._make_nutrient_step_stable(10)
    G = jnp.tile(jnp.asarray(ecology_jax.default_initial_state()), (nx * ny, 1))
    c = jnp.full((nx, ny), cfg.c_boundary)
    for _ in range(N_MACRO):
        G = reaction(G, c.reshape(nx * ny) * c_scale, p)
        phi = p2.diffusion_step_species_2d(G.reshape(nx, ny, 12)[:, :, :5],
                                           jnp.asarray(cfg.D_eff), cfg.dt_macro, cfg.dx, cfg.dy)
        G2 = G.reshape(nx, ny, 12).at[:, :, :5].set(phi)
        G2 = G2.at[:, :, 5].set(1.0 - jnp.sum(phi, axis=-1))
        G = G2.reshape(nx * ny, 12)
        c = nutrient(c, phi, cfg.D_c, cfg.k_monod, jnp.asarray(cfg.g_consumption),
                     cfg.c_boundary, cfg.dx, cfg.dy, cfg.dt_macro)
    return np.asarray(G), np.asarray(c)


def compare():
    """Rows: (label, node spread, c spread, max|PDE - 0D(C_STAR)|, max|PDE - 0D(c*=1)|)."""
    n = N_MACRO * N_SUB
    ref = {c: zero_d(c, n) for c in (C_STAR, 1.0)}
    out = {"1D (default)": run_1d(), "2D (default scale)": run_2d()}
    rows = []
    for label, (G, c) in out.items():
        spread = float(np.max(np.abs(G - G[0])))
        rows.append((label, spread, float(np.ptp(c)),
                     float(np.max(np.abs(G - ref[C_STAR]))), float(np.max(np.abs(G - ref[1.0])))))
    return ref, rows


if __name__ == "__main__":
    ref, rows = compare()
    print(f"{N_MACRO} macro x {N_SUB} reaction sub-steps, dt_h = {DT_H:g}, THETA_DEMO, "
          f"C_STAR = {C_STAR:g}, K_HILL = {K_HILL:g}\n")
    print(f"0D phi_tot after {N_MACRO * N_SUB} steps: c*={C_STAR:g} -> "
          f"{ecology_jax.living_fraction_total(ref[C_STAR]):.6f}, c*=1 -> "
          f"{ecology_jax.living_fraction_total(ref[1.0]):.6f}\n")
    print(f"{'PDE run':22s} {'node spread':>12s} {'c spread':>10s} {'max|PDE-0D(C_STAR)|':>21s} {'max|PDE-0D(c*=1)|':>19s}")
    for label, spread, cs, d, d1 in rows:
        print(f"{label:22s} {spread:12.3e} {cs:10.3e} {d:21.3e} {d1:19.3e}")

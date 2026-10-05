#!/usr/bin/env python3
"""Composition into growth, estimate in Python: species-weighted growth law.

The growth law of the material routine, alpha_dot = k_alpha phi (Klempt 2024
Eq. 36), becomes

    alpha_dot = k_alpha phi sum_i f_i phi_i / (phi_1 + phi_2),
    f_1 = 1 + d, f_2 = 1 - d,

so the composition from the point model reaches the growth and the stress
directly, in the material routine alone (no front term needed). d = 1/3 is the
assumption of two_way_step2.py (rates scale like 1/eta_i, eta = 1, 2);
d = 0 is the one-way law. Not from a paper.

Model: the partner's 8^3 model as in two_way_step3.py (seed of 32 elements,
growth 1.1e-3 = k_alpha phi_cap T* in the one-way law, E = 10 Pa, nu = 0.49,
B-bar); the composition per seed element from the coupled scheme with the
local nutrient (composition_local_nutrient.py, Thiele number Lambda), and the
growth factor is the time mean of sum_i f_i share_i over the coupling steps.

    python ansys_usermat/two_way_growth_rate.py [--s 0.15]

Result (5 Oct 2026, coupling step 0.025), seed averages against the one-way
law (vM 5.703e-5 Pa, p -3.768e-5 Pa):
  s = 0.15                       growth factor   seed vM     seed p
  case 3, uniform c, d = 1/3     1.085           +8.5 %      +8.5 %
  case 6, uniform c, d = 1/3     0.794           -21 %       -21 %
  case 3, Lambda = 4, d = 1/3    1.052-1.074     -46 %       +6 %
  case 6, Lambda = 4, d = 1/3    0.928-0.988     +194 %      -3 %
  case 6, Lambda = 4, d = 0.5    0.893-0.982     +294 %      -5 %
  s = 0.25: case 6, Lambda = 4, d = 1/3: factor 0.831-0.943, vM +414 %, p -11 %.
With a uniform composition the growth law only rescales the stress (by the
growth factor). Where the composition varies across the seed, the growth does
too, and the incompatible growth sets the von Mises stress: in case 6 it rises
by a factor of 3 to 5, while the mean stress changes by a few per cent. This is
the same mechanism as the alpha step at the seed surface in the ANSYS runs.
On the 8^3 grid only; the von Mises values of a growth gradient are mesh
dependent (see BETA002_THESIS_MAP.md), the mean stress is the robust number.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "apdl"))
import composition_local_nutrient as L  # noqa: E402
import mesh_study_seed as S  # noqa: E402
import two_way_step3 as W3  # noqa: E402

T = L.T


def mean_share(case, cvals, s):
    """time mean (over the coupling steps) and final value of phi_1/(phi_1+phi_2)."""
    eco = T.eco
    T.ms.set_case(case)
    th, hp = T.ms.ECOLOGY_CASE["theta"], dict(T.ms.ECOLOGY_CASE["hp"])
    T.ms.set_case(None)
    dt_pm = s * T.DT
    n_sub = max(1, math.ceil(dt_pm / 1e-4 - 1e-9))
    tht = jnp.asarray(th, dtype=jnp.float64)
    eta = jnp.asarray(hp["eta"], dtype=jnp.float64)
    run = jax.jit(jax.vmap(lambda g, c: eco._substep_scan(
        g, tht, dt_pm / n_sub, jnp.arange(n_sub), eco.active_mask(2), c,
        jnp.float64(hp["alpha"]), eta)[0]))
    cv = np.asarray(cvals, float)
    G = T.seed_state(cv.size)
    p = np.full(cv.size, T.PHI_CAP)
    cst = jnp.asarray(hp["c"] * cv)
    acc = np.zeros(cv.size); nst = int(round(1.0 / T.DT))
    for _ in range(nst):
        G = T.rescale(np.asarray(run(jnp.asarray(T.rescale(G, p)), cst)), p)
        acc += G[:, 0] / (G[:, 0] + G[:, 1])
    return acc / nst, G[:, 0] / (G[:, 0] + G[:, 1])


def stress(growth_e, ins):
    old = S.GROWTH
    S.GROWTH = growth_e                       # solve() uses np.where(ins, S.GROWTH, 0)
    try:
        Ee = np.where(ins, (1 + S.FLOOR) * S.E, S.FLOOR * S.E)
        return W3.solve(Ee, ins)
    finally:
        S.GROWTH = old


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--s", type=float, default=0.15, help="time link s (default 0.15)")
    ap.add_argument("--dt", type=float, default=0.025)
    a = ap.parse_args(argv)
    T.DT = a.dt
    L.N, L.H = 8, 0.25
    Sm, Dm = L.grid_masks()
    ins = Sm.ravel()
    g0 = 1.1e-3
    vm0, p0 = stress(np.full(ins.size, g0), ins)
    print(f"one-way: seed vM {vm0[ins].mean():.3e} Pa, p {p0[ins].mean():.3e} Pa")
    for case in ("2sp_case3", "2sp_case6"):
        for lam in (0.0, 4.0):
            c = np.ones((8,) * 3) if lam == 0 else L.nutrient(Sm, Dm, lam)
            mbar, fin = mean_share(case, c[Sm], a.s)
            for d in (1 / 3, 0.5):
                fac = np.ones((8,) * 3)
                fac[Sm] = (1 + d) * mbar + (1 - d) * (1 - mbar)
                gr = g0 * fac.ravel()
                vm, p = stress(gr, ins)
                print(f"{case} Lambda={lam:g} d={d:.3g}: final share {fin.min():.3f}-{fin.max():.3f}, "
                      f"growth factor {fac[Sm].min():.3f}-{fac[Sm].max():.3f}; seed vM {vm[ins].mean():.3e} "
                      f"({vm[ins].mean() / vm0[ins].mean() - 1:+.1%}), p {p[ins].mean():.3e} "
                      f"({p[ins].mean() / p0[ins].mean() - 1:+.1%}), vM spread {vm[ins].std() / vm[ins].mean():.2f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Is the Gauss-point ecology ODE the same model TMCMC was calibrated with?

    python ansys_usermat/coupling/tmcmc_ode_crosscheck.py --tmcmc <tmcmc202601 clone>

Runs TMCMC's BiofilmNewtonSolver5S (tmcmc/program2602/improved_5species_jit.py)
and this repository's ecology_jax.ecology_substeps from the same state, theta
and constants (c* = 25, Hill gate off, alpha* = 0) to T* = 2, and prints:

1. TMCMC at several dt -- which step the calibration solver needs;
2. ecology_jax against TMCMC at a fine step -- whether the two agree.

Findings (2026-10-01), see ODE_TMCMC_CROSSCHECK.md: the residuals are the same
equations, but "Hill gate off" zeroes species 5's interaction here and leaves
it unchanged in TMCMC; and TMCMC's own dt = 0.01 is unstable for THETA_DEMO.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmcmc", required=True, help="tmcmc202601 clone")
    ap.add_argument("--T", type=float, default=2.0)
    a = ap.parse_args()
    sys.path.insert(0, str(Path(a.tmcmc) / "tmcmc" / "program2602"))
    from improved_5species_jit import BiofilmNewtonSolver5S
    import ecology_jax as eco
    from ecology_constants import C_STAR, K_HILL, N_HILL
    from jax_hamilton_0d_5species_demo import THETA_DEMO
    th = np.array(THETA_DEMO, float)

    def tmcmc(dt):
        s = BiofilmNewtonSolver5S(dt=dt, maxtimestep=int(round(a.T / dt)),
                                  c_const=C_STAR, phi_init=0.01,
                                  K_hill=K_HILL, n_hill=N_HILL)
        return s.run_deterministic(th)[1]

    print(f"T* = {a.T}, c* = {C_STAR}, K_hill = {K_HILL}")
    for dt in (1e-2, 1e-3, 1e-4, 2.5e-5):
        g = tmcmc(dt)
        print(f"TMCMC dt={dt:<7g} phi={np.round(g[-1, :5], 4)} "
              f"min phi0={g[:, 5].min():+.3f}")
    ref = g[-1]
    g0 = np.zeros(12)
    g0[:5], g0[5], g0[6:11] = 0.01, 0.95, 0.999
    for dt in (1e-3, 1e-4):
        g, _ = eco.ecology_substeps(g0, th, a.T, int(round(a.T / dt)), 5)
        g = np.asarray(g)
        print(f"ours  dt={dt:<7g} phi={np.round(g[:5], 4)} "
              f"psi5={g[10]:.4f} max|diff|={np.max(np.abs(g[:11] - ref[:11])):.2e}")


if __name__ == "__main__":
    main()

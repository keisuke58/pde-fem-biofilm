#!/usr/bin/env python3
"""Sensitivity of the composition to s, the point model's clock (prop(28) = 7).

No paper links Klempt 2024's normalised T* (in [0, 1]) to the time unit of
Klempt et al. 2026, so the point model advances s * dt per substep and s is
an explicit assumption. This sweeps s = 0.05 / 0.15 / 0.5 over T* = 1 for
the amounts phi_3D(t) the 3D field produces -- seed, interior, a front
passing through, the void -- with the exact scheme the fragment runs
(composition_reference.reference, which the mock matches bit for bit).

Stress is not swept: in this mode growth is Eq. 36 with phi_3D only, so alpha
and the stress are identical for every s (pinned in
tests/test_composition_fragment.py). Only the composition moves.

    python ansys_usermat/composition_s_sweep.py [--dt 0.1]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "coupling"))
import composition_reference as cref                    # noqa: E402
import material_server as ms                            # noqa: E402

K_ALPHA, PHI_MIN = 1.0e-3, 1.0e-2
SCENARIOS = {
    "seed (phi_3D = 1)": lambda t: 1.0,
    "interior 0.9": lambda t: 0.9,
    "interior 0.4": lambda t: 0.4,
    "front 0.05 -> 1": lambda t: 0.05 + 0.95 * t,
    "void (K_LOCAL t)": lambda t: K_ALPHA * t,
}


def sweep(case, s_values, dt, T=1.0):
    ms.set_case(case)
    th, hp = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"]
    n = int(round(T / dt))
    out = {}
    for name, f in SCENARIOS.items():
        series = [f(k * dt) for k in range(1, n + 1)]
        for s in s_values:
            g = cref.reference(series, th, hp, dt, phi_min=PHI_MIN, s=s)[-1]
            out[(name, s)] = (g[0] / (g[0] + g[1]), g[6], g[7])
    ms.set_case(None)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--s", type=float, nargs="+", default=[0.05, 0.15, 0.5])
    a = ap.parse_args()
    for case in ("2sp_case3", "2sp_case6"):
        r = sweep(case, a.s, a.dt)
        print(f"\n{case}, deck dt = {a.dt}, T* = 1, phi_min = {PHI_MIN}: "
              "chi_1 (psi_1, psi_2) at T* = 1")
        print(f"{'':20s}" + "".join(f"{'s = ' + str(s):>26s}" for s in a.s))
        for name in SCENARIOS:
            print(f"{name:20s}" + "".join(
                f"{r[(name, s)][0]:>10.4f} ({r[(name, s)][1]:.3f}, "
                f"{r[(name, s)][2]:.3f})" for s in a.s))


if __name__ == "__main__":
    main()

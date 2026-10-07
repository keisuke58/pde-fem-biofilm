#!/usr/bin/env python3
"""Soleimani 2019 Eq. 32 against the exact solution of the ODE it was derived for.

`VISCOUS_UPDATE_SCHEME.md` specifies replacing this repository's explicit
viscous increment with an exponential update, and records why the change was
attempted and reverted. What it did not have was a measurement of the scheme
itself. This supplies one, against closed form, in about a second, without
ANSYS and without touching the verified UMAT.

------------------------------------------------------------------------------
The equation
------------------------------------------------------------------------------

Soleimani (2019) §2 puts viscosity in the free energy as an internal variable
rather than in a third multiplicative factor, and its evolution is the 1D
Maxwell spring-dashpot ODE of his Eq. 31:

    Q' + Q/tau = S'          S := dPsi_inf/dC   (the equilibrium driver)

His Eq. 32 integrates it by a mid-point scheme as a *recursive exponential*:

    Q_{n+1} = exp(-dt/tau) Q_n + exp(-dt/(2 tau)) (S_{n+1} - S_n)

Everything below is the scalar case, which is not a simplification chosen for
convenience -- it is the equation Eq. 32 is derived from. `SOLEIMANI2019_NOTES.md`
§2 has the reading.

For comparison, the forward-Euler update of the same ODE -- the family this
repository's `Fv` increment belongs to, with the driver taken at the old state:

    Q_{n+1} = Q_n + dt (S' - Q_n/tau)

------------------------------------------------------------------------------
The exact solutions
------------------------------------------------------------------------------

Both are standard and are written down here rather than taken from any code.

**Relaxation.** Hold the strain, so S' = 0, from Q(0) = Q0:

    Q(t) = Q0 exp(-t/tau)

**Ramp.** Load at a constant rate, S' = r, from Q(0) = 0:

    Q(t) = r tau (1 - exp(-t/tau))

The relaxation case is the interesting one for the stability question, because
Eq. 32 reduces on it to Q_{n+1} = exp(-dt/tau) Q_n, whose n-th iterate is
exp(-n dt/tau) Q0 -- the exact answer, **at every step size, with no error at
all**. Forward Euler gives (1 - dt/tau)^n Q0, which alternates in sign past
dt/tau = 1 and grows without bound past dt/tau = 2.

The ramp case separates the two on accuracy rather than stability: Eq. 32 is a
mid-point scheme and should show second order, forward Euler first.
"""
from __future__ import annotations

import numpy as np

TAU = 0.01          # Soleimani 2019 Table 2's relaxation time
RATE = 1.0          # S' for the ramp case
Q0 = 1.0            # initial overstress for the relaxation case


# ---------------------------------------------------------------------------
# exact solutions
# ---------------------------------------------------------------------------

def relax_exact(t, tau=TAU, q0=Q0):
    return q0 * np.exp(-t / tau)


def ramp_exact(t, tau=TAU, rate=RATE):
    return rate * tau * (1.0 - np.exp(-t / tau))


# ---------------------------------------------------------------------------
# the two integrators
# ---------------------------------------------------------------------------

def exponential_update(q, dS, dt, tau=TAU):
    """Soleimani 2019 Eq. 32."""
    return np.exp(-dt / tau) * q + np.exp(-dt / (2.0 * tau)) * dS


def forward_euler(q, dS, dt, tau=TAU):
    """The explicit family this repository's Fv increment belongs to."""
    return q + dS - dt * q / tau


def _march(step, n, dt, q_init, rate, tau=TAU):
    q = q_init
    for _ in range(n):
        q = step(q, rate * dt, dt, tau)
    return q


# ---------------------------------------------------------------------------
# study 1: relaxation -- Eq. 32 is exact at every step size
# ---------------------------------------------------------------------------

T_END = 5.0 * TAU


def relaxation_errors(ratios=(0.1, 0.5, 1.0, 2.0, 10.0, 100.0)):
    """Relative error at t = 5 tau, per dt/tau, for both integrators.

    `dt/tau >= 1` cannot be reached by taking whole steps to 5 tau, so each
    row marches to the nearest whole number of steps and compares against the
    exact solution *at that same time* -- the comparison stays honest at
    coarse steps instead of drifting off the target time.
    """
    rows = []
    for r in ratios:
        dt = r * TAU
        n = max(1, int(round(T_END / dt)))
        t = n * dt
        ex = relax_exact(t)
        qe = _march(exponential_update, n, dt, Q0, 0.0)
        qf = _march(forward_euler, n, dt, Q0, 0.0)
        rows.append((r, n, t / TAU, ex,
                     abs(qe - ex) / abs(ex), abs(qf - ex) / abs(ex)))
    return rows


# ---------------------------------------------------------------------------
# study 2: ramp -- second order against first
# ---------------------------------------------------------------------------

def ramp_orders(n_steps=(10, 20, 40, 80, 160)):
    """sup error at t = 5 tau under refinement, and the observed order."""
    rows = []
    for n in n_steps:
        dt = T_END / n
        ex = ramp_exact(T_END)
        qe = _march(exponential_update, n, dt, 0.0, RATE)
        qf = _march(forward_euler, n, dt, 0.0, RATE)
        rows.append((n, dt / TAU, abs(qe - ex), abs(qf - ex)))
    return rows


def _order(rows, col):
    return [np.log(a[col] / b[col]) / np.log(a[1] / b[1])
            for a, b in zip(rows[:-1], rows[1:])]


def main():
    print(f"tau = {TAU}  (Soleimani 2019 Table 2),  target t = 5 tau\n")
    print("=== relaxation, S' = 0: exact answer is Q0 exp(-t/tau) ===")
    print(f"{'dt/tau':>8} {'steps':>6} {'t/tau':>7} {'exact':>12} "
          f"{'Eq.32 rel err':>15} {'fwd Euler rel err':>19}")
    for r, n, tn, ex, ee, ef in relaxation_errors():
        print(f"{r:>8.1f} {n:>6} {tn:>7.1f} {ex:>12.4e} {ee:>15.2e} {ef:>19.2e}")
    print("\n  Eq. 32 reduces to Q_{n+1} = exp(-dt/tau) Q_n here, so it is the")
    print("  exact solution sampled at t_n -- no step restriction of any kind.")

    print("\n=== ramp, S' = r: exact answer is r tau (1 - exp(-t/tau)) ===")
    rows = ramp_orders()
    print(f"{'steps':>6} {'dt/tau':>8} {'Eq.32 err':>12} {'fwd Euler err':>15}")
    for n, r, ee, ef in rows:
        print(f"{n:>6} {r:>8.3f} {ee:>12.3e} {ef:>15.3e}")
    print("  Eq. 32   observed order: "
          + ", ".join(f"{p:.2f}" for p in _order(rows, 2)))
    print("  fwd Euler observed order: "
          + ", ".join(f"{p:.2f}" for p in _order(rows, 3)))


if __name__ == "__main__":
    main()

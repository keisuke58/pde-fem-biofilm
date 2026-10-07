"""Soleimani 2019 Eq. 32 against the exact solution of the ODE it came from.

Pins `ansys_usermat/viscous_integrator_verification.py`. See
`VISCOUS_UPDATE_SCHEME.md` for what the measurement is for: it substantiates
the specification of a change that was attempted, reverted, and left as a
specification, without touching the verified UMAT.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ansys_usermat"))

import viscous_integrator_verification as vi   # noqa: E402


def _order(rows, col):
    return [np.log(a[col] / b[col]) / np.log(a[1] / b[1])
            for a, b in zip(rows[:-1], rows[1:])]


def test_the_exponential_update_is_exact_on_relaxation_at_every_step_size():
    """Eq. 32 reduces to Q_{n+1} = exp(-dt/tau) Q_n when the strain is held,
    so its n-th iterate IS the exact solution sampled at t_n. Checked out to
    dt/tau = 100, which is the paper's unconditional-stability claim made
    quantitative."""
    for ratio, n, t_over_tau, exact, err_exp, _err_fe in vi.relaxation_errors():
        assert err_exp < 1.0e-14, (ratio, err_exp)


def test_forward_euler_fails_the_same_test_the_way_the_theory_says():
    """(1 - dt/tau)^n: already 23 % off at dt/tau = 0.1, useless at 1, and
    divergent past 2. This is the defect VISCOUS_UPDATE_SCHEME.md documents."""
    rows = {r[0]: r for r in vi.relaxation_errors()}
    assert rows[0.1][5] > 0.2          # 23 % at a step a tenth of tau
    assert rows[1.0][5] > 0.9          # no information left at dt = tau
    assert rows[2.0][5] > 10.0         # diverging
    assert rows[100.0][5] > 1.0e40     # spectacularly so


def test_the_exponential_update_is_second_order_on_a_ramp():
    rows = vi.ramp_orders()
    for p in _order(rows, 2)[-1:]:
        assert 1.9 < p < 2.1, p


def test_forward_euler_is_first_order_on_the_same_ramp():
    rows = vi.ramp_orders()
    for p in _order(rows, 3)[-1:]:
        assert 0.9 < p < 1.1, p


def test_the_exponential_update_is_not_uniformly_more_accurate():
    """Worth pinning so the case for the change stays honest: on the ramp at
    dt/tau = 0.5 forward Euler is the closer of the two. The exponential
    update's advantage is order and unconditional stability, not a smaller
    error at every step size."""
    coarse = vi.ramp_orders()[0]
    assert coarse[3] < coarse[2], coarse

#!/usr/bin/env python3
"""Closed-form Cauchy stress for the two growth verification cases.

This module exists to break a circularity. `reference_values.json` calls
itself "Closed-form Cauchy stress", but it is produced by `make_reference.py`
calling `stress_core` -- the implementation. It therefore records what the code
does rather than testing it, and no amount of agreement with it means anything.

Nothing here imports the implementation. Everything is derived from the
continuum statement of the model, so the comparison in
`tests/test_closed_form_reference.py` is a real one.

------------------------------------------------------------------------------
The two cases
------------------------------------------------------------------------------

Both use isotropic growth, Fg = (1+a) I, with Fv = I initially.

**Traction-free.** The body simply expands, F = (1+a) I, so Fe = I exactly:
J_e = 1 and the stress vanishes identically, for any material constants and any
viscosity. There is nothing to drive viscous flow, so Fv stays I.

**Fully constrained.** F = I, so Fe = I/(1+a), which is *isotropic*. For any
hyperelastic potential with an isochoric/volumetric split, the deviatoric
Kirchhoff stress of an isotropic elastic state is identically zero. Only the
volumetric term survives:

    J_e = (1+a)^-3           sigma = (2/D1) (J_e - 1) I

and since the deviatoric flow driver vanishes, Fv stays I here too -- so the
answer does not depend on eta. That last point is the sharp one: a correct
deviatoric dashpot must give the *same* stress at eta = 0 and eta > 0.

------------------------------------------------------------------------------
The known discrepancy
------------------------------------------------------------------------------

The implementation does not reproduce this, for the reason set out in
`DEVIATOR_SCALING_FINDING.md`: its isochoric split applies J^(-2/3) to the
subtracted trace but not to the tensor beside it. For the constrained case that
leaves a spurious, purely spherical term which `spurious_term` below predicts
in closed form. Predicting it rather than merely observing it is what makes the
test meaningful -- it pins the analysis, not just the number.
"""
import numpy as np

I3 = np.eye(3)


def elastic_jacobian(alpha: float) -> float:
    """J_e for the fully constrained case, F = I."""
    return 1.0 / (1.0 + alpha) ** 3


def constrained_stress(alpha: float, D1: float) -> np.ndarray:
    """Closed-form Cauchy stress, fully constrained isotropic growth.

    Independent of C10, C01, mtype and eta: the deviator of an isotropic
    elastic state vanishes, so only the volumetric term survives.
    Returned in ANSYS Voigt order (11,22,33,12,23,13).
    """
    J = elastic_jacobian(alpha)
    p = (2.0 / D1) * (J - 1.0)
    return np.array([p, p, p, 0.0, 0.0, 0.0])


def free_growth_stress() -> np.ndarray:
    """Closed-form Cauchy stress, traction-free isotropic growth: zero."""
    return np.zeros(6)


def spurious_term(alpha: float, C10: float) -> np.ndarray:
    """The purely spherical term the implementation adds, in closed form.

    From DEVIATOR_SCALING_FINDING.md §1, specialised to the constrained case.
    With b = (1+a)^-2 I and J^(-2/3) = (1+a)^2, the coded "deviator" reduces to

        tau = 2 C10 [ 1 - (1+a)^2 ] I     ->     sigma = tau / J_e

    which is zero only at a = 0.
    """
    s = 2.0 * C10 * (1.0 - (1.0 + alpha) ** 2) / elastic_jacobian(alpha)
    return np.array([s, s, s, 0.0, 0.0, 0.0])


if __name__ == "__main__":
    C10, D1 = 2.0e-4, 5000.0
    print(f"{'alpha':>7} {'Je':>10} {'closed form':>14} {'+ spurious':>14} "
          f"{'ratio':>7}")
    for a in (0.05, 0.20):
        cf = constrained_stress(a, D1)[0]
        sp = spurious_term(a, C10)[0]
        print(f"{a:7.2f} {elastic_jacobian(a):10.6f} {cf:14.6e} "
              f"{cf + sp:14.6e} {(cf + sp) / cf:7.3f}")

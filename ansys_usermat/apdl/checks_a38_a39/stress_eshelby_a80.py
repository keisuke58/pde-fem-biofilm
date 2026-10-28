"""A-80: the stress solver of the 3D videos and the coupled runs (wp3_stress_reference.solve_fields: hex8, B-bar,
E(phi) = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain (alpha - 1) I, free cube with rigid-body constraints) against
closed forms.

  1. uniform eigenstrain in a free cube: no stress;
  2. Eshelby: a sphere of radius a with the dilatational eigenstrain e* in a matrix of the same stiffness. Inside, the
     stress is uniform and hydrostatic, sigma = -2 E e* / (3 (1 - nu)) (infinite matrix); in a free sphere of radius b
     the magnitude is smaller by the factor (1 - (a/b)^3) (the cube of half width 1 lies between b = 1 and b = 1.24).

    python stress_eshelby_a80.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import wp3_stress_reference as S  # noqa: E402

E = S.E_BIO * (1.0 + S.FLOOR)
NU, ES, A = S.NU, 1.0e-3, 0.3
ref = -2.0 * E * ES / (3.0 * (1.0 - NU))
print(f"closed form, infinite matrix: mean stress {ref:.5e} Pa; free sphere b = 1: {ref * (1 - A**3):.5e}, "
      f"b = 1.24: {ref * (1 - (A / 1.24)**3):.5e}")
for n in (16, 32, 48):
    xc = -1 + (np.arange(n) + 0.5) * (2.0 / n)
    X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
    phi = np.ones((n,) * 3)
    vm0, p0 = S.solve_fields(n, phi, np.full((n,) * 3, ES))
    inside = X**2 + Y**2 + Z**2 < (0.6 * A)**2               # the inner part of the sphere (away from its stepped surface)
    am1 = np.where(X**2 + Y**2 + Z**2 < A**2, ES, 0.0)
    vm, p = S.solve_fields(n, phi, am1)
    print(f"{n}^3: uniform eigenstrain max |vM| {np.abs(vm0).max():.1e}, max |p| {np.abs(p0).max():.1e}; "
          f"sphere: mean stress inside {p[inside].mean():.5e} (spread {p[inside].std():.1e}), vM inside {vm[inside].mean():.2e}")

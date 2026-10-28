"""A-77: the Fortran point model of the 3D runs (ecology_native.f through ctypes, the file the coupled script writes for
--case k26c4) on its own, against the paper's Fig. 11 (Klempt et al. 2026, Table 2, case 4) and the repo's JAX
integrator (JAXFEM/klempt2026_reproduction.py 4sp_case4): 1500 steps of 1e-4, c* = 100, alpha* = 0 up to step 500
and 100 after it, start phi = 0.02 each, psi = 0.999.

    python k26_point_check_a77.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "ansys_usermat" / "coupling"))
import ecology_native_ctypes as FN  # noqa: E402

cst = FN.init(ROOT / "abaqus_composition" / "eco_klempt2026_t2c4_a100.txt")
g = np.zeros((1, 12))
g[0, :4] = 0.02
g[0, 5] = 1.0 - 0.08
g[0, 6:11] = 0.999
one, zero = np.ones(1), np.zeros(1)
for k in range(1, 1501):
    a0 = one if k - 1 > 500 else zero
    a1 = one if k > 500 else zero
    g, ok = FN.step_lin(g, 1.0e-4, 1, one, one, a0, a1)
    assert ok.all(), k
pb = g[0, :4] * g[0, 6:10]
print("Fortran point model, case 4, phibar at step 1500:", np.round(pb, 3))
print("JAX integrator (klempt2026_reproduction.py):      [0.    0.001 0.001 0.929]")
print("paper, Fig. 11 (read off):                        [0.    0.    0.    0.93 ]")

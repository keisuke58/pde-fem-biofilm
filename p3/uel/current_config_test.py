"""current_config_test.py -- the nutrient in the current configuration (JPROPS(3) = 1): amount J c per
reference volume, flux J D C^-1 Grad c. One element expanded uniformly between two states (u = 0.1 X, J =
1.331), no reaction, no flux: the amount is conserved, so c drops from 1 to 1/1.331 = 0.7513 (with the
option) and stays 1 without it. Then a bar stretched along x (F = diag(2, 1, 1)) with c fixed at both ends:
the steady flux per reference area is D J / F11^2 (c1 - c0) / L = D / 2 * (c1 - c0) / L, half the flux
without the option.

    python current_config_test.py
"""
import ctypes as ct
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
lib = ct.CDLL(str(HERE / "libunified_uel.so"))
D_P = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
lib.uel_element.argtypes = [D_P, D_P, D_P, D_P, D_P, ct.c_double, D_P, D_P]
lib.uel_set_finite.argtypes = [ct.c_int, ct.c_int]
lib.uel_set_current.argtypes = [ct.c_int]
F = lambda a: np.asfortranarray(a, dtype=np.float64)
X = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
P = F([1.0, 0.3, 1e-3, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.05, 0.0])   # L0 = 0, D = 1


def solve_c(U, Un, free_c, dt=0.1):
    sv = F(np.zeros((14, 8)))
    for it in range(20):
        R = F(np.zeros(40)); K = F(np.zeros((40, 40)))
        lib.uel_element(F(U), F(Un), F(X.T), sv, P, dt, R, K)
        if free_c.size == 0 or np.abs(R[free_c]).max() < 1e-13:
            return U, R
        U = U.copy(); U[free_c] -= np.linalg.solve(K[np.ix_(free_c, free_c)], R[free_c])
    return U, R


if __name__ == "__main__":
    lib.uel_set_finite(1, 0)
    cdof = np.arange(4, 40, 5)
    for cur in (0, 1):
        lib.uel_set_current(cur)
        Un = np.zeros(40); Un[3::5] = 0.5; Un[4::5] = 1.0
        U = Un.copy()
        for a in range(8):
            U[5 * a:5 * a + 3] = 0.1 * X[a]
        U, _ = solve_c(U, Un, cdof)
        print(f"current configuration {cur}: uniform expansion J = 1.331, c = {U[4]:.4f} (conservation: {1 / 1.331:.4f})")
    for cur in (0, 1):
        lib.uel_set_current(cur)
        U = np.zeros(40); U[3::5] = 0.5
        for a in range(8):
            U[5 * a] = X[a, 0]                      # u_x = X: F = diag(2, 1, 1), J = 2
        U[4::5] = np.where(X[:, 0] < 0.5, 0.0, 1.0)
        _, R = solve_c(U, U.copy(), np.array([], int), dt=1e30)
        flux = -R[cdof][X[:, 0] < 0.5].sum()      # reaction at x = 0: flux through the unit reference face
        print(f"current configuration {cur}: stretched bar, flux per reference area {flux:.4f} "
              f"(expected {0.5 if cur else 1.0:.4f})")

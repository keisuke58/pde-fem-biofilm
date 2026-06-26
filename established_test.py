import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import factorized
import tooth_pde3d as t3

nodes, tets = t3.parse_mesh("p23_biofilm_commensal_static.inp")
Ml, K = t3.assemble(nodes, tets)
outer = np.arange(4 * t3.NV, 5 * t3.NV)
N = len(Ml)
M = sp.diags(Ml)


def solve_established(phi_vec, dt=2e-3, ns=1000):
    """ESTABLISHED biofilm: phase-field phi=1 uniform over the full 0.2mm shell.
    Nutrient c diffuses from the oral (outer) surface and is consumed; alpha
    accumulates dalpha/dt = k_alpha*phi*monod(c)."""
    k_alpha = float(phi_vec @ t3.K_ALPHA)
    A_c = (M + dt * t3.D_C * K).tolil()
    for n in outer:
        A_c.rows[n] = [n]; A_c.data[n] = [1.0]
    sc = factorized(A_c.tocsc())
    phi = np.ones(N)
    c = np.ones(N)
    alpha = np.zeros(N)
    for _ in range(ns):
        monod = c / (t3.K_M + c)
        rhs = Ml * (c - dt * t3.GAMMA_C * phi * monod)
        rhs[outer] = 1.0
        c = sc(rhs)
        np.clip(c, 0, 1, out=c)
        alpha = alpha + dt * (k_alpha * phi * monod)
    return alpha


PHI = {"commensal_hobic": np.array([0.942, 0.012, 0.012, 0.011, 0.011]),
       "dysbiotic_hobic": np.array([0.097, 0.119, 0.474, 0.123, 0.093])}
for cond, phi in PHI.items():
    a = solve_established(phi)
    prof = [round(a[L * t3.NV:(L + 1) * t3.NV].mean(), 4) for L in range(5)]
    print("[%s] ESTABLISHED phi=1: alpha_max=%.4f mean=%.4f per-layer(in->out)=%s"
          % (cond, a.max(), a.mean(), prof))

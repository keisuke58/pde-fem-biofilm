"""Resolution check: fine 1D through-thickness (resolves the 0.045mm nutrient
penetration) vs the coarse 4-element 3D FE. Established biofilm phi=1.
Physical shell 0.2 mm; nutrient Dirichlet at outer (z=0.2), no-flux inner (z=0)."""
import numpy as np
import tooth_pde3d as t3

D_C, GAMMA_C, K_M = t3.D_C, t3.GAMMA_C, t3.K_M
K_ALPHA = t3.K_ALPHA
L = 0.2  # mm shell thickness
PHI = {"commensal_hobic": np.array([0.942, 0.012, 0.012, 0.011, 0.011]),
       "dysbiotic_hobic": np.array([0.097, 0.119, 0.474, 0.123, 0.093])}


def fine_1d(phi_vec, Nz, t_end=2.0):
    k_alpha = float(phi_vec @ K_ALPHA)
    z = np.linspace(0, L, Nz)
    dz = z[1] - z[0]
    # implicit diffusion (Thomas via numpy solve on tridiag), Dirichlet c=1 at z=L
    dt = 2e-3
    nst = int(t_end / dt)
    # build tridiagonal A = I - dt*D_C*Lap (interior), Dirichlet at last node
    import scipy.sparse as sp
    from scipy.sparse.linalg import factorized
    main = np.full(Nz, 1 + 2 * dt * D_C / dz**2)
    off = np.full(Nz - 1, -dt * D_C / dz**2)
    A = sp.diags([off, main, off], [-1, 0, 1]).tolil()
    # no-flux at z=0 (mirror): node0 couples 2x to node1
    A[0, 0] = 1 + 2 * dt * D_C / dz**2; A[0, 1] = -2 * dt * D_C / dz**2
    A[Nz - 1, :] = 0; A[Nz - 1, Nz - 1] = 1.0   # Dirichlet outer
    solve = factorized(A.tocsc())
    c = np.ones(Nz); alpha = np.zeros(Nz)
    for _ in range(nst):
        monod = c / (K_M + c)
        rhs = c - dt * GAMMA_C * monod
        rhs[-1] = 1.0
        c = solve(rhs); np.clip(c, 0, 1, out=c)
        alpha = alpha + dt * k_alpha * monod
    # sample at the 5 tooth layer depths (z/L = 0,0.25,0.5,0.75,1)
    layers = np.interp([0, 0.25, 0.5, 0.75, 1.0], z / L, alpha)
    return layers, alpha.max()


print("penetration sqrt(D_C/GAMMA_C) = %.4f mm  vs element size 0.05 mm\n"
      % np.sqrt(D_C / GAMMA_C))
for cond, phi in PHI.items():
    print("[%s]" % cond)
    for Nz in [5, 11, 21, 51, 101]:
        prof, amax = fine_1d(phi, Nz)
        print("  Nz=%3d  alpha@layers(in->out)=%s  amax=%.3f"
              % (Nz, [round(v, 3) for v in prof], amax))

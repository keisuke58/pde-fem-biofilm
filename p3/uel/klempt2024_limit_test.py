"""klempt2024_limit_test.py -- the unified user element (unified_uel.f90) in the limit of Klempt et al. 2024:
one species, living fraction frozen near 1, nutrient fixed at c = 1, phase fixed at phi = 1, lambda = 0.
The growth reaction then reduces to Eq. 36 of Klempt et al. 2024, alpha_dot = k_alpha phi, with
k_alpha = L0 * liv * (Y ln c - mu_b) = L0 * liv * ln 10 (c = 1, mu_b = ln 0.1).

The element is compared with references written independently of its Fortran source:
  1. free growth (supports 3-2-1): alpha = 1 + k_alpha t, the cube is stretched by alpha, no stress;
  2. fully constrained (all u = 0): hydrostatic Cauchy stress
       sigma = (phi^2 + f) K (alpha^-3 - 1),  K = E / (3 (1 - 2 nu)),
     from W = alpha^3 (phi^2 + f) [mu/2 (J_e^-2/3 tr C_e - 3) + K/2 (J_e - 1)^2], F_e = F / alpha;
  3. axially constrained, laterally free (u_x = 0 at x = 0 and x = 1): the axial stress and the lateral
     stretch from minimising the same W in numpy (an independent implementation of the energy).
Klempt et al. 2024, Table 2: E = 10 Pa, nu = 0.49, k_alpha = 1e-3; a large rate (alpha up to 1.5) as well.

    python klempt2024_limit_test.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import driver as dv  # noqa: E402
from driver import NDN, Model, props_base  # noqa: E402

E, NU, FV = 10.0, 0.49, 1e-3
MU, K = E / (2 * (1 + NU)), E / (3 * (1 - 2 * NU))


def cube():
    x = np.array([0.0, 1.0])
    nodes = np.array([[a, b, c] for a in x for b in x for c in x])
    nid = lambda i, j, k: 4 * i + 2 * j + k
    return nodes, np.array([[nid(0, 0, 0), nid(1, 0, 0), nid(1, 1, 0), nid(0, 1, 0),
                             nid(0, 0, 1), nid(1, 0, 1), nid(1, 1, 1), nid(0, 1, 1)]])


def run(k_alpha, T, steps, support):
    dv.setup(np.zeros(25), b_one=False, cstar=0.0)            # interactions off; the point model is frozen
    nodes, elems = cube()
    L0 = k_alpha / np.log(10.0)
    P = props_base(E=E, nu=NU, f=FV, lam=0.0, beta=0.0, L0=L0, Y=1.0, mu_b=np.log(0.1), s=1e-14, phi_min=2.0)
    m = Model(nodes, elems, P)
    for sv in m.sv:                                          # one species, living fraction ~ 1, frozen
        sv[:, :] = 0.0
        sv[0, :] = 0.999; sv[5, :] = 0.001; sv[6, :] = 0.999; sv[13, :] = 1.0
    liv = 0.999 * 0.999
    U = np.zeros(m.ndof); U[3::NDN] = 1.0; U[4::NDN] = 1.0
    m.fixed[3::NDN] = m.fixed[4::NDN] = True                 # phi = 1 and c = 1 everywhere
    X = nodes
    if support == "free":                                    # 3-2-1
        n0 = np.where((X == [0, 0, 0]).all(1))[0][0]; n1 = np.where((X == [1, 0, 0]).all(1))[0][0]
        n2 = np.where((X == [0, 1, 0]).all(1))[0][0]
        m.fixed[NDN * n0 + 0] = m.fixed[NDN * n0 + 1] = m.fixed[NDN * n0 + 2] = True
        m.fixed[NDN * n1 + 1] = m.fixed[NDN * n1 + 2] = True
        m.fixed[NDN * n2 + 2] = True
    elif support == "all":
        for d in range(3):
            m.fixed[d::NDN] = True
    elif support == "axial":                                 # u_x = 0 at x = 0 and 1; lateral 3-2-1 style
        m.fixed[0::NDN] = True
        n0 = np.where((X == [0, 0, 0]).all(1))[0][0]; n1 = np.where((X == [0, 1, 0]).all(1))[0][0]
        m.fixed[NDN * n0 + 1] = m.fixed[NDN * n0 + 2] = True; m.fixed[NDN * n1 + 2] = True
    dt = T / steps
    for _ in range(steps):
        U, it, r = m.step(U, dt)
    alpha = 1.0 + m.sv[0][12, :].mean()
    x = X + U.reshape(-1, NDN)[:, :3]
    return m, U, x, alpha, liv, P


def sigma_elem(m, U):
    """Cauchy stress of the element, from the residual of the mechanical dofs: the nodal forces of a single
    element with prescribed displacements are the boundary tractions, so the mean stress follows from
    the force on a face. Here: total force in x on the face x = 1 divided by its current area."""
    R, _ = m.assemble(U, U, 1e-30)
    X = m.nodes; x = X + U.reshape(-1, NDN)[:, :3]
    face = np.isclose(X[:, 0], 1.0)
    fx = R[NDN * np.where(face)[0] + 0].sum()
    ys = x[face][:, 1]; zs = x[face][:, 2]
    area = (ys.max() - ys.min()) * (zs.max() - zs.min())
    return fx / area


def W_np(lam_ax, lam_lat, alpha, phi=1.0):
    """independent energy per reference volume: F = diag(lam_ax, lam_lat, lam_lat), F_e = F / alpha"""
    Fe = np.diag([lam_ax, lam_lat, lam_lat]) / alpha
    Je = np.linalg.det(Fe); C = Fe.T @ Fe
    return alpha ** 3 * (phi ** 2 + FV) * (MU / 2 * (Je ** (-2 / 3) * np.trace(C) - 3) + K / 2 * (Je - 1) ** 2)


def axial_reference(alpha):
    h = 1e-6
    dW = lambda l, lt: (W_np(l, lt + h, alpha) - W_np(l, lt - h, alpha)) / (2 * h)
    lt = brentq(lambda s: dW(1.0, s), 0.5, 3.0)              # lateral faces traction free
    P11 = (W_np(1 + h, lt, alpha) - W_np(1 - h, lt, alpha)) / (2 * h)
    J = lt ** 2
    return lt, P11 / J * 1.0                                 # sigma_11 = P_11 F_11 / J, F_11 = 1


if __name__ == "__main__":
    print(f"E = {E}, nu = {NU}, f = {FV}; mu = {MU:.5f}, K = {K:.4f}")
    rows = []
    for k_alpha, T, steps, lab in ((1e-3, 1.1, 11, "Klempt 2024 rate"), (0.5, 1.0, 20, "large growth")):
        m, U, x, a, liv, P = run(k_alpha, T, steps, "free")
        a_ref = 1 + k_alpha * liv * T                        # Eq. 36 with phi = 1, explicit in time: exact
        stretch = x[:, 0].max() - x[:, 0].min()
        s1 = sigma_elem(m, U)
        print(f"[{lab}] free:        alpha-1 {a - 1:.6e} (Eq. 36: {a_ref - 1:.6e}), edge stretch {stretch:.8f} "
              f"(alpha {a:.8f}), sigma_xx {s1:.2e}")
        rows.append(("free", lab, abs(a - a_ref) / (a_ref - 1), abs(stretch - a) / a, abs(s1) / E))
        m, U, x, a, liv, P = run(k_alpha, T, steps, "all")
        s_ref = (1 + FV) * K * (a ** -3 - 1)
        s2 = sigma_elem(m, U)
        print(f"[{lab}] constrained: sigma_xx {s2:.6e} (closed form {s_ref:.6e}), rel. diff {abs(s2 / s_ref - 1):.1e}")
        rows.append(("all", lab, abs(s2 / s_ref - 1)))
        m, U, x, a, liv, P = run(k_alpha, T, steps, "axial")
        lt_ref, s_ref = axial_reference(a)
        lt = x[:, 1].max() - x[:, 1].min(); s3 = sigma_elem(m, U)
        print(f"[{lab}] axial:       lateral stretch {lt:.8f} (numpy {lt_ref:.8f}), sigma_xx {s3:.6e} "
              f"(numpy {s_ref:.6e}), rel. diff {abs(s3 / s_ref - 1):.1e}")
        rows.append(("axial", lab, abs(lt / lt_ref - 1), abs(s3 / s_ref - 1)))

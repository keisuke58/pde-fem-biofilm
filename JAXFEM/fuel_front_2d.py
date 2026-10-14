"""fuel_front_2d.py -- does the biofilm edge move towards the nutrient without a front term?
(UNIFIED_MODEL_CONCEPT.ja.md, section 9, the 2D check).

Reduced model (the composition is lumped: one living share sbar, constant, so that the 2D check is
about the kinematics only):
  growth as a nutrient reaction (Onsager, as unified_two_level_fuel_1d.py):
      xidot = L0 * phi * sbar * (Y ln c - mu_b),   alpha = 1 + xi,
      cdot  = D lap c - Y xidot,
  biofilm phase phi fixed in the reference configuration (no front term, no phase evolution);
  mechanics: plane strain, small strain, linear elastic with the growth eigenstrain (alpha - 1) I,
  stiffness E(phi) = (phi^2 + f) E_bio (Klempt et al. 2024's phi^2 weight, void floor f),
  substrate y = 0 fixed, other edges free.
Domain [0, 1]^2, colony = half disc of radius 0.25 at (0.5, 0), nutrient supplied on the right
edge x = 1 (c = 1), the other edges closed; the liquid starts at the stop concentration (no growth
until the supply arrives), so the nutrient comes only from the right.  If the edge moves towards the nutrient by the
kinematics alone, the colony's right flank moves out more than its left flank.

Assumptions of this check: lumped composition, small strain (alpha - 1 up to ~0.3), nu = 0.3
(compressible, to keep the bilinear elements free of locking), nutrient in the reference
configuration.  Output: assets/fig_fuel_front_2d.pdf.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))


def grid(n):
    x = np.linspace(0, 1, n + 1)
    X, Y = np.meshgrid(x, x, indexing="ij")
    return x, X, Y


def nutrient_growth(phi, n, *, L0=2.0, Y=1.0, c_stop=0.1, D=1.0, sbar=0.5, T=0.2375, c_init=None):
    h = 1.0 / n
    dt = 0.2 * h * h / D
    steps = int(np.ceil(T / dt)); dt = T / steps
    c = np.full_like(phi, c_stop if c_init is None else c_init); c[-1, :] = 1.0
    xi = np.zeros_like(phi); mu_b = Y * np.log(c_stop)
    for _ in range(steps):
        rate = L0 * phi * sbar * np.maximum(Y * np.log(np.maximum(c, 1e-12)) - mu_b, 0.0)
        cp = np.pad(c, 1, mode="edge")                       # closed edges (no flux)
        lap = (cp[2:, 1:-1] + cp[:-2, 1:-1] + cp[1:-1, 2:] + cp[1:-1, :-2] - 4 * c) / h ** 2
        c = c + dt * (D * lap - Y * rate)
        c[-1, :] = 1.0                                       # supply on the right edge x = 1
        c = np.maximum(c, 1e-12)
        xi = xi + dt * rate
    return 1.0 + xi, c


def elasticity(phi, alpha, n, *, f=1e-3, nu=0.3, E=1.0):
    """bilinear plane-strain elements on the (n+1)^2 grid, eigenstrain (alpha - 1) I per element"""
    h = 1.0 / n
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
    Cm = np.array([[lam + 2 * mu, lam, 0], [lam, lam + 2 * mu, 0], [0, 0, mu]])
    gp = np.array([-1, 1]) / np.sqrt(3)
    nodes = lambda i, j: [(i) * (n + 1) + j, (i + 1) * (n + 1) + j, (i + 1) * (n + 1) + j + 1, i * (n + 1) + j + 1]
    rows, cols, vals = [], [], []
    F = np.zeros(2 * (n + 1) ** 2)
    for i in range(n):
        for j in range(n):
            ph = 0.25 * (phi[i, j] + phi[i + 1, j] + phi[i + 1, j + 1] + phi[i, j + 1])
            ge = 0.25 * (alpha[i, j] + alpha[i + 1, j] + alpha[i + 1, j + 1] + alpha[i, j + 1]) - 1.0
            Ce = (ph ** 2 + f) * Cm
            Ke = np.zeros((8, 8)); fe = np.zeros(8)
            for a in gp:
                for b in gp:
                    dN = 0.25 * np.array([[-(1 - b), (1 - b), (1 + b), -(1 + b)], [-(1 - a), -(1 + a), (1 + a), (1 - a)]]) * 2 / h
                    Bm = np.zeros((3, 8))
                    Bm[0, 0::2] = dN[0]; Bm[1, 1::2] = dN[1]; Bm[2, 0::2] = dN[1]; Bm[2, 1::2] = dN[0]
                    w = (h / 2) ** 2
                    Ke += w * Bm.T @ Ce @ Bm
                    fe += w * Bm.T @ Ce @ np.array([ge, ge, 0.0])
            dofs = np.array([[2 * k, 2 * k + 1] for k in nodes(i, j)]).ravel()
            rows += list(np.repeat(dofs, 8)); cols += list(np.tile(dofs, 8)); vals += list(Ke.ravel())
            F[dofs] += fe
    K = sp.csr_matrix((vals, (rows, cols)), shape=(F.size, F.size))
    fixed = np.array([[2 * (i * (n + 1)), 2 * (i * (n + 1)) + 1] for i in range(n + 1)]).ravel()   # y = 0
    free = np.setdiff1d(np.arange(F.size), fixed)
    u = np.zeros(F.size)
    u[free] = spla.spsolve(K[free][:, free].tocsc(), F[free])
    return u.reshape(n + 1, n + 1, 2)


def edge_points(X, Y, U, R=0.25):
    """deformed position of the colony edge at angles 0..180 degrees (interpolated on the grid)"""
    from scipy.interpolate import RegularGridInterpolator
    x = X[:, 0]; y = Y[0, :]
    iu = RegularGridInterpolator((x, y), U[:, :, 0]); iv = RegularGridInterpolator((x, y), U[:, :, 1])
    th = np.linspace(0, np.pi, 91)
    P = np.stack([0.5 + R * np.cos(th), R * np.sin(th)], 1)
    return th, P, P + np.stack([iu(P), iv(P)], 1)


if __name__ == "__main__":
    n = 48
    x, X, Y = grid(n)
    r = np.hypot(X - 0.5, Y)
    phi = 0.5 * (1 - np.tanh((r - 0.25) / 0.02)) * 0.9 + 1e-4
    D = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05     # nutrient-limited: slow supply
    alpha, c = nutrient_growth(phi, n, D=D, L0=float(sys.argv[2]) if len(sys.argv) > 2 else 2.0)
    U = elasticity(phi, alpha, n)
    th, P, Q = edge_points(X, Y, U)
    disp = np.linalg.norm(Q - P, axis=1)
    radial = np.einsum("ij,ij->i", Q - P, (P - [0.5, 0]) / 0.25)
    i_r, i_l, i_t = 0, -1, 45
    print(f"alpha - 1 in the colony: right flank {alpha[int(0.74 * n), 1] - 1:.3f}, top {alpha[n // 2, int(0.24 * n)] - 1:.3f}, "
          f"left flank {alpha[int(0.26 * n), 1] - 1:.3f};  c at the colony centre {c[n // 2, 1]:.3f}")
    for deg in (30, 45, 60):
        print(f"outward edge displacement at {deg} deg (nutrient side) {radial[deg // 2]:.4f}, at {180 - deg} deg {radial[(180 - deg) // 2]:.4f}")
    print(f"top (90 deg) {radial[45]:.4f}")
    import figstyle
    figstyle.apply()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.4))
    m = ax[0].contourf(X, Y, c, 20, cmap="Blues_r"); fig.colorbar(m, ax=ax[0], label="nutrient $c$")
    ax[0].plot(P[:, 0], P[:, 1], "k--", lw=1)
    ax[0].set(title="nutrient (supplied at $x=1$)", xlabel="$x$", ylabel="$y$", aspect="equal")
    m = ax[1].contourf(X, Y, alpha - 1, 20, cmap="Oranges"); fig.colorbar(m, ax=ax[1], label=r"$\alpha-1$")
    ax[1].plot(P[:, 0], P[:, 1], "k--", lw=1, label="edge, reference")
    s = 1.0
    ax[1].plot(P[:, 0] + s * (Q[:, 0] - P[:, 0]), P[:, 1] + s * (Q[:, 1] - P[:, 1]), "k-", lw=1.5, label="edge, grown")
    ax[1].set(title="growth and the grown edge", xlabel="$x$", aspect="equal", xlim=(0.1, 0.9), ylim=(0, 0.5))
    ax[1].legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(ROOT / "assets/fig_fuel_front_2d.pdf")
    print(f"D = {D}")
    print("wrote assets/fig_fuel_front_2d.pdf")

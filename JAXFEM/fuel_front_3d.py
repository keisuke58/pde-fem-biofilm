"""fuel_front_3d.py -- the 3D version of fuel_front_2d.py: does the biofilm edge move towards the
nutrient without a front term?

Same reduced model: growth as a nutrient reaction (Onsager form, lumped composition with a constant
living share), xidot = L0 phi sbar (Y ln c - mu_b), cdot = D lap c - Y xidot; small-strain linear
elasticity with the growth eigenstrain (alpha - 1) I and stiffness (phi^2 + f) E; substrate z = 0
fixed, other faces free.  Domain [0, 1]^3, colony = half ball of radius 0.25 on the substrate at
(0.5, 0.5, 0); nutrient supplied on the face x = 1 only, the liquid starts at the stop
concentration.  Trilinear hexahedra on a uniform grid, so one element matrix scaled by the
element stiffness.  Output: assets/fig_fuel_front_3d.pdf.

Assumptions as in 2D: lumped composition, small strain, nu = 0.3, nutrient in the reference
configuration.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))


def nutrient_growth(phi, h, *, L0=2.0, Y=1.0, c_stop=0.1, D=0.2, sbar=0.5, T=0.2375):
    dt = 0.15 * h * h / D
    steps = int(np.ceil(T / dt)); dt = T / steps
    c = np.full_like(phi, c_stop); c[-1] = 1.0
    xi = np.zeros_like(phi); mu_b = Y * np.log(c_stop)
    for _ in range(steps):
        rate = L0 * phi * sbar * np.maximum(Y * np.log(np.maximum(c, 1e-12)) - mu_b, 0.0)
        cp = np.pad(c, 1, mode="edge")
        lap = (cp[2:, 1:-1, 1:-1] + cp[:-2, 1:-1, 1:-1] + cp[1:-1, 2:, 1:-1] + cp[1:-1, :-2, 1:-1]
               + cp[1:-1, 1:-1, 2:] + cp[1:-1, 1:-1, :-2] - 6 * c) / h ** 2
        c = np.maximum(c + dt * (D * lap - Y * rate), 1e-12); c[-1] = 1.0
        xi = xi + dt * rate
    return 1.0 + xi, c


def element_matrices(h, nu=0.3, E=1.0):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6)); C[:3, :3] = lam; C[np.arange(3), np.arange(3)] += 2 * mu; C[3:, 3:] = mu * np.eye(3)
    corners = np.array([[i, j, k] for k in (0, 1) for j in (0, 1) for i in (0, 1)]) * 2 - 1
    g = np.array([-1, 1]) / np.sqrt(3)
    K = np.zeros((24, 24)); f = np.zeros(24)
    for a in g:
        for b in g:
            for cc in g:
                dN = np.array([[corners[n, 0] * (1 + corners[n, 1] * b) * (1 + corners[n, 2] * cc),
                                corners[n, 1] * (1 + corners[n, 0] * a) * (1 + corners[n, 2] * cc),
                                corners[n, 2] * (1 + corners[n, 0] * a) * (1 + corners[n, 1] * b)] for n in range(8)]).T / 8 * 2 / h
                B = np.zeros((6, 24))
                for n in range(8):
                    B[0, 3 * n] = dN[0, n]; B[1, 3 * n + 1] = dN[1, n]; B[2, 3 * n + 2] = dN[2, n]
                    B[3, 3 * n] = dN[1, n]; B[3, 3 * n + 1] = dN[0, n]
                    B[4, 3 * n + 1] = dN[2, n]; B[4, 3 * n + 2] = dN[1, n]
                    B[5, 3 * n] = dN[2, n]; B[5, 3 * n + 2] = dN[0, n]
                w = (h / 2) ** 3
                K += w * B.T @ C @ B
                f += w * B.T @ C @ np.array([1, 1, 1, 0, 0, 0.0])
    return K, f, corners


def elasticity(phi, alpha, n, f_void=1e-3):
    h = 1.0 / n
    K0, f0, corners = element_matrices(h)
    nn = n + 1
    nid = lambda i, j, k: (i * nn + j) * nn + k
    I, J, Kk = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    I, J, Kk = I.ravel(), J.ravel(), Kk.ravel()
    off = (corners + 1) // 2
    enodes = np.stack([nid(I + o[0], J + o[1], Kk + o[2]) for o in off], 1)
    pe = np.mean([phi[I + o[0], J + o[1], Kk + o[2]] for o in off], 0)
    ge = np.mean([alpha[I + o[0], J + o[1], Kk + o[2]] for o in off], 0) - 1.0
    stiff = pe ** 2 + f_void
    dofs = np.stack([3 * enodes[:, m] + d for m in range(8) for d in range(3)], 1)
    rows = np.repeat(dofs, 24, axis=1).ravel(); cols = np.tile(dofs, (1, 24)).ravel()
    vals = (stiff[:, None] * K0.ravel()[None, :]).ravel()
    ndof = 3 * nn ** 3
    Kg = sp.csr_matrix((vals, (rows, cols)), shape=(ndof, ndof))
    F = np.zeros(ndof); np.add.at(F, dofs.ravel(), ((stiff * ge)[:, None] * f0[None, :]).ravel())
    fixed = np.array([3 * nid(i, j, 0) + d for i in range(nn) for j in range(nn) for d in range(3)])
    free = np.setdiff1d(np.arange(ndof), fixed)
    u = np.zeros(ndof)
    u[free] = spla.spsolve(Kg[free][:, free].tocsc(), F[free])
    return u.reshape(nn, nn, nn, 3)


if __name__ == "__main__":
    from scipy.interpolate import RegularGridInterpolator
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    D = float(sys.argv[2]) if len(sys.argv) > 2 else 0.2
    x = np.linspace(0, 1, n + 1)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt((X - 0.5) ** 2 + (Y - 0.5) ** 2 + Z ** 2)
    phi = 0.5 * (1 - np.tanh((r - 0.25) / 0.03)) * 0.9 + 1e-4
    alpha, c = nutrient_growth(phi, 1.0 / n, D=D)
    U = elasticity(phi, alpha, n)
    itp = [RegularGridInterpolator((x, x, x), U[..., d]) for d in range(3)]
    print(f"n = {n}, D = {D}: alpha - 1 at (0.7, 0.5, 0.05) {alpha[int(0.7 * n), n // 2, 1] - 1:.3f}, "
          f"(0.3, 0.5, 0.05) {alpha[int(0.3 * n), n // 2, 1] - 1:.3f}")
    for el in (30, 45):
        res = {}
        for name, az in (("+x (nutrient)", 0.0), ("-x", np.pi), ("+y (side)", np.pi / 2)):
            e = np.radians(el)
            dvec = np.array([np.cos(e) * np.cos(az), np.cos(e) * np.sin(az), np.sin(e)])
            p = np.array([0.5, 0.5, 0.0]) + 0.25 * dvec
            u = np.array([f(p)[0] for f in itp])
            res[name] = u @ dvec
        print(f"  elevation {el} deg, outward edge displacement: " + ", ".join(f"{k} {v:.4f}" for k, v in res.items()))
    import figstyle
    figstyle.apply()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    j = n // 2
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.4))
    m = ax[0].contourf(X[:, j, :], Z[:, j, :], c[:, j, :], 20, cmap="Blues_r"); fig.colorbar(m, ax=ax[0], label="nutrient $c$")
    th = np.linspace(0, np.pi, 91); P = np.stack([0.5 + 0.25 * np.cos(th), 0.5 + 0 * th, 0.25 * np.sin(th)], 1)
    ax[0].plot(P[:, 0], P[:, 2], "k--", lw=1)
    ax[0].set(title="nutrient, section $y=0.5$ (supplied at $x=1$)", xlabel="$x$", ylabel="$z$", aspect="equal")
    m = ax[1].contourf(X[:, j, :], Z[:, j, :], alpha[:, j, :] - 1, 20, cmap="Oranges"); fig.colorbar(m, ax=ax[1], label=r"$\alpha-1$")
    Q = P + np.stack([np.array([f(p)[0] for p in P]) for f in itp], 1)
    ax[1].plot(P[:, 0], P[:, 2], "k--", lw=1, label="edge, reference")
    ax[1].plot(Q[:, 0], Q[:, 2], "k-", lw=1.5, label="edge, grown")
    ax[1].set(title="growth and the grown edge", xlabel="$x$", aspect="equal", xlim=(0.1, 0.9), ylim=(0, 0.45))
    ax[1].legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(ROOT / "assets/fig_fuel_front_3d.pdf")
    print("wrote assets/fig_fuel_front_3d.pdf")

"""bilayer_test.py -- volumetric locking of the user element near incompressibility, and the F-bar
remedy: a strip L x b x h (1 x 0.05 x 0.1), clamped at x = 0, top half grown by xi (phi = 0.999, no
composition), bottom half not. Beam theory for two equal layers with an eigenstrain mismatch xi gives
the curvature 1.5 xi / h and the tip deflection 0.75 xi L^2 / h (Timoshenko 1925, equal thickness and
stiffness), for nu = 0 exactly and approximately otherwise (the strip is narrow).

    python bilayer_test.py
"""
import ctypes as ct
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
lib = ct.CDLL(str(HERE / "libunified_uel.so"))
D_P = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
lib.uel_element.argtypes = [D_P, D_P, D_P, D_P, D_P, ct.c_double, D_P, D_P]
lib.uel_set_finite.argtypes = [ct.c_int, ct.c_int]
NDN, NDOF, NGP, NSV = 5, 40, 8, 14
F = lambda a: np.asfortranarray(a, dtype=np.float64)


def strip(nx, nz, L=1.0, b=0.05, h=0.1):
    xs, ys, zs = np.linspace(0, L, nx + 1), np.array([0.0, b]), np.linspace(0, h, nz + 1)
    nodes = np.array([[x, y, z] for x in xs for y in ys for z in zs])
    nid = lambda i, j, k: (i * 2 + j) * (nz + 1) + k
    elems = np.array([[nid(i, 0, k), nid(i + 1, 0, k), nid(i + 1, 1, k), nid(i, 1, k),
                       nid(i, 0, k + 1), nid(i + 1, 0, k + 1), nid(i + 1, 1, k + 1), nid(i, 1, k + 1)]
                      for i in range(nx) for k in range(nz)])
    top = np.array([k >= nz // 2 for i in range(nx) for k in range(nz)])
    return nodes, elems, top


def solve(nu, xi, finite=1, fbar=0, nx=40, nz=4):
    lib.uel_set_finite(finite, fbar)
    nodes, elems, top = strip(nx, nz)
    props = F([1.0, nu, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.05, 0.0])
    U = np.zeros(NDN * len(nodes)); U[3::NDN] = 0.999; U[4::NDN] = 1.0          # phi below 1: the barrier term is 0/0 at 1
    fixed = np.zeros(len(U), bool); fixed[3::NDN] = fixed[4::NDN] = True
    clamp = np.where(np.isclose(nodes[:, 0], 0.0))[0]
    for d in range(3):
        fixed[NDN * clamp + d] = True
    free = np.where(~fixed)[0]
    EDS = [np.array([NDN * n + d for n in el for d in range(NDN)]) for el in elems]
    its = 0
    # growth raised in steps, as it accumulates over increments in a run: started from u = 0 with the
    # full xi, a nearly incompressible grown layer is in hydrostatic compression and the geometric
    # stiffness makes the tangent indefinite
    for lf in np.linspace(0.1, 1.0, 10):
        SV = []
        for t in top:
            sv = np.zeros((NSV, NGP)); sv[12] = lf * xi if t else 0.0; SV.append(F(sv))   # xi_n, not started
        for it in range(30):
            R = np.zeros(len(U)); rows, cols, vals = [], [], []
            for e, el in enumerate(elems):
                ed = EDS[e]; Re = F(np.zeros(NDOF)); Ke = F(np.zeros((NDOF, NDOF)))
                lib.uel_element(F(U[ed]), F(U[ed]), F(nodes[el].T), SV[e], props, 1.0, Re, Ke)
                R[ed] += Re; rows.append(np.repeat(ed, NDOF)); cols.append(np.tile(ed, NDOF)); vals.append(Ke.ravel())
            r = np.abs(R[free]).max()
            if r < 1e-12:
                break
            K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(U), len(U)))
            U[free] += spla.spsolve(K[free][:, free].tocsc(), -R[free])
        else:
            raise RuntimeError(f"no convergence at load factor {lf:.1f}")
        its = max(its, it)
    it = its
    tip = np.isclose(nodes[:, 0], 1.0)
    return U[NDN * np.where(tip)[0] + 2].mean(), it


if __name__ == "__main__":
    xi = 1e-3
    ref = 0.75 * xi * 1.0 ** 2 / 0.1
    print(f"beam theory tip deflection {ref:.4e} (curvature 1.5 xi/h)")
    for nu in (0.0, 0.3, 0.49):
        out = []
        for lab, fin, fb in (("small strain", 0, 0), ("finite", 1, 0), ("finite, mean dilatation", 1, 1)):
            w, it = solve(nu, xi, fin, fb)
            out.append(f"{lab} {abs(w) / ref:.3f} ({it} it)")
        print(f"nu = {nu}: tip deflection / beam theory: " + ", ".join(out), flush=True)

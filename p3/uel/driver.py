"""driver.py -- tests the Abaqus user element unified_uel.f90 without Abaqus: a small Newton solver
that assembles the element (through the C entry points of libunified_uel.so) on a column of
hexahedra, as Abaqus would.

Build:  gfortran -O2 -fPIC -c ../../ansys_usermat/coupling/ecology_native.f -o ecology_native.o
        gfortran -O2 -fPIC -ffree-line-length-none -c unified_uel.f90 -o unified_uel.o
        gfortran -shared -o libunified_uel.so unified_uel.o ecology_native.o
Run:    python driver.py

Checks (P3_STEP3_FEM_PLAN.ja.md, verification order):
  1. one element, uniform: the composition at every Gauss point equals the point model solved
     to convergence in Python (unified_1d_monolithic.point_model_converged);
  2. a column of elements along x, seed at x <= 0.25, nutrient supplied at x = 1: compared with
     the 1D two-level model with the growth reaction (unified_two_level_fuel_1d.py, finite
     differences, beta_s = 0.02 there, 0 here), and refined to show convergence.
"""
import ctypes as ct
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "JAXFEM", ROOT / "ansys_usermat/coupling"):
    sys.path.insert(0, str(p))
lib = ct.CDLL(str(HERE / "libunified_uel.so"))
D_P = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
lib.uel_set_point_model.argtypes = [D_P, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_int, ct.POINTER(ct.c_double)]
lib.uel_set_phi_init.argtypes = [D_P]
lib.uel_set_alstar.argtypes = [ct.c_double]
lib.uel_element.argtypes = [D_P, D_P, D_P, D_P, D_P, ct.c_double, D_P, D_P]
lib.uel_advance.argtypes = [D_P, D_P, D_P, D_P, D_P, ct.c_double, ct.POINTER(ct.c_int)]
lib.uel_store.argtypes = [D_P, D_P, D_P, D_P, ct.c_double]
NN, NDN, NDOF, NGP, NSV = 8, 5, 40, 8, 14
F = lambda a: np.asfortranarray(a, dtype=np.float64)


def setup(theta, b_one=True, alstar=0.0, cstar=25.0):
    th = np.array(theta, float)
    if b_one:
        th[[3, 4, 8, 9, 15]] = 1.0          # antibiotic susceptibility b_i = 1 (theta has none)
    lib.uel_set_point_model(F(th), ct.byref(ct.c_double(cstar)), ct.byref(ct.c_double(alstar)), 60, ct.byref(ct.c_double(1e-26)))
    return th


def column(n, L=1.0, h=None):
    """n hexahedra along x on [0, L], cross-section h x h; nodes (n+1) x 2 x 2"""
    h = L / n if h is None else h
    xs = np.linspace(0, L, n + 1)
    nodes = np.array([[x, y, z] for x in xs for y in (0, h) for z in (0, h)])
    nid = lambda i, j, k: 4 * i + 2 * j + k
    elems = [[nid(e, 0, 0), nid(e + 1, 0, 0), nid(e + 1, 1, 0), nid(e, 1, 0),
              nid(e, 0, 1), nid(e + 1, 0, 1), nid(e + 1, 1, 1), nid(e, 1, 1)] for e in range(n)]
    return nodes, np.array(elems), xs


class Model:
    def __init__(self, nodes, elems, props):
        self.nodes, self.elems, self.props = nodes, elems, F(props)
        self.nn = len(nodes); self.ndof = NDN * self.nn
        self.sv = [F(np.zeros((NSV, NGP))) for _ in elems]
        self.fixed = np.zeros(self.ndof, bool)

    def edofs(self, e):
        return np.array([NDN * n + d for n in self.elems[e] for d in range(NDN)])

    def assemble(self, U, Un, dt):
        R = np.zeros(self.ndof); K = np.zeros((self.ndof, self.ndof))
        for e in range(len(self.elems)):
            ed = self.edofs(e)
            Re = F(np.zeros(NDOF)); Ke = F(np.zeros((NDOF, NDOF)))
            lib.uel_element(F(U[ed]), F(Un[ed]), F(self.nodes[self.elems[e]].T), self.sv[e], self.props, dt, Re, Ke)
            R[ed] += Re; K[np.ix_(ed, ed)] += Ke
        return R, K

    def step(self, U, dt, tol=1e-10, maxit=30):
        Un = U.copy()
        for e in range(len(self.elems)):
            ed = self.edofs(e); ok = ct.c_int(0)
            lib.uel_advance(F(U[ed]), F(Un[ed]), F(self.nodes[self.elems[e]].T), self.sv[e], self.props, dt, ct.byref(ok))
            if not ok.value:
                raise RuntimeError("point model failed")
        free = ~self.fixed
        it, r = 0, 0.0
        for it in range(maxit if free.any() else 0):
            R, K = self.assemble(U, Un, dt)
            r = np.abs(R[free]).max()
            if r < tol:
                break
            d = np.linalg.solve(K[np.ix_(free, free)], -R[free])
            U = U.copy(); U[free] += d
        for e in range(len(self.elems)):
            ed = self.edofs(e)
            lib.uel_store(F(U[ed]), F(self.nodes[self.elems[e]].T), self.sv[e], self.props, dt)
        return U, it, r


def props_base(**kw):
    p = dict(E=1.0, nu=0.3, f=1e-3, beta=0.02, eta_phi=1.0, kphi=1e-5, lam=1.0, psiref=0.0, L0=0.7, Y=1.0,
             mu_b=np.log(0.1), D=1.0, s=0.2375, n_sub=1, phi_min=0.05, c_ref=0.0)
    p.update(kw)
    return np.array(list(p.values()), float)


if __name__ == "__main__":
    import unified_1d as u1
    import unified_1d_monolithic as m1
    g0, th = u1.ch_case()
    g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum(); g[5] = 0.1
    setup(th)
    lib.uel_set_phi_init(F(g[0:5]))
    # ---- check 1: one element, uniform phi = 0.8, c = 1 everywhere fixed
    nodes, elems, xs = column(1)
    P = props_base(lam=0.0, beta=0.0, n_sub=1)
    mod = Model(nodes, elems, P)
    U = np.zeros(mod.ndof); U[3::NDN] = 0.8; U[4::NDN] = 1.0
    mod.fixed[:] = True                                   # everything prescribed: only the composition moves
    steps, T = 200, 0.02 / 0.2375                         # point-model time 0.02 in 200 steps of 1e-4
    for k in range(steps):
        U, it, r = mod.step(U, T / steps)
    gref = m1.point_model_converged(np.r_[g[0:5], 0.1, np.clip(g0[6:11], 1e-10, 1 - 1e-10), 0.0], th, 1e-4, steps)
    print(f"check 1: max |Gauss-point composition - converged point model| = {np.abs(mod.sv[0][0:11, 0] - gref[0:11]).max():.2e}")
    # ---- check 2: column along x, seed x <= 0.25, nutrient supplied at x = 1; mechanics held fixed
    import time
    res = {}
    for n in (20, 40):
        nodes, elems, xs = column(n)
        # time convention of the 1D Python models: the field and the point model share one time,
        # T = 0.2375 (the calibration's Day 1 to Day 21), so s = 1
        P = props_base(s=1.0)
        mod = Model(nodes, elems, P)
        U = np.zeros(mod.ndof)
        xn = nodes[:, 0]
        U[3::NDN] = np.where(xn <= 0.25 + 1e-12, 0.9, 1e-4)
        U[4::NDN] = 1.0
        mod.fixed[0::NDN] = mod.fixed[1::NDN] = mod.fixed[2::NDN] = True     # no mechanics in this check
        mod.fixed[4::NDN] = np.isclose(xn, 1.0)                              # supply c = 1 at x = 1
        steps = 60; t0 = time.time(); worst = 0.0
        for k in range(steps):
            U, it, r = mod.step(U, 0.2375 / steps)
            worst = max(worst, r)
        phi = U[3::NDN].reshape(n + 1, 4)[:, 0]; c = U[4::NDN].reshape(n + 1, 4)[:, 0]
        xi0 = mod.sv[0][12].mean(); liv0 = (mod.sv[0][0:5] * mod.sv[0][6:11]).sum(0).mean()
        res[n] = (phi[0], c[0], xi0, liv0)
        print(f"check 2, {n} elements: alpha-1 near x = 0 {xi0:.4f}, phi(0) {phi[0]:.3f}, c(0) {c[0]:.3f}, "
              f"living near x = 0 {liv0:.3f}; Newton residual <= {worst:.1e}; {time.time() - t0:.0f} s", flush=True)
    print("reference (unified_two_level_fuel_1d.py, finite differences, 41 nodes): alpha-1(0) 0.137, phi(0) 0.77, "
          "c(0) 0.938, living(0) 0.53")

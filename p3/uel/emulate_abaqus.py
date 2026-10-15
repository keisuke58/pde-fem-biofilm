"""emulate_abaqus.py -- runs uel_column.inp's two steps by calling the Abaqus entry point UEL itself
(Abaqus's argument list, SVARS handed in at the start of every increment, the point model read from
the file named by BIOFILM_ECO_CASE), as a check of everything an Abaqus run uses except Abaqus.

    python make_column_inp.py 20 /tmp/x && BIOFILM_ECO_CASE=/tmp/x/uel_eco.txt python emulate_abaqus.py 20

Expected: the numbers of driver.py check 2 (20 elements: phi(0) 0.776, c(0) 0.936, alpha-1 0.1371).
"""
import ctypes as ct
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
assert os.environ.get("BIOFILM_ECO_CASE"), "set BIOFILM_ECO_CASE (make_column_inp.py writes uel_eco.txt)"
sys.path.insert(0, str(HERE))
from driver import NDN, NDOF, NGP, NSV, column, props_base  # noqa: E402

lib = ct.CDLL(str(HERE / "libunified_uel.so"))
F = lambda a: np.asfortranarray(a, dtype=np.float64)
I = lambda v: ct.byref(ct.c_int(v))
D = lambda v: ct.byref(ct.c_double(v))


def uel(U, DU, coords, svars, props, dt, kinc):
    rhs = F(np.zeros((NDOF, 1))); amat = F(np.zeros((NDOF, NDOF))); sv = F(svars.copy())
    z = F(np.zeros(NDOF)); energy = F(np.zeros(8)); time = F(np.zeros(2)); params = F(np.zeros(3))
    zz = F(np.zeros((1, 1))); predef = F(np.zeros((2, 1, 8))); lflags = np.asfortranarray([1, 1, 1, 0, 0], np.int32)
    jd = np.zeros((1, 1), np.int32, order="F"); jp = np.zeros(1, np.int32)
    pnewdt = ct.c_double(1.0)
    p = lambda a: a.ctypes.data_as(ct.c_void_p)
    lib.uel_(p(rhs), p(amat), p(sv), p(energy), I(NDOF), I(1), I(NSV * NGP), p(props), I(16), p(coords), I(3), I(8),
             p(U), p(DU), p(z), p(z), I(1), p(time), D(dt), I(1), I(kinc), I(1), p(params), I(0), p(jd), p(zz), p(predef),
             I(1), p(lflags), I(NDOF), p(zz), I(1), ct.byref(pnewdt), p(jp), I(0), D(0.0))
    return -rhs[:, 0], amat, sv, pnewdt.value


def increment(U0, U, fixed, nodes, elems, SV, props, dt, kinc, tol=1e-10):
    """one Abaqus increment: Newton on the free dofs, every element called with the start-of-increment SVARS"""
    for it in range(30):
        R = np.zeros(len(U)); K = np.zeros((len(U), len(U))); SVn = []
        for e, el in enumerate(elems):
            ed = np.array([NDN * n + d for n in el for d in range(NDN)])
            Re, Ke, sv, pn = uel(F(U[ed]), F((U[ed] - U0[ed])[:, None]), F(nodes[el].T), SV[e], props, dt, kinc)
            assert pn == 1.0, "UEL asked for a smaller increment"
            R[ed] += Re; K[np.ix_(ed, ed)] += Ke; SVn.append(sv)
        free = ~fixed
        r = np.abs(R[free]).max() if free.any() else 0.0
        if r < tol:
            return U, SVn, it
        U = U.copy(); U[free] += np.linalg.solve(K[np.ix_(free, free)], -R[free])
    raise RuntimeError("no convergence")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    nodes, elems, xs = column(n)
    props = F(props_base(s=1.0))
    xn = nodes[:, 0]
    SV = [F(np.zeros(NSV * NGP)) for _ in elems]
    # step INIT: everything prescribed (phi seed, c = 1), length 1e-8
    U0 = np.zeros(NDN * len(nodes)); U = U0.copy()
    U[3::NDN] = np.where(xn <= 0.25 + 1e-12, 0.9, 1e-4); U[4::NDN] = 1.0
    U, SV, _ = increment(U0, U, np.ones(len(U), bool), nodes, elems, SV, props, 1e-8, 1)
    # step GROW: 60 increments, c = 1 at x = 1, mechanics fixed
    fixed = np.zeros(len(U), bool); fixed[0::NDN] = fixed[1::NDN] = fixed[2::NDN] = True
    fixed[4::NDN] = np.isclose(xn, 1.0)
    for k in range(60):
        U, SV, it = increment(U.copy(), U, fixed, nodes, elems, SV, props, 0.2375 / 60, k + 1)
    phi = U[3::NDN].reshape(n + 1, 4)[:, 0]; c = U[4::NDN].reshape(n + 1, 4)[:, 0]
    sv0 = SV[0].reshape(NGP, NSV).T
    print(f"Abaqus emulation, {n} elements: alpha-1 near x = 0 {sv0[12].mean():.4f}, phi(0) {phi[0]:.3f}, "
          f"c(0) {c[0]:.3f}, living near x = 0 {(sv0[0:5] * sv0[6:11]).sum(0).mean():.3f}")

"""emulate_abaqus.py -- runs uel_column.inp's two steps by calling the Abaqus entry point UEL itself
(Abaqus's argument list, SVARS handed in at the start of every increment, the point model read from
the file named by BIOFILM_ECO_CASE), as a check of everything an Abaqus run uses except Abaqus.

    python make_column_inp.py 20 /tmp/x && BIOFILM_ECO_CASE=/tmp/x/uel_eco.txt python emulate_abaqus.py 20
    python make_seed_inp.py 8 /tmp/y && BIOFILM_ECO_CASE=/tmp/y/uel_eco.txt python emulate_abaqus.py --seed 8 [--antibiotic] [--small] [--nu 0.49] [--md]

Expected for the column: the numbers of driver.py check 2 (20 elements: phi(0) 0.776, c(0) 0.936,
alpha-1 0.1371). The seed (uel_seed.inp) is compared with JAXFEM/fuel_front_3d_coupled.py, lam = 0.
"""
import ctypes as ct
import os
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
assert os.environ.get("BIOFILM_ECO_CASE"), "set BIOFILM_ECO_CASE (make_column_inp.py writes uel_eco.txt)"
sys.path.insert(0, str(HERE))
from driver import NDN, NDOF, NGP, NSV, column, props_base  # noqa: E402

lib = ct.CDLL(str(HERE / "libunified_uel.so"))
NSVT = NSV * NGP + 2 * NGP          # + von Mises and pressure per Gauss point (SVARS 113:128)
F = lambda a: np.asfortranarray(a, dtype=np.float64)
I = lambda v: ct.byref(ct.c_int(v))
D = lambda v: ct.byref(ct.c_double(v))


def uel(U, DU, coords, svars, props, dt, kinc, ttot=0.0, finite=True, md=False, jelem=1):
    rhs = F(np.zeros((NDOF, 1))); amat = F(np.zeros((NDOF, NDOF))); sv = F(svars.copy())
    z = F(np.zeros(NDOF)); energy = F(np.zeros(8)); time = F(np.zeros(2)); params = F(np.zeros(3))
    time[:] = ttot
    zz = F(np.zeros((1, 1))); predef = F(np.zeros((2, 1, 8))); lflags = np.asfortranarray([1, 1, 1, 0, 0], np.int32)
    jd = np.zeros((1, 1), np.int32, order="F"); jp = np.array([1 if finite else 0, 1 if md else 0], np.int32)
    pnewdt = ct.c_double(1.0)
    p = lambda a: a.ctypes.data_as(ct.c_void_p)
    lib.uel_(p(rhs), p(amat), p(sv), p(energy), I(NDOF), I(1), I(len(svars)), p(props), I(len(props)), p(coords), I(3), I(8),
             p(U), p(DU), p(z), p(z), I(1), p(time), D(dt), I(1), I(kinc), I(jelem), p(params), I(0), p(jd), p(zz), p(predef),
             I(1), p(lflags), I(NDOF), p(zz), I(1), ct.byref(pnewdt), p(jp), I(2), D(0.0))
    return -rhs[:, 0], amat, sv, pnewdt.value


def increment(U0, U, fixed, nodes, elems, SV, props, dt, kinc, ttot=0.0, tol=1e-10, finite=True, md=False):
    """one Abaqus increment: Newton on the free dofs, every element called with the start-of-increment SVARS"""
    EDS = [np.array([NDN * n + d for n in el for d in range(NDN)]) for el in elems]
    for it in range(30):
        R = np.zeros(len(U)); rows, cols, vals = [], [], []; SVn = []
        for e, el in enumerate(elems):
            ed = EDS[e]
            Re, Ke, sv, pn = uel(F(U[ed]), F((U[ed] - U0[ed])[:, None]), F(nodes[el].T), SV[e], props, dt, kinc, ttot, finite, md, e + 1)
            assert pn == 1.0, "UEL asked for a smaller increment"
            R[ed] += Re; rows.append(np.repeat(ed, NDOF)); cols.append(np.tile(ed, NDOF)); vals.append(Ke.ravel())
            SVn.append(sv)
        free = np.where(~fixed)[0]
        r = np.abs(R[free]).max() if free.size else 0.0
        if r < tol:
            return U, SVn, it
        K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(U), len(U)))
        U = U.copy(); U[free] += spla.spsolve(K[free][:, free].tocsc(), -R[free])
    raise RuntimeError("no convergence")


def seed(n, antibiotic=False, finite=True, nu=0.3, md=False):
    import time as _t
    from make_seed_inp import cube, seed_phi
    sys.path.insert(0, str(HERE.parents[1] / "JAXFEM")); sys.path.insert(0, str(HERE.parents[1]))
    from fuel_front_3d_species import edge_disp
    nodes, elems = cube(n)
    T, nin, t0 = 0.2375, 60, 1e-8
    P = props_base(s=1.0, D=0.2, L0=2.0, lam=0.0, nu=nu)     # as make_seed_inp.deck
    if antibiotic:
        P = np.r_[P, t0 + 0.5 * T, 0.65]
    props = F(P)
    SV = [F(np.zeros(NSVT)) for _ in elems]
    U0 = np.zeros(NDN * len(nodes)); U = U0.copy(); U[3::NDN] = seed_phi(nodes); U[4::NDN] = 0.1
    U, SV, _ = increment(U0, U, np.ones(len(U), bool), nodes, elems, SV, props, t0, 1, finite=finite, md=md)
    fixed = np.zeros(len(U), bool)
    base = np.isclose(nodes[:, 2], 0.0)
    for d in range(3):
        fixed[NDN * np.where(base)[0] + d] = True
    fixed[NDN * np.where(np.isclose(nodes[:, 0], 1.0))[0] + 4] = True
    U[NDN * np.where(np.isclose(nodes[:, 0], 1.0))[0] + 4] = 1.0
    tt = time_start = _t.time()
    for k in range(nin):
        U, SV, it = increment(U.copy(), U, fixed, nodes, elems, SV, props, T / nin, k + 1, ttot=t0 + k * T / nin, finite=finite, md=md)
        if k % 10 == 9:
            print(f"  increment {k + 1}: {it} Newton iterations, {_t.time() - tt:.0f} s", flush=True)
    h = 1.0 / n; nn = n + 1
    phi = U[3::NDN].reshape(nn, nn, nn); c = U[4::NDN].reshape(nn, nn, nn)
    u = U.reshape(-1, NDN)[:, 0:3].reshape(nn, nn, nn, 3)
    # growth per element: mean over Gauss points, element index (i, j, k) as in cube()
    xi = np.array([sv[:NSV * NGP].reshape(NGP, NSV)[:, 12].mean() for sv in SV]).reshape(n, n, n)
    st = np.array([sv[NSV * NGP:].reshape(NGP, 2).mean(0) for sv in SV])          # element mean of vM, p
    pe = np.array([phi.ravel()[el].mean() for el in elems])
    inb = pe > 0.5
    if not inb.any():                       # coarse meshes: no element mean above 0.5
        inb = pe >= pe.max()
    j = n // 2
    print(f"seed {n}^3{', antibiotic' if antibiotic else ''}, {'finite' if finite else 'small'} strain, nu {nu}{', mean dilatation' if md else ''}: biofilm volume (reference) {(phi[:-1, :-1, :-1].sum() * h ** 3):.4f} "
          f"(node sum), phi at centre {phi[j, j, 1]:.3f}, alpha-1 nutrient side (element at x = 0.7, z near 0) "
          f"{xi[int(0.7 * n), j - 1:j + 1, 0].mean():.3f}, far side (x = 0.3) {xi[int(0.3 * n), j - 1:j + 1, 0].mean():.3f}; "
          f"edge 30 deg nutrient side / far side {edge_disp(np.linspace(0, 1, nn), u, 30, 0.0):.4f} / "
          f"{edge_disp(np.linspace(0, 1, nn), u, 30, np.pi):.4f}; biofilm (phi > 0.5) von Mises mean / max "
          f"{st[inb, 0].mean():.4e} / {st[inb, 0].max():.4e}, pressure mean {st[inb, 1].mean():.4e}; "
          f"{_t.time() - time_start:.0f} s")


if __name__ == "__main__":
    if "--seed" in sys.argv:
        nu = float(sys.argv[sys.argv.index("--nu") + 1]) if "--nu" in sys.argv else 0.3
        seed(int(sys.argv[sys.argv.index("--seed") + 1]), "--antibiotic" in sys.argv, "--small" not in sys.argv,
             nu, "--md" in sys.argv)
        sys.exit()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    nodes, elems, xs = column(n)
    props = F(props_base(s=1.0))
    xn = nodes[:, 0]
    SV = [F(np.zeros(NSVT)) for _ in elems]
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
    sv0 = SV[0][:NSV * NGP].reshape(NGP, NSV).T
    print(f"Abaqus emulation, {n} elements: alpha-1 near x = 0 {sv0[12].mean():.4f}, phi(0) {phi[0]:.3f}, "
          f"c(0) {c[0]:.3f}, living near x = 0 {(sv0[0:5] * sv0[6:11]).sum(0).mean():.3f}")

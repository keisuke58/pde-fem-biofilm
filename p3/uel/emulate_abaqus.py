"""emulate_abaqus.py -- runs uel_column.inp's two steps by calling the Abaqus entry point UEL itself
(Abaqus's argument list, SVARS handed in at the start of every increment, the point model read from
the file named by BIOFILM_ECO_CASE), as a check of everything an Abaqus run uses except Abaqus.

    python make_column_inp.py 20 /tmp/x && BIOFILM_ECO_CASE=/tmp/x/uel_eco.txt python emulate_abaqus.py 20
    python make_seed_inp.py 8 /tmp/y && BIOFILM_ECO_CASE=/tmp/y/uel_eco.txt python emulate_abaqus.py --seed 8 [--antibiotic] [--small] [--nu 0.49] [--md] [--L0 2] [--lam 1 --psiref 0.32] [--cur]

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


def uel(U, DU, coords, svars, props, dt, kinc, ttot=0.0, finite=True, md=False, jelem=1, cur=False):
    rhs = F(np.zeros((NDOF, 1))); amat = F(np.zeros((NDOF, NDOF))); sv = F(svars.copy())
    z = F(np.zeros(NDOF)); energy = F(np.zeros(8)); time = F(np.zeros(2)); params = F(np.zeros(3))
    time[:] = ttot
    zz = F(np.zeros((1, 1))); predef = F(np.zeros((2, 1, 8))); lflags = np.asfortranarray([1, 1, 1, 0, 0], np.int32)
    jd = np.zeros((1, 1), np.int32, order="F"); jp = np.array([1 if finite else 0, 1 if md else 0, 1 if cur else 0], np.int32)
    pnewdt = ct.c_double(1.0)
    p = lambda a: a.ctypes.data_as(ct.c_void_p)
    lib.uel_(p(rhs), p(amat), p(sv), p(energy), I(NDOF), I(1), I(len(svars)), p(props), I(len(props)), p(coords), I(3), I(8),
             p(U), p(DU), p(z), p(z), I(1), p(time), D(dt), I(1), I(kinc), I(jelem), p(params), I(0), p(jd), p(zz), p(predef),
             I(1), p(lflags), I(NDOF), p(zz), I(1), ct.byref(pnewdt), p(jp), I(3), D(0.0))
    return -rhs[:, 0], amat, sv, pnewdt.value


def increment(U0, U, fixed, nodes, elems, SV, props, dt, kinc, ttot=0.0, tol=1e-10, finite=True, md=False, cur=False,
              maxit=30, line_search=True, hist=None):
    """one Abaqus increment: Newton on the free dofs, every element called with the start-of-increment SVARS.
    line_search: when a full Newton step raises the residual by more than half, or the element returns NaN
    (an inverted element), the step is halved, up to 6 times (Abaqus: *CONTROLS, PARAMETERS=LINE SEARCH). The
    converged state is the same; without it the Newton of the implant sector with a seed of three element
    layers (6x18x24, 16 Oct 2026) diverged in the first growth increment at the stiff-soft corner of the top face.
    hist: a list; every iteration appends (iteration, largest residual, step factor, residual per dof kind u/phi/c).
    Raises RuntimeError when the residual is not finite or maxit is reached (NoConvergence); U is not changed."""
    EDS = [np.array([NDN * n + d for n in el for d in range(NDN)]) for el in elems]
    free = np.where(~fixed)[0]
    kind = free % NDN
    Up, dU, s, rp, cuts = None, None, 1.0, None, 0
    for it in range(maxit):
        R = np.zeros(len(U)); rows, cols, vals = [], [], []; SVn = []
        for e, el in enumerate(elems):
            ed = EDS[e]
            Re, Ke, sv, pn = uel(F(U[ed]), F((U[ed] - U0[ed])[:, None]), F(nodes[el].T), SV[e], props, dt, kinc, ttot, finite, md, e + 1, cur)
            if pn < 1.0:
                raise NoConvergence(f"increment {kinc}: element {e + 1} asked for a smaller increment (PNEWDT {pn:g})")
            R[ed] += Re; rows.append(np.repeat(ed, NDOF)); cols.append(np.tile(ed, NDOF)); vals.append(Ke.ravel())
            SVn.append(sv)
        r = np.abs(R[free]).max() if free.size else 0.0
        if hist is not None:
            Rf = np.abs(R[free])
            hist.append((it, float(r), s, tuple(float(Rf[m].max()) if m.any() else 0.0 for m in (kind < 3, kind == 3, kind == 4))))
        if line_search and dU is not None and cuts < 6 and (not np.isfinite(r) or r > 1.5 * rp):
            s *= 0.5; cuts += 1
            print(f"  increment {kinc}, iteration {it}: residual {r:.2e} after {rp:.2e}, step halved to {s:g}", file=sys.stderr, flush=True)
            U = Up.copy(); U[free] += s * dU
            continue
        if not np.isfinite(r):
            raise NoConvergence(f"increment {kinc}: residual not finite")
        if r < tol:
            return U, SVn, it
        K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(U), len(U)))
        dU = spla.spsolve(K[free][:, free].tocsc(), -R[free])
        Up, rp, s, cuts = U, r, 1.0, 0
        U = U.copy(); U[free] += dU
    raise NoConvergence(f"increment {kinc}: no convergence in {maxit} iterations")


class NoConvergence(RuntimeError):
    """Newton failed in an increment (not finite, or maxit reached); the caller may cut the increment back"""


def march(U, fixed, nodes, elems, SV, props, T, ninc, t0=0.0, kinc0=1, cutback=0.25, growth=1.5, max_cutbacks=5,
          easy=4, log=None, **kw):
    """a step of length T with Abaqus's automatic time incrementation (the GROW step of the implant and seed decks):
    dt starts at T/ninc, which is also the largest increment (DTMAX of the data line); when an increment fails
    (NoConvergence from `increment`: residual not finite, maxit reached, or PNEWDT < 1) it is retried from the same
    start state with cutback * dt, up to max_cutbacks times in one increment (Abaqus D_f = 0.25, I_A = 5), then
    RuntimeError; after two consecutive increments converged in at most `easy` iterations (I_G = 4) dt grows by
    `growth` (D_D = 1.5), never above T/ninc; the last increment ends at T. Without a failure the increments are the
    ninc fixed ones of the plain loop (same results). kw go to `increment` (finite, md, cur, tol, maxit).
    log: a list of (increment, time at its start, dt, iterations or 'cut') per attempt.
    Returns U, SV, info: increments, attempts, cutbacks, iterations, smallest dt."""
    dtmax = T / ninc
    dt, t, kinc, easy_run = dtmax, 0.0, kinc0, 0
    info = dict(increments=0, attempts=0, cutbacks=0, iterations=0, dt_min=dtmax)
    while t < T - 1e-12 * T:
        if T - t <= dt * (1 + 1e-9):
            dt = T - t
        cuts = 0
        while True:
            info["attempts"] += 1
            try:
                Un, SVn, it = increment(U.copy(), U, fixed, nodes, elems, SV, props, dt, kinc, ttot=t0 + t, **kw)
                break
            except NoConvergence as err:
                cuts += 1; info["cutbacks"] += 1
                if log is not None:
                    log.append((kinc, t, dt, "cut"))
                print(f"  {err}; increment {kinc} retried with dt {dt * cutback:.3e}", file=sys.stderr, flush=True)
                if cuts > max_cutbacks:
                    raise RuntimeError(f"increment {kinc}: {max_cutbacks} cutbacks did not converge (dt {dt:.3e})") from err
                dt *= cutback; easy_run = 0
        if log is not None:
            log.append((kinc, t, dt, it))
        U, SV = Un, SVn
        t += dt; kinc += 1
        info["increments"] += 1; info["iterations"] += it; info["dt_min"] = min(info["dt_min"], dt)
        easy_run = easy_run + 1 if it <= easy else 0
        if easy_run >= 2 and dt < dtmax:
            dt = min(dt * growth, dtmax)
    return U, SV, info


def seed(n, antibiotic=False, finite=True, nu=0.3, md=False, L0=2.0, lam=0.0, psiref=0.0, cur=False, Y=1.0):
    import time as _t
    from make_seed_inp import cube, seed_phi
    sys.path.insert(0, str(HERE.parents[1] / "JAXFEM")); sys.path.insert(0, str(HERE.parents[1]))
    from fuel_front_3d_species import edge_disp
    nodes, elems = cube(n)
    T, nin, t0 = 0.2375, 60, 1e-8
    P = props_base(s=1.0, D=0.2, L0=L0, lam=lam, psiref=psiref, nu=nu, Y=Y, mu_b=Y * np.log(0.1))  # as make_seed_inp.deck (c_stop = 0.1 for any Y)
    if antibiotic:
        P = np.r_[P, t0 + 0.5 * T, 0.65]
    props = F(P)
    SV = [F(np.zeros(NSVT)) for _ in elems]
    U0 = np.zeros(NDN * len(nodes)); U = U0.copy(); U[3::NDN] = seed_phi(nodes); U[4::NDN] = 0.1
    U, SV, _ = increment(U0, U, np.ones(len(U), bool), nodes, elems, SV, props, t0, 1, finite=finite, md=md, cur=cur)
    fixed = np.zeros(len(U), bool)
    base = np.isclose(nodes[:, 2], 0.0)
    for d in range(3):
        fixed[NDN * np.where(base)[0] + d] = True
    fixed[NDN * np.where(np.isclose(nodes[:, 0], 1.0))[0] + 4] = True
    U[NDN * np.where(np.isclose(nodes[:, 0], 1.0))[0] + 4] = 1.0
    tt = time_start = _t.time()
    for k in range(nin):
        U, SV, it = increment(U.copy(), U, fixed, nodes, elems, SV, props, T / nin, k + 1, ttot=t0 + k * T / nin, finite=finite, md=md, cur=cur)
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
    print(f"seed {n}^3, L0 {L0}, lam {lam}, psiref {psiref}{', antibiotic' if antibiotic else ''}, {'finite' if finite else 'small'} strain, nu {nu}{', mean dilatation' if md else ''}{', nutrient in the current configuration' if cur else ''}: biofilm volume (reference) {(phi[:-1, :-1, :-1].sum() * h ** 3):.4f} "
          f"(node sum), phi at centre {phi[j, j, 1]:.3f}, alpha-1 nutrient side (element at x = 0.7, z near 0) "
          f"{xi[int(0.7 * n), j - 1:j + 1, 0].mean():.3f}, far side (x = 0.3) {xi[int(0.3 * n), j - 1:j + 1, 0].mean():.3f}; "
          f"edge 30 deg nutrient side / far side {edge_disp(np.linspace(0, 1, nn), u, 30, 0.0):.4f} / "
          f"{edge_disp(np.linspace(0, 1, nn), u, 30, np.pi):.4f}; biofilm (phi > 0.5) von Mises mean / max "
          f"{st[inb, 0].mean():.4e} / {st[inb, 0].max():.4e}, pressure mean {st[inb, 1].mean():.4e}; "
          f"{_t.time() - time_start:.0f} s")
    return dict(alpha_seed=float(np.mean([sv[:NSV * NGP].reshape(NGP, NSV)[:, 12].mean() for sv, b in zip(SV, inb) if b])),
                vm_mean=float(st[inb, 0].mean()), vm_max=float(st[inb, 0].max()), p_mean=float(st[inb, 1].mean()),
                edge=float(edge_disp(np.linspace(0, 1, nn), u, 30, 0.0)), SV=SV, U=U, nodes=nodes, elems=elems)


if __name__ == "__main__":
    if "--seed" in sys.argv:
        get = lambda k, d: float(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
        seed(int(sys.argv[sys.argv.index("--seed") + 1]), "--antibiotic" in sys.argv, "--small" not in sys.argv,
             get("--nu", 0.3), "--md" in sys.argv, get("--L0", 2.0), get("--lam", 0.0), get("--psiref", 0.0),
             "--cur" in sys.argv)
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

"""coupled_living_consumption.py -- the coupling agreed with Oliver on 20 Oct 2026 (MEETING_2026-10-20.md), in Python,
with synthetic data only.

Macro (3D, Klempt et al. 2024), the Chapter 4 seed problem on the 2 mm cube, NF^3 cells (cell centres stand for the
quadrature points):
    d phi/dt   = beta lap phi + k_alpha alpha - penalty            Eq. 34; with --front F also the front term
                 - r c/(k + c) n_c . grad phi                        as printed (r = F x Table 2 converted, k = 1),
                                                                   upwinded
    d alpha/dt = k_alpha phi                                       Eq. 36, unchanged
    dc/dt = d lap c - g phi (sum_i phibar_i) c                     Eq. 35 (--transient-c d; without it quasi-static 0 = ...), consumption by the LIVING
                                                                   biomass only: phibar_i = phi_i psi_i of the point
                                                                   model (Klempt et al. 2026; Oliver's form, decided
                                                                   20 Oct, replacing the normalised share)
    0 = d_a lap a - g_a phi L a                                    the second nutrient field repurposed for the
                                                                   antibiotic, same equation
  c = 1 on the wall y = -1 mm (nutrient), a = 1 on the wall y = +1 mm from T* = T_AB on (antibiotic), zero flux
  elsewhere; Thiele numbers Lambda^2 = g L^2/d = 4 for both (L = 1 mm).

Meso (0D point model at every colonised cell): a black box, the equations of Klempt et al. 2026 unchanged, its own
constraint sum phi_i + phi_0 = 1 (no phi_0 = 1 - phi_3D). Synthetic data: Fritsch et al. 2025 Case I, two species,
true parameters theta* = [a11, a12, a22, b1, b2] = [1, 0.1, 1, 1, 2], eta = [1, 2], c* = 100, alpha* = 10, 1000 steps
of 1e-4 (so s = 0.1: T* = 1 spans their run). Inputs per cell, linearly interpolated in time over the inner steps of
each outer step:
    c*(t) = 100 c(x, t),   alpha*(t) = 10 a(x, t)
Output: phi_i, psi_i -> L = sum phi_i psi_i, fed back to the consumption of the next outer step (lagged one step).

Order in an outer step t_{i-1} -> t_i (dT* = 0.025): (1) phi, alpha (independent of the point model); (2) c_i, a_i
from phi(t_i) and L(t_{i-1}); (3) point models with c*, alpha* interpolated between (c_{i-1}, a_{i-1}) and
(c_i, a_i); (4) cells newly above phi_min = 0.01 start their point model from the phi-weighted mean state of their
colonised neighbours (an assumption: the composition is carried into new biofilm from where it grows from);
(5) L(t_i).

Examples (Oliver's sketch): "separate": species A (point-model species 1) in the centre seed (the 32 seed elements of
Chapter 4), species B in the four vertical edges of the cube (0.25 mm square columns); "mixed": A and B together in
the centre seed, Fritsch's start (0.25, 0.30). Start amounts: Fritsch's total 0.55; a mono-species point has the other
species at 1e-6 (the point model clips to 1e-10 and would revive a zero; 1e-6 makes that explicit).
Antibiotic from T* = 0.5; each example also without antibiotic.

    python ansys_usermat/apdl/coupled_living_consumption.py [--nf 32]
        -> ansys_usermat/apdl/results/2026-10-3d/coupled_living_consumption.json, figdata/coupled_living_<case>.npz,
           assets/fig_coupled_living_consumption.pdf/.png
"""
import json
import math
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

jax.config.update("jax_enable_x64", True)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (HERE, ROOT / "ansys_usermat", ROOT / "ansys_usermat/coupling", ROOT, ROOT / "JAXFEM"):
    sys.path.insert(0, str(p))
import ecology_jax as eco  # noqa: E402
import figstyle  # noqa: E402
import wp3_reference as W  # noqa: E402
import wp3_stress_reference as S  # noqa: E402
from fritsch2025_cases import CASES  # noqa: E402

NF = int(sys.argv[sys.argv.index("--nf") + 1]) if "--nf" in sys.argv else 32
FRONT = float(sys.argv[sys.argv.index("--front") + 1]) if "--front" in sys.argv else 0.0
FAST = "--fast" in sys.argv                # quarter domain x >= 0, z >= 0 (mirror symmetry) and AMG-preconditioned CG
ARTDIFF = "--artdiff" in sys.argv          # as the ANSYS element (addendum 27): central front term + beta + r h / 4
TVD = "--tvd" in sys.argv                  # second-order upwind (MUSCL, minmod) for the front term instead of first order
R_FRONT = FRONT * 100.0e-3 * (2.0 / 20.0e-3)  # Table 2 r = 100 um/T* on the 20 um cube -> 10 mm/T* on the 2 mm cube
K_MONOD = 1.0
def _opt(name, default, cast=float):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


NM, T_END = 16, 1.0   # stress mesh (NM = NF below 16^3, see after --nf)
T_END = float(sys.argv[sys.argv.index("--tend") + 1]) if "--tend" in sys.argv else T_END
DT_C = _opt("--dtc", 0.025)                 # outer (coupling) step; halved to test the one-step lag
PHI_MIN = _opt("--phi-min", 0.01)
ONLY = _opt("--only", None, str)            # run only this example, e.g. separate_ab
CHECK = "--check" in sys.argv               # c = 1 and a = 1 everywhere from the start: the seed point model must
                                            # equal Fritsch Case I run alone (klempt2026_reproduction)
LAM_C2 = LAM_A2 = 0.0 if CHECK else 4.0     # g L^2 / d, L = 1 mm
# example setting v2 (20 Oct, after the first runs): --v2 = B in four small corner blocks at the seed's height,
# antibiotic weakly taken up (Lambda_a^2 = 0.25), alpha*_0 = 50 (Fritsch's validation value), antibiotic from all
# outer faces from T* = 0.5; use with --phi-min 0.1 --front 0.1
V2 = "--v2" in sys.argv
if V2 and not CHECK:
    LAM_A2 = 0.25
ALPHA0 = 50.0 if V2 else None
# sensitivities (as the ANSYS fsens_/u1_ runs): --alpha0 A (alpha*_0 of the point model), --gab G (antibiotic uptake
# g_a L^2 / d_a; the ANSYS "uptake x16" run is 4)
if "--alpha0" in sys.argv:
    ALPHA0 = float(sys.argv[sys.argv.index("--alpha0") + 1])
if "--gab" in sys.argv:
    LAM_A2 = float(sys.argv[sys.argv.index("--gab") + 1])
# --gc G: nutrient consumption g L^2 / d (default 4, an assumption); --mix-share S: share of A in the mixed start
# (same total as Fritsch's; 0.5 = --mix-equal); --tab T: antibiotic from T* = T (default 0.5)
if "--gc" in sys.argv:
    LAM_C2 = float(sys.argv[sys.argv.index("--gc") + 1])
# --transient-c D: keep the time derivative of Eq. 35 (and of the antibiotic field), (1/D) du/dt = lap u - q phi L u,
# D in mm^2/T* (the element's MY_DIFF = 1), backward Euler with the field step; default: quasi-steady (0 = ...)
TRANSIENT_D = float(sys.argv[sys.argv.index("--transient-c") + 1]) if "--transient-c" in sys.argv else None
AB_TOP = "--ab-top" in sys.argv
C_SINE = "--c-sine" in sys.argv
CKPT_SEC = float(os.environ.get("CKPT_SEC", 120.0))   # checkpoint interval (seconds)
SNAP_EVERY = float(sys.argv[sys.argv.index("--snap-every") + 1]) if "--snap-every" in sys.argv else None
MIX_SHARE = float(sys.argv[sys.argv.index("--mix-share") + 1]) if "--mix-share" in sys.argv else None
# --fortran: the point model of the ANSYS runs itself (ecology_native.f through ctypes, ecology_native_ctypes.py),
# with the same point-model file as ANSYS; otherwise the JAX version (agrees to 1e-16, 21 Oct)
FORTRAN = "--fortran" in sys.argv
# --ansys-bc: the boundary conditions of the ANSYS decks (A-19/A-21): the nutrient held at c = 1 on the whole element
# layer y < -0.75 mm of the 8^3 deck, the antibiotic at a = 1 on the whole outer shell of elements (any |x_i| > 0.75 mm)
# from T* = 0.5; without it, c and a are held on the faces only
ANSYS_BC = "--ansys-bc" in sys.argv
# --start-position: newly colonised cells start from a fresh start state chosen by position, the rule of the ANSYS
# runs (A-19): in "separate" species B where |x| and |y| > 0.5625 mm, species A elsewhere; in "mixed" the mixed start.
# Without it they inherit the phi-weighted mean state of their colonised neighbours (psi included).
START_POSITION = "--start-position" in sys.argv
# --mix-equal: the mixed start with equal amounts of A and B (Oliver, 20 Oct: phi_A = phi_B at the start), keeping
# the file's total; without it the file's start (Fritsch case I: 0.25, 0.30)
MIX_EQUAL = "--mix-equal" in sys.argv
T_AB = -1.0 if CHECK else (float(sys.argv[sys.argv.index("--tab") + 1]) if "--tab" in sys.argv else 0.5)
CASE = CASES["fritsch_2sp"]
# --case k26c3: Klempt et al. 2026, Table 2, case 3 (four species): A = 1/2 [[1,5,5,5],[5,1,3,3],[5,3,1,2],[5,3,2,1]],
# eta = (0.8, 1, 1.5, 2), b = (0.4, 0.3, 0.2, 0.1), start phi = 0.02 each, c* = 50 + 50 sin(500 t) (use --c-sine),
# alpha* = 10 from the start (use --tab 0), dt = 1e-4, 1500 steps (the four-species runs; antibiotic of case 4 after
# one third = step 500). In the examples species 1, 2 play "A" (seed) and 3, 4 play "B" (corners).
K26C3 = "--case" in sys.argv and sys.argv[sys.argv.index("--case") + 1] == "k26c3"
if K26C3:
    CASE = dict(n=4, A=[[0.5, 2.5, 2.5, 2.5], [2.5, 0.5, 1.5, 1.5], [2.5, 1.5, 0.5, 1.0], [2.5, 1.5, 1.0, 0.5]],
                b=[0.4, 0.3, 0.2, 0.1], eta=[0.8, 1.0, 1.5, 2.0], phi_init=[0.02] * 4, c_star=100.0,
                alpha_star=10.0, dt=1.0e-4, steps=1500)
N_SP = CASE.get("n", 2)
GROUP_A = [0, 1] if N_SP == 4 else [0]
S_LINK = CASE["steps"] * CASE["dt"] / T_END  # 0.1: T* = 1 spans the 1000 steps (0.15 and 1500 for k26c3)
H = W.L / NF
NM = min(NM, NF)                              # the stress mesh cannot be finer than the field grid (8^3 runs)
SHAPE = (NF // 2, NF, NF // 2) if FAST else (NF,) * 3
QS = (slice(NF // 2, None), slice(None), slice(NF // 2, None)) if FAST else (slice(None),) * 3   # full -> computed


def unfold(u):
    """quarter field (x >= 0, z >= 0) -> full cube by mirroring; full fields pass through"""
    if not FAST or u.shape[0] == NF:
        return u
    u = np.concatenate([u[::-1], u], 0)
    return np.concatenate([u[:, :, ::-1], u], 2)
# --obstacle columns: four rigid bars after Klempt et al. 2024, Sect. 4.3 and Fig. 9 (growth against rigid obstacles):
# no biofilm and no nutrient inside, zero flux through their faces, frictionless contact with the biofilm in the
# mechanics (wp3_stress_reference.solve_fields, obstacles=...). Positions read off the schematic (an assumption): two
# bars along x beside the seed (y in [-0.125, 0.125] mm, |z| in [0.5, 0.75] mm) and two along z below and above it
# (x in [-0.125, 0.125] mm, |y| in [0.5, 0.75] mm, the lower one between the seed and the nutrient wall y = -1 mm);
# square 0.25 mm sections tangent to the seed, aligned to the 16^3 and 32^3 grids.
OBSTACLE = _opt("--obstacle", None, str)


def obstacle_boxes(kind):
    if kind != "columns":
        raise SystemExit(f"unknown obstacle {kind}")
    B = []
    for lo, hi in ((-0.75, -0.5), (0.5, 0.75)):
        B.append(((-1.0, -0.125, lo), (1.0, 0.125, hi)))      # along x, beside the seed
        B.append(((-0.125, lo, -1.0), (0.125, hi, 1.0)))      # along z, below / above the seed
    return B


def box_mask(n, boxes):
    """cells of the n^3 grid of the 2 mm cube whose centres lie in one of the boxes"""
    xc = -1 + (np.arange(n) + 0.5) * (2.0 / n)
    X_, Y_, Z_ = np.meshgrid(xc, xc, xc, indexing="ij")
    m = np.zeros((n,) * 3, bool)
    for lo, hi in boxes:
        m |= (X_ > lo[0]) & (X_ < hi[0]) & (Y_ > lo[1]) & (Y_ < hi[1]) & (Z_ > lo[2]) & (Z_ < hi[2])
    return m


BOXES = obstacle_boxes(OBSTACLE) if OBSTACLE else None
OBST = box_mask(NF, BOXES)[QS] if OBSTACLE else None
if OBSTACLE:
    _open = ~OBST
    NB = {}                                    # NB[a, s]: the neighbour in direction s along axis a is open fluid
    for _a in range(3):
        for _s in (1, -1):
            _m = np.roll(_open, -_s, _a) & _open
            _e = [slice(None)] * 3; _e[_a] = -1 if _s == 1 else 0; _m[tuple(_e)] = False
            NB[_a, _s] = _m
RES = HERE / "results" / ("2026-10-3d-case3" if (K26C3 or "--c-sine" in sys.argv) else "2026-10-3d")
RES.mkdir(exist_ok=True)
TAG = (("" if NF == 32 else f"_nf{NF}") + (f"_front{FRONT:g}" if FRONT else "")
       + ("" if DT_C == 0.025 else f"_dtc{DT_C:g}") + ("" if PHI_MIN == 0.01 else f"_phimin{PHI_MIN:g}")
       + ("_check" if CHECK else "") + ("_v2" if V2 else "") + ("_fortran" if FORTRAN else "")
       + ("_ansysbc" if ANSYS_BC else "") + ("_startpos" if START_POSITION else "")
       + ("_mixeq" if MIX_EQUAL else "") + ("" if T_END == 1.0 else f"_tend{T_END:g}")
       + ("" if "--alpha0" not in sys.argv else f"_a{ALPHA0:g}") + ("" if "--gab" not in sys.argv else f"_gab{LAM_A2:g}")
       + ("" if "--gc" not in sys.argv else f"_gc{LAM_C2:g}") + ("" if MIX_SHARE is None else f"_share{MIX_SHARE:g}")
       + ("" if "--tab" not in sys.argv else f"_tab{T_AB:g}") + ("" if SNAP_EVERY is None else "_snaps") + ("_abtop" if AB_TOP else "") + ("" if TRANSIENT_D is None else f"_transc{TRANSIENT_D:g}")
       + ("" if os.environ.get("VOID_FLOOR") is None else f"_void{float(os.environ['VOID_FLOOR']):g}")
       + ("_k26c3" if K26C3 else "") + ("_csine" if "--c-sine" in sys.argv else "")
       + ("_tvd" if TVD else "") + ("_artdiff" if ARTDIFF else "") + ("_fast" if FAST else "") + (f"_obst{OBSTACLE}" if OBSTACLE else "") + (f"_{ONLY}" if ONLY else ""))
if FORTRAN:
    import ecology_native_ctypes as FN
    _eco = ROOT / "abaqus_composition" / ("eco_fritsch_2sp_a50.txt" if V2 else "eco_fritsch_2sp.txt")
    if K26C3:                                   # point-model file of Klempt 2026 Table 2 case 3
        _A, _b = CASE["A"], CASE["b"]
        _th = [_A[0][0], _A[0][1], _A[1][1], _b[0], _b[1], _A[2][2], _A[2][3], _A[3][3], _b[2], _b[3],
               _A[0][2], _A[0][3], _A[1][2], _A[1][3], 0, 0, 0, 0, 0, 0]
        _eco = ROOT / "abaqus_composition" / "eco_klempt2026_t2c3.txt"
        _eco.write_text("4\n100.0 10.0\n" + " ".join(f"{v:g}" for v in CASE["eta"] + [1.0]) + "\n"
                        + " ".join(f"{v:g}" for v in _th) + "\nphi_init 0.02 0.02 0.02 0.02 0.0\nnewton 50 1e-20\n")
    if "--alpha0" in sys.argv:                  # the same file with alpha*_0 replaced (second line: c*_0 alpha*_0)
        _l = _eco.read_text().splitlines()
        _l[1] = f"{_l[1].split()[0]} {ALPHA0:g}"
        _eco = ROOT / "abaqus_composition" / f"{_eco.stem.replace('_a50', '')}_a{ALPHA0:g}.txt"
        _eco.write_text("\n".join(_l) + "\n")
    _cst = FN.init(_eco)
    assert _cst["nact"] == N_SP and _cst["c_star"] == 100.0
    assert ALPHA0 is None or _cst["alpha_star"] == ALPHA0


# ---------------------------------------------------------------- point model (black box)
def theta20():
    A, b = CASE["A"], CASE["b"]
    th = [0.0] * 20
    th[0], th[1], th[2], th[3], th[4] = A[0][0], A[0][1], A[1][1], b[0], b[1]
    if N_SP == 4:
        th[5:14] = [A[2][2], A[2][3], A[3][3], b[2], b[3], A[0][2], A[0][3], A[1][2], A[1][3]]
    return jnp.asarray(th, dtype=jnp.float64)


N_SUB = max(1, math.ceil(S_LINK * DT_C / CASE["dt"] - 1e-9))
DT_SUB = S_LINK * DT_C / N_SUB
_A, _B = eco.theta_to_matrices(theta20())
_P = eco.default_hparams(DT_SUB)
_P.update(A=_A, b_diag=_B, active_mask=eco.active_mask(N_SP),
          Eta=jnp.asarray(CASE["eta"] + [1.0] * (5 - N_SP), dtype=jnp.float64),
          EtaPhi=jnp.asarray(CASE["eta"] + [1.0] * (5 - N_SP), dtype=jnp.float64))
_P = jax.tree_util.tree_map(jnp.asarray, _P)
FRAC = (jnp.arange(N_SUB) + 1.0) / N_SUB    # interpolation weight at the end of each inner step


import jax_hamilton_0d_5species_demo as _demo  # noqa: E402
NW_MAX, NW_TOL = 50, 1e-20                  # as the Fortran point model (A-20): stop once |Q|^2 < tol


def _newton(g_prev, q):
    """Newton until |Q|^2 < NW_TOL (at most NW_MAX), as ecology_native.f; the first step from the start state
    needs ~40 iterations, the fixed 6 of the demo's newton_step leave it unconverged (A-20)."""
    mask = q["active_mask"]
    F = lambda gg: _demo.residual(gg, g_prev, q)  # noqa: E731

    def body(st):
        it, g, _ = st
        g = _demo.clip_state(g, mask)
        Q = F(g)
        done = jnp.sum(Q * Q) < NW_TOL
        step = _demo.clip_state(g + jnp.linalg.solve(jax.jacfwd(F)(g), -Q), mask)
        return it + 1, jnp.where(done, g, step), done

    return jax.lax.while_loop(lambda st: (st[0] < NW_MAX) & ~st[2], body, (0, g_prev, False))[1]


def _one_point(g, c0, c1, a0, a1):
    def body(gg, w):
        p = dict(_P)
        p["c"] = CASE["c_star"] * ((1 - w) * c0 + w * c1)
        p["alpha"] = (ALPHA0 or CASE["alpha_star"]) * ((1 - w) * a0 + w * a1)
        return _newton(gg, p), None
    return jax.lax.scan(body, g, FRAC)[0]


STEP = jax.jit(jax.vmap(_one_point))


def start_state(*phis):
    g = np.zeros(12)
    g[:len(phis)] = phis
    g[5] = 1.0 - sum(phis)
    g[6:11] = 0.999
    return g


def group_state(tot, group):
    """start state with the amount tot shared by the species in group, 1e-6 for the others"""
    other = [i for i in range(N_SP) if i not in group]
    ph = [1e-6] * N_SP
    for i in group:
        ph[i] = (tot - 1e-6 * len(other)) / len(group)
    return start_state(*ph)


def living(G):
    """living volume fraction of the point model, sum_i phibar_i = sum_i phi_i psi_i (Klempt et al. 2026 notation)"""
    return (G[:, :N_SP] * G[:, 6:6 + N_SP]).sum(1)


def grad(u):
    """central differences, one-sided at the walls (and at obstacle faces)"""
    if OBSTACLE:
        out = []
        for a in range(3):
            mp, mm = NB[a, 1], NB[a, -1]
            fw = np.where(mp, (np.roll(u, -1, a) - u) / H, 0.0); bw = np.where(mm, (u - np.roll(u, 1, a)) / H, 0.0)
            if FAST and a in (0, 2):             # mirror cell on the symmetry plane: a neighbour with difference 0
                mm = mm.copy(); e = [slice(None)] * 3; e[a] = 0
                mm[tuple(e)] = ~OBST[tuple(e)]
            out.append((fw + bw) / np.maximum(mp.astype(float) + mm, 1.0))
        return out
    if FAST:                                 # mirror layer on the symmetry planes x = 0, z = 0
        up = np.concatenate([u[:1], u], 0); up = np.concatenate([up[:, :, :1], up], 2)
        return [np.gradient(up, H, axis=a)[1:, :, 1:] for a in range(3)]
    return [np.gradient(u, H, axis=a) for a in range(3)]


def lap_obst(u, out):
    """7-point Laplacian with zero flux through the walls and the obstacle faces"""
    out[...] = 0.0
    for a in range(3):
        for s_ in (1, -1):
            out += np.where(NB[a, s_], np.roll(u, -s_, a) - u, 0.0)
    out /= H * H
    out[OBST] = 0.0
    return out




def _nb(a, s_):
    """neighbour in direction s_ along axis a exists and is open (walls, obstacles)"""
    if OBSTACLE:
        return NB[a, s_]
    m = np.ones(SHAPE, bool)
    e = [slice(None)] * 3; e[a] = -1 if s_ == 1 else 0; m[tuple(e)] = False
    return m


def upwind_adv_tvd(u, v):
    """v . grad u with MUSCL face values (minmod slopes), upwind; zero derivative into walls as the first-order form"""
    out = np.zeros_like(u)
    for a in range(3):
        mp, mm = _nb(a, 1), _nb(a, -1)
        up_, um_ = np.roll(u, -1, a), np.roll(u, 1, a)
        dp = np.where(mp, up_ - u, 0.0); dm = np.where(mm, u - um_, 0.0)
        sl = np.where(dp * dm > 0, np.sign(dp) * np.minimum(np.abs(dp), np.abs(dm)), 0.0)
        # v > 0: (u_{i+1/2}^- - u_{i-1/2}^-) / H, faces from the left cells
        fr = u + 0.5 * sl; fl = np.roll(u + 0.5 * sl, 1, a)
        dP = np.where(mm, (fr - fl) / H, 0.0)
        # v < 0: faces from the right cells
        fl2 = u - 0.5 * sl; fr2 = np.roll(u - 0.5 * sl, -1, a)
        dN = np.where(mp, (fr2 - fl2) / H, 0.0)
        out += np.maximum(v[a], 0) * dP + np.minimum(v[a], 0) * dN
    return out


def upwind_adv(u, v):
    if TVD:
        return upwind_adv_tvd(u, v)
    out = np.zeros_like(u)
    for a in range(3):
        fw = (np.roll(u, -1, a) - u) / H; bw = (u - np.roll(u, 1, a)) / H
        if OBSTACLE:
            fw = np.where(NB[a, 1], fw, 0.0); bw = np.where(NB[a, -1], bw, 0.0)
        else:
            e = [slice(None)] * 3; e[a] = -1; fw[tuple(e)] = 0.0
            e[a] = 0; bw[tuple(e)] = 0.0
        out += np.maximum(v[a], 0) * bw + np.minimum(v[a], 0) * fw
    return out


def front_central(phi, c):
    """the front term with central differences for grad phi (the element's NEM gradient is of central type)"""
    g = grad(c)
    nrm = np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2 + 1e-6)
    sp_ = R_FRONT * c / (K_MONOD + c)
    gp = grad(phi)
    return -sum(sp_ * g[a] / nrm * gp[a] for a in range(3))


def front(phi, c):
    """-r c/(k + c) n_c . grad phi (Klempt 2024 Eq. 34 as printed), upwinded"""
    if ARTDIFF:
        return front_central(phi, c)
    g = grad(c)
    nrm = np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2 + 1e-6)
    sp_ = R_FRONT * c / (K_MONOD + c)
    return -upwind_adv(phi, [sp_ * x / nrm for x in g])


# ---------------------------------------------------------------- macro fields
class Reaction:
    """0 = lap u - q phi L u, u = value on the face y = side (ghost cell), zero flux elsewhere."""

    def __init__(self, q, side):
        n3 = int(np.prod(SHAPE))
        idx = np.arange(n3).reshape(SHAPE)
        rows, cols, diag = [], [], np.zeros(SHAPE)
        for ax in range(3):
            for sgn in (1, -1):
                ok = np.ones(SHAPE, bool)
                e = [slice(None)] * 3; e[ax] = -1 if sgn == 1 else 0; ok[tuple(e)] = False
                rows.append(idx[ok]); cols.append(np.roll(idx, -sgn, ax)[ok]); diag += ok
        faces = [(1, 0 if side < 0 else -1)] if side != 0 else [(a, e) for a in range(3) for e in (0, -1)
                                                                 if not (FAST and a in (0, 2) and e == 0)]
        for ax_, j in faces:
            sl = [slice(None)] * 3; sl[ax_] = j
            diag[tuple(sl)] += 2.0
        if OBSTACLE:                          # links only between open cells; obstacle rows u = 0
            o = (~OBST).ravel()
            keep = [o[r] & o[c_] for r, c_ in zip(rows, cols)]
            rows = [r[k] for r, k in zip(rows, keep)]; cols = [c_[k] for c_, k in zip(cols, keep)]
            diag = np.zeros(n3)
            for r in rows:
                np.add.at(diag, r, 1.0)
            diag = diag.reshape(SHAPE)
            for ax_, j in faces:
                sl = [slice(None)] * 3; sl[ax_] = j
                diag[tuple(sl)] += 2.0
            diag[OBST] = 1.0
        self.A0 = sp.csr_matrix((-np.ones(sum(r.size for r in rows)) / H ** 2,
                                 (np.concatenate(rows), np.concatenate(cols))), shape=(n3, n3))
        self.d0 = diag.ravel() / H ** 2
        b = np.zeros(SHAPE)
        for ax_, j in faces:
            sl = [slice(None)] * 3; sl[ax_] = j
            b[tuple(sl)] += 2.0 / H ** 2
        if OBSTACLE:
            b[OBST] = 0.0
        self.b, self.q = b.ravel(), q

    def solve(self, phi, L, value=1.0, u_prev=None, dtau=None):
        if value == 0.0:
            return np.zeros(SHAPE)
        m = 0.0 if (TRANSIENT_D is None or u_prev is None) else 1.0 / (TRANSIENT_D * dtau)
        A = (self.A0 + sp.diags(self.d0 + m + self.q * (np.clip(phi, 0, 1) * L).ravel())).tocsr()
        b = self.b.copy() + (m * u_prev.ravel() / value if m else 0.0)
        if getattr(self, "fixed", None) is not None:      # cells held at the value (the ANSYS element layers)
            f = self.fixed.ravel()
            D = sp.diags((~f).astype(float))
            A = (D @ A + sp.diags(f.astype(float))).tocsr()
            b = np.where(f, 1.0, b)
        if FAST:
            return (value * self._amg(A, b)).reshape(SHAPE)
        return (value * spla.spsolve(A.tocsc(), b)).reshape(SHAPE)

    def _amg(self, A, b):
        """CG on the unknown cells (symmetric), preconditioned by an AMG hierarchy that is rebuilt every 25 solves;
        the matrix changes only on the diagonal between solves. Warm start from the previous solution."""
        import pyamg
        f = getattr(self, "fixed", None)
        f = np.zeros(A.shape[0], bool) if f is None else f.ravel()
        fr = ~f
        x = np.where(f, 1.0, 0.0)
        Aff = A[fr][:, fr].tocsr()
        rhs = b[fr] - A[fr][:, f] @ x[f] if f.any() else b[fr]
        self._n = getattr(self, "_n", 0) + 1
        if getattr(self, "_ml", None) is None or self._n % 25 == 1:
            self._ml = pyamg.smoothed_aggregation_solver(Aff, symmetry="symmetric")
        x0 = getattr(self, "_last", None)
        sol, info = spla.cg(Aff, rhs, x0=x0, rtol=1e-12, maxiter=500, M=self._ml.aspreconditioner())
        if info != 0:                              # fall back to the direct solve
            sol = spla.spsolve(Aff.tocsc(), rhs)
        self._last = sol
        x[fr] = sol
        return x


def c_supply(T):
    """held nutrient value at T*: 1, or with --c-sine 0.5 + 0.5 sin(500 t), t = S_LINK T* the point-model time
    (Klempt et al. 2026, Table 2, case 3: c* = 50 + 50 sin(500 t) with c* = c*_0 c, c*_0 = 100)"""
    return 0.5 + 0.5 * math.sin(500.0 * S_LINK * T) if C_SINE else 1.0


def neighbours_mean(G, w):
    """phi-weighted mean of the 6-neighbour states (w = weight per cell, 0 where not colonised)."""
    num = np.zeros_like(G); den = np.zeros(G.shape[:3])
    for ax in range(3):
        for sgn in (1, -1):
            wn = np.roll(w, sgn, ax); gn = np.roll(G, sgn, ax)
            e = [slice(None)] * 3; e[ax] = 0 if sgn == 1 else -1
            wn[tuple(e)] = 0.0
            num += wn[..., None] * gn; den += wn
    return num, den


def run(example, antibiotic):
    cells8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
    seed = W.seed_mask(NF, cells8)
    ax = -1 + (np.arange(NF) + 0.5) * H
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    edges = (np.abs(X) > 0.75) & (np.abs(Y) > 0.75)
    if V2:                                        # four small corner blocks at the seed's height
        zmax = np.abs(Z[seed]).max() + 1e-9
        edges = edges & (np.abs(Z) <= zmax)
    seed, edges, X, Y, Z = seed[QS], edges[QS], X[QS], Y[QS], Z[QS]
    phi = (seed | (edges if example == "separate" else False)).astype(float)
    alpha = np.ones_like(phi)
    G = np.zeros(SHAPE + (12,))
    tot = sum(CASE["phi_init"])
    GROUP_B = [i for i in range(N_SP) if i not in GROUP_A]
    if MIX_SHARE is not None:
        MIX_START = tuple(tot * MIX_SHARE / len(GROUP_A) if i in GROUP_A else tot * (1 - MIX_SHARE) / len(GROUP_B)
                          for i in range(N_SP))
    else:
        MIX_START = tuple([tot / N_SP] * N_SP) if MIX_EQUAL else tuple(CASE["phi_init"])
    if example == "separate":
        G[seed] = group_state(tot, GROUP_A)
        G[edges] = group_state(tot, GROUP_B)
    else:
        G[seed] = start_state(*MIX_START)
    on = phi >= PHI_MIN
    nut, ab = Reaction(LAM_C2, -1), Reaction(LAM_A2, 0 if V2 else +1)
    if ANSYS_BC:
        nut.fixed = Y < -0.75
        ab.fixed = ((np.abs(X) > 0.75) | (np.abs(Y) > 0.75) | (np.abs(Z) > 0.75)) if (V2 and not AB_TOP) else (Y > 0.75)
        if OBSTACLE:
            nut.fixed = nut.fixed & ~OBST; ab.fixed = ab.fixed & ~OBST
    L = np.where(on, living(G.reshape(-1, 12)).reshape(SHAPE), 0.0)
    c = nut.solve(phi, L, c_supply(0.0))
    a = ab.solve(phi, L, 1.0) if (antibiotic and T_AB < 0) else np.zeros_like(phi)
    lap = np.empty_like(phi)
    dt = 0.15 * H * H / (W.BETA + (R_FRONT * H / 4.0 if ARTDIFF else 0.0))
    if FRONT:
        dt = min(dt, 0.4 * H / (R_FRONT / (1 + K_MONOD)))
    per = int(np.ceil(DT_C / dt)); dt = DT_C / per
    hist = dict(t=[0.0], living_seed=[float(L[seed].mean())], share_a_seed=[], share_a_edges=[], living_edges=[],
                c_seed=[float(c[seed].mean())], a_seed=[0.0], colonised=[int(on.sum()) * (4 if FAST else 1)])

    def shares(m):
        ph = G[m][:, :N_SP]
        return float((ph[:, GROUP_A].sum(1) / ph.sum(1)).mean())
    hist["share_a_seed"].append(shares(seed)); hist["share_a_edges"].append(shares(edges) if example == "separate" else np.nan)
    hist["living_edges"].append(float(L[edges].mean()) if example == "separate" else np.nan)
    snaps = {}
    t0 = time.time()
    # checkpoint (a cloud container stops when idle): the state is saved every CKPT_SEC seconds and a rerun with the
    # same options continues from it; the file is removed when the run has finished
    ckpt = ROOT / "figdata" / "ckpt" / f"{example}{'_ab' if antibiotic else ''}{TAG}.npz"
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    k0 = 0
    if ckpt.exists():
        z = np.load(ckpt, allow_pickle=False)
        phi, alpha, c, a, L, on, G = (z[n].copy() for n in ("phi", "alpha", "c", "a", "L", "on", "G"))
        hist = json.loads(str(z["hist"]))
        snaps = {}
        for name in z.files:
            if name.startswith("snap|"):
                _, T, f = name.split("|")
                snaps.setdefault(float(T), {})[f] = z[name].copy()
        k0 = int(z["k_done"]) + 1
        print(f"[{example}, antibiotic {antibiotic}] resumed from the checkpoint after step {k0}", flush=True)
    t_ck = time.time()
    for k in range(k0, int(round(T_END / DT_C))):
        T1 = (k + 1) * DT_C
        # (1) macro phi, alpha
        cf = c
        for _ in range(per):
            if FRONT:
                cf = nut.solve(phi, L, c_supply(T1), u_prev=cf, dtau=dt)  # the front follows the nutrient
            def f(ph, al):
                (lap_obst(ph, lap) if OBSTACLE else W.laplacian(ph, H, lap))
                d = (W.BETA + (R_FRONT * H / 4.0 if ARTDIFF else 0.0)) * lap + W.K_ALPHA * al - W.PENALTY * (np.maximum(ph - 1, 0) + np.minimum(ph, 0))
                if FRONT:
                    d = d + front(ph, cf)
                if OBSTACLE:                    # no biofilm inside the obstacles
                    d = np.where(OBST, 0.0, d)
                    return d, np.where(OBST, 0.0, W.K_ALPHA * ph)
                return d, W.K_ALPHA * ph
            k1 = f(phi, alpha); k2 = f(phi + 0.5 * dt * k1[0], alpha + 0.5 * dt * k1[1])
            k3 = f(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1]); k4 = f(phi + dt * k3[0], alpha + dt * k3[1])
            phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            alpha = alpha + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        # (2) nutrient and antibiotic at t_i, with L(t_{i-1})
        if TRANSIENT_D is not None and FRONT:
            c_new = cf                        # already advanced over the outer step in the inner steps
        else:
            c_new = nut.solve(phi, L, c_supply(T1), u_prev=c, dtau=DT_C)
        a_new = ab.solve(phi, L, 1.0 if (antibiotic and T1 > T_AB + 1e-12) else 0.0, u_prev=a, dtau=DT_C)
        # (3) point models of the colonised cells
        idx = np.flatnonzero(on.ravel())
        Gf = G.reshape(-1, 12)
        if FORTRAN:
            Gf[idx], okf = FN.step_lin(Gf[idx], S_LINK * DT_C, N_SUB, c.ravel()[idx], c_new.ravel()[idx],
                                       a.ravel()[idx], a_new.ravel()[idx])
            assert okf.all(), "Fortran point model failed"
        else:
            Gf[idx] = np.asarray(STEP(jnp.asarray(Gf[idx]), jnp.asarray(c.ravel()[idx]),
                                      jnp.asarray(c_new.ravel()[idx]), jnp.asarray(a.ravel()[idx]),
                                      jnp.asarray(a_new.ravel()[idx])))
        G = Gf.reshape(SHAPE + (12,))
        c, a = c_new, a_new
        # (4) newly colonised cells inherit the neighbours' state
        new = (phi >= PHI_MIN) & ~on
        if new.any() and START_POSITION:
            if example == "separate":
                isb = (np.abs(X) > 0.5625) & (np.abs(Y) > 0.5625)
                G[new & isb] = group_state(tot, GROUP_B)
                G[new & ~isb] = group_state(tot, GROUP_A)
            else:
                G[new] = start_state(*MIX_START)
            on = on | new
        elif new.any():
            num, den = neighbours_mean(G, np.where(on, phi, 0.0))
            fill = new & (den > 0)
            gs = num[fill] / den[fill][:, None]
            gs[:, :N_SP] = np.clip(gs[:, :N_SP], 1e-10, None)
            gs[:, 5] = 1.0 - gs[:, :N_SP].sum(1)               # the point model's own constraint
            G[fill] = gs
            on = on | fill
        # (5) living biomass for the next step
        L = np.where(on, living(G.reshape(-1, 12)).reshape(SHAPE), 0.0)
        hist.setdefault("phi_max", []).append(float(phi.max())); hist.setdefault("phi_min", []).append(float(phi.min()))
        hist["t"].append(T1); hist["living_seed"].append(float(L[seed].mean()))
        hist["share_a_seed"].append(shares(seed))
        hist["share_a_edges"].append(shares(edges) if example == "separate" else np.nan)
        hist["living_edges"].append(float(L[edges].mean()) if example == "separate" else np.nan)
        hist["c_seed"].append(float(c[seed].mean())); hist["a_seed"].append(float(a[seed].mean()))
        hist["colonised"].append(int(on.sum()) * (4 if FAST else 1))
        if (abs(T1 - T_AB) < 1e-9 or abs(T1 - T_END) < 1e-9
                or (SNAP_EVERY and abs(T1 / SNAP_EVERY - round(T1 / SNAP_EVERY)) < 1e-6)):
            ph = G[..., :N_SP]
            shA = np.where(on, ph[..., 0] / np.maximum(ph.sum(-1), 1e-30), np.nan)
            snaps[round(T1, 3)] = dict(phi=unfold(phi).copy(), c=unfold(c).copy(), a=unfold(a).copy(),
                                       L=unfold(L).copy(), shareA=unfold(shA).copy())
        if time.time() - t_ck > CKPT_SEC:
            tmp = ckpt.with_suffix(".tmp.npz")
            np.savez(tmp, k_done=k, phi=phi, alpha=alpha, c=c, a=a, L=L, on=on, G=G,
                     hist=json.dumps(hist, default=float),
                     **{f"snap|{T}|{f}": v for T, dd in snaps.items() for f, v in dd.items()})
            os.replace(tmp, ckpt)
            t_ck = time.time()
    print(f"[{example}, antibiotic {antibiotic}] {time.time() - t0:.0f} s; colonised {int(on.sum())} cells; "
          f"seed: living {hist['living_seed'][-1]:.3f}, share A {hist['share_a_seed'][-1]:.3f}, "
          f"c {hist['c_seed'][-1]:.3f}", flush=True)
    if ckpt.exists():
        ckpt.unlink()
    hist["seed_state_end"] = G[seed].mean(0).tolist()
    hist["seed_state_spread"] = float(np.abs(G[seed] - G[seed].mean(0)).max())
    seed16 = S.seed_mask(NM, cells8)
    vm, p = S.solve_fields(NM, W.block_mean(unfold(phi), NM), W.block_mean(unfold(alpha) - 1, NM), obstacles=BOXES)
    if OBSTACLE:
        snaps.setdefault(round(T_END, 3), {}).update(vm16=vm, p16=p, obst=unfold(OBST).astype(np.int8))
        hist["contact"] = getattr(S, "LAST_CONTACT", None)
    hist["seed_am1"] = float((alpha - 1)[seed].mean()); hist["seed_vm_Pa"] = float(vm[seed16].mean())
    return hist, snaps


def mid(u):
    k = NF // 2
    return 0.5 * (u[:, :, k - 1] + u[:, :, k])


def main():
    R, SN = {}, {}
    keys = [f"{ex}{'_ab' if abx else ''}" for ex in ("separate", "mixed") for abx in (False, True)]
    if CHECK:
        keys = ["mixed_ab"]
    elif ONLY:
        keys = [ONLY]
    for key in keys:
            ex, abx = key.replace("_ab", ""), key.endswith("_ab")
            R[key], SN[key] = run(ex, abx)
            np.savez_compressed(ROOT / "figdata" / f"coupled_living_{key}{TAG}.npz",
                                **{f"{k}_{T}": v for T, d in SN[key].items() for k, v in d.items()})
    res = {"setting": {"field": f"{NF}^3", "point_model": "Fritsch et al. 2025 Case I, true parameters",
                       "s": S_LINK, "inner_steps_per_outer": N_SUB, "Lambda2_c": LAM_C2, "Lambda2_a": LAM_A2,
                       "antibiotic_from": T_AB, "phi_min": PHI_MIN}, "runs": R}
    if CHECK:
        import klempt2026_reproduction as K
        K.DT = CASE["dt"]
        ref = K.simulate(CASE, variant="newton30")   # converged Newton (6 fixed iterations do not converge for 4 species)
        got = np.array(R["mixed_ab"]["seed_state_end"])
        d = dict(phi=float(np.abs(got[:N_SP] - ref["phi"][-1]).max()),
                 psi=float(np.abs(got[6:6 + N_SP] - ref["psi"][-1]).max()),
                 phi0=float(abs(got[5] - ref["phi0"][-1])))
        res["check_vs_point_model_alone"] = dict(max_abs_diff=d, coupled_end=got.tolist(),
                                                 alone_phi=ref["phi"][-1].tolist(), alone_psi=ref["psi"][-1].tolist(),
                                                 spread_over_seed_cells=R["mixed_ab"]["seed_state_spread"])
        print("CHECK: coupled seed point model vs Fritsch Case I alone, max |diff|:", d,
              "spread over seed cells", R["mixed_ab"]["seed_state_spread"])
    (RES / f"coupled_living_consumption{TAG}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    if CHECK or ONLY:
        return
    print("\n| run | seed living T*=1 | seed share A | edges living | edges share A | seed c | seed a | colonised |")
    for k, h in R.items():
        print(f"| {k} | {h['living_seed'][-1]:.3f} | {h['share_a_seed'][-1]:.3f} | {h['living_edges'][-1]:.3f} | "
              f"{h['share_a_edges'][-1]:.3f} | {h['c_seed'][-1]:.3f} | {h['a_seed'][-1]:.3f} | {h['colonised'][-1]} |")
    print("seed alpha-1 and vM:",
          {k: (round(h['seed_am1'], 8), round(h['seed_vm_Pa'], 8)) for k, h in R.items()})

    figstyle.apply(size=8)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6.3, 4.6))
    gs = fig.add_gridspec(2, 4, hspace=0.45, wspace=0.55)
    e = np.linspace(-1, 1, NF + 1)
    panels = [("separate_ab", "shareA", r"share of A, $\phi_A/(\phi_A+\phi_B)$", "coolwarm_r", (0, 1)),
              ("separate_ab", "L", r"living $\Sigma\phi_i\psi_i$", "Greens", (0, None)),
              ("separate_ab", "c", "nutrient $c$", "Blues", (0, 1)),
              ("separate_ab", "a", "antibiotic $a$", "Reds", (0, 1))]
    for j, (key, fld, lab, cm, lim) in enumerate(panels):
        a_ = fig.add_subplot(gs[0, j])
        u = mid(SN[key][1.0][fld])
        im = a_.pcolormesh(e, e, np.ma.masked_invalid(u).T, cmap=cm, vmin=lim[0], vmax=lim[1], shading="flat",
                           rasterized=True)
        cb = fig.colorbar(im, ax=a_, fraction=0.046, pad=0.03, orientation="horizontal")
        cb.ax.tick_params(labelsize=6)
        a_.set_title(f"({'abcd'[j]}) {lab}", fontsize=7.5, pad=2)
        a_.set_aspect("equal"); a_.set_xticks([-1, 0, 1]); a_.set_yticks([-1, 0, 1]); a_.grid(False)
        a_.tick_params(labelsize=6)
    sty = {"separate": dict(color="C0"), "separate_ab": dict(color="C0", ls="--"),
           "mixed": dict(color="C2"), "mixed_ab": dict(color="C2", ls="--")}
    ax1 = fig.add_subplot(gs[1, 0:2]); ax2 = fig.add_subplot(gs[1, 2:4])
    for k, h in R.items():
        ax1.plot(h["t"], h["living_seed"], lw=1.1, label=k.replace("_ab", ", antibiotic"), **sty[k])
        ax2.plot(h["t"], h["share_a_seed"], lw=1.1, **sty[k])
    for a_ in (ax1, ax2):
        a_.axvline(T_AB, color="0.6", lw=0.6, ls=":")
        a_.set_xlabel("$T^*$"); a_.set_xlim(0, 1)
    ax1.set_ylabel(r"seed mean $\Sigma_i\phi_i\psi_i$"); ax2.set_ylabel(r"seed mean share of A")
    ax1.set_title("(e)", fontsize=7.5, loc="left"); ax2.set_title("(f)", fontsize=7.5, loc="left")
    ax1.legend(fontsize=6.5, frameon=False, ncol=2)
    out = ROOT / "assets" / f"fig_coupled_living_consumption{TAG.replace('.', 'p')}"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), dpi=300,
                                                                          bbox_inches="tight")
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

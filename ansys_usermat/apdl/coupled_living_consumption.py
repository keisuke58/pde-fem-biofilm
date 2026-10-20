"""coupled_living_consumption.py -- the coupling agreed with Oliver on 20 Oct 2026 (MEETING_2026-10-20.md), in Python,
with synthetic data only.

Macro (3D, Klempt et al. 2024), the Chapter 4 seed problem on the 2 mm cube, NF^3 cells (cell centres stand for the
quadrature points):
    d phi/dt   = beta lap phi + k_alpha alpha - penalty            Eq. 34; with --front F also the front term
                 - r c/(k + c) n_c . grad phi                        as printed (r = F x Table 2 converted, k = 1),
                                                                   upwinded
    d alpha/dt = k_alpha phi                                       Eq. 36, unchanged
    0 = d lap c - g phi (sum_i phibar_i) c                         Eq. 35, quasi-static, consumption by the LIVING
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
R_FRONT = FRONT * 100.0e-3 * (2.0 / 20.0e-3)  # Table 2 r = 100 um/T* on the 20 um cube -> 10 mm/T* on the 2 mm cube
K_MONOD = 1.0
def _opt(name, default, cast=float):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


NM, T_END = 16, 1.0
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
T_AB = -1.0 if CHECK else 0.5
CASE = CASES["fritsch_2sp"]
N_SP = 2
S_LINK = CASE["steps"] * CASE["dt"] / T_END  # 0.1: T* = 1 spans the 1000 steps
H = W.L / NF
RES = HERE / "results" / "2026-10-3d"
TAG = (("" if NF == 32 else f"_nf{NF}") + (f"_front{FRONT:g}" if FRONT else "")
       + ("" if DT_C == 0.025 else f"_dtc{DT_C:g}") + ("" if PHI_MIN == 0.01 else f"_phimin{PHI_MIN:g}")
       + ("_check" if CHECK else "") + ("_v2" if V2 else ""))


# ---------------------------------------------------------------- point model (black box)
def theta20():
    A, b = CASE["A"], CASE["b"]
    th = [0.0] * 20
    th[0], th[1], th[2], th[3], th[4] = A[0][0], A[0][1], A[1][1], b[0], b[1]
    return jnp.asarray(th, dtype=jnp.float64)


N_SUB = max(1, math.ceil(S_LINK * DT_C / CASE["dt"] - 1e-9))
DT_SUB = S_LINK * DT_C / N_SUB
_A, _B = eco.theta_to_matrices(theta20())
_P = eco.default_hparams(DT_SUB)
_P.update(A=_A, b_diag=_B, active_mask=eco.active_mask(N_SP),
          Eta=jnp.asarray(CASE["eta"] + [1.0] * 3, dtype=jnp.float64),
          EtaPhi=jnp.asarray(CASE["eta"] + [1.0] * 3, dtype=jnp.float64))
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


def start_state(phi_a, phi_b):
    g = np.zeros(12)
    g[0], g[1] = phi_a, phi_b
    g[5] = 1.0 - phi_a - phi_b
    g[6:11] = 0.999
    return g


def living(G):
    """living volume fraction of the point model, sum_i phibar_i = sum_i phi_i psi_i (Klempt et al. 2026 notation)"""
    return (G[:, :N_SP] * G[:, 6:6 + N_SP]).sum(1)


def grad(u):
    """central differences, one-sided at the walls"""
    return [np.gradient(u, H, axis=a) for a in range(3)]


def upwind_adv(u, v):
    out = np.zeros_like(u)
    for a in range(3):
        fw = (np.roll(u, -1, a) - u) / H; bw = (u - np.roll(u, 1, a)) / H
        e = [slice(None)] * 3; e[a] = -1; fw[tuple(e)] = 0.0
        e[a] = 0; bw[tuple(e)] = 0.0
        out += np.maximum(v[a], 0) * bw + np.minimum(v[a], 0) * fw
    return out


def front(phi, c):
    """-r c/(k + c) n_c . grad phi (Klempt 2024 Eq. 34 as printed), upwinded"""
    g = grad(c)
    nrm = np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2 + 1e-6)
    sp_ = R_FRONT * c / (K_MONOD + c)
    return -upwind_adv(phi, [sp_ * x / nrm for x in g])


# ---------------------------------------------------------------- macro fields
class Reaction:
    """0 = lap u - q phi L u, u = value on the face y = side (ghost cell), zero flux elsewhere."""

    def __init__(self, q, side):
        n3 = NF ** 3
        idx = np.arange(n3).reshape((NF,) * 3)
        rows, cols, diag = [], [], np.zeros((NF,) * 3)
        for ax in range(3):
            for sgn in (1, -1):
                ok = np.ones((NF,) * 3, bool)
                e = [slice(None)] * 3; e[ax] = -1 if sgn == 1 else 0; ok[tuple(e)] = False
                rows.append(idx[ok]); cols.append(np.roll(idx, -sgn, ax)[ok]); diag += ok
        faces = [(1, 0 if side < 0 else -1)] if side != 0 else [(a, e) for a in range(3) for e in (0, -1)]
        for ax_, j in faces:
            sl = [slice(None)] * 3; sl[ax_] = j
            diag[tuple(sl)] += 2.0
        self.A0 = sp.csr_matrix((-np.ones(sum(r.size for r in rows)) / H ** 2,
                                 (np.concatenate(rows), np.concatenate(cols))), shape=(n3, n3))
        self.d0 = diag.ravel() / H ** 2
        b = np.zeros((NF,) * 3)
        for ax_, j in faces:
            sl = [slice(None)] * 3; sl[ax_] = j
            b[tuple(sl)] += 2.0 / H ** 2
        self.b, self.q = b.ravel(), q

    def solve(self, phi, L, value=1.0):
        if value == 0.0:
            return np.zeros((NF,) * 3)
        A = (self.A0 + sp.diags(self.d0 + self.q * (np.clip(phi, 0, 1) * L).ravel())).tocsc()
        return (value * spla.spsolve(A, self.b)).reshape((NF,) * 3)


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
    phi = (seed | (edges if example == "separate" else False)).astype(float)
    alpha = np.ones_like(phi)
    G = np.zeros((NF, NF, NF, 12))
    tot = sum(CASE["phi_init"])
    if example == "separate":
        G[seed] = start_state(tot - 1e-6, 1e-6)
        G[edges] = start_state(1e-6, tot - 1e-6)
    else:
        G[seed] = start_state(*CASE["phi_init"])
    on = phi >= PHI_MIN
    nut, ab = Reaction(LAM_C2, -1), Reaction(LAM_A2, 0 if V2 else +1)
    L = np.where(on, living(G.reshape(-1, 12)).reshape((NF,) * 3), 0.0)
    c = nut.solve(phi, L)
    a = ab.solve(phi, L, 1.0) if (antibiotic and T_AB < 0) else np.zeros_like(phi)
    lap = np.empty_like(phi)
    dt = 0.15 * H * H / W.BETA
    if FRONT:
        dt = min(dt, 0.4 * H / (R_FRONT / (1 + K_MONOD)))
    per = int(np.ceil(DT_C / dt)); dt = DT_C / per
    hist = dict(t=[0.0], living_seed=[float(L[seed].mean())], share_a_seed=[], share_a_edges=[], living_edges=[],
                c_seed=[float(c[seed].mean())], a_seed=[0.0], colonised=[int(on.sum())])

    def shares(m):
        ph = G[m][:, :N_SP]
        return float((ph[:, 0] / ph.sum(1)).mean())
    hist["share_a_seed"].append(shares(seed)); hist["share_a_edges"].append(shares(edges) if example == "separate" else np.nan)
    hist["living_edges"].append(float(L[edges].mean()) if example == "separate" else np.nan)
    snaps = {}
    t0 = time.time()
    for k in range(int(round(T_END / DT_C))):
        T1 = (k + 1) * DT_C
        # (1) macro phi, alpha
        cf = c
        for _ in range(per):
            if FRONT:
                cf = nut.solve(phi, L)        # the front follows the nutrient within the outer step
            def f(ph, al):
                W.laplacian(ph, H, lap)
                d = W.BETA * lap + W.K_ALPHA * al - W.PENALTY * (np.maximum(ph - 1, 0) + np.minimum(ph, 0))
                if FRONT:
                    d = d + front(ph, cf)
                return d, W.K_ALPHA * ph
            k1 = f(phi, alpha); k2 = f(phi + 0.5 * dt * k1[0], alpha + 0.5 * dt * k1[1])
            k3 = f(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1]); k4 = f(phi + dt * k3[0], alpha + dt * k3[1])
            phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            alpha = alpha + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        # (2) nutrient and antibiotic at t_i, with L(t_{i-1})
        c_new = nut.solve(phi, L)
        a_new = ab.solve(phi, L, 1.0 if (antibiotic and T1 > T_AB + 1e-12) else 0.0)
        # (3) point models of the colonised cells
        idx = np.flatnonzero(on.ravel())
        Gf = G.reshape(-1, 12)
        Gf[idx] = np.asarray(STEP(jnp.asarray(Gf[idx]), jnp.asarray(c.ravel()[idx]), jnp.asarray(c_new.ravel()[idx]),
                                  jnp.asarray(a.ravel()[idx]), jnp.asarray(a_new.ravel()[idx])))
        G = Gf.reshape((NF,) * 3 + (12,))
        c, a = c_new, a_new
        # (4) newly colonised cells inherit the neighbours' state
        new = (phi >= PHI_MIN) & ~on
        if new.any():
            num, den = neighbours_mean(G, np.where(on, phi, 0.0))
            fill = new & (den > 0)
            gs = num[fill] / den[fill][:, None]
            gs[:, :N_SP] = np.clip(gs[:, :N_SP], 1e-10, None)
            gs[:, 5] = 1.0 - gs[:, :N_SP].sum(1)               # the point model's own constraint
            G[fill] = gs
            on = on | fill
        # (5) living biomass for the next step
        L = np.where(on, living(G.reshape(-1, 12)).reshape((NF,) * 3), 0.0)
        hist.setdefault("phi_max", []).append(float(phi.max())); hist.setdefault("phi_min", []).append(float(phi.min()))
        hist["t"].append(T1); hist["living_seed"].append(float(L[seed].mean()))
        hist["share_a_seed"].append(shares(seed))
        hist["share_a_edges"].append(shares(edges) if example == "separate" else np.nan)
        hist["living_edges"].append(float(L[edges].mean()) if example == "separate" else np.nan)
        hist["c_seed"].append(float(c[seed].mean())); hist["a_seed"].append(float(a[seed].mean()))
        hist["colonised"].append(int(on.sum()))
        if abs(T1 - T_AB) < 1e-9 or abs(T1 - T_END) < 1e-9:
            ph = G[..., :N_SP]
            shA = np.where(on, ph[..., 0] / np.maximum(ph.sum(-1), 1e-30), np.nan)
            snaps[round(T1, 3)] = dict(phi=phi.copy(), c=c.copy(), a=a.copy(), L=L.copy(), shareA=shA)
    print(f"[{example}, antibiotic {antibiotic}] {time.time() - t0:.0f} s; colonised {int(on.sum())} cells; "
          f"seed: living {hist['living_seed'][-1]:.3f}, share A {hist['share_a_seed'][-1]:.3f}, "
          f"c {hist['c_seed'][-1]:.3f}", flush=True)
    hist["seed_state_end"] = G[seed].mean(0).tolist()
    hist["seed_state_spread"] = float(np.abs(G[seed] - G[seed].mean(0)).max())
    seed16 = S.seed_mask(NM, cells8)
    vm, p = S.solve_fields(NM, W.block_mean(phi, NM), W.block_mean(alpha - 1, NM))
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
        ref = K.simulate(CASE)
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

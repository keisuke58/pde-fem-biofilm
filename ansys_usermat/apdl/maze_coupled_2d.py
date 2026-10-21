#!/usr/bin/env python3
"""maze_coupled_2d.py -- the maze of Klempt et al. 2024 (Sect. 4.4, Fig. 11, Table 6) with the agreed coupling.

Geometry (an assumption: the paper's maze is not given as data): a circular maze of the paper's size, 560 um across,
walls 10 um wide, generated as a perfect maze on polar cells with a fixed random seed (one path between any two
cells, dead ends as in Fig. 11). Nutrient c = 1 held in a central disk, the colony (phi = 1) at the entrance in the
outer wall at the top. Semi-2D in the paper (one element of 4 um): here 2D, cell size H.

Equations: the coupled model of coupled_living_consumption.py (Klempt 2024 Eqs. 34-36 with the front term, Eq. 35
with g phi_3D (sum_i phibar_i) c, the point model of Klempt et al. 2026 used as published at every colonised cell,
c* = c*_0 c and alpha* = alpha*_0 a linear in time over the inner steps), zero flux through all walls.
Constants: Table 6 of Klempt 2024 for the maze (d = 1e13 um^2/T*, beta = 10 um^2/T*, k_alpha = 1e-6 1/T*,
r = 200 um/T*), all others Table 2 (g = 1e8 1/T*, k = 1); point model Fritsch case I (c*_0 = 100, s = 0.1),
alpha*_0 = 100, phi_min = 0.001. With d = 1e13 the fields adjust within ~1e-8 T*: quasi-static.
Antibiotic (--ab): a = 1 held at the entrance from T* = T_END/3 on, same diffusivity and uptake as the nutrient
(as the cube examples: Klempt 2026 case 4 inserts it after one third).

    python ansys_usermat/apdl/maze_coupled_2d.py [--h 2.5] [--tend 1] [--ab]
        -> figdata/maze_<tag>.npz, results/2026-10-maze/maze_<tag>.json
"""
import json
import math
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (HERE, ROOT / "ansys_usermat/coupling"):
    sys.path.insert(0, str(p))
import ecology_native_ctypes as FN  # noqa: E402


def _opt(name, default, cast=float):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


H = _opt("--h", 2.5)                  # cell size [um]
T_END = _opt("--tend", 1.0)
AB = "--ab" in sys.argv
SEED = _opt("--seed", 1, int)
# Klempt 2024 Table 6 (maze), the rest Table 2
D_NUT, BETA, K_ALPHA, R_FRONT, K_MONOD, G_CONS = 1e13, 10.0, 1e-6, 200.0, 1.0, 1e8
PENALTY = 5.0                         # as wp3_reference.PENALTY (the cube runs)
PHI_MIN, DT_C, S_LINK, N_SUB = 0.001, 0.025, 0.1, None
ALPHA0 = 100.0
T_AB = T_END / 3.0
START = (0.25, 0.30)                  # Fritsch case I, both species in the colony
R_OUT, WALL, R_HUB, R_NUT = 280.0, 10.0, 30.0, 20.0
TAG = f"h{H:g}_tend{T_END:g}" + ("_ab" if AB else "") + ("" if SEED == 1 else f"_seed{SEED}")


# ---------------------------------------------------------------- maze
def make_maze(seed=SEED, period=40.0):
    """polar perfect maze: rings of width `period` (wall WALL + path), sectors doubled when the arc gets long"""
    rng = np.random.default_rng(seed)
    K = int((R_OUT - WALL - R_HUB) // period)
    nsec = []
    n = 6
    for k in range(K):
        rm = R_HUB + k * period + WALL + (period - WALL) / 2
        if 2 * math.pi * rm / n > 110.0:
            n *= 2
        nsec.append(n)
    nodes = [("hub",)] + [(k, j) for k in range(K) for j in range(nsec[k])]
    idx = {nd: i for i, nd in enumerate(nodes)}

    def parent(k, j):
        return ("hub",) if k == 0 else (k - 1, j * nsec[k - 1] // nsec[k])
    adj = {nd: [] for nd in nodes}
    for k in range(K):
        for j in range(nsec[k]):
            a, b = (k, j), (k, (j + 1) % nsec[k])
            adj[a].append(b); adj[b].append(a)
            p = parent(k, j)
            adj[(k, j)].append(p); adj[p].append((k, j))
    start = (K - 1, 0)                         # entrance cell at the top (sector 0 is centred on +y)
    seen, stack, links = {start}, [start], set()
    while stack:
        cur = stack[-1]
        nb = [x for x in adj[cur] if x not in seen]
        if not nb:
            stack.pop(); continue
        nx = nb[rng.integers(len(nb))]
        links.add(frozenset((cur, nx))); seen.add(nx); stack.append(nx)
    return K, nsec, links, period


def rasterise(K, nsec, links, period):
    """open (path) mask, nutrient disk, entrance colony and antibiotic source on the cell grid"""
    n = int(round(2 * (R_OUT + 5 * H) / H))
    xc = -n * H / 2 + (np.arange(n) + 0.5) * H
    X, Y = np.meshgrid(xc, xc, indexing="ij")
    R = np.hypot(X, Y)
    TH = np.mod(np.arctan2(X, Y), 2 * math.pi)            # 0 at +y, clockwise
    open_ = R < R_HUB
    for k in range(K):
        r0 = R_HUB + k * period
        ring = (R >= r0 + WALL) & (R < r0 + period)
        n_k = nsec[k]
        dth = 2 * math.pi / n_k
        sec = np.floor(np.mod(TH + dth / 2, 2 * math.pi) / dth).astype(int) % n_k   # sector j centred on j*dth
        # radial walls on the sector borders unless linked
        rel = np.mod(TH + dth / 2, 2 * math.pi) - sec * dth                          # angle from the border
        near = np.minimum(rel, dth - rel) * R < WALL / 2
        bord = np.where(rel < dth / 2, sec, (sec + 1) % n_k)                         # border index b between b-1, b
        blocked = np.zeros_like(ring)
        for b in range(n_k):
            if frozenset(((k, (b - 1) % n_k), (k, b))) not in links:
                blocked |= near & (bord == b)
        open_ |= ring & ~blocked
        # openings in the circular wall below the ring towards the parent cell
        cw = (R >= r0) & (R < r0 + WALL)
        for j in range(n_k):
            p = ("hub",) if k == 0 else (k - 1, j * nsec[k - 1] // n_k)
            if frozenset(((k, j), p)) in links:
                thj = j * dth
                dang = np.abs(np.mod(TH - thj + math.pi, 2 * math.pi) - math.pi)
                open_ |= cw & (dang * R < 10.0)        # 20 um openings
    r_last = R_HUB + K * period
    gate = (R >= r_last) & (np.abs(X) < (period - WALL) / 2) & (Y > 0) & (R < R_OUT + 3 * H)
    open_ |= gate
    nut = R < R_NUT
    colony = gate & (R > R_OUT - 2 * H)
    ab_src = gate & (R > R_OUT + H)
    return X, Y, open_, nut, colony, ab_src


# ---------------------------------------------------------------- operators on the open cells
class Ops:
    def __init__(self, open_):
        self.open = open_
        self.NB = {}
        for a in range(2):
            for s in (1, -1):
                m = np.roll(open_, -s, a) & open_
                e = [slice(None)] * 2; e[a] = -1 if s == 1 else 0; m[tuple(e)] = False
                self.NB[a, s] = m
        n2 = open_.size
        idx = np.arange(n2).reshape(open_.shape)
        rows, cols = [], []
        for (a, s), m in self.NB.items():
            rows.append(idx[m]); cols.append(np.roll(idx, -s, a)[m])
        rows, cols = np.concatenate(rows), np.concatenate(cols)
        deg = np.bincount(rows, minlength=n2).astype(float)
        self.L = sp.csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n2, n2)) - sp.diags(deg)  # H^2 lap
        self.n2 = n2

    def lap(self, u):
        out = np.zeros_like(u)
        for (a, s), m in self.NB.items():
            out += np.where(m, np.roll(u, -s, a) - u, 0.0)
        return out / H ** 2

    def grad(self, u):
        out = []
        for a in range(2):
            mp, mm = self.NB[a, 1], self.NB[a, -1]
            fw = np.where(mp, (np.roll(u, -1, a) - u) / H, 0.0); bw = np.where(mm, (u - np.roll(u, 1, a)) / H, 0.0)
            out.append((fw + bw) / np.maximum(mp.astype(float) + mm, 1.0))
        return out

    def upwind(self, u, v):
        out = np.zeros_like(u)
        for a in range(2):
            fw = np.where(self.NB[a, 1], (np.roll(u, -1, a) - u) / H, 0.0)
            bw = np.where(self.NB[a, -1], (u - np.roll(u, 1, a)) / H, 0.0)
            out += np.maximum(v[a], 0) * bw + np.minimum(v[a], 0) * fw
        return out

    def solve(self, phi, Lv, fixed, value):
        """0 = d lap u - g phi L u on the open cells, u = value on `fixed`, u = 0 in the walls (quasi-static)"""
        if value == 0.0:
            return np.zeros(self.open.shape)
        q = (G_CONS / D_NUT) * H ** 2 * (np.clip(phi, 0, 1) * Lv)
        A = (-self.L + sp.diags(q.ravel())).tolil() if False else (-self.L + sp.diags(q.ravel()))
        keep = (self.open & ~fixed).ravel().astype(float)
        A = (sp.diags(keep) @ A + sp.diags(1.0 - keep)).tocsc()
        b = np.where(fixed.ravel(), 1.0, 0.0)
        return value * spla.spsolve(A, b).reshape(self.open.shape)


def start_state(*phis):
    g = np.zeros(12)
    g[:len(phis)] = phis
    g[5] = 1.0 - sum(phis)
    g[6:11] = 0.999
    return g


def living(G):
    return (G[:, :2] * G[:, 6:8]).sum(1)


def main():
    eco = ROOT / "abaqus_composition" / "eco_fritsch_2sp_a100.txt"
    cst = FN.init(eco)
    assert cst["nact"] == 2 and cst["c_star"] == 100.0 and cst["alpha_star"] == ALPHA0
    n_sub = max(1, math.ceil(S_LINK * DT_C / 1e-4 - 1e-9))
    K, nsec, links, period = make_maze()
    X, Y, open_, nut, colony, ab_src = rasterise(K, nsec, links, period)
    O = Ops(open_)
    shp = open_.shape
    phi = np.where(colony, 1.0, 0.0)
    alpha = np.ones(shp)
    G = np.zeros(shp + (12,))
    G[colony] = start_state(*START)
    on = phi >= PHI_MIN
    Lv = np.where(on, living(G.reshape(-1, 12)).reshape(shp), 0.0)
    c = O.solve(phi, Lv, nut, 1.0)
    a = np.zeros(shp)
    dt = min(0.15 * H * H / BETA, 0.4 * H / (R_FRONT / (1 + K_MONOD)))
    per = int(np.ceil(DT_C / dt)); dt = DT_C / per
    eps2 = (1e-9) ** 2                          # regularisation of |grad c| [1/um]^2

    def front(ph, cf):
        g = O.grad(cf)
        nrm = np.sqrt(g[0] ** 2 + g[1] ** 2 + eps2)
        sp_ = R_FRONT * cf / (K_MONOD + cf)
        return -O.upwind(ph, [sp_ * x / nrm for x in g])

    def f(ph, al, cf):
        d = BETA * O.lap(ph) + K_ALPHA * al - PENALTY * (np.maximum(ph - 1, 0) + np.minimum(ph, 0)) + front(ph, cf)
        return np.where(open_, d, 0.0), np.where(open_, K_ALPHA * ph, 0.0)
    snaps, hist = {}, dict(t=[0.0], colonised=[int(on.sum())], reach_um=[], living_mean=[], c_mean=[])
    path_r = np.hypot(X, Y)
    nsteps = int(round(T_END / DT_C))
    every = max(1, nsteps // 10)
    t0 = time.time()
    print(f"maze: {K} rings, sectors {nsec}, grid {shp}, open cells {int(open_.sum())}, dt {dt:.2e} x {per}",
          flush=True)
    for k in range(nsteps):
        T1 = (k + 1) * DT_C
        cf = c
        for _ in range(per):
            cf = O.solve(phi, Lv, nut, 1.0)
            k1 = f(phi, alpha, cf); k2 = f(phi + 0.5 * dt * k1[0], alpha + 0.5 * dt * k1[1], cf)
            k3 = f(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1], cf); k4 = f(phi + dt * k3[0], alpha + dt * k3[1], cf)
            phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            alpha = alpha + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        c_new = O.solve(phi, Lv, nut, 1.0)
        a_new = O.solve(phi, Lv, ab_src, 1.0 if (AB and T1 > T_AB + 1e-12) else 0.0)
        idx = np.flatnonzero(on.ravel())
        Gf = G.reshape(-1, 12)
        if idx.size:
            Gf[idx], okf = FN.step_lin(Gf[idx], S_LINK * DT_C, n_sub, c.ravel()[idx], c_new.ravel()[idx],
                                       a.ravel()[idx], a_new.ravel()[idx])
            assert okf.all(), "point model failed"
        G = Gf.reshape(shp + (12,))
        c, a = c_new, a_new
        new = (phi >= PHI_MIN) & ~on & open_
        if new.any():                          # neighbours' phi-weighted mean state (as the cube examples)
            w = np.where(on, phi, 0.0)
            num = np.zeros_like(G); den = np.zeros(shp)
            for (ax, s), m in O.NB.items():
                wn = np.where(m, np.roll(w, -s, ax), 0.0)
                num += wn[..., None] * np.roll(G, -s, ax); den += wn
            fill = new & (den > 0)
            gs = num[fill] / den[fill][:, None]
            gs[:, :2] = np.clip(gs[:, :2], 1e-10, None)
            gs[:, 5] = 1.0 - gs[:, :2].sum(1)
            G[fill] = gs
            on = on | fill
        Lv = np.where(on, living(G.reshape(-1, 12)).reshape(shp), 0.0)
        hist["t"].append(T1); hist["colonised"].append(int(on.sum()))
        hist["reach_um"].append(float(path_r[phi >= 0.5].min()) if (phi >= 0.5).any() else None)
        hist["living_mean"].append(float(Lv[on].mean()) if on.any() else 0.0)
        hist["c_mean"].append(float(c[open_].mean()))
        if (k + 1) % every == 0 or k == nsteps - 1:
            ph = G[..., :2]
            snaps[round(T1, 4)] = dict(phi=phi.copy(), c=c.copy(), a=a.copy(), L=Lv.copy(),
                                       shareA=np.where(on, ph[..., 0] / np.maximum(ph.sum(-1), 1e-30), np.nan))
            print(f"  T* {T1:.3f}: colonised {int(on.sum())}, closest phi>=0.5 to the centre "
                  f"{hist['reach_um'][-1]} um, mean c {hist['c_mean'][-1]:.3f}, {time.time() - t0:.0f} s", flush=True)
    (ROOT / "figdata").mkdir(exist_ok=True)
    np.savez_compressed(ROOT / "figdata" / f"maze_{TAG}.npz", open=open_, nut=nut, colony=colony, ab_src=ab_src,
                        X=X, Y=Y, **{f"{k}_{T}": v for T, d in snaps.items() for k, v in d.items()})
    res_dir = HERE / "results" / "2026-10-maze"
    res_dir.mkdir(exist_ok=True)
    (res_dir / f"maze_{TAG}.json").write_text(json.dumps(dict(
        setting=dict(H_um=H, T_end=T_END, antibiotic=AB, T_ab=T_AB if AB else None, maze_seed=SEED, rings=K,
                     sectors=nsec, d=D_NUT, beta=BETA, k_alpha=K_ALPHA, r=R_FRONT, g=G_CONS, k=K_MONOD,
                     c_star0=100.0, alpha_star0=ALPHA0, s=S_LINK, phi_min=PHI_MIN, start=START),
        hist=hist), indent=1) + "\n")
    print(f"done in {time.time() - t0:.0f} s -> figdata/maze_{TAG}.npz")


if __name__ == "__main__":
    main()

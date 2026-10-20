#!/usr/bin/env python3
"""Mesh convergence of the biofilm stress on the real molar (anim_real_tooth_3d.py).

The demonstration on the molar FDI 46 (results/2026-10-3d/REAL_TOOTH.md) computes its growth field on 96^3 cells and
its stress on 32^3 hex8 elements, where the 0.5 mm film is one element thick. This script repeats the final state on
coarser and finer grids and separates the numbers that converge with the mesh from those that do not, as
wall_stress_convergence.py does for the 2 mm cube.

Field (same model and geometry as the demonstration: d phi/dt = beta lap phi + k_alpha alpha - penalty, d alpha/dt =
k_alpha phi, beta = 0.08 mm^2/T*, k_alpha = 1e-3 /T*, penalty 5, T* = 1.1, RK4, zero flux at the walls and at the
voxelised tooth; seed phi = 1 within 0.5 mm of the tooth on the -x flank, z = -3 .. -1 mm): final state only, on
48^3, 96^3, 144^3 and 192^3 cells (h = 0.333, 0.167, 0.111, 0.083 mm), cached in figdata/real_tooth_conv_field_<n>.npz.
The tooth and the seed are re-voxelised on every grid (cell centre inside the STL / within the seed region), so the
field convergence includes the geometry. dt = min(0.15 h^2/beta, T/66), the time step of the demonstration on 96^3
(a run with dt/2 on 96^3 checks that the time error is below the mesh differences).

Stress (mesh_study_seed.element_matrices, hex8 B-bar, small strain, E = 10 Pa (phi^2 + 1e-3) in film and space,
E = 1e4 Pa in the tooth without growth, nu = 0.49, eigenstrain (alpha - 1) I, bottom face clamped) on 16^3, 24^3,
32^3 and 48^3 elements (64^3 optional), each with the fields block-averaged from the 96^3 grid and from the 192^3
grid; both divide every stress mesh, so between two stress meshes only the mechanics changes, and between the two
field grids only the field. pyamg smoothed aggregation + CG, tol 1e-10, residual and time in the JSON.

Element sets in physical space, the same on every mesh: they are defined on the field grid (film = cells with
phi > 0.3 that are not tooth; cervical band = film with z in [-2.1, -1.1] mm, 0.5 mm either side of the CEJ at
z = -1.6; core = the rest of the film) and every stress element is weighted by the fraction of its field cells in
the set, so the means are volume means over the same region on every mesh (a 1 mm element holding half a
millimetre of film counts with weight 0.5; with the block-averaged phi > 0.3 instead, the 16^3 mesh would hold
13 film elements and no band). Maximum, location and counts use the majority-film elements (fraction > 0.5), the
percentiles are volume-weighted over every element with film in it; first layer = majority-film elements with a
tooth element among their six neighbours, outer = the rest (a layer one element thick, so mesh-dependent by
construction). Reported per mesh: film mean von Mises, band mean, core mean, band/core, first-layer and outer
means, maximum with its element centre, 95th and 99th percentiles, maximum and mean pressure p = -tr(sigma)/3.
Fits as in wall_stress_convergence.py: v = v0 + c h^p for the means (observed order and h -> 0 value), v = a h^-s
for the maxima and percentiles; Richardson order from the 48/96/192 field triple.

    python ansys_usermat/apdl/real_tooth_convergence.py --field-only            # the four field grids (about 5 min)
    python ansys_usermat/apdl/real_tooth_convergence.py                         # fields + stress meshes (about 15 min)
    python ansys_usermat/apdl/real_tooth_convergence.py --meshes 64 --stress-fields 192   # optional 64^3 (slow)
    python ansys_usermat/apdl/real_tooth_convergence.py --figure
        -> ansys_usermat/apdl/results/2026-10-3d/real_tooth_convergence.json, assets/fig_real_tooth_convergence.pdf/.png
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")        # two agents share four cores; threaded BLAS only slows the CG down

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import mesh_study_seed as S  # noqa: E402

# the demonstration's settings (anim_real_tooth_3d.py)
L = 16.0
BETA0, KA, PEN, TEND = 0.02, 1e-3, 5.0, 1.1
SEED_THICK = 0.5
BETA = BETA0 * (SEED_THICK / 0.25) ** 2        # 0.08 mm^2/T*
E_BIO, E_TOOTH, FLOOR, NU = 10.0, 1e4, 1e-3, 0.49
Z_TOP = 5.5
SEED_Z, SEED_HALF_ANGLE = (-3.0, -1.0), 60.0
NFRAMES = 66                                   # dt of the demonstration = TEND / NFRAMES
STL = ROOT / "assets" / "geometry" / "mendeley_xjsx7nfhj8_D46.stl"
CACHE_DEMO = ROOT / "figdata" / "anim_real_tooth_3d.npz"
OUT_JSON = HERE / "results" / "2026-10-3d" / "real_tooth_convergence.json"
FIG = ROOT / "assets" / "fig_real_tooth_convergence"
Z_CEJ, BAND_HALF = -1.6, 0.5
FILM_PHI = 0.3
FIELD_GRIDS = (48, 96, 144, 192)
STRESS_MESHES = (16, 24, 32, 48)
STRESS_FIELDS = (96, 192)


def down(a, n):
    f = a.shape[0] // n
    assert f * n == a.shape[0], (a.shape, n)
    return a.reshape(n, f, n, f, n, f).mean((1, 3, 5))


# ---------------------------------------------------------------- geometry (copied from anim_real_tooth_3d.py)
def load_tooth():
    """The molar, centred on its crown in x and y, crown top at Z_TOP."""
    import pyvista as pv
    m = pv.read(STL).triangulate().clean()
    P = np.asarray(m.points)
    top = P[:, 2].max()
    crown = P[P[:, 2] > top - 7.0]
    shift = np.array([crown[:, 0].mean(), crown[:, 1].mean(), top - Z_TOP])
    m.points = P - shift
    return m


def geometry(nf):
    """Signed distance, tooth mask, open mask and seed mask on the nf^3 cell grid."""
    import pyvista as pv
    h = L / nf
    ax = -L / 2 + (np.arange(nf) + 0.5) * h
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    cen = pv.ImageData(dimensions=(nf,) * 3, spacing=(h,) * 3, origin=(ax[0],) * 3)
    sdf = np.asarray(cen.compute_implicit_distance(load_tooth())["implicit_distance"]).reshape((nf,) * 3, order="F")
    tooth = sdf < 0
    ang = np.degrees(np.arctan2(Y, -X))
    seed = ~tooth & (sdf < SEED_THICK) & (Z > SEED_Z[0]) & (Z < SEED_Z[1]) & (np.abs(ang) < SEED_HALF_ANGLE)
    return sdf.astype(np.float32), tooth, ~tooth, seed


def make_lap(OPEN, h):
    """7-point Laplacian with zero flux at the walls and at the tooth (masks precomputed)."""
    oks = []
    for axn in range(3):
        for sgn in (1, -1):
            on = np.roll(OPEN, -sgn, axn)
            edge = [slice(None)] * 3
            edge[axn] = -1 if sgn == 1 else 0
            ok = OPEN & on
            ok[tuple(edge)] = False
            oks.append((axn, sgn, ok))

    def lap(u):
        out = np.zeros_like(u)
        for axn, sgn, ok in oks:
            out += np.where(ok, np.roll(u, -sgn, axn) - u, 0.0)
        return out / h ** 2
    return lap


def field(nf, dt_factor=1.0, force=False):
    """Final state (phi, alpha - 1) on the nf^3 grid at T* (cached) and the field measures."""
    tag = f"{nf}" if dt_factor == 1.0 else f"{nf}_dt{dt_factor:g}"
    cache = ROOT / "figdata" / f"real_tooth_conv_field_{tag}.npz"
    if cache.exists() and not force:
        d = np.load(cache)
        return d["phi"], d["am1"], d["tooth"], d["seed"], json.loads(str(d["info"]))
    h = L / nf
    t0 = time.time()
    sdf, tooth, OPEN, seed = geometry(nf)
    t_sdf = time.time() - t0
    dt = min(0.15 * h * h / BETA, TEND / NFRAMES) * dt_factor
    nstep = int(np.ceil(TEND / dt))
    dt = TEND / nstep
    print(f"field {nf}^3: h = {h:.4f} mm, tooth {tooth.sum()} cells = {tooth.sum() * h ** 3:.1f} mm^3, seed "
          f"{seed.sum()} cells = {seed.sum() * h ** 3:.2f} mm^3, sdf {t_sdf:.0f} s; dt = {dt:.4f}, {nstep} steps",
          flush=True)
    lap = make_lap(OPEN, h)
    phi = np.where(seed, 1.0, 0.0)
    alpha = np.ones_like(phi)

    def rhs(phi, alpha):
        dphi = BETA * lap(phi) + KA * alpha - PEN * (np.maximum(phi - 1, 0) + np.minimum(phi, 0))
        return np.where(OPEN, dphi, 0.0), np.where(OPEN, KA * phi, 0.0)

    t0 = time.time()
    for s in range(nstep):
        k1 = rhs(phi, alpha)
        k2 = rhs(phi + 0.5 * dt * k1[0], alpha + 0.5 * dt * k1[1])
        k3 = rhs(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1])
        k4 = rhs(phi + dt * k3[0], alpha + dt * k3[1])
        phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        alpha = alpha + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if (s + 1) % max(1, nstep // 4) == 0:
            print(f"  step {s + 1}/{nstep}  phi volume {phi.sum() * h ** 3:.3f} mm^3  {time.time() - t0:.0f} s",
                  flush=True)
    am1 = np.where(OPEN, alpha - 1.0, 0.0)
    info = {"h": h, "cells": nf ** 3, "dt": dt, "steps": nstep, "dt_factor": dt_factor,
            "tooth cells": int(tooth.sum()), "tooth volume mm^3": float(tooth.sum() * h ** 3),
            "seed cells": int(seed.sum()), "seed volume mm^3": float(seed.sum() * h ** 3),
            "seconds sdf": t_sdf, "seconds field": time.time() - t0}
    info.update(field_measures(phi, am1, OPEN, seed, h))
    print("  " + ", ".join(f"{k} {v:.4g}" if isinstance(v, float) else f"{k} {v}" for k, v in info.items()
                           if k not in ("h", "cells", "dt_factor")), flush=True)
    cache.parent.mkdir(exist_ok=True)
    np.savez_compressed(cache, phi=phi.astype(np.float32), am1=am1.astype(np.float32), tooth=tooth, seed=seed,
                        sdf=sdf, info=json.dumps(info))
    return phi.astype(np.float32), am1.astype(np.float32), tooth, seed, info


def field_measures(phi, am1, OPEN, seed, h):
    return {"film volume mm^3": float(phi[OPEN].sum() * h ** 3),
            "volume phi>0.5 mm^3": float((phi > 0.5).sum() * h ** 3),
            "volume phi>0.3 mm^3": float((phi > FILM_PHI).sum() * h ** 3),
            "seed phi mean": float(phi[seed].mean()),
            "seed am1 mean": float(am1[seed].mean()), "am1 max": float(am1[OPEN].max()),
            "am1 mean phi>0.3": float(am1[phi > FILM_PHI].mean())}


# ---------------------------------------------------------------- mechanics (Mech of anim_real_tooth_3d.py, tol 1e-10)
def mechanics(n, phi_m, am1_m, tooth_m, tol=1e-10):
    import pyamg
    h = L / n
    Ke0, fe0, B0, D, m = S.element_matrices(h, NU, False)
    idx = np.arange((n + 1) ** 3, dtype=np.int64).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1).astype(np.int32)
    ndof = 3 * (n + 1) ** 3
    Ee = np.where(tooth_m.ravel(), E_TOOTH, E_BIO * (np.clip(phi_m.ravel(), 0, 1) ** 2 + FLOOR))
    eps0 = np.where(tooth_m.ravel(), 0.0, am1_m.ravel())
    F = np.zeros(ndof)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    fix = (3 * idx[:, :, 0].ravel()[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(ndof), fix)
    fmap = np.full(ndof, -1, np.int32)                       # dof -> index among the free dofs
    fmap[free] = np.arange(len(free), dtype=np.int32)
    # assembled in chunks of elements straight into the free-dof matrix: the 96^3 mesh (2.7M dof, 510M entries
    # before summing) does not fit as one COO triplet list
    Kff = None
    fd = fmap[dofs]
    for c0 in range(0, len(dofs), 65536):
        d = fd[c0:c0 + 65536]
        rows = np.repeat(d, 24, 1).ravel()
        cols = np.tile(d, (1, 24)).ravel()
        val = (Ke0[None] * Ee[c0:c0 + 65536, None, None]).ravel()
        keep = (rows >= 0) & (cols >= 0)
        Kc = sp.csr_matrix((val[keep], (rows[keep], cols[keep])), shape=(len(free),) * 2)
        Kff = Kc if Kff is None else Kff + Kc
        del rows, cols, val, keep, Kc
    Kff.sum_duplicates()
    xn = np.linspace(-L / 2, L / 2, n + 1)
    Xn, Yn, Zn = np.meshgrid(xn, xn, xn, indexing="ij")
    P = np.stack([Xn.ravel(), Yn.ravel(), Zn.ravel()], 1)
    B = np.zeros((ndof, 6))                                  # rigid body modes for the AMG near-nullspace
    B[0::3, 0] = 1; B[1::3, 1] = 1; B[2::3, 2] = 1
    B[0::3, 3] = -P[:, 1]; B[1::3, 3] = P[:, 0]
    B[1::3, 4] = -P[:, 2]; B[2::3, 4] = P[:, 1]
    B[0::3, 5] = P[:, 2]; B[2::3, 5] = -P[:, 0]
    t0 = time.time()
    ml = pyamg.smoothed_aggregation_solver(Kff, B=B[free], max_coarse=1000)
    t_setup = time.time() - t0
    res = []
    u = np.zeros(ndof)
    u[free] = ml.solve(F[free], tol=tol, accel="cg", maxiter=2000, residuals=res)
    rel = float(np.linalg.norm(Kff @ u[free] - F[free]) / np.linalg.norm(F[free]))
    info = {"solver": "pyamg smoothed aggregation + cg", "dof": int(len(free)), "iterations": len(res),
            "relative residual": rel, "seconds setup": t_setup, "seconds solve": time.time() - t0 - t_setup}
    print(f"    amg-cg {len(free)} dof: {len(res)} iterations, relative residual {rel:.1e}, "
          f"setup {t_setup:.0f} s, solve {info['seconds solve']:.0f} s", flush=True)
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    sig = Ee[:, None] * (eps @ D.T)
    p = -sig[:, :3].mean(1)
    s = sig[:, :3] + p[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm.reshape((n,) * 3), p.reshape((n,) * 3), info


def fine_sets(nf, phi_f, tooth_f):
    """Film, band and core on the field grid: the physical regions that every stress mesh is measured on."""
    h = L / nf
    z = -L / 2 + (np.arange(nf) + 0.5) * h
    Z = z[None, None, :]
    film = (phi_f > FILM_PHI) & ~tooth_f
    band = film & (Z > Z_CEJ - BAND_HALF) & (Z < Z_CEJ + BAND_HALF)
    return dict(film=film, band=band, core=film & ~band)


def element_sets(n, tooth_m, fine):
    """Weights per stress element = fraction of its field cells in the film / band / core (the same physical
    region on every mesh, so a 1 mm element that holds half a millimetre of film counts with weight 0.5);
    majority-film elements (weight > 0.5) for counts, maximum and location; first layer = majority-film
    elements with a tooth element among the six neighbours, outer = the other majority-film elements."""
    h = L / n
    ax = -L / 2 + (np.arange(n) + 0.5) * h
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    w = {k: down(v.astype(float), n) for k, v in fine.items()}
    for k in w:
        w[k][tooth_m] = 0.0
    film = w["film"] > 0.5
    touch = np.zeros_like(film)
    for axn in range(3):
        for sgn in (1, -1):
            nb = np.roll(tooth_m, sgn, axn)
            edge = [slice(None)] * 3
            edge[axn] = 0 if sgn == 1 else -1
            nb[tuple(edge)] = False
            touch |= nb
    return w, dict(film=film, first=film & touch, outer=film & ~touch), (X, Y, Z)


def wmean(v, w):
    return float((v * w).sum() / w.sum()) if w.sum() > 0 else None


def wpercentile(v, w, q):
    """Volume-weighted percentile: the value below which q % of the weight lies."""
    i = np.argsort(v)
    v, w = v[i], w[i]
    c = np.cumsum(w) - 0.5 * w
    return float(np.interp(q / 100 * w.sum(), c, v))


def measures(n, phi_m, tooth_m, vm, p, fine):
    w, sets, (X, Y, Z) = element_sets(n, tooth_m, fine)
    film, first, outer = sets["film"], sets["first"], sets["outer"]
    any_film = w["film"] > 0
    imax = np.unravel_index(np.argmax(np.where(film, vm, -np.inf)), vm.shape)
    ipmax = np.unravel_index(np.argmax(np.where(film, p, -np.inf)), p.shape)
    vmf = vm[film]
    h = L / n
    r = {"h": h, "elements": n ** 3,
         "film elements": int(film.sum()), "film elements any": int(any_film.sum()),
         "film elements block phi>0.3": int(((phi_m > FILM_PHI) & ~tooth_m).sum()),
         "film volume mm^3": float(w["film"].sum() * h ** 3), "band volume mm^3": float(w["band"].sum() * h ** 3),
         "core volume mm^3": float(w["core"].sum() * h ** 3),
         "first layer elements": int(first.sum()), "outer elements": int(outer.sum()),
         "film vM mean": wmean(vm, w["film"]),
         "band vM mean": wmean(vm, w["band"]), "core vM mean": wmean(vm, w["core"]),
         "first layer vM mean": float(vm[first].mean()) if first.any() else None,
         "outer vM mean": float(vm[outer].mean()) if outer.any() else None,
         "vM max": float(vmf.max()), "vM max at": [float(X[imax]), float(Y[imax]), float(Z[imax])],
         "vM max in band": bool(Z_CEJ - BAND_HALF < Z[imax] < Z_CEJ + BAND_HALF),
         "vM max in first layer": bool(first[imax]),
         "vM max any film": float(vm[any_film].max()),
         "vM p95": wpercentile(vm[any_film], w["film"][any_film], 95),
         "vM p99": wpercentile(vm[any_film], w["film"][any_film], 99),
         "vM median": wpercentile(vm[any_film], w["film"][any_film], 50),
         "p max": float(p[film].max()), "p max at": [float(X[ipmax]), float(Y[ipmax]), float(Z[ipmax])],
         "p mean": wmean(p, w["film"]), "p band mean": wmean(p, w["band"]),
         "vM max all": float(vm[~tooth_m].max())}
    r["band/core"] = r["band vM mean"] / r["core vM mean"] if r["core vM mean"] else None
    r["first/outer"] = r["first layer vM mean"] / r["outer vM mean"] if r["outer vM mean"] else None
    return r


# ---------------------------------------------------------------- fits (as wall_stress_convergence.py)
def fit_converging(hs, v):
    """Convergence of a mean: change per refinement, first- and second-order extrapolations (linear least squares
    on all meshes), and a free fit v = v0 + c h^p (order bounded to [0.3, 4]). The free fit is marked reliable only
    when the changes are monotone and the last one is above the solver's resolution (0.2 %); otherwise the first-
    order value is the h -> 0 estimate."""
    from scipy.optimize import least_squares
    hs, v = np.asarray(hs, float), np.asarray(v, float)
    change = [float((v[i + 1] / v[i] - 1) * 100) for i in range(len(v) - 1)]
    ext = {}
    for pp in (1, 2):
        A = np.stack([np.ones_like(hs), hs ** pp], 1)
        ext[pp] = float(np.linalg.lstsq(A, v, rcond=None)[0][0])
    f = lambda q: q[0] + q[1] * hs ** q[2] - v                      # noqa: E731
    best = None
    for p0 in (1.0, 2.0):
        q = least_squares(f, [v[-1], v[0] - v[-1], p0], bounds=([-np.inf, -np.inf, 0.3], [np.inf, np.inf, 4.0]))
        if best is None or q.cost < best.cost:
            best = q
    order = float(best.x[2])
    mono = all(c > 0 for c in change) or all(c < 0 for c in change)
    reliable = bool(mono and abs(change[-1]) > 0.2 and 0.31 < order < 3.99)
    h0 = float(best.x[0]) if reliable else ext[1]
    return {"change % per refinement": change, "h->0 p=1": ext[1], "h->0 p=2": ext[2],
            "free fit order": order, "free fit h->0": float(best.x[0]), "free fit reliable": reliable,
            "h->0": h0, "error % at h->0": [float((x / h0 - 1) * 100) for x in v],
            "spread % over meshes": float((v.max() / v.min() - 1) * 100)}


def fit_power(hs, v):
    """v = a h^-s (log-log least squares) and the slope of each successive pair."""
    hs, v = np.asarray(hs, float), np.asarray(v, float)
    s, loga = np.polyfit(np.log(1 / hs), np.log(v), 1)
    pair = [float(np.log(v[i + 1] / v[i]) / np.log(hs[i] / hs[i + 1])) for i in range(len(v) - 1)]
    return {"s": float(s), "a": float(np.exp(loga)), "pair slopes": pair,
            "change % per refinement": [float((v[i + 1] / v[i] - 1) * 100) for i in range(len(v) - 1)],
            "spread % over meshes": float((v.max() / v.min() - 1) * 100)}


def richardson(hs, v):
    """Observed order and extrapolated value from three grids with a constant refinement ratio."""
    hs, v = np.asarray(hs, float), np.asarray(v, float)
    r = hs[0] / hs[1]
    d1, d2 = v[1] - v[0], v[2] - v[1]
    if d1 == 0 or d2 == 0 or d1 / d2 <= 0:
        return {"order": None, "h->0": float(v[2]), "note": "differences not monotone"}
    p = float(np.log(d1 / d2) / np.log(r))
    v0 = float(v[2] + d2 / (r ** p - 1))
    return {"order": p, "h->0": v0, "ratio": float(r), "error % finest": float((v[2] / v0 - 1) * 100)}


# ---------------------------------------------------------------- driver
def load_json():
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text())
    return {"model": {"L mm": L, "beta mm^2/T*": BETA, "k_alpha 1/T*": KA, "penalty": PEN, "T*": TEND,
                      "E_bio Pa": E_BIO, "E_tooth Pa": E_TOOTH, "nu": NU, "floor": FLOOR,
                      "seed": {"thickness mm": SEED_THICK, "z mm": list(SEED_Z), "half angle deg": SEED_HALF_ANGLE},
                      "film": f"field cells with phi > {FILM_PHI}, not tooth; each stress element weighted by "
                              "its film volume fraction (means, percentiles); majority-film elements "
                              "(fraction > 0.5) for counts, maximum, location",
                      "band": f"film, z in [{Z_CEJ - BAND_HALF}, {Z_CEJ + BAND_HALF}] mm (CEJ at {Z_CEJ})",
                      "core": "film outside the band",
                      "first layer": "majority-film elements with a tooth element among the six neighbours",
                      "p": "-tr(sigma)/3, compression positive", "units": "Pa, mm"},
            "field": {}, "stress": {}}


def save_json(out):
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=1))


def run_fields(grids, dt_check=True):
    out = load_json()
    for nf in grids:
        _, _, _, _, info = field(nf)
        out["field"][str(nf)] = info
        save_json(out)
    if dt_check:
        _, _, _, _, info = field(96, dt_factor=0.5)
        out["field"]["96 dt/2"] = info
    if CACHE_DEMO.exists() and "96" in out["field"]:
        d = np.load(CACHE_DEMO)
        phi_c, am1_c = d["phi"][-1].astype(float), d["am1"][-1].astype(float)
        phi, am1, tooth, seed, _ = field(96)
        out["field"]["96 demo cache (float16)"] = dict(
            field_measures(phi_c, am1_c, ~tooth, seed, L / 96),
            **{"max |phi - phi_demo|": float(np.abs(phi - phi_c).max()),
               "max |am1 - am1_demo|": float(np.abs(am1 - am1_c).max())})
    out["field fits"] = field_fits(out["field"])
    save_json(out)
    print("wrote", OUT_JSON, flush=True)


def field_fits(F):
    ns = sorted(int(k) for k in F if k.isdigit())
    hs = [L / n for n in ns]
    fits = {}
    keys = ["film volume mm^3", "seed am1 mean", "am1 max", "volume phi>0.5 mm^3", "seed volume mm^3",
            "tooth volume mm^3", "seed phi mean"]
    for k in keys:
        v = [F[str(n)][k] for n in ns]
        fits[k] = {"grids": ns, "values": v}
        if len(ns) >= 3:
            fits[k]["fit"] = fit_converging(hs, v)
        tri = [n for n in (48, 96, 192) if n in ns]
        if len(tri) == 3:
            fits[k]["richardson 48/96/192"] = richardson([L / n for n in tri], [F[str(n)][k] for n in tri])
    if "96 dt/2" in F:
        fits["dt check 96"] = {k: {"dt": F["96"][k], "dt/2": F["96 dt/2"][k],
                                   "change %": float((F["96 dt/2"][k] / F["96"][k] - 1) * 100)}
                               for k in keys if F["96"][k]}
    return fits


def run_stress(meshes, field_grids, force=False):
    """One solve per (field grid, stress mesh); each result goes to the JSON and a small npz as soon as it is
    finished, and a pair already in the JSON is skipped (unless force), so a restart loses at most one solve."""
    for nf in field_grids:
        phi_f, am1_f, tooth_f, seed_f, _ = field(nf)
        fine = fine_sets(nf, phi_f, tooth_f)
        for n in meshes:
            if nf % n:
                print(f"skip stress {n}^3 on field {nf}^3 ({nf} not a multiple of {n})")
                continue
            if not force and str(n) in load_json()["stress"].get(str(nf), {}):
                print(f"stress {n}^3 on field {nf}^3 already in {OUT_JSON.name}, skipped", flush=True)
                continue
            t0 = time.time()
            phi_m, am1_m = down(phi_f.astype(float), n), down(am1_f.astype(float), n)
            tooth_m = down(tooth_f.astype(float), n) > 0.5
            print(f"stress {n}^3 on field {nf}^3: tooth {tooth_m.sum()} elements, film "
                  f"{(down(fine['film'].astype(float), n) > 0.5).sum()} majority elements", flush=True)
            vm, p, info = mechanics(n, phi_m, am1_m, tooth_m)
            r = measures(n, phi_m, tooth_m, vm, p, fine)
            r["tooth elements"] = int(tooth_m.sum())
            r["solver"] = info
            r["seconds"] = time.time() - t0
            fmt = lambda x: "n/a" if x is None else f"{x:.4g}"                                  # noqa: E731
            print(f"  film {fmt(r['film vM mean'])} band {fmt(r['band vM mean'])} core {fmt(r['core vM mean'])} "
                  f"ratio {fmt(r['band/core'])} first {fmt(r['first layer vM mean'])} outer "
                  f"{fmt(r['outer vM mean'])} max {r['vM max']:.4g} at {r['vM max at']} p95 {r['vM p95']:.4g} "
                  f"p99 {r['vM p99']:.4g} pmax {r['p max']:.4g} pmean {fmt(r['p mean'])} ({r['seconds']:.0f} s)",
                  flush=True)
            np.savez_compressed(ROOT / "figdata" / f"real_tooth_conv_stress_f{nf}_m{n}.npz",
                                vm=vm.astype(np.float32), p=p.astype(np.float32), phi=phi_m.astype(np.float32),
                                tooth=tooth_m)
            out = load_json()
            out["stress"].setdefault(str(nf), {})[str(n)] = r
            out["stress fits"] = stress_fits(out["stress"])
            save_json(out)
    print("wrote", OUT_JSON, flush=True)


CONV_KEYS = ["film vM mean", "band vM mean", "core vM mean", "band/core", "first layer vM mean",
             "outer vM mean", "first/outer", "p mean", "p band mean", "vM median"]
POWER_KEYS = ["vM max", "vM p99", "vM p95", "p max"]


def stress_fits(Sd):
    fits = {}
    for nf, rows in Sd.items():
        ns = sorted(int(k) for k in rows)
        if len(ns) < 3:
            continue
        hs = [L / n for n in ns]
        f = {"meshes": ns}
        for k in CONV_KEYS:
            v = [rows[str(n)][k] for n in ns]
            if any(x is None for x in v):
                continue
            f[k] = fit_converging(hs, v)
        for k in POWER_KEYS:
            f[k] = fit_power(hs, [rows[str(n)][k] for n in ns])
        fits[nf] = f
    # field influence at fixed stress mesh
    if len(Sd) >= 2:
        a, b = sorted(Sd, key=int)[0], sorted(Sd, key=int)[-1]
        common = sorted(set(Sd[a]) & set(Sd[b]), key=int)
        infl = {"fields": [a, b]}
        for k in CONV_KEYS + POWER_KEYS:
            infl[k] = {n: float((Sd[b][n][k] / Sd[a][n][k] - 1) * 100) for n in common
                       if Sd[a][n][k] not in (None, 0) and Sd[b][n][k] is not None}
        fits["field influence"] = infl
    return fits


def refit():
    out = load_json()
    out["field fits"] = field_fits(out["field"])
    out["stress fits"] = stress_fits(out["stress"])
    save_json(out)
    print_summary(out)


def print_summary(out):
    for nf, f in out.get("stress fits", {}).items():
        if nf == "field influence":
            continue
        print(f"== stress on field {nf}^3, meshes {f['meshes']}")
        for k in CONV_KEYS:
            if k not in f:
                continue
            g = f[k]
            tag = f"order {g['free fit order']:.2f}" if g["free fit reliable"] else "order n/a (p = 1 used)"
            print(f"  {k:20s} change {' / '.join(f'{e:+.2f}' for e in g['change % per refinement'])} %  {tag}"
                  f"  h->0 {g['h->0']:.4e}  error % {' / '.join(f'{e:+.1f}' for e in g['error % at h->0'])}")
        for k in POWER_KEYS:
            g = f[k]
            print(f"  {k:20s} ~ h^({-g['s']:+.2f}), pair slopes {' / '.join(f'{e:.2f}' for e in g['pair slopes'])},"
                  f" change {' / '.join(f'{e:+.1f}' for e in g['change % per refinement'])} %")
    for k, g in out.get("field fits", {}).items():
        if "fit" in g:
            print(f"  field {k:20s} {g['values']}  order {g['fit']['free fit order']:.2f} "
                  f"({'ok' if g['fit']['free fit reliable'] else 'n/a'}) h->0 {g['fit']['h->0']:.5g}"
                  + (f"  richardson order {g['richardson 48/96/192']['order']}"
                     if "richardson 48/96/192" in g and g["richardson 48/96/192"]["order"] is not None else ""))


# ---------------------------------------------------------------- figure
def figure():
    """(a) the stress measures against the element size on log axes, (b) the ratios between them; the 96^3 field
    solid, the finest field dashed (same symbols)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, NullLocator
    sys.path.insert(0, str(ROOT / "ansys_usermat"))
    import figstyle
    figstyle.apply(size=9)
    out = load_json()
    Sd = out["stress"]
    fields = sorted(Sd, key=int)
    ls_of = {nf: ls for nf, ls in zip(fields, ("-", "--", ":"))}
    col = {"max": "#2a78d6", "p99": "#6250d6", "p95": "#1baf7a", "band": "#eb6834", "core": "#8c6a00",
           "film": "#444444", "first": "#b03060", "outer": "#777777"}
    mk = {"max": "o", "p99": "^", "p95": "v", "band": "s", "core": "D", "film": "x", "first": "p", "outer": "h"}

    def rows_of(nf):
        rows = Sd[nf]
        ns = sorted(int(x) for x in rows)
        return rows, ns, [L / n for n in ns]

    def axes_h(ax, hs, ns):
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(FixedLocator(hs)); ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xticklabels([f"{h:.2g}\n${n}^3$" for h, n in zip(hs, ns)], fontsize=7.5)
        ax.set_xlim(hs[0] * 1.15, hs[-1] / 1.15)
        ax.set_xlabel("stress element $h$ (mm), mesh")

    fig, axs = plt.subplots(1, 2, figsize=(6.3, 3.0))
    # (a) measures
    ax = axs[0]
    series = (("vM max", "maximum", "max"), ("vM p99", "99th percentile", "p99"),
              ("vM p95", "95th percentile", "p95"), ("band vM mean", "band mean", "band"),
              ("film vM mean", "film mean", "film"), ("core vM mean", "core mean", "core"))
    for nf in fields:
        rows, ns, hs = rows_of(nf)
        for key, lab, c in series:
            v = np.array([rows[str(n)][key] for n in ns], float) * 1e3
            ax.plot(hs, v, ls=ls_of[nf], marker=mk[c], color=col[c], lw=1.3, ms=4.5,
                    mfc="white" if c != "film" else col[c], label=lab if nf == fields[0] else None)
    rows, ns, hs = rows_of(fields[0])
    axes_h(ax, hs, ns)
    ax.set_yscale("log")
    ax.set_ylabel("von Mises in the film (mPa)")
    ax.set_title("(a) measures", fontsize=9, loc="left")
    f = out.get("stress fits", {}).get(fields[-1], {})
    if "vM max" in f:
        ax.text(0.03, 0.03, f"field ${fields[-1]}^3$: maximum $\\sim h^{{-{f['vM max']['s']:.2f}}}$, "
                            f"film mean $\\sim h^{{-{f['film vM mean']['free fit order'] * 0 + fit_power(hs if False else [L / n for n in f['meshes']], [Sd[fields[-1]][str(n)]['film vM mean'] for n in f['meshes']])['s']:.2f}}}$",
                transform=ax.transAxes, va="bottom", ha="left", fontsize=7)
    # (b) ratios
    ax = axs[1]
    ratios = (("vM max", "film vM mean", "maximum / film mean", "max"),
              ("vM p99", "film vM mean", "99th percentile / film mean", "p99"),
              ("vM p95", "film vM mean", "95th percentile / film mean", "p95"),
              ("band vM mean", "core vM mean", "band mean / core mean", "band"),
              ("first layer vM mean", "outer vM mean", "first layer / outer", "first"))
    for nf in fields:
        rows, ns, hs = rows_of(nf)
        for a, b, lab, c in ratios:
            pts = [(h, rows[str(n)][a] / rows[str(n)][b]) for h, n in zip(hs, ns)
                   if rows[str(n)][a] is not None and rows[str(n)][b]]
            if pts:
                ax.plot(*zip(*pts), ls=ls_of[nf], marker=mk[c], color=col[c], lw=1.3, ms=4.5, mfc="white",
                        label=lab if nf == fields[0] else None)
    ax.axhline(1.0, color="#999999", lw=0.8)
    rows, ns, hs = rows_of(fields[0])
    axes_h(ax, hs, ns)
    ax.set_ylim(0, None)
    ax.set_ylabel("ratio")
    ax.set_title("(b) ratios", fontsize=9, loc="left")
    ax.legend(fontsize=7, frameon=False, loc="upper right", handlelength=1.8)
    hl = axs[0].get_legend_handles_labels()
    h2 = [Line2D([], [], color="#333333", ls=ls_of[nf], lw=1.3) for nf in fields]
    l2 = [f"field ${nf}^3$" for nf in fields]
    fig.tight_layout(w_pad=1.0, rect=[0, 0.11, 1, 1])
    fig.legend(hl[0] + h2, hl[1] + l2, loc="lower center", ncol=4, fontsize=7, frameon=False,
               bbox_to_anchor=(0.5, 0.0), handlelength=1.8, columnspacing=1.0)
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG.with_suffix(".pdf"))
    fig.savefig(FIG.with_suffix(".png"), dpi=300)
    print("wrote", FIG.with_suffix(".pdf"), FIG.with_suffix(".png"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", type=int, nargs="+", default=list(FIELD_GRIDS))
    ap.add_argument("--meshes", type=int, nargs="+", default=list(STRESS_MESHES))
    ap.add_argument("--stress-fields", type=int, nargs="+", default=list(STRESS_FIELDS))
    ap.add_argument("--field-only", action="store_true")
    ap.add_argument("--stress-only", action="store_true", help="fields from the cache")
    ap.add_argument("--no-dt-check", action="store_true")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--refit", action="store_true", help="recompute the fits and print the summary")
    ap.add_argument("--force-stress", action="store_true", help="recompute stress meshes already in the JSON")
    a = ap.parse_args(argv)
    if a.refit:
        refit()
    if a.figure:
        figure()
    if a.refit or a.figure:
        return
    if not a.stress_only:
        run_fields(a.fields, dt_check=not a.no_dt_check)
    if not a.field_only:
        run_stress(a.meshes, a.stress_fields, force=a.force_stress)
        print_summary(load_json())


if __name__ == "__main__":
    main()

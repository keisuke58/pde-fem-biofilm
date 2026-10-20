#!/usr/bin/env python3
"""Growth stress in the measured 3D structure of an oral biofilm (CLSM/FISH stacks).

The stacks are the five-species in-vitro biofilms of N. Heine et al. (MHH/NIFE,
HOBIC flow chamber on titanium), decoded per voxel into five species by the
nife_intern pipeline (`results/fish_3d/ic/phi_xyz_<CH|DH>_d<day>.npy`). They are
read from that repository at run time and are not copied here: this script
writes only computed fields (stress, alpha) and aggregated numbers.

Processing of a stack (``load_stack``):
  - S = sum of the five decoded channels (each channel was divided by its own
    99.5th percentile in the decode, so S is a dimensionless intensity, 0 to ~15);
  - Gaussian smoothing, sigma 1.44 um laterally and 1 um in z, against the
    per-voxel decode speckle;
  - phi = clip((S - S_10) / (S_95 - S_10), 0, 1): the 10th percentile of the
    stack is the background (fluid or substratum, present in every stack), the
    95th percentile is phi = 1 (the densest twentieth of the box). The map is an
    assumption (intensity is not a calibrated volume fraction); its effect is
    checked with (0, S_95) and (S_10, S_99) on one stack;
  - the substratum: in every stack the film has one sharp, flat edge (the Ti
    disc) and a ragged side with protrusions (the fluid side). The sharp edge is
    at the low z index (assumption: the Leica stack starts at the substratum).
    The base z_b(x, y) = first slice with phi > 0.3 in the column; a plane is
    fitted to z_b (least squares, outliers beyond two slices dropped), the stack
    is tilted by up to 11 degrees in DH day 10; each column is shifted by an
    integer number of slices so that the plane is z = 0, slices below it are
    dropped, and the box is cut above the film (H_top = the top of the last
    slice where the film covers 5 % of the area, plus a margin of 30 %). The mirror assumption (substratum at the far side) is
    run as a check on one stack;
  - voxel size: 0.18 um pixels (nife_intern notes), 1024 px per field of view,
    stride 8 in the saved arrays: 1.44 um laterally (0.72 um for the 256-wide
    stack of DH day 6). The z step of a stack is not stored with the arrays; it
    is taken as the pooled stack depth of that condition and day in
    `results/diffusion_fit/zprofiles_all_ti.csv` divided by the number of slices
    (0.8 to 3.5 um, nominal 2 um). The Leica metadata would settle it.

Model (the thesis growth model, `wp3_reference.py`): on the box above the
substratum,
    d phi / dt   = beta lap phi + k_alpha alpha - P (max(0, phi - 1) + min(0, phi)),
    d alpha / dt = k_alpha phi,
with k_alpha = 1e-3 /T*, P = 5, T* = 1.1, zero flux on every face, RK4,
phi(0) = the CLSM phi, alpha(0) = 1. beta = 0.02 mm^2/T* (H / 0.25 mm)^2 with
H the film thickness of the stack (biovolume per substratum area, phi > 0.3), so
the spreading relative to the thickness is that of the thesis runs on the
0.25 mm film (assumption; the molar runs scale beta the same way). beta = 0 is
run as well: the stress of the structure as measured, without spreading.
Mechanics at T*: small strain, hex8 with B-bar, E = 10 Pa (phi^2 + 1e-3),
nu = 0.49, eigenstrain (alpha - 1) I, substratum clamped, top free, lateral
faces on rollers (normal displacement zero: the mirror symmetry of the zero-flux
growth field; free and periodic faces are run as checks), pyamg AMG-CG.

Measures over the film (phi(0) > 0.3): mean von Mises, by height (bottom quarter
on the substratum, middle half, top quarter of H_top), 95th percentile, mean pressure (compression positive), the
Pearson and Spearman correlation of the local von Mises with the local phi and
with the local thickness of the column, and the same for a laterally uniform
film with the stack's own mean profile phi(z), weighted with the film's voxel
count per slice.

    python ansys_usermat/apdl/clsm_microscale.py run                # all stacks, coarse mesh
    python ansys_usermat/apdl/clsm_microscale.py checks             # mesh, sub-region, faces, phi map, flip
    python ansys_usermat/apdl/clsm_microscale.py figure
    python ansys_usermat/apdl/clsm_microscale.py video --case DH_d21
Results: results/2026-10-3d/clsm_microscale.json; fields: figdata/clsm_microscale_*.npz.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import warnings
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NIFE = Path(os.environ.get("NIFE_INTERN", "/home/user/keisuke58/nife_intern"))
IC_DIR = NIFE / "results/fish_3d/ic"
ZPROF_CSV = NIFE / "results/diffusion_fit/zprofiles_all_ti.csv"
OUT_JSON = HERE / "results/2026-10-3d/clsm_microscale.json"
FIGDATA = ROOT / "figdata"
ASSETS = ROOT / "assets"

CASES = [f"{c}_d{d}" for c in ("CH", "DH") for d in (1, 6, 10, 15, 21)]
PIX_UM_FULL = 0.18          # Leica pixel, nife_intern docs
PX_FOV = 1024               # pixels per field of view (stride 8 -> 128)
SIG_XY_UM, SIG_Z_UM = 1.44, 1.0
P_BG, P_REF = 10.0, 95.0
THR = 0.3                   # film = phi > THR
COV_TOP = 0.05              # H_top: top of the last slice where the film covers this fraction of the area
E_BIO, NU, FLOOR = 10.0, 0.49, 1e-3
K_ALPHA, PENALTY, T_END = 1e-3, 5.0, 1.1
BETA0, H_REF = 0.02, 0.25   # mm^2/T*, mm


# ----------------------------------------------------------------------------- data
def depth_table():
    """Pooled stack depth (um) per (condition, day) from the nife_intern z profiles."""
    d = {}
    with open(ZPROF_CSV) as f:
        for r in csv.DictReader(f):
            k = (r["condition"], int(r["day"]))
            d[k] = max(d.get(k, 0.0), float(r["z_um"]))
    return d


def fit_plane(zb, max_res=2.0, iters=4):
    """Least-squares plane through the base heights zb(x, y) (nan = no film), outliers dropped."""
    nx, ny = zb.shape
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    ok = np.isfinite(zb)
    A = np.stack([np.ones(ok.sum()), X[ok], Y[ok]], 1)
    b = zb[ok]
    keep = np.ones(len(b), bool)
    for _ in range(iters):
        c, *_ = np.linalg.lstsq(A[keep], b[keep], rcond=None)
        r = A @ c - b
        keep_new = np.abs(r) <= max_res
        if keep_new.sum() < 50 or np.array_equal(keep_new, keep):
            break
        keep = keep_new
    plane = c[0] + c[1] * X + c[2] * Y
    return plane, c, float(np.sqrt(np.mean(r[keep] ** 2))), int(keep.sum())


def load_stack(case, flip=False, p_bg=P_BG, p_ref=P_REF, s0=False, verbose=True):
    """Processed stack: dict(phi (nx, ny, nz) above the substratum, hx, hz [mm], thick [um] per column, meta)."""
    cond, day = case.split("_d")
    day = int(day)
    path = NIFE / "results/fish_3d/phi_xyz_DH_d6_s0.npy" if s0 else IC_DIR / f"phi_xyz_{case}.npy"
    a = np.load(path, mmap_mode="r")
    S = np.asarray(a, np.float32).sum(-1)                    # (nx, ny, nz), dimensionless intensity
    nx, ny, nz = S.shape
    pix_um = PIX_UM_FULL * PX_FOV / nx
    dz_um = depth_table()[(cond, day)] / nz
    if flip:
        S = S[:, :, ::-1]
    Ss = gaussian_filter(S, (SIG_XY_UM / pix_um, SIG_XY_UM / pix_um, SIG_Z_UM / dz_um), mode="nearest")
    lo, hi = np.percentile(Ss, [p_bg, p_ref])
    phi = np.clip((Ss - lo) / (hi - lo), 0.0, 1.0)
    film = phi > THR
    has = film.any(2)
    zb = np.where(has, film.argmax(2), np.nan).astype(float)      # first film slice of the column
    plane, c, rms, nfit = fit_plane(zb)
    tilt = np.degrees(np.arctan(np.hypot(c[1] * dz_um / pix_um, c[2] * dz_um / pix_um)))
    shift = np.rint(plane).astype(int)
    # flatten: column (x, y) slice z goes to z - shift; the box keeps z' = 0 .. nz_box - 1
    nz_all = nz - int(np.min(shift))
    out = np.zeros((nx, ny, nz_all), np.float32)
    for sft in np.unique(shift):
        m = shift == sft
        src = phi[m]                                          # (ncol, nz)
        k0, k1 = max(sft, 0), min(nz, nz_all + sft)
        if k1 > k0:
            out[m, (k0 - sft):(k1 - sft)] = src[:, k0:k1]
    cov = (out > THR).mean((0, 1))                            # film coverage per slice
    k_top = int(np.where(cov >= COV_TOP)[0].max()) + 1       # slices up to the top of the film
    H_top = k_top * dz_um
    nz_box = max(2, min(k_top + max(3, int(round(0.3 * k_top))), nz_all))
    out = np.ascontiguousarray(out[:, :, :nz_box])
    filmb = out > THR
    thick = filmb.sum(2) * dz_um                               # biovolume per column, um
    ztop = np.where(filmb.any(2), (nz_box - filmb[:, :, ::-1].argmax(2)) * dz_um, np.nan)
    H = float(thick.mean())
    meta = {"file": path.name, "shape": [nx, ny, nz], "pixel um": pix_um, "dz um": dz_um,
            "S percentiles": {"bg": float(lo), "ref": float(hi)}, "p_bg": p_bg, "p_ref": p_ref, "flip": flip,
            "plane slices": [float(v) for v in c], "plane rms slices": rms, "plane columns": nfit,
            "tilt deg": float(tilt), "shift range": [int(shift.min()), int(shift.max())], "nz box": int(nz_box),
            "box um": [nx * pix_um, ny * pix_um, nz_box * dz_um], "H um": H, "H_top um": H_top,
            "film fraction of box": float(filmb.mean()), "columns with film": float(filmb.any(2).mean()),
            "mean phi film": float(out[filmb].mean()) if filmb.any() else 0.0}
    if verbose:
        print(f"  {case}{' flipped' if flip else ''}{' s0' if s0 else ''}: {nx}x{ny}x{nz} voxels of "
              f"{pix_um:.2f} x {pix_um:.2f} x {dz_um:.2f} um, S bg/ref {lo:.2f}/{hi:.2f}, substratum plane tilt "
              f"{tilt:.1f} deg (rms {rms:.2f} slices), box {nz_box} slices = {nz_box * dz_um:.1f} um, "
              f"H {H:.1f} um, H_top {H_top:.1f} um, film {filmb.mean() * 100:.1f} % of the box", flush=True)
    return {"phi": out, "hx": pix_um * 1e-3, "hz": dz_um * 1e-3, "thick": thick, "ztop": ztop, "meta": meta}


def coarsen(st, f):
    """Lateral block mean by the factor f (z kept)."""
    if f == 1:
        return st
    phi = st["phi"]
    nx, ny, nz = phi.shape
    nx2, ny2 = nx // f, ny // f
    phi2 = phi[:nx2 * f, :ny2 * f].reshape(nx2, f, ny2, f, nz).mean((1, 3))
    thick = st["thick"][:nx2 * f, :ny2 * f].reshape(nx2, f, ny2, f).mean((1, 3))
    zt = st["ztop"][:nx2 * f, :ny2 * f].reshape(nx2, f, ny2, f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ztop = np.nanmean(zt, (1, 3))
    return {"phi": phi2.astype(np.float32), "hx": st["hx"] * f, "hz": st["hz"], "thick": thick, "ztop": ztop,
            "meta": dict(st["meta"], coarsen=f)}


def subregion(st, i0, i1, j0, j1):
    return {"phi": st["phi"][i0:i1, j0:j1], "hx": st["hx"], "hz": st["hz"], "thick": st["thick"][i0:i1, j0:j1],
            "ztop": st["ztop"][i0:i1, j0:j1], "meta": dict(st["meta"], subregion=[i0, i1, j0, j1])}


def uniform_stack(st):
    """Laterally uniform film with the stack's mean profile phi(z) on a 2 x 2 column (rollers: exact 1D)."""
    prof = st["phi"].mean((0, 1))
    phi = np.broadcast_to(prof[None, None, :], (2, 2, len(prof))).astype(np.float32).copy()
    return {"phi": phi, "hx": st["hx"], "hz": st["hz"], "thick": np.full((2, 2), st["thick"].mean()),
            "ztop": st["ztop"], "meta": dict(st["meta"], uniform=True), "profile": prof}


# ----------------------------------------------------------------------------- growth
def laplacian(u, h):
    hx, hy, hz = h
    up = np.pad(u, 1, mode="edge")                              # mirror cells: zero flux on every face
    c = up[1:-1, 1:-1, 1:-1]
    return ((up[2:, 1:-1, 1:-1] + up[:-2, 1:-1, 1:-1] - 2 * c) / hx ** 2
            + (up[1:-1, 2:, 1:-1] + up[1:-1, :-2, 1:-1] - 2 * c) / hy ** 2
            + (up[1:-1, 1:-1, 2:] + up[1:-1, 1:-1, :-2] - 2 * c) / hz ** 2)


def grow(phi0, h, beta, T=T_END, lam=0.45, verbose=False):
    """RK4 of the growth system to T; returns phi(T), alpha(T) - 1."""
    phi = np.array(phi0, np.float64)
    alpha = np.ones_like(phi)
    if beta > 0:
        dt = lam / (beta * sum(1.0 / hh ** 2 for hh in h))   # RK4 stability for the 7-point stencil
    else:
        dt = T / 10
    nstep = max(10, int(np.ceil(T / dt)))
    dt = T / nstep

    def rhs(p, a):
        lap = beta * laplacian(p, h) if beta > 0 else 0.0
        dp = lap + K_ALPHA * a - PENALTY * (np.maximum(p - 1.0, 0.0) + np.minimum(p, 0.0))
        return dp, K_ALPHA * p
    t0 = time.time()
    for s in range(nstep):
        k1p, k1a = rhs(phi, alpha)
        k2p, k2a = rhs(phi + 0.5 * dt * k1p, alpha + 0.5 * dt * k1a)
        k3p, k3a = rhs(phi + 0.5 * dt * k2p, alpha + 0.5 * dt * k2a)
        k4p, k4a = rhs(phi + dt * k3p, alpha + dt * k3a)
        phi += dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        alpha += dt / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a)
    if verbose:
        print(f"    growth: beta {beta:.3e} mm^2/T*, {nstep} RK4 steps, phi in [{phi.min():.4f}, {phi.max():.4f}], "
              f"alpha - 1 max {alpha.max() - 1:.3e} ({time.time() - t0:.0f} s)", flush=True)
    return phi, alpha - 1.0, nstep


# ----------------------------------------------------------------------------- mechanics
def element_matrices(h, nu=NU):
    """hex8 of edges (hx, hy, hz): deviatoric part at 2x2x2 points, volumetric at the centre (B-bar)."""
    hx, hy, hz = h
    g = np.array([-1, 1]) / np.sqrt(3)
    corners = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                        [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)

    def bmat(xi):
        dN = np.array([[cx * (1 + xi[1] * cy) * (1 + xi[2] * cz) * 2 / hx,
                        cy * (1 + xi[0] * cx) * (1 + xi[2] * cz) * 2 / hy,
                        cz * (1 + xi[0] * cx) * (1 + xi[1] * cy) * 2 / hz] for cx, cy, cz in corners]) / 8
        b = np.zeros((6, 24))
        for a in range(8):
            dx, dy, dz = dN[a]
            i = 3 * a
            b[0, i] = dx; b[1, i + 1] = dy; b[2, i + 2] = dz
            b[3, i] = dy; b[3, i + 1] = dx
            b[4, i + 1] = dz; b[4, i + 2] = dy
            b[5, i] = dz; b[5, i + 2] = dx
        return b
    lam = nu / ((1 + nu) * (1 - 2 * nu))
    mu = 1 / (2 * (1 + nu))
    Dd = np.diag([2 * mu] * 3 + [mu] * 3)
    Dv = np.zeros((6, 6)); Dv[:3, :3] = lam
    m = np.array([1, 1, 1, 0, 0, 0.])
    vol = hx * hy * hz
    Kd = np.zeros((24, 24)); fd = np.zeros(24)
    for a in g:
        for b_ in g:
            for c in g:
                B = bmat((a, b_, c))
                Kd += B.T @ Dd @ B * vol / 8
                fd += B.T @ Dd @ m * vol / 8
    B0 = bmat((0, 0, 0))
    return Kd + B0.T @ Dv @ B0 * vol, fd + B0.T @ Dv @ m * vol, B0, Dd + Dv, m


def solve_mech(phi_T, am1, h, lateral="roller", tol=1e-9, chunk=32768, verbose=True):
    """Stress at T* on the voxel mesh; substratum (z = 0) clamped, top free, lateral faces roller | free | periodic.
    Returns dict(vm, p [Pa] per element, solver)."""
    import pyamg
    t_start = time.time()
    n1, n2, n3 = phi_T.shape
    hx, hy, hz = h
    per = lateral == "periodic"
    ng = (n1 if per else n1 + 1, n2 if per else n2 + 1, n3 + 1)
    idx = np.arange(np.prod(ng), dtype=np.int64).reshape(ng)
    idxp = idx[np.arange(n1 + 1) % ng[0]][:, np.arange(n2 + 1) % ng[1]]
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idxp[i:i + n1, j:j + n2, k:k + n3].ravel() for i, j, k in off], 1)
    ne = len(conn)
    nnod = int(np.prod(ng))
    Ee = E_BIO * (np.clip(phi_T, 0, 1).ravel() ** 2 + FLOOR)
    eps0 = am1.ravel()
    fixed = np.zeros((nnod, 3), bool)
    fixed[idxp[:, :, 0].ravel(), :] = True                     # substratum clamped
    if lateral == "roller":
        fixed[idxp[0, :, :].ravel(), 0] = True; fixed[idxp[-1, :, :].ravel(), 0] = True
        fixed[idxp[:, 0, :].ravel(), 1] = True; fixed[idxp[:, -1, :].ravel(), 1] = True
    elif lateral not in ("free", "periodic"):
        raise ValueError(lateral)
    free = ~fixed.ravel()
    nfree = int(free.sum())
    dmap = np.full(3 * nnod, -1, np.int64); dmap[free] = np.arange(nfree)
    dofs = (3 * conn[:, :, None] + np.arange(3)[None, None, :]).reshape(ne, 24)
    fd = dmap[dofs].astype(np.int32)
    Ke0, fe0, B0, D, mvec = element_matrices(h)
    t0 = time.time()
    F = np.zeros(nfree)
    fe = fe0[None] * (Ee * eps0)[:, None]
    ok = fd >= 0
    np.add.at(F, fd[ok], fe[ok])
    Kff = None
    for c0 in range(0, ne, chunk):
        d = fd[c0:c0 + chunk]
        rows = np.repeat(d, 24, 1).ravel()
        cols = np.tile(d, (1, 24)).ravel()
        val = (Ke0[None] * Ee[c0:c0 + chunk, None, None]).ravel()
        keep = (rows >= 0) & (cols >= 0)
        Kc = sp.csr_matrix((val[keep], (rows[keep], cols[keep])), shape=(nfree, nfree))
        Kff = Kc if Kff is None else Kff + Kc
        del rows, cols, val, keep, Kc
    Kff.sum_duplicates()
    t_asm = time.time() - t0
    # near-nullspace: rigid-body modes on the free dofs
    t0 = time.time()
    gx, gy, gz = np.meshgrid(np.arange(ng[0]) * hx, np.arange(ng[1]) * hy, np.arange(ng[2]) * hz, indexing="ij")
    P = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)
    fdof = np.where(free)[0]
    node, comp = fdof // 3, fdof % 3
    Pn = P[node]
    B = np.zeros((nfree, 6))
    for k in range(3):
        B[comp == k, k] = 1.0
    B[comp == 0, 3] = -Pn[comp == 0, 1]; B[comp == 1, 3] = Pn[comp == 1, 0]
    B[comp == 1, 4] = -Pn[comp == 1, 2]; B[comp == 2, 4] = Pn[comp == 2, 1]
    B[comp == 0, 5] = Pn[comp == 0, 2]; B[comp == 2, 5] = -Pn[comp == 2, 0]
    ml = pyamg.smoothed_aggregation_solver(Kff, B=B, max_coarse=1000, coarse_solver="splu")
    t_setup = time.time() - t0
    res = []
    uf = ml.solve(F, tol=tol, accel="cg", maxiter=4000, residuals=res)
    rel = float(np.linalg.norm(Kff @ uf - F) / np.linalg.norm(F))
    t_solve = time.time() - t0 - t_setup
    del ml, B, Kff
    u = np.zeros(3 * nnod)
    u[fdof] = uf
    vm = np.empty(ne); pr = np.empty(ne)
    for c0 in range(0, ne, 4 * chunk):
        sl = slice(c0, c0 + 4 * chunk)
        eps = (B0 @ u[dofs[sl]].T).T - eps0[sl, None] * mvec
        sig = Ee[sl, None] * (eps @ D.T)
        p = -sig[:, :3].mean(1)
        s = sig[:, :3] + p[:, None]
        vm[sl] = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
        pr[sl] = p
    solver = {"lateral": lateral, "dof": nfree, "elements": ne, "iterations": len(res), "relative residual": rel,
              "seconds assembly": t_asm, "seconds setup": t_setup, "seconds solve": t_solve,
              "seconds": time.time() - t_start}
    if verbose:
        print(f"    mechanics ({lateral}): {ne} elements, {nfree} dof, {len(res)} it, residual {rel:.1e}, "
              f"assembly {t_asm:.0f} s, setup {t_setup:.0f} s, solve {t_solve:.0f} s", flush=True)
    return {"vm": vm.reshape(phi_T.shape), "p": pr.reshape(phi_T.shape), "solver": solver}


# ----------------------------------------------------------------------------- measures
def wpercentile(v, w, q):
    o = np.argsort(v)
    cw = np.cumsum(w[o]) - 0.5 * w[o]
    return float(np.interp(q / 100 * w.sum(), cw, v[o]))


def measures(st, vm, p, ref_film=None):
    """Means and percentiles over the film (phi(0) > THR) and by height; ref_film: the film mask and thickness of the
    real stack when st is its uniform counterpart (weights per slice)."""
    from scipy.stats import pearsonr, spearmanr
    dz_um = st["hz"] * 1e3
    nz = st["phi"].shape[2]
    z = (np.arange(nz) + 0.5) * dz_um
    if ref_film is None:
        film = st["phi"] > THR
        H_top = st["meta"]["H_top um"]
        thick = st["thick"]
    else:
        film, H_top, thick = ref_film["film"], ref_film["H_top um"], ref_film["thick"]
    out = {"film voxels": int(film.sum()), "H_top um": H_top}
    zb = max(0.25 * H_top, dz_um)                              # the bottom band holds at least the first slice
    bands = {"all": np.ones(nz, bool), "bottom": z < zb, "middle": (z >= zb) & (z < 0.75 * H_top),
             "top": (z >= 0.75 * H_top) & (z >= zb)}
    if ref_film is None:
        vmf, pf = vm[film], p[film]
        zk = np.broadcast_to(z[None, None, :], film.shape)[film]
        for name, bm in bands.items():
            sel = bm[np.searchsorted(z, zk)] if name != "all" else np.ones(len(vmf), bool)
            out[f"vM mean {name}"] = float(vmf[sel].mean()) if sel.any() else None
            out[f"p mean {name}"] = float(pf[sel].mean()) if sel.any() else None
            out[f"voxels {name}"] = int(sel.sum())
        out["vM p95"] = float(np.percentile(vmf, 95))
        out["vM p99"] = float(np.percentile(vmf, 99))
        out["vM max"] = float(vmf.max())
        out["p max"] = float(pf.max())
        phif = st["phi"][film]
        tk = np.broadcast_to(thick[:, :, None], film.shape)[film]
        out["corr vM phi"] = {"pearson": float(pearsonr(vmf, phif)[0]), "spearman": float(spearmanr(vmf, phif)[0])}
        out["corr vM thickness"] = {"pearson": float(pearsonr(vmf, tk)[0]), "spearman": float(spearmanr(vmf, tk)[0])}
        out["corr vM height"] = {"pearson": float(pearsonr(vmf, zk)[0]), "spearman": float(spearmanr(vmf, zk)[0])}
        cnt = film.sum((0, 1))
        prof = np.array([vm[:, :, k][film[:, :, k]].mean() if cnt[k] else np.nan for k in range(nz)])
        pprof = np.array([p[:, :, k][film[:, :, k]].mean() if cnt[k] else np.nan for k in range(nz)])
        # stress against the thickness of the column: binned means (10 bins of thickness over the film columns)
        tcol = thick[film.any((0, 1)) if False else slice(None)]
        edges = np.percentile(tk, np.linspace(0, 100, 11))
        bi = np.clip(np.searchsorted(edges, tk, side="right") - 1, 0, 9)
        out["thickness bins um"] = [float(0.5 * (edges[i] + edges[i + 1])) for i in range(10)]
        out["vM by thickness bin"] = [float(vmf[bi == i].mean()) if (bi == i).any() else None for i in range(10)]
    else:
        cnt = film.sum((0, 1)).astype(float)
        vmz, pz = vm[0, 0, :], p[0, 0, :]
        for name, bm in bands.items():
            w = cnt * bm
            out[f"vM mean {name}"] = float((vmz * w).sum() / w.sum()) if w.sum() else None
            out[f"p mean {name}"] = float((pz * w).sum() / w.sum()) if w.sum() else None
            out[f"voxels {name}"] = int(w.sum())
        out["vM p95"] = wpercentile(vmz, cnt, 95)
        out["vM p99"] = wpercentile(vmz, cnt, 99)
        out["vM max"] = float(vmz[cnt > 0].max())
        out["p max"] = float(pz[cnt > 0].max())
        prof, pprof = vmz, pz
    out["z um"] = [float(v) for v in z]
    out["vM profile"] = [None if np.isnan(v) else float(v) for v in prof]
    out["p profile"] = [None if np.isnan(v) else float(v) for v in pprof]
    out["film count per slice"] = [int(v) for v in film.sum((0, 1))]
    return out


# ----------------------------------------------------------------------------- runs
def load_json():
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text())
    return {"model": {
        "growth": "d phi/dt = beta lap phi + k_alpha alpha - P(max(0, phi - 1) + min(0, phi)); d alpha/dt = k_alpha phi; "
                  "RK4, zero flux on every face of the box, phi(0) = CLSM phi, alpha(0) = 1",
        "k_alpha 1/T*": K_ALPHA, "penalty": PENALTY, "T*": T_END,
        "beta": f"{BETA0} mm^2/T* (H / {H_REF} mm)^2, H = biovolume per substratum area of the stack (phi > {THR})",
        "mechanics": f"hex8 B-bar, E = {E_BIO} Pa (phi^2 + {FLOOR}), nu = {NU}, eigenstrain (alpha - 1) I, substratum "
                     "clamped, top free, lateral faces on rollers; pyamg smoothed aggregation + CG",
        "phi": f"clip((S - S_p{P_BG:g}) / (S_p{P_REF:g} - S_p{P_BG:g}), 0, 1) of the smoothed channel sum S "
               f"(sigma {SIG_XY_UM} um lateral, {SIG_Z_UM} um in z); film = phi > {THR}",
        "substratum": "plane fitted to the first film slice of every column (the sharp edge at the low z index), "
                      "columns shifted to z = 0 (integer slices), box cut at H_top + 30 % (at least 3 slices)",
        "voxel": f"{PIX_UM_FULL} um pixels, {PX_FOV} px per field, lateral stride of the arrays; z step = pooled "
                 "stack depth of the day (zprofiles_all_ti.csv) / slices of the stack",
        "bands": "bottom: z < H_top/4 (at least the first slice), middle: to 3 H_top/4, top: above; H_top = top of "
                 "the last slice where the film covers at least 5 % of the area", "units": "stress Pa, lengths mm in the solver, um in the tables"},
        "cases": {}, "checks": {}}


def save_json(out):
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=1))


def beta_of(st):
    return BETA0 * (st["meta"]["H um"] * 1e-3 / H_REF) ** 2


def run_one(st, beta, lateral="roller", tag=None, cache=True, verbose=True):
    """Growth + mechanics on the stack st; returns (measures, fields)."""
    h = (st["hx"], st["hx"], st["hz"])
    phi_T, am1, nstep = grow(st["phi"], h, beta, verbose=verbose)
    mech = solve_mech(phi_T, am1, h, lateral, verbose=verbose)
    m = measures(st, mech["vm"], mech["p"])
    m.update({"beta mm2/T*": beta, "rk4 steps": nstep, "phi(T) max": float(phi_T.max()),
              "alpha-1 max": float(am1.max()), "alpha-1 film mean": float(am1[st["phi"] > THR].mean()),
              "mesh": list(st["phi"].shape), "h um": [v * 1e3 for v in h], "solver": mech["solver"]})
    if cache and tag:
        FIGDATA.mkdir(exist_ok=True)
        np.savez_compressed(FIGDATA / f"clsm_microscale_{tag}.npz", vm=mech["vm"].astype(np.float32),
                            p=mech["p"].astype(np.float32), am1=am1.astype(np.float32),
                            h_mm=np.array(h), shape=np.array(st["phi"].shape))
    return m, {"phi_T": phi_T, "am1": am1, "vm": mech["vm"], "p": mech["p"]}


def run_uniform(st, beta, verbose=True):
    u = uniform_stack(st)
    h = (u["hx"], u["hx"], u["hz"])
    phi_T, am1, nstep = grow(u["phi"], h, beta, verbose=False)
    mech = solve_mech(phi_T, am1, h, "roller", verbose=False)
    ref = {"film": st["phi"] > THR, "H_top um": st["meta"]["H_top um"], "thick": st["thick"]}
    m = measures(u, mech["vm"], mech["p"], ref_film=ref)
    m.update({"beta mm2/T*": beta, "rk4 steps": nstep, "mean profile phi": [float(v) for v in u["profile"]],
              "alpha-1 max": float(am1.max())})
    if verbose:
        print_m("uniform film", m)
    return m


def print_m(label, m):
    f = lambda k: (f"{m[k] * 1e3:.4f}" if m.get(k) is not None else "   -  ")
    print(f"    {label}: vM mean {f('vM mean all')} mPa, bottom {f('vM mean bottom')}, middle {f('vM mean middle')}, "
          f"top {f('vM mean top')}, p95 {f('vM p95')}, p mean {f('p mean all')}, max {f('vM max')}"
          + (f"; corr phi {m['corr vM phi']['spearman']:+.2f} thick {m['corr vM thickness']['spearman']:+.2f} (Spearman)"
             if "corr vM phi" in m else ""), flush=True)


def cmd_run(a):
    out = load_json()
    f = {"coarse": 2, "full": 1}[a.level]
    for case in a.cases:
        rec = out["cases"].setdefault(case, {"runs": {}})
        st0 = load_stack(case)
        rec["data"] = st0["meta"]
        st = coarsen(st0, f)
        beta = beta_of(st0)
        for bname, b in (("beta", beta), ("beta0", 0.0)):
            if a.only_beta and bname != a.only_beta:
                continue
            key = f"{a.level}_{bname}"
            if key in rec["runs"] and not a.force:
                print(f"  {case} {key}: done, skipped", flush=True)
                continue
            t0 = time.time()
            m, _ = run_one(st, b, tag=f"{case}_{key}")
            print_m(f"{case} {key}", m)
            rec["runs"][key] = m
            rec["runs"][f"uniform_{bname}"] = run_uniform(st, b)
            rec["runs"][key]["seconds total"] = time.time() - t0
            save_json(out)
    print("wrote", OUT_JSON)


def cmd_checks(a):
    out = load_json()
    C = out["checks"]
    want = set(a.which) if a.which else {"faces", "phimap", "flip", "mesh", "subregion", "mesh0", "subregion0"}
    # 1. lateral faces, phi map, flip on DH_d6 coarse
    case = "DH_d6"
    st0 = load_stack(case)
    st = coarsen(st0, 2)
    beta = beta_of(st0)
    if "faces" in want and ("faces" not in C or a.force):
        C["faces"] = {"case": case, "level": "coarse", "beta": beta, "runs": {}}
        for lat in ("roller", "free", "periodic"):
            m, _ = run_one(st, beta, lateral=lat, cache=False)
            # interior region: columns more than a quarter of the width away from the lateral faces
            n1, n2 = st["phi"].shape[:2]
            print_m(f"faces {lat}", m)
            C["faces"]["runs"][lat] = m
            save_json(out)
    if "phimap" in want and ("phimap" not in C or a.force):
        C["phimap"] = {"case": case, "level": "coarse", "runs": {}}
        for pb, pr in ((P_BG, P_REF), (0.0, P_REF), (P_BG, 99.0), (P_BG, 90.0)):
            s0 = load_stack(case, p_bg=pb, p_ref=pr)
            s = coarsen(s0, 2)
            b = beta_of(s0)
            m, _ = run_one(s, b, cache=False)
            m["data"] = s0["meta"]
            print_m(f"phi map p{pb:g}-p{pr:g}", m)
            C["phimap"]["runs"][f"p{pb:g}_p{pr:g}"] = m
            save_json(out)
    if "flip" in want and ("flip" not in C or a.force):
        s0 = load_stack(case, flip=True)
        s = coarsen(s0, 2)
        b = beta_of(s0)
        m, _ = run_one(s, b, cache=False)
        m["data"] = s0["meta"]
        m["uniform"] = run_uniform(s, b)
        print_m("flipped", m)
        C["flip"] = {"case": case, "level": "coarse", "run": m}
        save_json(out)
    # 2. mesh: full lateral resolution (1.44 um) against the coarse mesh (2.88 um) on two stacks
    for key in ("mesh", "mesh0"):
        if key not in want:
            continue
        C.setdefault(key, {})
        for case in a.mesh_cases:
            if case in C[key] and not a.force:
                continue
            s0 = load_stack(case)
            b = beta_of(s0) if key == "mesh" else 0.0
            rows = {}
            for lev, f in (("coarse", 2), ("full", 1)):
                s = coarsen(s0, f)
                m, _ = run_one(s, b, tag=f"{case}_{lev}_{'beta' if key == 'mesh' else 'beta0'}")
                m["uniform"] = run_uniform(s, b)
                print_m(f"{key} {case} {lev}", m)
                rows[lev] = m
                C[key][case] = rows
                save_json(out)
    # 3. sub-region of the 256-wide stack of DH day 6 at 0.72, 1.44 and 2.88 um
    for key in ("subregion", "subregion0"):
        if key not in want or (key in C and not a.force):
            continue
        s0 = load_stack("DH_d6", s0=True)
        sub = subregion(s0, 64, 192, 64, 192)                 # central 92 x 92 um
        b = beta_of(s0) if key == "subregion" else 0.0
        C[key] = {"case": "DH_d6 (256-wide stack)", "region": "central 128 x 128 voxels = 92 x 92 um",
                  "beta": b, "data": s0["meta"], "runs": {}}
        for f in (4, 2, 1):
            s = coarsen(sub, f)
            m, _ = run_one(s, b, tag=f"DH_d6_sub_f{f}_{'beta' if key == 'subregion' else 'beta0'}")
            print_m(f"{key} h = {s['hx'] * 1e3:.2f} um", m)
            C[key]["runs"][f"h{s['hx'] * 1e3:.2f}"] = m
            save_json(out)
    print("wrote", OUT_JSON)


# ----------------------------------------------------------------------------- figure
def cmd_figure(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "ansys_usermat"))
    import figstyle
    figstyle.apply(9)
    out = load_json()
    days = (1, 6, 10, 15, 21)
    cm = plt.get_cmap("viridis")
    fig, axs = plt.subplots(2, 2, figsize=(6.3, 5.6))
    for ax, cond in zip(axs[0], ("CH", "DH")):
        for i, d in enumerate(days):
            rec = out["cases"].get(f"{cond}_d{d}")
            if not rec or "coarse_beta" not in rec["runs"]:
                continue
            r, u = rec["runs"]["coarse_beta"], rec["runs"]["uniform_beta"]
            z = np.array(r["z um"])
            vm = np.array([np.nan if v is None else v for v in r["vM profile"]]) * 1e3
            vu = np.array(u["vM profile"]) * 1e3
            cnt = np.array(r["film count per slice"])
            okk = (cnt > 0) & (z <= rec["data"]["H_top um"])
            col = cm(i / (len(days) - 1))
            ax.plot(vm[okk], z[okk], color=col, lw=1.3, label=f"day {d}")
            ax.plot(vu[okk], z[okk], color=col, lw=1.0, ls="--")
        ax.set_xlabel("mean von Mises in the film (mPa)")
        ax.set_ylabel("height above the substratum ($\\mu$m)")
        ax.set_title(f"({'ab'[cond == 'DH']}) {cond}: real structure (solid), uniform film (dashed)", loc="left",
                     fontsize=9)
        ax.set_ylim(0, None)
        ax.legend(frameon=False, fontsize=7)
    # (c) section of the von Mises field of one stack
    ax = axs[1, 0]
    case = a.case
    npz = FIGDATA / f"clsm_microscale_{case}_coarse_beta.npz"
    if npz.exists():
        d = np.load(npz)
        vm = d["vm"] * 1e3
        h = d["h_mm"] * 1e3
        st = coarsen(load_stack(case, verbose=False), 2)
        j = vm.shape[1] // 2
        sec = vm[:, j, :].T
        film = (st["phi"] > THR)[:, j, :].T
        vmax = np.percentile(vm[st["phi"] > THR], 95)
        im = ax.imshow(np.where(film, sec, np.nan), origin="lower", cmap="viridis", vmin=0, vmax=vmax, aspect="auto",
                       extent=[0, vm.shape[0] * h[0], 0, vm.shape[2] * h[2]], interpolation="nearest")
        ax.imshow(np.where(film, np.nan, 0.0), origin="lower", cmap="Greys", vmin=-1, vmax=1, alpha=0.25, aspect="auto",
                  extent=[0, vm.shape[0] * h[0], 0, vm.shape[2] * h[2]], interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label("von Mises (mPa)")
        ax.set_xlabel("x ($\\mu$m)")
        ax.set_ylabel("z ($\\mu$m)")
        ax.set_title(f"(c) {case.replace('_d', ' day ')}: section y = {j * h[1]:.0f} $\\mu$m, film only", loc="left",
                     fontsize=9)
        ax.grid(False)
    # (d) stress against the local thickness
    ax = axs[1, 1]
    for cond, base in (("CH", "Blues"), ("DH", "Reds")):
        cmc = plt.get_cmap(base)
        for i, d in enumerate(days):
            rec = out["cases"].get(f"{cond}_d{d}")
            if not rec or "coarse_beta" not in rec["runs"]:
                continue
            r = rec["runs"]["coarse_beta"]
            tb = np.array(r["thickness bins um"])
            vb = np.array([np.nan if v is None else v for v in r["vM by thickness bin"]]) * 1e3
            ax.plot(tb, vb / r["vM mean all"] / 1e3, color=cmc(0.35 + 0.6 * i / (len(days) - 1)), lw=0.9,
                    marker="o", ms=2, label=f"{cond} day {d}")
    ax.axhline(1, color="grey", lw=0.6)
    ax.set_xlabel("local film thickness of the column ($\\mu$m)")
    ax.set_ylabel("von Mises / film mean")
    ax.set_title("(d) stress against the local thickness (deciles)", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=6, ncol=2)
    fig.tight_layout()
    ASSETS.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(ASSETS / f"fig_clsm_microscale.{ext}", dpi=400 if ext == "png" else None)
    print("wrote", ASSETS / "fig_clsm_microscale.pdf")


# ----------------------------------------------------------------------------- video
def cmd_video(a):
    import pyvista as pv
    pv.OFF_SCREEN = True
    case = a.case
    st = coarsen(load_stack(case), 2)
    beta = beta_of(load_stack(case, verbose=False))
    npz = FIGDATA / f"clsm_microscale_{case}_coarse_{a.run}.npz"
    d = np.load(npz)
    vm = d["vm"] * 1e3
    h = (st["hx"], st["hx"], st["hz"])
    phi_T, am1, _ = grow(st["phi"], h, beta if a.run == "beta" else 0.0)
    nx, ny, nz = vm.shape
    sp_um = tuple(v * 1e3 for v in h)
    g = pv.ImageData(dimensions=(nx + 1, ny + 1, nz + 1), spacing=sp_um, origin=(0, 0, 0))
    g.cell_data["phi"] = phi_T.ravel(order="F")
    g.cell_data["vm"] = vm.ravel(order="F")
    gp = g.cell_data_to_point_data()
    iso = gp.contour([THR], scalars="phi")
    vmax = float(np.percentile(iso.point_data["vm"], 95))
    Lx, Ly, Lz = nx * sp_um[0], ny * sp_um[1], nz * sp_um[2]
    pl = pv.Plotter(off_screen=True, window_size=(1600, 1000))
    pl.set_background("white")
    pl.add_mesh(iso, scalars="vm", cmap="viridis", clim=(0, vmax), smooth_shading=True, ambient=0.35, diffuse=0.7,
                specular=0.1, scalar_bar_args=dict(title=" ", vertical=False, position_x=0.25,
                                                   position_y=0.03, width=0.5, height=0.06, label_font_size=22,
                                                   fmt="%.1f", n_labels=5, font_family="times", color="black"))
    pl.add_text("von Mises (mPa)", position=(0.44, 0.115), viewport=True, font_size=13, font="times", color="black")
    pl.add_mesh(pv.Plane(center=(Lx / 2, Ly / 2, -0.3), direction=(0, 0, 1), i_size=1.15 * Lx, j_size=1.15 * Ly),
                color="#b8b8b8", ambient=0.4)
    pl.add_mesh(pv.Box((0, Lx, 0, Ly, 0, Lz)).extract_feature_edges(), color="#8a8a8a", line_width=1.0)
    what = "grown with spreading" if a.run == "beta" else "structure as measured, no spreading (beta = 0)"
    pl.add_text(f"{case.replace('_d', ' day ')}: film surface phi = {THR}, {what}, "
                f"{Lx:.0f} x {Ly:.0f} x {Lz:.0f} um", position="upper_left", font_size=13, font="times",
                color="black")
    pl.add_text("colour: von Mises at T* on the surface (its 95th percentile at the top of the scale)",
                position="lower_right", font_size=9, font="times", color="#555555")
    pl.camera_position = [(Lx / 2 - 1.6 * Lx, Ly / 2 - 1.4 * Ly, 2.2 * Lx * Lz / Lx if False else 1.4 * Lx),
                          (Lx / 2, Ly / 2, Lz / 3), (0, 0, 1)]
    pl.camera.zoom(1.15)
    out = ASSETS / f"anim_clsm_microscale_{case}_3d"
    ASSETS.mkdir(exist_ok=True)
    nfr = 120
    pl.open_movie(str(out.with_suffix(".mp4")), framerate=20, quality=8)
    az0 = pl.camera.azimuth
    for i in range(nfr):
        pl.camera.azimuth = az0 + 360.0 * i / nfr
        pl.write_frame()
    pl.camera.azimuth = az0 + 30.0
    pl.screenshot(str(out.with_suffix(".png")))
    pl.close()
    print("wrote", out.with_suffix(".mp4"), "iso points", iso.n_points, "vmax mPa", vmax)


def cmd_table(a):
    """Markdown tables of the results for the report."""
    out = load_json()
    mp = lambda v: "-" if v is None else f"{v * 1e3:.2f}"
    print("| stack | dz (um) | tilt (deg) | H (um) | H_top (um) | film (% of box) | beta (mm^2/T*) |")
    print("|---|---|---|---|---|---|---|")
    for case in CASES:
        rec = out["cases"].get(case)
        if not rec:
            continue
        d = rec["data"]
        b = rec["runs"].get("coarse_beta", {}).get("beta mm2/T*", float("nan"))
        print(f"| {case.replace('_d', ' day ')} | {d['dz um']:.2f} | {d['tilt deg']:.1f} | {d['H um']:.1f} | "
              f"{d['H_top um']:.1f} | {d['film fraction of box'] * 100:.0f} | {b:.2e} |")
    for bname in ("beta", "beta0"):
        print(f"\n{bname}: von Mises in mPa, real structure / uniform film with the same mean profile")
        print("| stack | mean | bottom | middle | top | p95 | p mean | max | r_S(phi) | r_S(thickness) |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for case in CASES:
            rec = out["cases"].get(case)
            if not rec or f"coarse_{bname}" not in rec["runs"]:
                continue
            r, u = rec["runs"][f"coarse_{bname}"], rec["runs"][f"uniform_{bname}"]
            cells = [f"{mp(r[k])} / {mp(u[k])}" for k in ("vM mean all", "vM mean bottom", "vM mean middle",
                                                            "vM mean top", "vM p95", "p mean all", "vM max")]
            print(f"| {case.replace('_d', ' day ')} | " + " | ".join(cells)
                  + f" | {r['corr vM phi']['spearman']:+.2f} | {r['corr vM thickness']['spearman']:+.2f} |")
    C = out.get("checks", {})
    for key in ("mesh", "mesh0"):
      if key not in C:
        continue
      print(f"\n{key}: coarse 2.88 um / full 1.44 um lateral, mPa")
      for case, rows in C[key].items():
            for lev, m in rows.items():
                print(f"| {case} {lev} | {mp(m['vM mean all'])} | {mp(m['vM mean bottom'])} | {mp(m['vM mean middle'])} "
                      f"| {mp(m['vM mean top'])} | {mp(m['vM p95'])} | {mp(m['p mean all'])} | {mp(m['vM max'])} "
                      f"| uniform {mp(m['uniform']['vM mean all'])} |")
    for key in ("subregion", "subregion0"):
      if key not in C:
        continue
      print(f"\n{key} of the 256-wide stack, mPa")
      for k, m in C[key]["runs"].items():
            print(f"| {k} | {mp(m['vM mean all'])} | {mp(m['vM mean bottom'])} | {mp(m['vM mean middle'])} "
                  f"| {mp(m['vM mean top'])} | {mp(m['vM p95'])} | {mp(m['p mean all'])} | {mp(m['vM max'])} |")
    for key in ("faces", "phimap"):
        if key in C:
            print(f"\n{key}, mPa")
            for k, m in C[key]["runs"].items():
                print(f"| {k} | {mp(m['vM mean all'])} | {mp(m['vM mean bottom'])} | {mp(m['vM mean middle'])} "
                      f"| {mp(m['vM mean top'])} | {mp(m['vM p95'])} | {mp(m['p mean all'])} | {mp(m['vM max'])} "
                      f"| {m['corr vM phi']['spearman']:+.2f} | {m['corr vM thickness']['spearman']:+.2f} |")
    if "flip" in C:
        m = C["flip"]["run"]
        print(f"\nflip: | {mp(m['vM mean all'])} | {mp(m['vM mean bottom'])} | {mp(m['vM mean middle'])} "
              f"| {mp(m['vM mean top'])} | {mp(m['vM p95'])} | {mp(m['p mean all'])} | {mp(m['vM max'])} "
              f"| uniform {mp(m['uniform']['vM mean all'])} | H {m['data']['H um']:.1f} H_top {m['data']['H_top um']:.1f}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--cases", nargs="+", default=CASES)
    r.add_argument("--level", default="coarse", choices=("coarse", "full"))
    r.add_argument("--only-beta", choices=("beta", "beta0"))
    r.add_argument("--force", action="store_true")
    c = sub.add_parser("checks")
    c.add_argument("--which", nargs="+", choices=("faces", "phimap", "flip", "mesh", "subregion", "mesh0", "subregion0"))
    c.add_argument("--mesh-cases", nargs="+", default=["DH_d6", "CH_d6"])
    c.add_argument("--force", action="store_true")
    f = sub.add_parser("figure")
    f.add_argument("--case", default="DH_d21")
    v = sub.add_parser("video")
    v.add_argument("--case", default="DH_d21")
    v.add_argument("--run", default="beta0", choices=("beta", "beta0"))
    sub.add_parser("table")
    a = ap.parse_args(argv)
    {"run": cmd_run, "checks": cmd_checks, "figure": cmd_figure, "video": cmd_video, "table": cmd_table}[a.cmd](a)


if __name__ == "__main__":
    main()

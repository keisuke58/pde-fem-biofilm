#!/usr/bin/env python3
"""Submodel of the biofilm on the real molar: the stress across the 0.5 mm film on meshes fine enough to resolve it.

The full-box meshes of real_tooth_convergence.py (16^3 to 64^3 on the 16 mm cube, h = 1 to 0.25 mm) do not converge
the stress of the film, because the film is only half an element to two elements thick on them, and a full-box mesh
with h = 0.05 mm (320^3) is out of reach. Here only the film is meshed:

  - the tooth (E = 1e4 Pa, 1000 times the film, clamped at the root) is taken as rigid: u = 0 on every node of a tooth
    element; checked once on the full box against the elastic tooth;
  - a box around the film (the bounding box of the field cells with phi > 0.01 at the final state, 1 mm margin; the
    background phi is k_alpha T = 1.1e-3 everywhere, so 1e-3 is no threshold) is cut out of the cube; outside the film
    the material is the void of the demonstration (E = 1e-3 x 10 Pa), so the crop faces are left free; checked with a
    2 mm margin and against the full box;
  - on the finest meshes the void elements (phi below a threshold) are dropped altogether and the film surface is free;
    checked against the crop with the void on a coarser mesh;
  - uniform hex8 B-bar elements (mesh_study_seed.element_matrices, small strain, nu = 0.49, eigenstrain (alpha - 1) I,
    E = 10 Pa (phi^2 + 1e-3)) aligned with the field grid: h = m h_f (block means of the field) or h = h_f / r (the
    field is piecewise constant per field cell, or trilinear from the cell centres; both are run on one mesh);
  - the tooth is classified per stress element from the signed distance to the STL at the element centre (so the tooth
    surface gets finer with the mesh, while the field keeps the 0.083 mm staircase it was computed on; the cells of
    the field's own staircase tooth that are open on the fine mesh take the field of the nearest open cell), or from
    the field's voxel tooth (check);
  - pyamg smoothed aggregation + CG (relative residual 1e-10), chunked assembly, tooth and unused dofs eliminated.

Element sets as in real_tooth_convergence.py: film = field cells with phi > 0.3 that are not tooth, band = film with
z in [-2.1, -1.1] mm (CEJ at -1.6), core = the rest, each stress element weighted by its fraction of field cells in
the set; counts, maximum and location from the majority-film elements; first layer = majority-film elements next to a
tooth element (one element thick, mesh-dependent). Added here: the layer within 0.1 mm of the tooth surface (signed
distance at the element centre), a physical region that is the same on every mesh.

    python ansys_usermat/apdl/real_tooth_submodel.py --checks            # rigid tooth, margin, void, interpolation
    python ansys_usermat/apdl/real_tooth_submodel.py --levels m4 m2 r1 r2  # the convergence levels (h = 0.333 .. 0.042)
    python ansys_usermat/apdl/real_tooth_submodel.py --figure
        -> ansys_usermat/apdl/results/2026-10-3d/real_tooth_submodel.json, assets/fig_real_tooth_submodel.pdf/.png,
           figdata/real_tooth_submodel_<level>.npz

From another script:

    from real_tooth_submodel import load_field, solve_submodel
    fld = load_field(192)             # figdata/real_tooth_conv_field_192.npz
    fld = load_field("gate")          # figdata/real_tooth_gate.npz of real_tooth_model.py (final frame, 96^3)
    r = solve_submodel(fld["phi"], fld["am1"], h=fld["h_f"], tooth=fld["tooth"], sdf=fld["sdf"], h_f=fld["h_f"])
    r["measures"]["film vM mean"], r["seconds"], r["vm"] (crop grid, nan outside the elements)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")        # other jobs share the four cores; threaded BLAS only slows the CG down

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import mesh_study_seed as S  # noqa: E402
import real_tooth_convergence as RC  # noqa: E402  (constants, load_tooth, the fits)

L, E_BIO, E_TOOTH, FLOOR, NU = RC.L, RC.E_BIO, RC.E_TOOTH, RC.FLOOR, RC.NU
Z_CEJ, BAND_HALF, FILM_PHI = RC.Z_CEJ, RC.BAND_HALF, RC.FILM_PHI
CROP_PHI, MARGIN, VOID_PHI, LAYER = 0.01, 1.0, 0.01, 0.1
FIELD_DEFAULT = 192
OUT_JSON = HERE / "results" / "2026-10-3d" / "real_tooth_submodel.json"
FIG = ROOT / "assets" / "fig_real_tooth_submodel"
CONV_JSON = RC.OUT_JSON
# the convergence levels: name -> (h / h_f, void, void threshold). The void is kept where the crop is affordable
# (h >= h_f) and dropped on the refined meshes (phi <= 0.05 left out, free film surface; checked: 0.02 % on the film mean)
LEVELS = {"m4": (4.0, "keep", VOID_PHI), "m2": (2.0, "keep", VOID_PHI), "r1": (1.0, "keep", VOID_PHI),
          "r2": (0.5, "drop", 0.05), "r3": (1 / 3, "drop", 0.1)}
H_F192 = L / 192


# ---------------------------------------------------------------- fields
def load_field(src):
    """The final state of a field run. src: a grid size (figdata/real_tooth_conv_field_<n>.npz of
    real_tooth_convergence.py), a run tag (figdata/real_tooth_<tag>.npz of real_tooth_model.py, or "base" for the
    cache of anim_real_tooth_3d.py; last frame, float16, 96^3) or a path to an npz with phi, am1 [, tooth, sdf].
    Returns dict(phi, am1, tooth, sdf, nf, h_f, source)."""
    if isinstance(src, (int, np.integer)):
        path = ROOT / "figdata" / f"real_tooth_conv_field_{int(src)}.npz"
    elif isinstance(src, str) and not src.endswith(".npz"):
        path = RC.CACHE_DEMO if src == "base" else ROOT / "figdata" / f"real_tooth_{src}.npz"
    else:
        path = Path(src)
    d = np.load(path)
    phi, am1 = d["phi"], d["am1"]
    if phi.ndim == 4:                                  # a run cache: frames x grid
        phi, am1 = phi[-1], am1[-1]
    phi, am1 = phi.astype(np.float64), am1.astype(np.float64)
    nf = phi.shape[0]
    if "sdf" in d.files:
        sdf = d["sdf"].astype(np.float64)
    else:                                              # the run caches of real_tooth_model.py have no sdf: 96^3 of the demo
        sdf = np.load(RC.CACHE_DEMO)["sdf"].astype(np.float64)
        assert sdf.shape == phi.shape, (sdf.shape, phi.shape)
    tooth = d["tooth"] if "tooth" in d.files else sdf < 0
    return dict(phi=phi, am1=am1, tooth=np.asarray(tooth, bool), sdf=sdf, nf=nf, h_f=L / nf, source=str(path))


def fill_tooth(a, tooth):
    """The field with every tooth cell replaced by the value of the nearest open cell (for the fine-mesh elements that
    lie in a staircase tooth cell but outside the STL, and for the trilinear interpolation across the staircase)."""
    from scipy.ndimage import distance_transform_edt
    ind = distance_transform_edt(tooth, return_distances=False, return_indices=True)
    return a[tuple(ind)]


def crop_box(phi, tooth, h_f, thr=CROP_PHI, margin=MARGIN, align=4):
    """Index box [i0, i1) per axis on the field grid around the cells with phi > thr (not tooth), widened by the
    margin and rounded outwards to multiples of `align` cells (so that the coarsened meshes h = 2 h_f, 4 h_f align),
    clipped to the grid."""
    nf = phi.shape[0]
    sel = (phi > thr) & ~tooth
    box = []
    for k in range(3):
        idx = np.where(sel.any(axis=tuple(j for j in range(3) if j != k)))[0]
        m = int(np.ceil(margin / h_f))
        i0 = max(0, idx.min() - m)
        i1 = min(nf, idx.max() + 1 + m)
        i0 = (i0 // align) * align
        i1 = min(nf, int(np.ceil(i1 / align)) * align)
        box.append((int(i0), int(i1)))
    return box


def box_mm(box, h_f):
    return [[-L / 2 + i0 * h_f, -L / 2 + i1 * h_f] for i0, i1 in box]


def flat_film_reference(fld):
    """Mesh-free reference from the field alone: for a flat film bonded to a rigid wall with a free surface the
    stress is local, sigma_xx = sigma_yy = -E (alpha - 1) / (1 - nu), sigma_zz = 0, so von Mises = E (alpha - 1) /
    (1 - nu) and p = 2/3 of it in every cell. Means over the film, band and core cells of the field grid (the
    curvature of the flank, the film ends and the staircase of the wall are what the FEM adds to this)."""
    phi, am1, tooth, sdf, nf = fld["phi"], fld["am1"], fld["tooth"], fld["sdf"], fld["nf"]
    sets = RC.fine_sets(nf, phi, tooth)
    s = E_BIO * (np.clip(phi, 0, 1) ** 2 + FLOOR) * am1 / (1 - NU)
    lay = (sdf >= 0) & (sdf < LAYER)
    out = {"note": "von Mises = E (alpha - 1) / (1 - nu), p = 2/3 of it, per field cell; flat film on a rigid wall"}
    for k, m in sets.items():
        out[f"{k} vM mean"] = float(s[m].mean())
        out[f"{k} p mean"] = float(2 / 3 * s[m].mean())
    out["layer vM mean"] = float(s[sets["film"] & lay].mean())
    out["layer band vM mean"] = float(s[sets["band"] & lay].mean())
    for q in (50, 95, 99):
        out[f"vM p{q}"] = float(np.percentile(s[sets["film"]], q))
    out["vM max"] = float(s[sets["film"]].max())
    return out


# ---------------------------------------------------------------- the stress mesh on the crop
class Mesh:
    """Uniform hex8 mesh on the field box [i0, i1) x ... with element size h = ratio * h_f (ratio = m or 1/r)."""

    def __init__(self, box, h_f, ratio):
        self.box, self.h_f, self.ratio = box, h_f, ratio
        self.h = ratio * h_f
        self.origin = np.array([-L / 2 + i0 * h_f for i0, _ in box])
        self.shape = tuple(int(round((i1 - i0) / ratio)) for i0, i1 in box)
        for (i0, i1), n in zip(box, self.shape):
            assert abs(n * ratio - (i1 - i0)) < 1e-9, (box, ratio)
        self.axes = [self.origin[k] + (np.arange(self.shape[k]) + 0.5) * self.h for k in range(3)]

    def centres(self):
        X, Y, Z = np.meshgrid(*self.axes, indexing="ij")
        return X, Y, Z

    def sample(self, a, mode="pc", fill=None):
        """Per-element values of a field-grid array: block mean when h >= h_f; the containing cell ("pc") or the
        trilinear interpolation from the cell centres ("tri", the tooth cells filled first) when h < h_f."""
        (i0, i1), (j0, j1), (k0, k1) = self.box
        sub = a[i0:i1, j0:j1, k0:k1]
        if self.ratio >= 1:
            m = int(round(self.ratio))
            n = self.shape
            return sub.reshape(n[0], m, n[1], m, n[2], m).mean((1, 3, 5))
        r = int(round(1 / self.ratio))
        if mode == "pc":
            return np.repeat(np.repeat(np.repeat(sub, r, 0), r, 1), r, 2)
        from scipy.interpolate import RegularGridInterpolator
        if fill is not None:
            a = fill_tooth(a, fill)
        nf = a.shape[0]
        axf = -L / 2 + (np.arange(nf) + 0.5) * self.h_f
        f = RegularGridInterpolator((axf,) * 3, a, bounds_error=False, fill_value=None)
        X, Y, Z = self.centres()
        return f(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(self.shape)

    def sdf(self):
        """Signed distance to the STL at the element centres (pyvista)."""
        import pyvista as pv
        g = pv.ImageData(dimensions=self.shape, spacing=(self.h,) * 3, origin=tuple(self.axes[k][0] for k in range(3)))
        d = np.asarray(g.compute_implicit_distance(RC.load_tooth())["implicit_distance"])
        return d.reshape(self.shape, order="F").astype(np.float64)


def _neighbours(mask):
    """Elements with a True element among their six neighbours."""
    out = np.zeros_like(mask)
    for axn in range(3):
        for sgn in (1, -1):
            nb = np.roll(mask, sgn, axn)
            edge = [slice(None)] * 3
            edge[axn] = 0 if sgn == 1 else -1
            nb[tuple(edge)] = False
            out |= nb
    return out


# ---------------------------------------------------------------- the solve
def solve_submodel(phi, am1, h, *, tooth=None, sdf=None, h_f=None, box=None, margin=MARGIN, crop_thr=CROP_PHI,
                   void="keep", void_thr=VOID_PHI, interp="pc", tooth_mode="sdf", rigid=True, tol=1e-10,
                   chunk=32768, keep_fields=False, verbose=True):
    """Stress of the film on a uniform hex8 mesh of element size h (a multiple or an integer fraction of the field
    cell h_f) on the crop box around the film (box = None: the bounding box of phi > crop_thr plus the margin;
    box = "full": the whole cube, h must then divide L). void: "keep" (E = 10 Pa (phi^2 + 1e-3) everywhere outside
    the tooth) or "drop" (elements with phi <= void_thr are left out, the film surface is free). tooth_mode: "sdf"
    (signed distance to the STL at the element centre; needs sdf of the field grid only for the layer measure) or
    "voxel" (the field's tooth, block majority or containing cell). rigid: u = 0 on every node of a tooth element;
    False: the tooth is elastic (E_TOOTH, no growth) and the bottom face of the mesh is clamped (the full-box model).
    Returns dict(measures, solver, seconds, mesh, vm, p, tooth, active, w_film, ...)."""
    import pyamg
    t_start = time.time()
    phi, am1 = np.asarray(phi, float), np.asarray(am1, float)
    nf = phi.shape[0]
    h_f = L / nf if h_f is None else h_f
    if tooth is None:
        tooth = sdf < 0 if sdf is not None else np.zeros_like(phi, bool)
    open_f = ~tooth
    ratio = h / h_f
    if box is None:
        box = crop_box(phi, tooth, h_f, crop_thr, margin, align=max(4, int(round(ratio)) if ratio >= 1 else 1))
    elif box == "full":
        box = [(0, nf)] * 3
    M = Mesh(box, h_f, ratio)
    n1, n2, n3 = M.shape
    # --- per-element data
    t0 = time.time()
    if tooth_mode == "sdf":
        d_el = M.sdf()
        tooth_e = d_el < 0
    else:
        tooth_e = M.sample(tooth.astype(float)) > 0.5
        d_el = M.sample(sdf) if sdf is not None else None
    phi_e = M.sample(phi, interp, fill=tooth)
    am1_e = M.sample(am1, interp, fill=tooth)
    if ratio < 1 and interp == "pc":                 # fine elements in a staircase tooth cell but outside the STL
        gap = ~tooth_e & (M.sample(tooth.astype(float)) > 0.5)
        if gap.any():
            phi_e = np.where(gap, M.sample(fill_tooth(phi, tooth)), phi_e)
            am1_e = np.where(gap, M.sample(fill_tooth(am1, tooth)), am1_e)
    else:
        gap = np.zeros_like(tooth_e)
    phi_e = np.clip(phi_e, 0, 1)
    active = ~tooth_e if rigid else np.ones_like(tooth_e)
    if void == "drop":
        active &= phi_e > void_thr
    # element sets in physical space (field grid), as real_tooth_convergence.py
    fine = RC.fine_sets(nf, phi, tooth)
    w = {k: M.sample(v.astype(float)) for k, v in fine.items()}
    for k in w:
        w[k][~active | tooth_e] = 0.0
    t_prep = time.time() - t0
    # --- connectivity of the active elements
    t0 = time.time()
    nn = (n1 + 1, n2 + 1, n3 + 1)
    idx = np.arange(np.prod(nn), dtype=np.int64).reshape(nn)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    ea = np.where(active.ravel())[0]
    conn = np.stack([idx[i:i + n1, j:j + n2, k:k + n3].ravel()[ea] for i, j, k in off], 1)
    ne = len(ea)
    Ee = E_BIO * (phi_e.ravel()[ea] ** 2 + FLOOR)
    eps0 = am1_e.ravel()[ea]
    if not rigid:
        te = tooth_e.ravel()[ea]
        Ee = np.where(te, E_TOOTH, Ee)
        eps0 = np.where(te, 0.0, eps0)
    used = np.zeros(np.prod(nn), bool)
    used[conn.ravel()] = True
    fixed = np.zeros(np.prod(nn), bool)
    if rigid:
        tconn = np.stack([idx[i:i + n1, j:j + n2, k:k + n3].ravel()[tooth_e.ravel()] for i, j, k in off], 1)
        fixed[tconn.ravel()] = True
        del tconn
    else:
        fixed[idx[:, :, 0].ravel()] = True             # the bottom face of the mesh clamped
    # components of active nodes without a fixed node would float: drop them (none expected)
    freen = used & ~fixed
    if freen.sum() == 0:
        raise RuntimeError("no free node")
    ne_node = np.repeat(np.arange(ne), 8)
    A = sp.coo_matrix((np.ones(8 * ne), (ne_node, conn.ravel())), shape=(ne, np.prod(nn))).tocsr()
    An = A[:, used]
    G = An.T @ An
    ncomp, lab = sp.csgraph.connected_components(G, directed=False)
    lab_full = np.full(np.prod(nn), -1); lab_full[used] = lab
    held = np.zeros(ncomp, bool); held[lab_full[fixed & used]] = True
    floating = 0
    if not held.all():
        bad_nodes = used & ~held[np.where(used, lab_full, 0)] & ~fixed
        bad_el = bad_nodes[conn].any(1)
        floating = int(bad_el.sum())
        keep = ~bad_el
        ea, conn, Ee, eps0 = ea[keep], conn[keep], Ee[keep], eps0[keep]
        ne = len(ea)
        active_r = np.zeros(active.size, bool); active_r[ea] = True
        active = active_r.reshape(active.shape)
        used[:] = False; used[conn.ravel()] = True
        freen = used & ~fixed
        for k in w:
            w[k][~active] = 0.0
    del A, An, G
    # --- assembly on the free dofs, in chunks
    Ke0, fe0, B0, D, mvec = S.element_matrices(h, NU, False)
    nfree = int(freen.sum())
    nmap = np.full(np.prod(nn), -1, np.int64); nmap[freen] = np.arange(nfree)
    dofs = (3 * conn[:, :, None] + np.arange(3)[None, None, :]).reshape(ne, 24)
    nm = nmap[conn]
    fd = np.where(nm[:, :, None] >= 0, 3 * nm[:, :, None] + np.arange(3), -1).reshape(ne, 24).astype(np.int32)
    del nm
    F = np.zeros(3 * nfree)
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
        Kc = sp.csr_matrix((val[keep], (rows[keep], cols[keep])), shape=(3 * nfree,) * 2)
        Kff = Kc if Kff is None else Kff + Kc
        del rows, cols, val, keep, Kc
    Kff.sum_duplicates()
    t_asm = time.time() - t0
    # --- solve
    t0 = time.time()
    Xn, Yn, Zn = np.meshgrid(*[M.origin[k] + np.arange(nn[k]) * h for k in range(3)], indexing="ij")
    P = np.stack([Xn.ravel()[freen], Yn.ravel()[freen], Zn.ravel()[freen]], 1)
    B = np.zeros((3 * nfree, 6))
    B[0::3, 0] = 1; B[1::3, 1] = 1; B[2::3, 2] = 1
    B[0::3, 3] = -P[:, 1]; B[1::3, 3] = P[:, 0]
    B[1::3, 4] = -P[:, 2]; B[2::3, 4] = P[:, 1]
    B[0::3, 5] = P[:, 2]; B[2::3, 5] = -P[:, 0]
    ml = pyamg.smoothed_aggregation_solver(Kff, B=B, max_coarse=1000, coarse_solver="splu")
    t_setup = time.time() - t0
    res = []
    uf = ml.solve(F, tol=tol, accel="cg", maxiter=3000, residuals=res)
    rel = float(np.linalg.norm(Kff @ uf - F) / np.linalg.norm(F))
    t_solve = time.time() - t0 - t_setup
    del ml, B, P
    # --- stresses
    u = np.zeros(3 * np.prod(nn))
    fr = np.where(freen)[0]
    u[(3 * fr[:, None] + np.arange(3)).ravel()] = uf
    vm = np.full(ne, np.nan); pr = np.full(ne, np.nan)
    for c0 in range(0, ne, 4 * chunk):
        sl = slice(c0, c0 + 4 * chunk)
        eps = (B0 @ u[dofs[sl]].T).T - eps0[sl, None] * mvec
        sig = Ee[sl, None] * (eps @ D.T)
        p = -sig[:, :3].mean(1)
        s = sig[:, :3] + p[:, None]
        vm[sl] = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
        pr[sl] = p
    VM = np.full(M.shape, np.nan); VM.ravel()[ea] = vm
    PR = np.full(M.shape, np.nan); PR.ravel()[ea] = pr
    solver = {"solver": "pyamg smoothed aggregation + cg", "dof": int(3 * nfree), "iterations": len(res),
              "relative residual": rel, "seconds prepare": t_prep, "seconds assembly": t_asm,
              "seconds setup": t_setup, "seconds solve": t_solve, "elements": int(ne), "floating dropped": floating}
    meas = measures(M, VM, PR, tooth_e, active, w, d_el)
    meas.update({"h": h, "h/h_f": ratio, "box cells": [list(b) for b in box], "box mm": box_mm(box, h_f),
                 "mesh": list(M.shape), "elements in box": int(np.prod(M.shape)), "tooth elements": int(tooth_e.sum()),
                 "active elements": int(ne), "gap elements": int(gap.sum()), "void": void, "void thr": void_thr,
                 "interp": interp, "tooth mode": tooth_mode, "rigid": rigid})
    out = {"measures": meas, "solver": solver, "seconds": time.time() - t_start, "mesh": M, "vm": VM, "p": PR,
           "tooth": tooth_e, "active": active, "w_film": w["film"], "sdf_el": d_el}
    if keep_fields:
        out.update(phi_e=phi_e, am1_e=am1_e, u=u.reshape(nn + (3,)))
    if verbose:
        r = meas
        print(f"    h = {h:.4f} mm ({ratio:g} h_f), box {M.shape}, {ne} active / {tooth_e.sum()} tooth elements, "
              f"{solver['dof']} dof: {len(res)} it, residual {rel:.1e}, prepare {t_prep:.0f} s, assembly {t_asm:.0f} s, "
              f"setup {t_setup:.0f} s, solve {t_solve:.0f} s, total {out['seconds']:.0f} s", flush=True)
        print(f"    film {r['film vM mean'] * 1e3:.4f} band {r['band vM mean'] * 1e3:.4f} core {r['core vM mean'] * 1e3:.4f} "
              f"ratio {r['band/core']:.3f} layer0.1 {(r['layer vM mean'] or 0) * 1e3:.4f} p95 {r['vM p95'] * 1e3:.3f} "
              f"p99 {r['vM p99'] * 1e3:.3f} max {r['vM max'] * 1e3:.3f} at {np.round(r['vM max at'], 2).tolist()} "
              f"pmax {r['p max'] * 1e3:.3f} pmean {r['p mean'] * 1e3:.4f} mPa", flush=True)
    return out


def measures(M, vm, p, tooth_e, active, w, d_el):
    X, Y, Z = M.centres()
    film = w["film"] > 0.5
    any_film = w["film"] > 0
    touch = _neighbours(tooth_e)
    first, outer = film & touch, film & ~touch
    h = M.h
    vmf = np.where(active, vm, np.nan)
    pf = np.where(active, p, np.nan)
    imax = np.unravel_index(np.nanargmax(np.where(film, vmf, -np.inf)), vm.shape)
    ipmax = np.unravel_index(np.nanargmax(np.where(film, pf, -np.inf)), vm.shape)
    r = {"film elements": int(film.sum()), "film elements any": int(any_film.sum()),
         "film volume mm^3": float(w["film"].sum() * h ** 3), "band volume mm^3": float(w["band"].sum() * h ** 3),
         "core volume mm^3": float(w["core"].sum() * h ** 3),
         "first layer elements": int(first.sum()), "outer elements": int(outer.sum()),
         "film vM mean": RC.wmean(np.nan_to_num(vmf), w["film"]),
         "band vM mean": RC.wmean(np.nan_to_num(vmf), w["band"]), "core vM mean": RC.wmean(np.nan_to_num(vmf), w["core"]),
         "first layer vM mean": float(vmf[first].mean()) if first.any() else None,
         "outer vM mean": float(vmf[outer].mean()) if outer.any() else None,
         "vM max": float(vmf[imax]), "vM max at": [float(X[imax]), float(Y[imax]), float(Z[imax])],
         "vM max in band": bool(Z_CEJ - BAND_HALF < Z[imax] < Z_CEJ + BAND_HALF),
         "vM max in first layer": bool(first[imax]),
         "vM p95": RC.wpercentile(vmf[any_film], w["film"][any_film], 95),
         "vM p99": RC.wpercentile(vmf[any_film], w["film"][any_film], 99),
         "vM median": RC.wpercentile(vmf[any_film], w["film"][any_film], 50),
         "p max": float(pf[ipmax]), "p max at": [float(X[ipmax]), float(Y[ipmax]), float(Z[ipmax])],
         "p mean": RC.wmean(np.nan_to_num(pf), w["film"]), "p band mean": RC.wmean(np.nan_to_num(pf), w["band"]),
         "vM max all": float(np.nanmax(vmf))}
    r["band/core"] = r["band vM mean"] / r["core vM mean"] if r["core vM mean"] else None
    r["first/outer"] = r["first layer vM mean"] / r["outer vM mean"] if r["outer vM mean"] else None
    if d_el is not None:                                # the 0.1 mm next to the tooth surface, a physical layer
        lay = (d_el >= 0) & (d_el < LAYER)
        r["layer vM mean"] = RC.wmean(np.nan_to_num(vmf), w["film"] * lay)
        r["layer band vM mean"] = RC.wmean(np.nan_to_num(vmf), w["band"] * lay)
        r["layer volume mm^3"] = float((w["film"] * lay).sum() * h ** 3)
    else:
        r["layer vM mean"] = r["layer band vM mean"] = r["layer volume mm^3"] = None
    return r


# ---------------------------------------------------------------- driver
def load_json():
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text())
    return {"model": dict(RC.load_json()["model"], **{
        "tooth": "rigid: u = 0 on every node of a tooth element (checked against the elastic tooth on the full box)",
        "crop": f"bounding box of field cells with phi > {CROP_PHI} (not tooth) + {MARGIN} mm margin, free faces",
        "layer": f"film within {LAYER} mm of the tooth surface (signed distance at the element centre)"}),
        "field": {}, "checks": {}, "levels": {}}


def save_json(out):
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=1))


def record(res, keys_extra=None):
    r = dict(res["measures"])
    r["solver"] = res["solver"]
    r["seconds"] = res["seconds"]
    return r


def save_npz(name, res):
    M = res["mesh"]
    np.savez_compressed(ROOT / "figdata" / f"real_tooth_submodel_{name}.npz",
                        vm=res["vm"].astype(np.float32), p=res["p"].astype(np.float32), tooth=res["tooth"],
                        active=res["active"], w_film=res["w_film"].astype(np.float32),
                        sdf=(res["sdf_el"] if res["sdf_el"] is not None else np.zeros(0)).astype(np.float32),
                        origin=M.origin, h=M.h, shape=np.array(M.shape))


def run_checks(field_src=FIELD_DEFAULT, force=False):
    """The assumptions, each on a coarse mesh: (1) the code against real_tooth_convergence.py (full box, elastic
    tooth, voxel tooth, block means: must agree to the digit); (2) rigid against elastic tooth on the full box;
    (3) the crop with 1 mm margin against the full box and the 2 mm margin, rigid tooth; (4) void kept against
    dropped; (5) tooth from the signed distance against the voxel tooth; (6) piecewise-constant against trilinear
    field on a refined mesh."""
    fld = load_field(field_src)
    phi, am1, tooth, sdf, h_f = fld["phi"], fld["am1"], fld["tooth"], fld["sdf"], fld["h_f"]
    out = load_json()
    out["field"] = {"source": fld["source"], "nf": fld["nf"], "h_f": h_f, "crop thr": CROP_PHI, "margin mm": MARGIN,
                    "crop box cells": crop_box(phi, tooth, h_f), "crop box mm": box_mm(crop_box(phi, tooth, h_f), h_f)}
    save_json(out)
    C = out["checks"]

    def do(name, label, fld=fld, **k):
        if name in C and not force:
            print(f"check {name} already in {OUT_JSON.name}, skipped", flush=True)
            return
        print(f"check {name}: {label}", flush=True)
        res = solve_submodel(fld["phi"], fld["am1"], tooth=fld["tooth"], sdf=fld["sdf"], h_f=fld["h_f"], **k)
        C[name] = dict(record(res), label=label)
        out = load_json()                              # reload: a level run may have written meanwhile
        out["checks"] = C
        save_json(out)
    h4 = 4 * h_f                                       # 48^3 on the full box, h = 0.333 mm
    do("full_elastic_voxel_m4", "full box, elastic tooth, voxel tooth, block means (= real_tooth_convergence 48^3)",
       h=h4, box="full", rigid=False, tooth_mode="voxel")
    do("full_rigid_voxel_m4", "full box, rigid tooth, voxel tooth", h=h4, box="full", rigid=True, tooth_mode="voxel")
    do("full_rigid_sdf_m4", "full box, rigid tooth, sdf tooth", h=h4, box="full", rigid=True, tooth_mode="sdf")
    do("crop1_rigid_sdf_m4", "crop 1 mm, rigid, sdf tooth, void kept", h=h4, margin=1.0)
    do("crop2_rigid_sdf_m4", "crop 2 mm, rigid, sdf tooth, void kept", h=h4, margin=2.0)
    do("crop1_drop_m4", "crop 1 mm, rigid, sdf tooth, void dropped (phi <= 0.01)", h=h4, margin=1.0, void="drop")
    h2 = 2 * h_f
    do("crop1_rigid_sdf_m2", "crop 1 mm, rigid, sdf tooth, void kept, h = 2 h_f", h=h2, margin=1.0)
    do("crop2_rigid_sdf_m2", "crop 2 mm, rigid, sdf tooth, void kept, h = 2 h_f", h=h2, margin=2.0)
    do("crop1_drop_m2", "crop 1 mm, rigid, void dropped (phi <= 0.01), h = 2 h_f", h=h2, margin=1.0, void="drop")
    do("crop1_drop05_m2", "crop 1 mm, rigid, void dropped (phi <= 0.05), h = 2 h_f", h=h2, margin=1.0, void="drop",
       void_thr=0.05)
    do("crop1_voxel_m2", "crop 1 mm, rigid, voxel tooth, void kept, h = 2 h_f", h=h2, margin=1.0, tooth_mode="voxel")
    do("crop1_drop_r1_pc", "crop 1 mm, rigid, sdf tooth, void dropped (phi <= 0.01), h = h_f", h=h_f, void="drop")
    # the interpolation matters only on refined meshes (h < h_f): checked at ratio 2 on the 96^3 field (h = 0.083),
    # which also gives the influence of the field grid at fixed h against the 192^3 field
    f96 = load_field(96)
    do("field96_r2_pc", "96^3 field, crop 1 mm, rigid, sdf tooth, void dropped (phi <= 0.01), h = h_f/2 = 0.083, "
       "piecewise constant", fld=f96, h=f96["h_f"] / 2, void="drop")
    do("field96_r2_tri", "96^3 field, as above, trilinear", fld=f96, h=f96["h_f"] / 2, void="drop", interp="tri")
    if "r2" in out.get("levels", {}):                  # the same on the 192^3 field, against the r2 level (2M dof)
        do("crop1_r2_tri", "192^3 field, crop 1 mm, rigid, sdf tooth, void dropped (phi <= 0.05), h = h_f/2 = 0.042, "
           "trilinear (level r2 is piecewise constant)", h=h_f / 2, void="drop", void_thr=0.05, interp="tri")
    out = load_json()
    out["checks summary"] = checks_summary(C)
    save_json(out)
    print(json.dumps(out["checks summary"], indent=1))


CMP_KEYS = ["film vM mean", "band vM mean", "core vM mean", "band/core", "layer vM mean", "vM p95", "vM p99",
            "vM max", "p mean", "p max"]


def pct(a, b):
    return None if a in (None, 0) or b is None else float((b / a - 1) * 100)


def checks_summary(C):
    pairs = [("code vs real_tooth_convergence 48^3 (field 192)", "conv:192:48", "full_elastic_voxel_m4"),
             ("rigid vs elastic tooth (full box, h = 0.333)", "full_elastic_voxel_m4", "full_rigid_voxel_m4"),
             ("sdf vs voxel tooth (full box, rigid, h = 0.333)", "full_rigid_voxel_m4", "full_rigid_sdf_m4"),
             ("crop 1 mm vs full box (rigid, sdf, h = 0.333)", "full_rigid_sdf_m4", "crop1_rigid_sdf_m4"),
             ("crop 2 mm vs 1 mm (h = 0.333)", "crop1_rigid_sdf_m4", "crop2_rigid_sdf_m4"),
             ("void dropped vs kept (h = 0.333)", "crop1_rigid_sdf_m4", "crop1_drop_m4"),
             ("crop 2 mm vs 1 mm (h = 0.167)", "crop1_rigid_sdf_m2", "crop2_rigid_sdf_m2"),
             ("void dropped (phi <= 0.01) vs kept (h = 0.167)", "crop1_rigid_sdf_m2", "crop1_drop_m2"),
             ("void dropped (phi <= 0.05) vs kept (h = 0.167)", "crop1_rigid_sdf_m2", "crop1_drop05_m2"),
             ("voxel vs sdf tooth (crop, h = 0.167)", "crop1_rigid_sdf_m2", "crop1_voxel_m2"),
             ("trilinear vs piecewise-constant field (96^3 field, h = h_f/2 = 0.083)", "field96_r2_pc",
              "field96_r2_tri"),
             ("field 96^3 vs 192^3 at h = 0.083 (void dropped, phi <= 0.01)", "crop1_drop_r1_pc", "field96_r2_pc"),
             ("void dropped (phi <= 0.01) vs kept (h = 0.083)", "level:r1", "crop1_drop_r1_pc"),
             ("trilinear vs piecewise-constant field (192^3 field, h = h_f/2 = 0.042)", "level:r2", "crop1_r2_tri")]
    conv = RC.load_json()["stress"].get("192", {}).get("48")
    levels = load_json().get("levels", {})
    S = {}
    for label, a, b in pairs:
        ra = conv if a == "conv:192:48" else levels.get(a[6:]) if a.startswith("level:") else C.get(a)
        rb = C.get(b)
        if ra is None or rb is None:
            continue
        S[label] = {k: pct(ra.get(k), rb.get(k)) for k in CMP_KEYS if k in ra and k in rb}
    return S


def run_levels(levels, field_src=FIELD_DEFAULT, void=None, force=False):
    """The convergence levels on the 1 mm crop with the rigid tooth from the signed distance; each level goes to the
    JSON and an npz as soon as it is solved and is skipped on a rerun."""
    fld = load_field(field_src)
    phi, am1, tooth, sdf, h_f = fld["phi"], fld["am1"], fld["tooth"], fld["sdf"], fld["h_f"]
    out = load_json()
    out["field"]["flat film"] = flat_film_reference(fld)
    save_json(out)
    for lev in levels:
        out = load_json()
        if lev in out["levels"] and not force:
            print(f"level {lev} already in {OUT_JSON.name}, skipped", flush=True)
            continue
        ratio, v, thr = LEVELS[lev]
        v = void or v
        print(f"level {lev}: h = {ratio * h_f:.4f} mm, void {v} (phi <= {thr})", flush=True)
        res = solve_submodel(phi, am1, ratio * h_f, tooth=tooth, sdf=sdf, h_f=h_f, void=v, void_thr=thr)
        save_npz(lev, res)
        out = load_json()
        out["levels"][lev] = record(res)
        out["fits"] = fits(out["levels"])
        save_json(out)
    print_summary(load_json())


CONV_KEYS = ["film vM mean", "band vM mean", "core vM mean", "band/core", "layer vM mean", "layer band vM mean",
             "p mean", "p band mean", "vM median", "vM p95", "vM p99"]
POWER_KEYS = ["vM max", "vM p99", "vM p95", "p max", "first layer vM mean"]


def fits(levels):
    rows = sorted(levels.values(), key=lambda r: -r["h"])
    hs = [r["h"] for r in rows]
    f = {"h": hs, "levels": [k for k, _ in sorted(levels.items(), key=lambda kv: -kv[1]["h"])]}
    if len(rows) < 3:
        return f
    for k in CONV_KEYS:
        v = [r.get(k) for r in rows]
        if any(x is None for x in v):
            continue
        f[k] = RC.fit_converging(hs, v)
        # Richardson on the last three (ratio 2 only when the last three are m2, r1, r2 or r1, r2, r4)
        h3, v3 = hs[-3:], v[-3:]
        if abs(h3[0] / h3[1] - h3[1] / h3[2]) < 1e-6:
            f[k]["richardson last 3"] = RC.richardson(h3, v3)
    for k in POWER_KEYS:
        v = [r.get(k) for r in rows]
        if any(x is None for x in v):
            continue
        f[k + " power"] = RC.fit_power(hs, v)
    return f


def print_summary(out):
    f = out.get("fits", {})
    if "h" not in f:
        return
    print("== submodel levels", f["levels"], "h", [f"{h:.4f}" for h in f["h"]])
    flat = out.get("field", {}).get("flat film")
    if flat:
        print("   flat-film reference: " + ", ".join(f"{k} {v * 1e3:.4g} mPa" for k, v in flat.items()
                                                      if isinstance(v, float)))
    for k in CONV_KEYS:
        if k not in f:
            continue
        g = f[k]
        tag = f"order {g['free fit order']:.2f}" if g["free fit reliable"] else "order n/a (p = 1 used)"
        ri = g.get("richardson last 3", {})
        rtag = f"  richardson order {ri['order']:.2f} h->0 {ri['h->0']:.4e}" if ri.get("order") else ""
        print(f"  {k:20s} change {' / '.join(f'{e:+.2f}' for e in g['change % per refinement'])} %  {tag}"
              f"  h->0 {g['h->0']:.4e}  error % {' / '.join(f'{e:+.1f}' for e in g['error % at h->0'])}{rtag}")
    for k in POWER_KEYS:
        if k + " power" not in f:
            continue
        g = f[k + " power"]
        print(f"  {k:20s} ~ h^({-g['s']:+.2f}), pair slopes {' / '.join(f'{e:.2f}' for e in g['pair slopes'])},"
              f" change {' / '.join(f'{e:+.1f}' for e in g['change % per refinement'])} %")


# ---------------------------------------------------------------- figure
def figure():
    """(a) the measures against h on log axes, submodel solid, the full-box values of real_tooth_convergence.py
    (192^3 field, elastic tooth) dotted; (b) the von Mises stress on the finest mesh in the section z = -1.6 mm
    (the middle of the cervical band), tooth grey; (c) the vertical section y = 0 through the film thickness."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, NullLocator
    sys.path.insert(0, str(ROOT / "ansys_usermat"))
    import figstyle
    figstyle.apply(size=9)
    out = load_json()
    lev = out["levels"]
    names = sorted(lev, key=lambda k: -lev[k]["h"])
    hs = np.array([lev[k]["h"] for k in names])
    conv = RC.load_json()["stress"].get("192", {})
    cn = sorted(conv, key=int)
    ch = np.array([L / int(n) for n in cn])
    col = {"max": "#2a78d6", "p99": "#6250d6", "p95": "#1baf7a", "band": "#eb6834", "core": "#8c6a00",
           "film": "#444444", "layer": "#b03060"}
    mk = {"max": "o", "p99": "^", "p95": "v", "band": "s", "core": "D", "film": "x", "layer": "p"}
    series = (("vM max", "maximum", "max"), ("vM p99", "99th percentile", "p99"), ("vM p95", "95th percentile", "p95"),
              ("layer vM mean", "0.1 mm at the tooth, mean", "layer"), ("band vM mean", "band mean", "band"),
              ("film vM mean", "film mean", "film"), ("core vM mean", "core mean", "core"))
    fig = plt.figure(figsize=(6.3, 6.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.0], hspace=0.75, wspace=0.3,
                          left=0.1, right=0.98, top=0.96, bottom=0.08)
    ax = fig.add_subplot(gs[0, :])
    for key, lab, c in series:
        v = np.array([lev[k][key] for k in names], float) * 1e3
        ax.plot(hs, v, "-", marker=mk[c], color=col[c], lw=1.3, ms=4.5, mfc="white" if c != "film" else col[c],
                label=lab)
        if key in conv.get(cn[0], {}):
            cv = np.array([conv[n][key] for n in cn], float) * 1e3
            ax.plot(ch, cv, ":", marker=mk[c], color=col[c], lw=1.0, ms=3.5, mfc="white", alpha=0.6)
    flat = out.get("field", {}).get("flat film")
    if flat:
        ax.axhline(flat["film vM mean"] * 1e3, color=col["film"], lw=0.8, ls="--")
        ax.text(hs.min() / 1.1, flat["film vM mean"] * 1e3 * 1.06, "film mean, flat-film\nestimate from the field",
                fontsize=6.5, color=col["film"], ha="right", va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    allh = np.unique(np.concatenate([hs, ch]))
    ax.xaxis.set_major_locator(FixedLocator(allh)); ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f"{h:.3g}" for h in allh], fontsize=7.5)
    ax.set_xlim(allh.max() * 1.2, allh.min() / 1.2)
    ax.set_xlabel("stress element $h$ (mm)")
    ax.set_ylabel("von Mises in the film (mPa)")
    ax.axvspan(0.5 / 5, allh.max() * 1.2, color="#000000", alpha=0.04, lw=0)
    ax.text(np.sqrt(0.1 * allh.max()), ax.get_ylim()[0] * 1.1, "fewer than 5 elements across the film", fontsize=7,
            ha="center", va="bottom", color="#555555")
    ax.set_title("(a) measures against the element size", fontsize=9, loc="left")
    hl = ax.get_legend_handles_labels()
    h2 = [Line2D([], [], color="#333333", ls="-", lw=1.3), Line2D([], [], color="#333333", ls=":", lw=1.0, alpha=0.6)]
    l2 = ["film submodel, rigid tooth", "full box, elastic tooth"]
    ax.legend(hl[0] + h2, hl[1] + l2, fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=3, handlelength=1.8, columnspacing=1.2)
    # sections on the finest level
    fine = names[-1]
    d = np.load(ROOT / "figdata" / f"real_tooth_submodel_{fine}.npz")
    vm, tooth, wfilm = d["vm"] * 1e3, d["tooth"], d["w_film"]
    org, h, shp = d["origin"], float(d["h"]), d["shape"]
    axes = [org[k] + (np.arange(shp[k]) + 0.5) * h for k in range(3)]
    vmax = lev[fine]["vM p99"] * 1e3
    cmap = plt.get_cmap("viridis")
    grey = ListedColormap(["#d0ccc4"])

    def section(axp, img, tooth2, film2, ext, ax1, ax2, xl, yl, title, pad=0.4):
        axp.imshow(np.ma.masked_where(~tooth2.T, np.ones_like(img.T)), cmap=grey, extent=ext, origin="lower",
                   vmin=0, vmax=1, interpolation="nearest")
        im = axp.imshow(np.ma.masked_invalid(img.T), cmap=cmap, extent=ext, origin="lower", vmin=0, vmax=vmax,
                        interpolation="nearest")
        if film2.any():                                   # zoom on the film in this section
            ii, kk = np.where(film2)
            axp.set_xlim(ax1[ii.min()] - pad, ax1[ii.max()] + pad)
            axp.set_ylim(ax2[kk.min()] - pad, ax2[kk.max()] + pad)
        axp.set_xlabel(xl); axp.set_ylabel(yl)
        axp.set_title(title, fontsize=9, loc="left")
        axp.grid(False)
        axp.set_aspect("equal")
        return im
    kz = int(np.argmin(np.abs(axes[2] - Z_CEJ)))
    ext = [org[0], org[0] + shp[0] * h, org[1], org[1] + shp[1] * h]
    axb = fig.add_subplot(gs[1, 0])
    section(axb, vm[:, :, kz], tooth[:, :, kz], wfilm[:, :, kz] > 0, ext, axes[0], axes[1], "$x$ (mm)", "$y$ (mm)",
            f"(b) section $z$ = {axes[2][kz]:.2f} mm (band), $h$ = {h:.3f} mm")
    ky = int(np.argmin(np.abs(axes[1] - 0.0)))
    extc = [org[0], org[0] + shp[0] * h, org[2], org[2] + shp[2] * h]
    axc = fig.add_subplot(gs[1, 1])
    im = section(axc, vm[:, ky, :], tooth[:, ky, :], wfilm[:, ky, :] > 0, extc, axes[0], axes[2], "$x$ (mm)",
                 "$z$ (mm)", f"(c) section $y$ = {axes[1][ky]:.2f} mm")
    for zz in (Z_CEJ - BAND_HALF, Z_CEJ + BAND_HALF):
        axc.axhline(zz, color="#eb6834", lw=0.7, ls="--")
    axc.text(axc.get_xlim()[0] + 0.05, Z_CEJ, "band", color="#eb6834", fontsize=7, va="center")
    cb = fig.colorbar(im, ax=[axb, axc], orientation="vertical", fraction=0.04, pad=0.03, aspect=30)
    cb.set_label("von Mises (mPa), scale: 99th percentile", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG.with_suffix(".pdf"))
    fig.savefig(FIG.with_suffix(".png"), dpi=300)
    print("wrote", FIG.with_suffix(".pdf"), FIG.with_suffix(".png"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--checks", action="store_true")
    ap.add_argument("--levels", nargs="+", choices=list(LEVELS), help="e.g. m4 m2 r1 r2 (h = 4, 2, 1, 1/2 h_f)")
    ap.add_argument("--field", default=str(FIELD_DEFAULT), help="field grid size or run tag")
    ap.add_argument("--void", choices=["keep", "drop"], help="override the default (keep for h >= h_f, drop below)")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--refit", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    src = int(a.field) if a.field.isdigit() else a.field
    if a.checks:
        run_checks(src, force=a.force)
    if a.levels:
        run_levels(a.levels, src, void=a.void, force=a.force)
    if a.refit:
        out = load_json()
        out["field"]["flat film"] = flat_film_reference(load_field(src))
        out["fits"] = fits(out["levels"])
        out["checks summary"] = checks_summary(out["checks"])
        save_json(out)
        print_summary(out)
    if a.figure:
        figure()


if __name__ == "__main__":
    main()

"""coupled_amr.py -- prototype of a two-level block-structured adaptive grid for the coupled living-consumption
solver (coupled_living_consumption.py), 2 Nov 2026. Numerics only; the model equations are those of the main
script (Klempt et al. 2024 Eqs. 34-36 with the front term, the point model of Klempt et al. 2026 per cell, the
coupling agreed on 20 Oct 2026) and are not changed.

Question it is meant to answer: the seed von Mises stress of the "mixed" example converges slowly with the grid
because the front term keeps the biofilm edge about beta/v = 0.004 mm wide, thinner than a cell up to 192^3. Does a
base grid N^3 with a factor-2 refinement around the moving edge give the seed stress of the uniform (2N)^3 run at a
lower cost?

Design (the simplest consistent choice at every point; everything marked "prototype" is what a full version
would do better):
  grid      quarter domain x >= 0, z >= 0 of the 2 mm cube (mirror symmetry, as --fast), base cells H = 2/N,
            cell centred finite volumes. Level 1: blocks of B^3 base cells, each refined to (2B)^3 cells of H/2.
  tagging   base cells whose phi_3D differs from a neighbour by more than --tag-tol, dilated by --dilate base
            cells (box), every --retag inner steps; a block is active when it holds a tagged cell. New blocks take
            phi_3D and alpha by limited linear prolongation (minmod slopes: conservative, monotone; exact for the
            seed's step); dropped blocks leave their restricted (8-cell mean) values on the base grid.
  transport Eq. 34 (beta Laplacian, k_alpha alpha, penalty, front term with TVD MUSCL/minmod upwind, exactly the
            stencils of the main script) and Eq. 36 on both levels with classical RK4 and ONE time step for both
            levels, the fine step (no subcycling; the base grid is cheap, and it keeps the two levels in step at every
            RK stage). The fine blocks take two ghost layers per face: from the neighbouring active block, by
            mirroring at a wall or symmetry plane, or by limited linear prolongation of the base-grid stage values
            where the neighbour is unrefined. After every step the fine values are restricted (mean of the 8
            children) onto the base grid, so the base grid always carries the fine solution's cell means.
            Prototype: no flux correction (refluxing) at the coarse-fine faces: the base cell next to a fine block uses
            the restricted mean of the block, the block uses the prolonged base value. The interface is --dilate cells
            away from the edge, where phi_3D is smooth, so the inconsistency is small; Eq. 34 is not conservative
            anyway (source, penalty).
  nutrient  c and the antibiotic a are solved on the base grid only (d = 1e8 mm^2/T*, quasi-static, smooth); the
            front velocity v = r c/(k + c) n_c is computed on the base grid and interpolated (linear, central slopes)
            to the fine cells. c is solved at every inner (fine) step as in the main script.
            --c-fine (diagnostic): c on the whole (2N)^3 grid, v from it at the fine spacing. --c-composite (4 Nov):
            c and a on the composite grid (fine cells in the active blocks, base cells elsewhere; class Composite),
            one sparse matrix, conservative coarse-fine flux matching with ghost values interpolated from the coarse
            neighbour, BiCGSTAB warm started with a two-grid preconditioner (Gauss-Seidel on the composite grid, one
            AMG V-cycle on the base grid, so nothing expensive is rebuilt when the blocks change); v on the fine
            cells from the composite c (ghosts as for phi, central slopes). The base grid carries the restricted
            means of the fine c, so the point model (still per base cell), the coarse velocity and the histories see
            the composite solution. Quasi-static (as the main runs); --transient-c keeps dc/dt there too (backward
            Euler on the composite vector). Not done: no flux correction of phi at the coarse-fine faces, no
            subcycling, point model per base cell.
  point     the point model runs per BASE cell on the restricted phi_3D (colonisation at phi_3D >= phi_min of the
  model     base cell, inheritance from base neighbours as --inherit-living), so the living fraction, the
            consumption and the composition have the base resolution. Matching the composition was not the aim.
  stress    at T* = 1 phi_3D and alpha - 1 are written on the equivalent uniform (2N)^3 grid (fine blocks where
            active, limited linear prolongation of the base grid elsewhere), unfolded to the full cube, averaged per
            element of the 16^3 stress mesh and solved with wp3_stress_reference.solve_fields, as the main script.
            Because N is a multiple of 16 and the prolongation is conservative, these element means equal those of
            the base grid.
  check     --no-refine gives the uniform N^3 scheme of the main script (same stencils, same time step) and must
            reproduce its results.

    python ansys_usermat/apdl/coupled_amr.py --nf 32 --example mixed_ab --mix-equal [--no-refine]
        [--block 4 --dilate 2 --retag 1 --tag-tol 0.01] [--dtc 0.025] [--transient-c 1e8] [--c-fine | --c-composite]
    -> results/2026-11-amr/amr_<tag>.json (+ .npz), python coupled_amr.py --table prints the comparison
Model options have the values of the main runs (--front 1 --gc 1 --gab 1 --alpha0 100 --tab 0.333333 --phi-min
0.001 --tvd --wall-grad mirror --v2 --ansys-bc --inherit-living --fortran) and can be changed as there.
"""
import argparse
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
from scipy.ndimage import maximum_filter  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (HERE, ROOT / "ansys_usermat", ROOT / "ansys_usermat/coupling"):
    sys.path.insert(0, str(p))
import wp3_reference as W  # noqa: E402
import wp3_stress_reference as S  # noqa: E402

RES = HERE / "results" / "2026-11-amr"


def parse(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nf", type=int, default=32, help="base grid N (multiple of 16)")
    ap.add_argument("--example", default="mixed_ab", choices=["mixed", "mixed_ab", "separate", "separate_ab"])
    ap.add_argument("--mix-equal", action="store_true", help="mixed start with equal A and B (0.275 each)")
    ap.add_argument("--no-refine", action="store_true", help="uniform N^3 (check against the main script)")
    ap.add_argument("--block", type=int, default=4, help="base cells per block side (divides N/2)")
    ap.add_argument("--dilate", type=int, default=2, help="dilation of the tagged cells, base cells")
    ap.add_argument("--retag", type=int, default=1, help="inner steps between retaggings")
    ap.add_argument("--tag-tol", type=float, default=0.01, help="tag a base cell when a neighbour differs by more")
    ap.add_argument("--dtc", type=float, default=0.025, help="outer (coupling) step")
    ap.add_argument("--tend", type=float, default=1.0)
    ap.add_argument("--phi-min", type=float, default=0.001)
    ap.add_argument("--front", type=float, default=1.0)
    ap.add_argument("--gc", type=float, default=1.0)
    ap.add_argument("--gab", type=float, default=1.0)
    ap.add_argument("--alpha0", type=float, default=100.0)
    ap.add_argument("--tab", type=float, default=0.333333)
    ap.add_argument("--transient-c", type=float, default=None, help="keep dc/dt (D in mm^2/T*), backward Euler")
    ap.add_argument("--no-fortran", action="store_true", help="(not implemented: the Fortran point model is used)")
    ap.add_argument("--tag", default="", help="extra tag for the output name")
    ap.add_argument("--snap-every", type=float, default=None, help="store the base phi_3D (and the tag mask) every T")
    ap.add_argument("--c-fine-basephi", action="store_true", help="diagnostic: as --c-fine, but c from the base-grid "
                    "phi_3D (prolonged): only the gradient of c is taken at the fine spacing")
    ap.add_argument("--c-fine", action="store_true", help="diagnostic: nutrient and front velocity on the full (2N)^3 "
                    "grid (phi_3D from the two levels, living fraction of the base cell), point model still per base cell")
    ap.add_argument("--c-composite", action="store_true", help="nutrient and antibiotic on the composite grid (fine "
                    "cells in the active blocks, base cells elsewhere, class Composite), front velocity on the fine "
                    "cells from the composite c; point model still per base cell")
    ap.add_argument("--c-composite-nolateral", action="store_true", help="diagnostic: --c-composite without the "
                    "lateral interpolation of the ghost values (symmetric matrix, CG)")
    ap.add_argument("--table", action="store_true", help="print the comparison table of results/2026-11-amr")
    return ap.parse_args(argv)


# ---------------------------------------------------------------- stencils on a stack of blocks (nb, n0, n1, n2)
def along(u, a, sl):
    """slice sl along grid axis a (stack axis a + 1)"""
    idx = [slice(None)] * u.ndim
    idx[a + 1] = sl
    return u[tuple(idx)]


def pad_mirror(u, a, g=2):
    """zero-flux (mirror) ghost cells, g layers along grid axis a"""
    w = [(0, 0)] * u.ndim
    w[a + 1] = (g, g)
    return np.pad(u, w, mode="symmetric")


def lap_axis(up, a, h):
    """second difference along a of an array padded by 2 along a"""
    return (along(up, a, slice(3, -1)) - 2.0 * along(up, a, slice(2, -2)) + along(up, a, slice(1, -3))) / (h * h)


def adv_axis(up, va, a, h):
    """v_a d u/dx_a, MUSCL face values with minmod slopes, upwind (coupled_living_consumption.upwind_adv_tvd);
    up padded by 2 along a; with mirror ghosts this equals the main script's wall treatment"""
    u1 = along(up, a, slice(1, -1))
    dp = along(up, a, slice(2, None)) - u1
    dm = u1 - along(up, a, slice(0, -2))
    sl = np.where(dp * dm > 0, np.sign(dp) * np.minimum(np.abs(dp), np.abs(dm)), 0.0)
    fr = u1 + 0.5 * sl
    dP = (along(fr, a, slice(1, -1)) - along(fr, a, slice(0, -2))) / h
    fl = u1 - 0.5 * sl
    dN = (along(fl, a, slice(2, None)) - along(fl, a, slice(1, -1))) / h
    return np.maximum(va, 0.0) * dP + np.minimum(va, 0.0) * dN


def slopes(u, limited):
    """cell slopes (differences, not divided by h) per axis of a stack array, mirror ghosts;
    limited: minmod (prolongation of phi, alpha); else central (smooth fields)"""
    out = []
    for a in range(3):
        up = pad_mirror(u, a, 1)
        dp = along(up, a, slice(2, None)) - along(up, a, slice(1, -1))
        dm = along(up, a, slice(1, -1)) - along(up, a, slice(0, -2))
        if limited:
            out.append(np.where(dp * dm > 0, np.sign(dp) * np.minimum(np.abs(dp), np.abs(dm)), 0.0))
        else:
            out.append(0.5 * (dp + dm))
    return out


def children(u, s0, s1, s2):
    """(…, n0, n1, n2) parent values and slopes -> (…, 2n0, 2n1, 2n2) children at +-h/4: conservative"""
    sh = u.shape
    out = np.empty(sh[:-3] + (sh[-3], 2, sh[-2], 2, sh[-1], 2))
    for p in (0, 1):
        for q in (0, 1):
            for r in (0, 1):
                out[..., :, p, :, q, :, r] = u + 0.25 * ((2 * p - 1) * s0 + (2 * q - 1) * s1 + (2 * r - 1) * s2)
    return out.reshape(sh[:-3] + (2 * sh[-3], 2 * sh[-2], 2 * sh[-1]))


def region_index(origins, off, size):
    """index arrays (n, size0, size1, size2) of the base cells at origins + off"""
    I = origins[:, 0, None, None, None] + off[0] + np.arange(size[0])[None, :, None, None]
    J = origins[:, 1, None, None, None] + off[1] + np.arange(size[1])[None, None, :, None]
    K = origins[:, 2, None, None, None] + off[2] + np.arange(size[2])[None, None, None, :]
    I, J, K = np.broadcast_arrays(I, J, K)
    return I, J, K


class Level1:
    """the refined blocks: fine data Fphi, Fal (nb, 2B, 2B, 2B), block bookkeeping and ghost tables"""

    def __init__(self, shape, B):
        assert all(s % B == 0 for s in shape), (shape, B)
        self.B = B
        self.NB = tuple(s // B for s in shape)
        self.active = np.zeros(self.NB, bool)
        self.bid = -np.ones(self.NB, int)
        self.nb = 0
        self.Fphi = np.zeros((0, 2 * B, 2 * B, 2 * B))
        self.Fal = np.zeros((0, 2 * B, 2 * B, 2 * B))

    def set_active(self, active, phc, alc, smooth=None):
        """new active set: keep the fine data of blocks that stay, prolong (limited) the base grid into new ones;
        smooth: {name: base field} of smooth fine fields (c, a of --c-composite) carried in self.extra[name] the
        same way, new blocks by central-slope prolongation"""
        B = self.B
        old_bid, old_phi, old_al = self.bid, self.Fphi, self.Fal
        old_extra = getattr(self, "extra", {})
        self.active = active
        self.nb = int(active.sum())
        self.bid = -np.ones(self.NB, int)
        self.bid[active] = np.arange(self.nb)
        self.bc = np.argwhere(active)                               # block coordinates (nb, 3)
        self.origins = self.bc * B
        self.CI, self.CJ, self.CK = region_index(self.origins, (0, 0, 0), (B, B, B))
        self.FI, self.FJ, self.FK = region_index(2 * self.origins, (0, 0, 0), (2 * B,) * 3)   # on the (2N)^3 grid
        self.Fphi = np.empty((self.nb, 2 * B, 2 * B, 2 * B))
        self.Fal = np.empty_like(self.Fphi)
        keep = old_bid[active]                                       # old id per new block (-1: new)
        kept = keep >= 0
        self.Fphi[kept] = old_phi[keep[kept]]
        self.Fal[kept] = old_al[keep[kept]]
        if (~kept).any():
            sel = np.flatnonzero(~kept)
            I, J, K = self.CI[sel], self.CJ[sel], self.CK[sel]
            for src, dst in ((phc, self.Fphi), (alc, self.Fal)):
                s = slopes(src[None], True)
                dst[sel] = children(src[I, J, K], s[0][0][I, J, K], s[1][0][I, J, K], s[2][0][I, J, K])
        self.extra = {}
        for name, src in (smooth or {}).items():
            dst = np.empty_like(self.Fphi)
            if kept.any() and name in old_extra:
                dst[kept] = old_extra[name][keep[kept]]
            if (~kept).any():
                sel = np.flatnonzero(~kept)
                I, J, K = self.CI[sel], self.CJ[sel], self.CK[sel]
                s = slopes(src[None], False)
                dst[sel] = children(src[I, J, K], s[0][0][I, J, K], s[1][0][I, J, K], s[2][0][I, J, K])
            self.extra[name] = dst
        # ghost tables per (axis, side): neighbour block, wall or base grid
        self.gh = {}
        for a in range(3):
            for side in (0, 1):
                nbc = self.bc.copy()
                nbc[:, a] += 1 if side else -1
                inside = (nbc[:, a] >= 0) & (nbc[:, a] < self.NB[a])
                nid = np.where(inside, self.bid[tuple(np.clip(nbc, 0, np.array(self.NB) - 1).T)], -1)
                nbr = np.flatnonzero(inside & (nid >= 0))
                wall = np.flatnonzero(~inside)
                crs = np.flatnonzero(inside & (nid < 0))
                off = [0, 0, 0]; off[a] = B if side else -1
                size = [B, B, B]; size[a] = 1
                self.gh[a, side] = dict(nbr=nbr, nid=nid[nbr], wall=wall, crs=crs,
                                        cidx=region_index(self.origins[crs], off, size))

    def ghosts(self, F, a, side, uc, sc):
        """two ghost layers of the fine stage array F on face (a, side); uc, sc: base stage array and its slopes"""
        B = self.B
        sh = [self.nb, 2 * B, 2 * B, 2 * B]; sh[a + 1] = 2
        G = np.empty(sh)
        t = self.gh[a, side]
        if t["nbr"].size:
            G[t["nbr"]] = along(F[t["nid"]], a, slice(0, 2) if side else slice(-2, None))
        if t["wall"].size:
            G[t["wall"]] = np.flip(along(F[t["wall"]], a, slice(-2, None) if side else slice(0, 2)), axis=a + 1)
        if t["crs"].size:
            I, J, K = t["cidx"]
            G[t["crs"]] = children(uc[I, J, K], sc[0][I, J, K], sc[1][I, J, K], sc[2][I, J, K])
        return G

    def padded(self, F, a, uc, sc):
        return np.concatenate([self.ghosts(F, a, 0, uc, sc), F, self.ghosts(F, a, 1, uc, sc)], axis=a + 1)

    def restrict(self, F):
        B = self.B
        return F.reshape(self.nb, B, 2, B, 2, B, 2).mean(axis=(2, 4, 6))

    def to_base(self, F, uc):
        """write the restricted block means into the base array uc (in place)"""
        if self.nb:
            uc[self.CI, self.CJ, self.CK] = self.restrict(F)

    def interp(self, uc):
        """base field -> fine blocks, linear with central slopes (smooth fields: the front velocity)"""
        s = slopes(uc[None], False)
        I, J, K = self.CI, self.CJ, self.CK
        return children(uc[I, J, K], s[0][0][I, J, K], s[1][0][I, J, K], s[2][0][I, J, K])

    def fine_cells(self):
        return self.nb * (2 * self.B) ** 3


# ---------------------------------------------------------------- macro fields on the base grid (as the main script)
def held_faces(side):
    """faces (axis, index) held at the value: side -1 the y = -1 face (nutrient), 0 every face except the symmetry
    planes x = 0, z = 0 (antibiotic)"""
    if side != 0:
        return [(1, 0 if side < 0 else -1)]
    return [(a, e) for a in range(3) for e in (0, -1) if not (a in (0, 2) and e == 0)]


class Reaction:
    """0 = lap u - q phi L u (quasi-static) or (1/D) du/dt = ..., u held on a face (ghost) or in fixed cells,
    zero flux elsewhere; quarter domain (symmetry planes x = 0, z = 0 are zero-flux)"""

    def __init__(self, shape, h, q, side, transient):
        self.shape, self.h, self.q, self.transient = shape, h, q, transient
        n3 = int(np.prod(shape))
        idx = np.arange(n3).reshape(shape)
        rows, cols, diag = [], [], np.zeros(shape)
        for ax in range(3):
            for sgn in (1, -1):
                ok = np.ones(shape, bool)
                e = [slice(None)] * 3; e[ax] = -1 if sgn == 1 else 0; ok[tuple(e)] = False
                rows.append(idx[ok]); cols.append(np.roll(idx, -sgn, ax)[ok]); diag += ok
        b = np.zeros(shape)
        for ax_, j in held_faces(side):
            sl = [slice(None)] * 3; sl[ax_] = j
            diag[tuple(sl)] += 2.0; b[tuple(sl)] += 2.0 / h ** 2
        self.A0 = sp.csr_matrix((-np.ones(sum(r.size for r in rows)) / h ** 2,
                                 (np.concatenate(rows), np.concatenate(cols))), shape=(n3, n3))
        self.d0 = diag.ravel() / h ** 2
        self.b = b.ravel()
        self.fixed = None

    def solve(self, phi, L, value=1.0, u_prev=None, dtau=None):
        if value == 0.0:
            return np.zeros(self.shape)
        m = 0.0 if (self.transient is None or u_prev is None) else 1.0 / (self.transient * dtau)
        A = (self.A0 + sp.diags(self.d0 + m + self.q * (np.clip(phi, 0, 1) * L).ravel())).tocsr()
        b = self.b.copy() + (m * u_prev.ravel() / value if m else 0.0)
        if self.fixed is not None:
            f = self.fixed.ravel()
            A = (sp.diags((~f).astype(float)) @ A + sp.diags(f.astype(float))).tocsr()
            b = np.where(f, 1.0, b)
        return (value * self._amg(A, b)).reshape(self.shape)

    def _amg(self, A, b):
        import pyamg
        f = np.zeros(A.shape[0], bool) if self.fixed is None else self.fixed.ravel()
        fr = ~f
        x = np.where(f, 1.0, 0.0)
        Aff = A[fr][:, fr].tocsr()
        rhs = b[fr] - A[fr][:, f] @ x[f] if f.any() else b[fr]
        self._n = getattr(self, "_n", 0) + 1
        if getattr(self, "_ml", None) is None or self._n % 25 == 1:
            self._ml = pyamg.smoothed_aggregation_solver(Aff, symmetry="symmetric")
        sol, info = spla.cg(Aff, rhs, x0=getattr(self, "_last", None), rtol=1e-12, maxiter=500,
                            M=self._ml.aspreconditioner())
        if info != 0:
            sol = spla.spsolve(Aff.tocsc(), rhs)
        self._last = sol
        x[fr] = sol
        return x


class Composite:
    """the problem of Reaction on the composite grid of one Level1 (option --c-composite): unknowns on the coarse
    cells outside the active blocks and on the fine cells inside, one sparse matrix over both. Flux form divided by
    H^3: sum_faces T (u_i - u_j) + q phi L V_i u_i / H^3 = b_i, T = face area / centre distance / H^3, so coarse-coarse
    1/H^2 (= Reaction), fine-fine 1/(2 H^2). Coarse-fine face: each of the four fine cells exchanges T = 1/(3 H^2)
    (centre distance 3H/4) with the coarse value interpolated to its lateral position (ghost value: central slopes
    from the two lateral coarse neighbours, dropped where one of them is not a coarse cell), and the coarse face
    flux is the sum of the four fine fluxes (the lateral terms cancel in the sum: conservative). Walls, symmetry
    planes, held faces and fixed cells as in Reaction (base: the Reaction of the base grid, whose fixed mask and
    matrix are used). Solver: Krylov (BiCGSTAB; CG when the lateral terms are off and the matrix symmetric), warm
    started from the previous solution, with a geometric two-grid preconditioner whose setup does not depend on the
    block pattern: one Gauss-Seidel sweep on the composite matrix, the residual summed onto the base grid
    (restriction = sum over the 8 children, prolongation = piecewise constant), one AMG V-cycle of the base-grid
    matrix (Reaction's operator with the restricted reaction term; hierarchy rebuilt every 25 solves), one backward
    sweep. (An AMG hierarchy of the composite matrix itself would have to be rebuilt whenever the blocks change, at
    nearly every inner step.) The composite pattern is rebuilt when the active blocks change."""

    def __init__(self, base, side, lev, lateral=True):
        self.base, self.lev, self.lateral = base, lev, lateral
        self.shape, self.H, self.q, self.transient = base.shape, base.h, base.q, base.transient
        self.faces, self.fixed_base = held_faces(side), base.fixed
        self.key, self._ml, self._n, self.nbuild = None, None, 0, 0
        self.t = dict(build=0.0, amg=0.0, solve=0.0)
        self.iters = []
        fb = self.fixed_base.ravel()
        self.frb = ~fb
        self.A0bff = base.A0[self.frb][:, self.frb].tocsr()        # the base operator on its free cells

    @staticmethod
    def up8(u):
        """(nb, B, B, B) base values of the blocks -> (nb, 2B, 2B, 2B), each repeated on its 8 children"""
        return np.repeat(np.repeat(np.repeat(u, 2, 1), 2, 2), 2, 3)

    def _build(self):
        t0 = time.time()
        lev, B, sh, H = self.lev, self.lev.B, self.shape, self.H
        T_cc, T_ff, T_cf = 1 / H ** 2, 1 / (2 * H ** 2), 1 / (3 * H ** 2)
        cov = np.repeat(np.repeat(np.repeat(lev.active, B, 0), B, 1), B, 2)    # base cells covered by a block
        self.Mc = ~cov
        nc = int(self.Mc.sum())
        idxC = -np.ones(sh, int)
        idxC[self.Mc] = np.arange(nc)
        nF = lev.fine_cells()
        idxF = (nc + np.arange(nF)).reshape((lev.nb,) + (2 * B,) * 3)
        self.nc, self.n = nc, nc + nF
        S, Lt = [], []                        # (rows, cols, vals): symmetric part, lateral interpolation terms

        def face(i, j, T):
            """flux T (u_i - u_j) into the rows of i and j"""
            i, j = i.ravel(), j.ravel()
            S.append((np.concatenate([i, j, i, j]), np.concatenate([i, j, j, i]),
                      np.concatenate([np.full(2 * i.size, T), np.full(2 * i.size, -T)])))

        for a in range(3):
            i, j = along(idxC[None], a, slice(None, -1))[0], along(idxC[None], a, slice(1, None))[0]
            ok = (i >= 0) & (j >= 0)
            face(i[ok], j[ok], T_cc)
            if lev.nb:
                face(along(idxF, a, slice(None, -1)), along(idxF, a, slice(1, None)), T_ff)
                t = lev.gh[a, 1]
                if t["nbr"].size:             # fine-fine across the face to the next active block
                    face(along(idxF[t["nbr"]], a, slice(-1, None)), along(idxF[t["nid"]], a, slice(0, 1)), T_ff)
            lat = [l for l in range(3) if l != a]
            for side in (0, 1):
                t = lev.gh[a, side]
                if not t["crs"].size:
                    continue
                I, J, K = t["cidx"]           # the coarse cells across the face (n, ., ., .), size 1 along a

                def up2(u):                   # coarse face layer -> the fine lateral resolution
                    return np.repeat(np.repeat(u, 2, lat[0] + 1), 2, lat[1] + 1)
                fl = along(idxF[t["crs"]], a, slice(0, 1) if side == 0 else slice(-1, None))
                face(up2(idxC[I, J, K]), fl, T_cf)
                if not self.lateral:
                    continue
                for l in lat:                 # ghost = u_c + sgn/8 (u_{c,l+} - u_{c,l-}), sgn = -1 / +1 for the
                    sgn = np.ones(fl.shape)   # fine cell at -H/4 / +H/4 from the coarse centre along l
                    ix = [slice(None)] * 4
                    ix[l + 1] = slice(0, None, 2)
                    sgn[tuple(ix)] = -1.0
                    ids = []
                    for d in (1, -1):
                        P = [I, J, K]
                        P[l] = P[l] + d
                        inside = (P[l] >= 0) & (P[l] < sh[l])
                        Pc = tuple(np.clip(P[m], 0, sh[m] - 1) for m in range(3))
                        ids.append(up2(np.where(inside, idxC[Pc], -1)))
                    ok = (ids[0] >= 0) & (ids[1] >= 0)
                    for idn, s in zip(ids, (-1.0, 1.0)):      # row f holds -T ghost
                        Lt.append((fl[ok], idn[ok], s * T_cf * sgn[ok] / 8))
        n = self.n
        self.S0 = sp.csr_matrix((np.concatenate([x[2] for x in S]),
                                 (np.concatenate([x[0] for x in S]), np.concatenate([x[1] for x in S]))), (n, n))
        self.Lat = sp.csr_matrix((n, n)) if not Lt else sp.csr_matrix(
            (np.concatenate([x[2] for x in Lt]), (np.concatenate([x[0] for x in Lt]), np.concatenate([x[1] for x in Lt]))),
            (n, n))
        # held faces: ghost = 2 value - u, T = area / (half spacing) / H^3: coarse 2/H^2, fine 1/H^2
        self.d0, self.b0 = np.zeros(n), np.zeros(n)
        for ax_, j in self.faces:
            ids = along(idxC[None], ax_, slice(0, 1) if j == 0 else slice(-1, None))[0].ravel()
            ids = ids[ids >= 0]
            self.d0[ids] += 2 * T_cc
            self.b0[ids] += 2 * T_cc
            t = lev.gh[ax_, 0 if j == 0 else 1]
            if t["wall"].size:
                ids = along(idxF[t["wall"]], ax_, slice(0, 1) if j == 0 else slice(-1, None)).ravel()
                self.d0[ids] += T_cc
                self.b0[ids] += T_cc
        self.vol = np.concatenate([np.ones(nc), np.full(nF, 0.125)])
        self.fixed = np.concatenate([self.fixed_base[self.Mc],
                                     self.up8(self.fixed_base[lev.CI, lev.CJ, lev.CK]).ravel()])
        # the free unknowns (fixed cells eliminated, value 1), and the transfer to the base grid
        f, fr = self.fixed, ~self.fixed
        self.fr = fr
        self.S0ff = self.S0[fr][:, fr].tocsr()
        self.Latff = self.Lat[fr][:, fr].tocsr()
        self.Afx = (self.S0 + self.Lat)[fr][:, f] @ np.ones(int(f.sum()))      # A_free,fixed * 1
        idxb = np.arange(int(np.prod(sh))).reshape(sh)
        parent = np.concatenate([idxb[self.Mc], self.up8(idxb[lev.CI, lev.CJ, lev.CK]).ravel()])
        P = sp.csr_matrix((np.ones(n), (np.arange(n), parent)), (n, idxb.size))
        self.R = P.T.tocsr()                                             # sum over the children
        self.Pf = P[fr][:, self.frb].tocsr()
        self.Rf = self.Pf.T.tocsr()
        self.key, self.nbuild = lev.active.copy(), self.nbuild + 1
        self.t["build"] += time.time() - t0

    def vector(self, ub, Ff):
        """composite vector from a base field (restricted means under the blocks) and the fine field"""
        return np.concatenate([ub[self.Mc], Ff.ravel()])

    def solve(self, phc, Fp, Lc, value=1.0, u_prev=None, dtau=None):
        """phc, Lc: base phi_3D and living fraction; Fp: fine phi_3D; u_prev = (base field, fine field) of the last
        solution (warm start; previous value of the transient option). Returns (base field with the restricted
        fine means, fine field)."""
        import pyamg
        from pyamg.relaxation.relaxation import gauss_seidel
        lev = self.lev
        if self.key is None or not np.array_equal(self.key, lev.active):
            self._build()
        if value == 0.0:
            return np.zeros(self.shape), np.zeros_like(Fp)
        t0 = time.time()
        nc, fr = self.nc, self.fr
        Lf = self.up8(Lc[lev.CI, lev.CJ, lev.CK])             # the living fraction of the base cell
        r = self.q * np.concatenate([(np.clip(phc, 0, 1) * Lc)[self.Mc], (np.clip(Fp, 0, 1) * Lf).ravel() / 8])
        x0 = None if u_prev is None else self.vector(*u_prev) / value
        m = 0.0 if (self.transient is None or u_prev is None) else 1.0 / (self.transient * dtau)
        d = r + m * self.vol                                   # reaction (and transient) diagonal, fine cells / 8
        b = self.b0 + (m * self.vol * x0 if m else 0.0)
        Aff = (self.S0ff + sp.diags((self.d0 + d)[fr]) + self.Latff).tocsr()
        rhs = b[fr] - self.Afx
        self._n += 1
        if self._ml is None or self._n % 25 == 1:             # base-grid AMG, reaction term = sum of the children
            t1 = time.time()
            Abff = (self.A0bff + sp.diags((self.base.d0 + self.R @ d)[self.frb])).tocsr()
            self._ml = pyamg.smoothed_aggregation_solver(Abff, symmetry="symmetric")
            self.t["amg"] += time.time() - t1
        Mb, Pf, Rf = self._ml.aspreconditioner(), self.Pf, self.Rf

        def two_grid(rr):
            rr = np.ascontiguousarray(rr, dtype=np.float64).ravel()
            xx = np.zeros_like(rr)
            gauss_seidel(Aff, xx, rr, iterations=1, sweep="forward")
            xx += Pf @ Mb.matvec(Rf @ (rr - Aff @ xx))
            gauss_seidel(Aff, xx, rr, iterations=1, sweep="backward")
            return xx
        it = [0]

        def cb(*_):
            it[0] += 1
        kw = dict(x0=None if x0 is None else x0[fr], rtol=1e-12, maxiter=500, callback=cb,
                  M=spla.LinearOperator(Aff.shape, matvec=two_grid, dtype=np.float64))
        sol, info = (spla.bicgstab if self.Latff.nnz else spla.cg)(Aff, rhs, **kw)
        if info != 0:                                          # stagnation near round-off: finish with GMRES
            sol, info = spla.gmres(Aff, rhs, x0=sol, rtol=1e-10, restart=50, maxiter=20, M=kw["M"])
            assert info == 0, "composite nutrient solve did not converge"
        self.iters.append(it[0])
        x = np.where(self.fixed, 1.0, 0.0)
        x[fr] = sol
        x *= value
        ub = np.zeros(self.shape)
        ub[self.Mc] = x[:nc]
        Ff = x[nc:].reshape(Fp.shape)
        lev.to_base(Ff, ub)
        self.t["solve"] += time.time() - t0
        return ub, Ff


def grad_mirror(u, h):
    """central differences with mirror ghosts (= the main script's grad with --wall-grad mirror and --fast)"""
    out = []
    for a in range(3):
        up = pad_mirror(u[None], a, 1)
        out.append((along(up, a, slice(2, None)) - along(up, a, slice(0, -2)))[0] / (2 * h))
    return out


# ---------------------------------------------------------------- point model helpers (copied from the main script)
N_SP = 2


def start_state(*phis):
    g = np.zeros(12)
    g[:len(phis)] = phis
    g[5] = 1.0 - sum(phis)
    g[6:11] = 0.999
    return g


def group_state(tot, group):
    other = [i for i in range(N_SP) if i not in group]
    ph = [1e-6] * N_SP
    for i in group:
        ph[i] = (tot - 1e-6 * len(other)) / len(group)
    return start_state(*ph)


def living(G):
    return (G[:, :N_SP] * G[:, 6:6 + N_SP]).sum(1)


def living_start(gs):
    ph = np.clip(gs[:, :N_SP], 0.0, None)
    lv = ph * gs[:, 6:6 + N_SP]
    tot = ph.sum(1, keepdims=True)
    sl = lv.sum(1, keepdims=True)
    share = np.where(sl > 1e-12, lv / np.maximum(sl, 1e-30), ph / np.maximum(tot, 1e-30))
    out = np.zeros_like(gs)
    out[:, :N_SP] = tot * share
    out[:, 6:11] = 0.999
    out[:, 11] = 0.0
    return out


def neighbours_mean(G, w):
    num = np.zeros_like(G); den = np.zeros(G.shape[:3])
    for ax in range(3):
        for sgn in (1, -1):
            wn = np.roll(w, sgn, ax); gn = np.roll(G, sgn, ax)
            e = [slice(None)] * 3; e[ax] = 0 if sgn == 1 else -1
            wn[tuple(e)] = 0.0
            num += wn[..., None] * gn; den += wn
    return num, den


# ---------------------------------------------------------------- the run
def run(o):
    import ecology_native_ctypes as FN
    NF, B = o.nf, o.block
    assert NF % 16 == 0, "N must be a multiple of 16 (seed and stress mesh)"
    refine = not o.no_refine
    BOX = W.L
    H = BOX / NF
    hf = H / 2 if refine else H
    SHAPE = (NF // 2, NF, NF // 2)
    QS = (slice(NF // 2, None), slice(None), slice(NF // 2, None))
    R_FRONT = o.front * 100.0e-3 * (2.0 / 20.0e-3)
    K_MONOD, EPS2 = 1.0, 1e-6
    T_AB = o.tab
    DT_C = o.dtc
    example, antibiotic = o.example.replace("_ab", ""), o.example.endswith("_ab")
    sep = example == "separate"
    # point model (Fortran, the file of the ANSYS runs with alpha*_0 = --alpha0)
    eco = ROOT / "abaqus_composition" / "eco_fritsch_2sp_a50.txt"
    lines = eco.read_text().splitlines()
    lines[1] = f"{lines[1].split()[0]} {o.alpha0:g}"
    eco = ROOT / "abaqus_composition" / f"eco_fritsch_2sp_a{o.alpha0:g}.txt"
    eco.write_text("\n".join(lines) + "\n")
    cst = FN.init(eco)
    assert cst["nact"] == N_SP and cst["alpha_star"] == o.alpha0
    S_LINK = 1000 * 1.0e-4 / 1.0                     # Fritsch case I: 1000 steps of 1e-4 span T* = 1
    N_SUB = max(1, math.ceil(S_LINK * DT_C / 1.0e-4 - 1e-9))
    phi_init = (0.25, 0.30)
    tot = sum(phi_init)
    MIX_START = (tot / 2, tot / 2) if o.mix_equal else phi_init
    # geometry (as the main script with --v2)
    cells8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
    seed = W.seed_mask(NF, cells8)
    ax = -BOX / 2 + (np.arange(NF) + 0.5) * H
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    zmax = np.abs(Z[seed]).max() + 1e-9
    edges = (np.abs(X) > 0.75) & (np.abs(Y) > 0.75) & (np.abs(Z) <= zmax)
    seed, edges, X, Y, Z = seed[QS], edges[QS], X[QS], Y[QS], Z[QS]
    phi = (seed | (edges if sep else False)).astype(float)
    alpha = np.ones_like(phi)
    G = np.zeros(SHAPE + (12,))
    if sep:
        G[seed] = group_state(tot, [0]); G[edges] = group_state(tot, [1])
    else:
        G[seed] = start_state(*MIX_START)
    on = phi >= o.phi_min
    nut = Reaction(SHAPE, H, o.gc, -1, o.transient_c)
    ab = Reaction(SHAPE, H, o.gab, 0, o.transient_c)
    nut.fixed = Y < -0.75                                        # --ansys-bc, --v2
    ab.fixed = (np.abs(X) > 0.75) | (np.abs(Y) > 0.75) | (np.abs(Z) > 0.75)
    L = np.where(on, living(G.reshape(-1, 12)).reshape(SHAPE), 0.0)
    c = nut.solve(phi, L, 1.0)
    a = np.zeros_like(phi)
    # time step of the main script for the spacing hf (the fine level when refining)
    dt = 0.15 * hf * hf / W.BETA
    dt = min(dt, 0.4 * hf / (R_FRONT / (1 + K_MONOD)))
    dt = min(dt, 2.2 / (12 * W.BETA / hf ** 2 + 2 * R_FRONT / (1 + K_MONOD) / hf))
    per = int(np.ceil(DT_C / dt)); dt = DT_C / per
    lev = Level1(SHAPE, B) if refine else None
    o.c_fine = o.c_fine or o.c_fine_basephi
    o.c_composite = o.c_composite or o.c_composite_nolateral
    assert not (o.c_fine and o.c_composite), "--c-fine and --c-composite exclude each other"
    if o.c_composite:                                            # nutrient and antibiotic on the composite grid
        if not refine:                                           # (an empty block level: the base grid alone)
            lev = Level1(SHAPE, B)
            lev.set_active(np.zeros(lev.NB, bool), phi, alpha, smooth=dict(c=c, a=a))
        comp = Composite(nut, -1, lev, lateral=not o.c_composite_nolateral)
        comp_a = Composite(ab, 0, lev, lateral=not o.c_composite_nolateral)
    if o.c_fine:                                                 # diagnostic: the nutrient on the (2N)^3 grid
        assert refine
        SH2 = tuple(2 * s_ for s_ in SHAPE)
        nut2 = Reaction(SH2, hf, o.gc, -1, o.transient_c)
        ax2 = -BOX / 2 + (np.arange(2 * NF) + 0.5) * hf
        Y2 = np.meshgrid(ax2, ax2, ax2, indexing="ij")[1][tuple(slice(s_, None) if i_ != 1 else slice(None)
                                                                   for i_, s_ in enumerate((NF, 0, NF)))]
        nut2.fixed = Y2 < -0.75

    def rep2N(u, limited=True):
        """the equivalent uniform (2N)^3 field: fine blocks where active, prolonged base grid elsewhere"""
        s = slopes(u[None], limited)
        u2 = children(u, s[0][0], s[1][0], s[2][0])
        return u2

    def restrict8(u2):
        n0, n1, n2 = SHAPE
        return u2.reshape(n0, 2, n1, 2, n2, 2).mean(axis=(1, 3, 5))

    def tag_blocks(phc):
        d = np.zeros(SHAPE)
        for a_ in range(3):
            up = pad_mirror(phc[None], a_, 1)
            d = np.maximum(d, np.abs(along(up, a_, slice(2, None))[0] - phc))
            d = np.maximum(d, np.abs(along(up, a_, slice(0, -2))[0] - phc))
        t = d > o.tag_tol
        if o.dilate > 0:
            t = maximum_filter(t.astype(np.uint8), size=2 * o.dilate + 1, mode="constant", cval=0) > 0
        return t.reshape(lev.NB[0], B, lev.NB[1], B, lev.NB[2], B).any(axis=(1, 3, 5))

    def rhs_base(ph, al, v):
        """Eq. 34 and 36 on the base grid (stack of one block, mirror ghosts)"""
        d = W.K_ALPHA * al - W.PENALTY * (np.maximum(ph - 1, 0) + np.minimum(ph, 0))
        for a_ in range(3):
            up = pad_mirror(ph[None], a_, 2)
            d = d + (W.BETA * lap_axis(up, a_, H) - adv_axis(up, v[a_][None], a_, H))[0]
        return d, W.K_ALPHA * ph

    def rhs_fine(Fp, Fa, phc, vf):
        d = W.K_ALPHA * Fa - W.PENALTY * (np.maximum(Fp - 1, 0) + np.minimum(Fp, 0))
        sc = [s_[0] for s_ in slopes(phc[None], True)]
        for a_ in range(3):
            up = lev.padded(Fp, a_, phc, sc)
            d = d + W.BETA * lap_axis(up, a_, hf) - adv_axis(up, vf[a_], a_, hf)
        return d, W.K_ALPHA * Fp

    def velocity(cc, h, g=None):
        """v = r c/(k + c) n_c, n_c = grad c/|grad c|; g: the gradient, else central differences at spacing h"""
        g = grad_mirror(cc, h) if g is None else g
        nrm = np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2 + EPS2)
        sp_ = R_FRONT * cc / (K_MONOD + cc)
        return [sp_ * x / nrm for x in g]

    def grad_fine(Fc, cb):
        """gradient of the composite c on the fine cells, central differences at hf; ghosts from the neighbouring
        block, by mirroring at a wall, or by central-slope prolongation of the base field cb (restricted means
        under the blocks) where the neighbour is coarse"""
        sc = [s_[0] for s_ in slopes(cb[None], False)]
        return [(along(up, a_, slice(3, -1)) - along(up, a_, slice(1, -3))) / (2 * hf)
                for a_, up in ((a_, lev.padded(Fc, a_, cb, sc)) for a_ in range(3))]

    def phi_2N():
        p2 = rep2N(phi)
        if lev.nb and not o.c_fine_basephi:
            p2[lev.FI, lev.FJ, lev.FK] = lev.Fphi
        return p2

    def shares(m):
        ph = G[m][:, :N_SP]
        return float((ph[:, 0] / ph.sum(1)).mean())

    hist = dict(t=[0.0], living_seed=[float(L[seed].mean())], share_a_seed=[shares(seed)],
                c_seed=[float(c[seed].mean())], a_seed=[0.0], colonised=[int(on.sum()) * 4],
                fine_cells=[], active_blocks=[], phi_max=[], phi_min=[])
    tm = dict(base=0.0, fine=0.0, nutrient=0.0, point=0.0, retag=0.0, restrict=0.0)
    snaps = {}
    t0 = time.time()
    nstep_total = 0
    if refine:
        lev.set_active(tag_blocks(phi), phi, alpha, smooth=dict(c=c, a=a) if o.c_composite else None)
    for k in range(int(round(o.tend / DT_C))):
        T1 = (k + 1) * DT_C
        cf, Lf = c, L
        for i_in in range(per):
            # (0) retag: the edge moves, the blocks follow (base grid carries the restricted fine values)
            if refine and nstep_total % o.retag == 0 and nstep_total > 0:
                t_ = time.time()
                lev.set_active(tag_blocks(phi), phi, alpha, smooth=dict(c=cf, a=a) if o.c_composite else None)
                tm["retag"] += time.time() - t_
            # nutrient of this inner step on the base grid, the front velocity on both levels
            t_ = time.time()
            if o.c_fine:                                          # diagnostic: c and v on the (2N)^3 grid
                L2 = np.repeat(np.repeat(np.repeat(Lf, 2, 0), 2, 1), 2, 2)     # living fraction of the base cell
                cf2 = nut2.solve(phi_2N(), L2, 1.0, u_prev=cf2 if i_in or k else None, dtau=dt)
                v2 = velocity(cf2, hf)
                cf = restrict8(cf2)
                v = [restrict8(x) for x in v2]
                vf = [x[lev.FI, lev.FJ, lev.FK] for x in v2] if lev.nb else None
            elif o.c_composite:                                   # c on the composite grid, v on both levels
                cf, lev.extra["c"] = comp.solve(phi, lev.Fphi, Lf, 1.0, u_prev=(cf, lev.extra["c"]), dtau=dt)
                v = velocity(cf, H)
                vf = velocity(lev.extra["c"], hf, grad_fine(lev.extra["c"], cf)) if (refine and lev.nb) else None
            else:
                cf = nut.solve(phi, Lf, 1.0, u_prev=cf, dtau=dt)
                v = velocity(cf, H)
                vf = [lev.interp(x) for x in v] if (refine and lev.nb) else None
            tm["nutrient"] += time.time() - t_

            def f(ph, al, Fp, Fa):
                t1 = time.time()
                dc = rhs_base(ph, al, v)
                tm["base"] += time.time() - t1
                if vf is None:
                    return dc, (None, None)
                t1 = time.time()
                df = rhs_fine(Fp, Fa, ph, vf)
                tm["fine"] += time.time() - t1
                return dc, df

            Fp, Fa = (lev.Fphi, lev.Fal) if vf is not None else (None, None)
            k1c, k1f = f(phi, alpha, Fp, Fa)
            st = lambda kc, kf, w: (phi + w * dt * kc[0], alpha + w * dt * kc[1],  # noqa: E731
                                    None if vf is None else Fp + w * dt * kf[0],
                                    None if vf is None else Fa + w * dt * kf[1])
            k2c, k2f = f(*st(k1c, k1f, 0.5))
            k3c, k3f = f(*st(k2c, k2f, 0.5))
            k4c, k4f = f(*st(k3c, k3f, 1.0))
            phi = phi + dt / 6 * (k1c[0] + 2 * k2c[0] + 2 * k3c[0] + k4c[0])
            alpha = alpha + dt / 6 * (k1c[1] + 2 * k2c[1] + 2 * k3c[1] + k4c[1])
            if vf is not None:
                t_ = time.time()
                lev.Fphi = Fp + dt / 6 * (k1f[0] + 2 * k2f[0] + 2 * k3f[0] + k4f[0])
                lev.Fal = Fa + dt / 6 * (k1f[1] + 2 * k2f[1] + 2 * k3f[1] + k4f[1])
                lev.to_base(lev.Fphi, phi); lev.to_base(lev.Fal, alpha)
                tm["restrict"] += time.time() - t_
            hist["fine_cells"].append(lev.fine_cells() if refine else 0)
            hist["active_blocks"].append(lev.nb if refine else 0)
            nstep_total += 1
        # (2) nutrient and antibiotic at t_i with L(t_{i-1}); (3) point models; (4) inheritance; (5) living
        t_ = time.time()
        if o.transient_c is not None:
            c_new = cf
        elif o.c_fine:
            c_new = restrict8(nut2.solve(phi_2N(), L2, 1.0))
        elif o.c_composite:
            c_new, lev.extra["c"] = comp.solve(phi, lev.Fphi, Lf, 1.0, u_prev=(cf, lev.extra["c"]))
        else:
            c_new = nut.solve(phi, Lf, 1.0, u_prev=c, dtau=DT_C)
        a_val = (1.0 if T1 > T_AB + 1e-12 else 0.0) if antibiotic else 0.0
        if o.c_composite:
            a_new, lev.extra["a"] = comp_a.solve(phi, lev.Fphi, Lf, a_val, u_prev=(a, lev.extra["a"]), dtau=DT_C)
        else:
            a_new = ab.solve(phi, Lf, a_val, u_prev=a, dtau=DT_C)
        tm["nutrient"] += time.time() - t_
        t_ = time.time()
        idx = np.flatnonzero(on.ravel())
        Gf = G.reshape(-1, 12)
        Gf[idx], okf = FN.step_lin(Gf[idx], S_LINK * DT_C, N_SUB, c.ravel()[idx], c_new.ravel()[idx],
                                   a.ravel()[idx], a_new.ravel()[idx])
        assert okf.all(), "Fortran point model failed"
        G = Gf.reshape(SHAPE + (12,))
        c, a = c_new, a_new
        while True:
            new = (phi >= o.phi_min) & ~on
            num, den = neighbours_mean(G, np.where(on, phi, 0.0))
            fill = new & (den > 0)
            if not fill.any():
                break
            gs = living_start(num[fill] / den[fill][:, None])
            gs[:, :N_SP] = np.clip(gs[:, :N_SP], 1e-10, None)
            gs[:, 5] = 1.0 - gs[:, :N_SP].sum(1)
            G[fill] = gs
            on = on | fill
        L = np.where(on, living(G.reshape(-1, 12)).reshape(SHAPE), 0.0)
        tm["point"] += time.time() - t_
        hist["t"].append(T1); hist["living_seed"].append(float(L[seed].mean()))
        hist["share_a_seed"].append(shares(seed))
        hist["c_seed"].append(float(c[seed].mean())); hist["a_seed"].append(float(a[seed].mean()))
        hist["colonised"].append(int(on.sum()) * 4)
        hist["phi_max"].append(float(phi.max())); hist["phi_min"].append(float(phi.min()))
        if o.snap_every and abs(T1 / o.snap_every - round(T1 / o.snap_every)) < 1e-6:
            snaps[f"phi_{T1:.3f}"] = phi.copy()
            if refine:
                snaps[f"active_{T1:.3f}"] = lev.active.copy()
        if (k + 1) % 4 == 0 or k == 0:
            nfc = hist["fine_cells"][-1]
            print(f"  T* {T1:.3f}  {time.time() - t0:6.0f} s  colonised {int(on.sum()) * 4}  "
                  f"living seed {L[seed].mean():.4f}  fine cells {nfc} ({lev.nb if refine else 0} blocks, "
                    f"{100 * nfc / (8 * phi.size):.0f} % of (2N)^3)", flush=True)
    wall = time.time() - t0
    # the equivalent uniform (2N)^3 representation and the stress on the 16^3 mesh
    if refine:
        s = slopes(phi[None], True); phi2 = children(phi, s[0][0], s[1][0], s[2][0])
        s = slopes(alpha[None], True); al2 = children(alpha, s[0][0], s[1][0], s[2][0])
        if lev.nb:
            I2, J2, K2 = region_index(2 * lev.origins, (0, 0, 0), (2 * B,) * 3)
            phi2[I2, J2, K2] = lev.Fphi; al2[I2, J2, K2] = lev.Fal
    else:
        phi2, al2 = phi, alpha

    def unfold(u):
        u = np.concatenate([u[::-1], u], 0)
        return np.concatenate([u[:, :, ::-1], u], 2)
    NM = 16
    vm, p = S.solve_fields(NM, W.block_mean(unfold(phi2), NM), W.block_mean(unfold(al2) - 1, NM))
    seed16 = S.seed_mask(NM, cells8)
    hist["seed_am1"] = float((alpha - 1)[seed].mean()); hist["seed_vm_Pa"] = float(vm[seed16].mean())
    hist["seed_vm_max_Pa"] = float(vm[seed16].max()); hist["seed_p_Pa"] = float(p[seed16].mean())
    hist["seed_am1_2N"] = float((al2 - 1)[W.seed_mask(2 * NF if refine else NF, cells8)[
        (slice(phi2.shape[0], None), slice(None), slice(phi2.shape[2], None))]].mean())
    hist["seed_state_end"] = G[seed].mean(0).tolist()
    fc = np.array(hist["fine_cells"], float)
    setting = dict(base=f"{NF}^3", refine=refine, block=B, dilate=o.dilate, retag=o.retag, tag_tol=o.tag_tol,
                   dt_inner=dt, inner_per_outer=per, dtc=DT_C, phi_min=o.phi_min, front=o.front, r=R_FRONT,
                   gc=o.gc, gab=o.gab, alpha0=o.alpha0, tab=T_AB, transient_c=o.transient_c, mix_equal=o.mix_equal,
                   example=o.example, wall_s=wall, timing_s=tm, fine_cells_mean=float(fc.mean()),
                   fine_cells_max=int(fc.max()), fine_fraction_of_2N_mean=float(fc.mean() / (8 * phi.size)),
                   base_cells_quarter=int(phi.size), uniform_2N_cells_quarter=int(8 * phi.size))
    if o.c_composite:
        setting["composite"] = dict(lateral=not o.c_composite_nolateral, timing_s=comp.t, timing_a_s=comp_a.t,
                                    iters_mean=float(np.mean(comp.iters)), iters_max=int(np.max(comp.iters)),
                                    unknowns_last=comp.n, pattern_builds=comp.nbuild, solves=comp._n)
    print(f"[{o.example}] {wall:.0f} s; seed: living {hist['living_seed'][-1]:.4f}, alpha-1 {hist['seed_am1']:.4e}, "
          f"vM {hist['seed_vm_Pa']:.4e} Pa; fine cells mean {fc.mean():.0f} (max {int(fc.max())}, "
          f"{100 * fc.mean() / (8 * phi.size):.1f} % of (2N)^3 quarter); timing {{"
          + ", ".join(f"{k_}: {v_:.0f}" for k_, v_ in tm.items()) + "}", flush=True)
    RES.mkdir(exist_ok=True)
    name = (f"amr_nf{NF}" + ("_norefine" if not refine else f"_B{B}D{o.dilate}R{o.retag}")
            + ("" if o.tag_tol == 0.01 else f"_tol{o.tag_tol:g}") + ("" if DT_C == 0.025 else f"_dtc{DT_C:g}")
            + ("" if o.transient_c is None else f"_transc{o.transient_c:g}") + ("_cfinebasephi" if o.c_fine_basephi else "_cfine" if o.c_fine else "")
            + ("_ccompnolat" if o.c_composite_nolateral else "_ccomp" if o.c_composite else "")
            + ("_mixeq" if o.mix_equal else "")
            + (f"_{o.tag}" if o.tag else "") + f"_{o.example}")
    (RES / f"{name}.json").write_text(json.dumps(dict(setting=setting, hist=hist), indent=1, default=float) + "\n")
    np.savez_compressed(RES / f"{name}.npz", phi=phi, alpha=alpha, c=c, a=a, L=L, vm16=vm, p16=p,
                        fine_cells=fc, active_blocks=np.array(hist["active_blocks"]), **snaps)
    print("wrote", RES / f"{name}.json")
    return hist


def table():
    """comparison table of the runs in results/2026-11-amr with the uniform runs of results/2026-10-3d"""
    rows = []
    for f in sorted(RES.glob("amr_*.json")):
        d = json.loads(f.read_text()); s, h = d["setting"], d["hist"]
        rows.append((f.stem, s["base"], "uniform" if not s["refine"] else f"B{s['block']} D{s['dilate']}",
                     h["living_seed"][-1], h["seed_am1"], h["seed_vm_Pa"], s["wall_s"], s["fine_cells_mean"],
                     s["fine_cells_max"]))
    print("| run | base | grid | seed living | seed alpha-1 | seed vM [Pa] | wall [s] | fine cells mean | max |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]:.4f} | {r[4]:.4e} | {r[5]:.4e} | {r[6]:.0f} | {r[7]:.0f} | {r[8]} |")
    ref = HERE / "results" / "2026-10-3d"
    pat = "coupled_living_consumption{nf}_front1_phimin0.001_v2_fortran_ansysbc_inhliv{mx}_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_{k}.json"
    print("\nuniform runs of the main script (dtc 0.025):")
    for nf in ("", "_nf64", "_nf128"):
        for mx, k in (("_mixeq", "mixed_ab"), ("", "separate_ab")):
            p = ref / pat.format(nf=nf, mx=mx, k=k)
            if p.exists():
                h = json.loads(p.read_text())["runs"][k]
                print(f"| {nf or '_nf32'} {k} | living {h['living_seed'][-1]:.4f} | alpha-1 {h['seed_am1']:.4e} | "
                      f"vM {h['seed_vm_Pa']:.4e} |")


if __name__ == "__main__":
    o = parse()
    if o.table:
        table()
    else:
        run(o)

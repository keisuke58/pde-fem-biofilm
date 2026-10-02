#!/usr/bin/env python3
"""Klempt 2024 test case 1 (Fig. 3/4): two readings of the set-up not yet tested.

klempt2024_quantitative.py solves Eq. 34-36 with the seed as an INITIAL value
and the nutrient held on an EDGE of the cube. The paper (sec. 4.1, Fig. 2,
references/Klempt2024_Hamilton_biofilm_growth_BMMB.txt) says:

  "The biofilm state variable phi is applied as a constant value of phi = 1
   to all nodes within a radius of 5 um at the center of the cube. In one of
   the corners, the value of the nutrients is fixed to be c = 1 at all times
   ... These boundary conditions are visualized in Fig. 2."

and Fig. 2 draws both node sets as boundary conditions. So two readings are
tested here, one at a time and together, against the paper's Fig. 4:

  seed = "ic"         phi = 1 on the seed nodes at t = 0 only (as before)
  seed = "dirichlet"  phi = 1 held on the seed nodes at all times
  nutrient = "edge"   c = 1 on the edge x = y = L (as before)
  nutrient = "corner" c = 1 on the corner node (L, L, L)

Everything else is klempt2024_quantitative.run: Eq. 34 as printed (transport
form), Eq. 35 quasi-static, Eq. 36, Table 2, dt = 1e-3, phi clipped to [0, 1].
Why the seed matters: as printed, the front term of Eq. 34 is
-r c/(k+c) grad(phi) . n_c, an advection of phi up the nutrient gradient. With
the seed only an initial value it moves the colony and erodes its back, so
the mean of phi barely grows; with the seed held, the back cannot erode.

    python JAXFEM/klempt2024_case1_bc.py            # four combinations
    python JAXFEM/klempt2024_case1_bc.py --seed dirichlet --nutrient band3 --growth abs

Result (2026-10-02), mean phi / mean c against Fig. 4 (paper phi 0.13 /
0.325 / 0.74 and c 0.33 / 0.165 / 0.085 at T* = 0.05 / 0.2 / 1). All runs:
Table 2, dt = 1e-3, 21^3 nodes. Largest |difference| over the curve:

  set-up                                                     phi     c
  Eq. 34/35 as printed, nutrient on one edge line (before)   0.63    0.24
    + seed held at phi = 1                                   0.60    0.26
    nutrient on one corner node                              0.67    0.33  (c -> 0)
  Eq. 34/35 as printed, nutrient on the Fig. 2 strip w=4     0.34    0.31
    no clip of phi / of c, or 2nd-order ENO advection        no change (< 0.01)
  first-order consumption g phi c (not Eq. 24/35), edge      0.30    0.10
  + growth on both faces |grad phi . n_c| (not Eq. 34), edge 0.16    0.10
  + the same, strip w = 2                                    0.08    0.30
  + the same, strip w = 1                                    0.05    0.23
  DIAGNOSTIC (not the paper): w=1, g x4                      0.27    0.03
  DIAGNOSTIC (not the paper): w=1, g x4, k = 0.25            0.07    0.09

What this establishes:
  1. Where the nutrient is held matters most. Fig. 2 draws a strip along the
     edge, not one line of nodes; with the strip the nutrient stays at the
     paper's level early on and the early phi curve follows.
  2. Eq. 34 as printed is a non-conservative advection of phi up the nutrient
     gradient: the colony is moved, not grown, so mean phi cannot rise from
     0.065 to 0.74 under any of the numerical variants (clipping, scheme,
     held seed). Growth on both faces, |grad phi . n_c|, does; it matches the
     paper's text ("growth is also enhanced" toward the opposite corner, none
     on the faces perpendicular to grad c) but not Eq. 34's sign.
  3. The paper's phi and c curves are not consistent with each other under
     Table 2's Monod constant k = 1: with the paper's c ~ 0.1 the front factor
     c/(k+c) is ~0.09, too slow for the paper's phi. Matching both needs a
     stronger consumption and a smaller k (last line), i.e. two values not in
     Table 2. The same direction was found for Fig. 7 (klempt2024_sensitivity:
     k = 0.01). A question for the authors, not something to tune here.
  4. Table 1's weak form differs from Eq. 34 (sign of the front term, no
     k_alpha alpha) and Eq. 35 (consumption g, not g phi); which form was run
     is not recoverable from the text.
So phi alone can be reproduced to 0.05 with two departures from the printed
equations (both-face growth, first-order consumption) and a strip width read
off a sketch; phi and c together are not reproduced with Table 2.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parent))
import klempt2024_quantitative as K  # noqa: E402

OUT = Path(__file__).resolve().parent / "klempt2024_results" / "case1_bc.json"


def nutrient_mask(kind):
    if kind == "edge":
        return np.isclose(K.X, K.L) & np.isclose(K.Y, K.L)
    if kind.startswith("band"):
        # Fig. 2 draws c = 1 as a strip along the edge, a few um wide on both
        # faces: nodes within w um of the edge x = y = L ("band4" -> w = 4)
        w = float(kind[4:])
        return (K.X >= K.L - w - 1e-9) & (K.Y >= K.L - w - 1e-9)
    if kind == "corner":
        return np.isclose(K.X, K.L) & np.isclose(K.Y, K.L) & np.isclose(K.Z, K.L)
    raise ValueError(kind)


class QuasiStaticC:
    """Eq. 35 quasi-static, printed (zero-order) consumption: -d lap c = -g phi,
    c = 1 on the mask. The matrix does not change, so factorise once."""

    def __init__(self, mask, g, clip_c=True):
        self.clip_c = clip_c
        self.free = ~mask.ravel()
        A = (-K.LAP).tocsr()
        self.Aff = spla.splu(A[self.free][:, self.free].tocsc())
        self.bfix = A[self.free][:, ~self.free] @ np.ones((~self.free).sum())
        self.g = g

    def __call__(self, phi):
        rhs = -(self.g / K.D) * phi.ravel()[self.free] - self.bfix
        c = np.ones(K.N ** 3)
        c[self.free] = self.Aff.solve(rhs)
        if self.clip_c:
            c = np.clip(c, 0.0, 1.0)
        return c.reshape(K.N, K.N, K.N)


class FirstOrderC:
    """Eq. 35 with consumption g phi c (first order in c), quasi-static:
    -d lap c + g phi c = 0, c = 1 on the mask. Not the paper's Eq. 24/35."""

    def __init__(self, mask, g):
        self.free = ~mask.ravel()
        A = (-K.LAP).tocsr()
        self.Lff = A[self.free][:, self.free].tocsr()
        self.Lfx = A[self.free][:, ~self.free].tocsr()
        self.nfix = int((~self.free).sum())
        self.g = g
        self.x = None

    def __call__(self, phi):
        q = (self.g / K.D) * phi.ravel()[self.free]
        Aff = self.Lff + sp.diags(q)
        rhs = -(self.Lfx @ np.ones(self.nfix))
        # the reflected-wall Laplacian is not symmetric, so no CG; BiCGSTAB from
        # the previous step's solution, checked against the residual
        x, info = spla.bicgstab(Aff, rhs, x0=self.x, rtol=1e-12, maxiter=20000)
        if info != 0:
            x = spla.spsolve(Aff.tocsc(), rhs)
        self.x = x
        c = np.ones(K.N ** 3)
        c[self.free] = x
        return np.clip(c, 0.0, 1.0).reshape(K.N, K.N, K.N)


def _minmod(a, b):
    return np.where(a * b > 0, np.where(np.abs(a) < np.abs(b), a, b), 0.0)


def eno2_dot(phi, v):
    """v . grad(phi) with second-order ENO upwind differences per axis
    (Shu & Osher): far less numerical diffusion than first-order upwind,
    whose v h / 2 ~ 25 um^2/T* swamps beta = 2 here."""
    h = K.H
    p = np.pad(phi, 2, mode="reflect")
    out = np.zeros_like(phi)
    n = phi.shape[0]
    for ax in range(3):
        def sh(k):
            sl = [slice(2, 2 + n)] * 3
            sl[ax] = slice(2 + k, 2 + k + n)
            return p[tuple(sl)]
        pm2, pm1, p0, pp1, pp2 = sh(-2), sh(-1), sh(0), sh(1), sh(2)
        dm = (p0 - pm1) / h + 0.5 * _minmod(p0 - 2 * pm1 + pm2, pp1 - 2 * p0 + pm1) / h
        dp = (pp1 - p0) / h - 0.5 * _minmod(pp1 - 2 * p0 + pm1, pp2 - 2 * pp1 + p0) / h
        out += np.where(v[ax] > 0, v[ax] * dm, v[ax] * dp)
    return out


def run(seed, nutrient, consumption="printed", growth="printed", clip="both",
        clip_c=True, scheme="upwind1", dt=K.DT, t_end=K.T_END):
    phi0, _, g = K.setup("fig4_edge")
    held = phi0 > 0.5
    mask = nutrient_mask(nutrient)
    solve = QuasiStaticC(mask, g, clip_c) if consumption == "printed" else FirstOrderC(mask, g)
    phi = phi0.copy()
    alpha = np.ones_like(phi)
    c = solve(phi)
    n = int(round(t_end / dt))
    rec = {"t": [0.0], "phi": [K.avg(phi)], "c": [K.avg(c)]}
    for s in range(1, n + 1):
        gx, gy, gz = K.grad_c(c)
        mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
        speed = np.where(mag > 1e-14, K.R * c / (K.K_M + c) / np.maximum(mag, 1e-14), 0.0)
        v = (speed * gx, speed * gy, speed * gz)

        def rhs(f):
            adv = (eno2_dot if scheme == "eno2" else K.upwind_dot)(f, v)
            src = np.abs(adv) if growth == "abs" else -adv
            return K.BETA * K.lap(f) + K.K_A * alpha + src

        if scheme == "eno2":           # Heun / SSP-RK2, nutrient frozen over the step
            f1 = phi + dt * rhs(phi)
            phi = 0.5 * (phi + f1 + dt * rhs(f1))
        else:
            phi = phi + dt * rhs(phi)
        if clip == "both":
            phi = np.clip(phi, 0.0, 1.0)
        elif clip == "lower":
            phi = np.maximum(phi, 0.0)
        if seed == "dirichlet":
            phi[held] = 1.0
        alpha = alpha + dt * K.K_A * phi
        c = solve(phi)
        if s % max(1, n // 100) == 0:
            rec["t"].append(s * dt); rec["phi"].append(K.avg(phi)); rec["c"].append(K.avg(c))
    return rec, phi, c


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", nargs="+", default=["ic", "dirichlet"])
    ap.add_argument("--nutrient", nargs="+", default=["edge", "corner"])
    ap.add_argument("--consumption", default="printed", choices=("printed", "first_order"))
    ap.add_argument("--growth", default="printed", choices=("printed", "abs"))
    ap.add_argument("--clip", default="both", choices=("both", "lower", "none"),
                    help="phi clipped to [0,1] (as before), >= 0 only, or not at all (the paper)")
    ap.add_argument("--scheme", default="upwind1", choices=("upwind1", "eno2"))
    ap.add_argument("--g-scale", type=float, default=1.0,
                    help="DIAGNOSTIC ONLY: multiply Table 2's g (not the paper)")
    ap.add_argument("--k-monod", type=float, default=1.0,
                    help="DIAGNOSTIC ONLY: Monod constant k (Table 2: 1)")
    ap.add_argument("--no-clip-c", action="store_true",
                    help="leave c unclipped (zero-order consumption can drive it below 0)")
    a = ap.parse_args(argv)
    paper = K.PAPER["fig4_edge"]
    K.G_TABLE = 1e8 * a.g_scale
    K.K_M = a.k_monod
    try:
        results = json.loads(OUT.read_text())
    except (OSError, ValueError):
        results = {}
    for nut in a.nutrient:
        for sd in a.seed:
            t0 = time.time()
            rec, phi, c = run(sd, nut, a.consumption, a.growth, a.clip, not a.no_clip_c,
                              a.scheme)
            cmp = K.compare(rec, paper)
            key = (f"seed={sd}/nutrient={nut}/consumption={a.consumption}/growth={a.growth}"
                   f"/clip={a.clip}" + ("/c=unclipped" if a.no_clip_c else "")
                   + f"/scheme={a.scheme}" + (f"/g_x{a.g_scale:g}" if a.g_scale != 1 else "")
                   + (f"/k={a.k_monod:g}" if a.k_monod != 1 else ""))
            print(f"   phi range at T*=1: [{phi.min():.3f}, {phi.max():.3f}], "
                  f"c range [{c.min():.3f}, {c.max():.3f}]")
            results[key] = {**cmp, "curve": rec}
            print(f"\n== {key}  ({time.time() - t0:.0f}s)")
            print("   t      " + "  ".join(f"{t:5.2f}" for t in paper["t"]))
            for k in ("phi", "c"):
                print(f"   {k:3s} sim " + "  ".join(f"{v:5.3f}" for v in cmp[k]["sim"]))
                print(f"   {k:3s} pap " + "  ".join(f"{v:5.3f}" for v in cmp[k]["paper"]))
                print(f"   max |diff| {k}: {cmp[k]['max_abs_diff']}", flush=True)
            OUT.write_text(json.dumps(results, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()

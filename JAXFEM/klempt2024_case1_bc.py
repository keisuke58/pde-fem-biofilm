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

Result (2026-10-02), mean phi against Fig. 4 (paper 0.13 / 0.325 / 0.74 at
T* = 0.05 / 0.2 / 1), Eq. 34 as printed unless noted:
  - seed held vs initial value: small (edge: 0.109 -> 0.141 at T* = 1).
  - nutrient on one edge line (as before) or one corner node: the nutrient
    runs out (mean c 0.07 / 0.00 at T* = 0.05 against 0.33), the front
    r c/(k+c) barely moves, mean phi stays near 0.1.
  - nutrient on the strip Fig. 2 draws (nodes within w of the edge): the
    early curve is reproduced. w = 4 um, seed held: phi 0.157 / 0.315 at
    T* = 0.05 / 0.2 (paper 0.13 / 0.325); w = 2 um gives c closest to the
    paper (0.43 -> 0.067 against 0.33 -> 0.085).
  - after T* ~ 0.3 every variant levels off (0.30-0.45 at T* = 1, paper
    0.74). Growth on both faces (|grad phi . n_c|, growth="abs", not Eq. 34
    as printed) helps late (w = 3 um: 0.446) but does not close it.
  - Remaining candidates, untested: the strip width is read off a sketch;
    the paper's implicit Galerkin FEM (no upwinding, no clip of phi) against
    explicit upwind FD here; Table 1's weak form differs from Eq. 34 in the
    sign of the front term and drops k_alpha alpha, and from Eq. 35 in the
    consumption (g, not g phi).
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

    def __init__(self, mask, g):
        self.free = ~mask.ravel()
        A = (-K.LAP).tocsr()
        self.Aff = spla.splu(A[self.free][:, self.free].tocsc())
        self.bfix = A[self.free][:, ~self.free] @ np.ones((~self.free).sum())
        self.g = g

    def __call__(self, phi):
        rhs = -(self.g / K.D) * phi.ravel()[self.free] - self.bfix
        c = np.ones(K.N ** 3)
        c[self.free] = self.Aff.solve(rhs)
        return np.clip(c, 0.0, 1.0).reshape(K.N, K.N, K.N)


class FirstOrderC:
    """Eq. 35 with consumption g phi c (first order in c), quasi-static:
    -d lap c + g phi c = 0, c = 1 on the mask. Not the paper's Eq. 24/35."""

    def __init__(self, mask, g):
        self.mask, self.g = mask, g

    def __call__(self, phi):
        return K.solve_c(phi, self.g, self.mask, "first_order")


def run(seed, nutrient, consumption="printed", growth="printed", dt=K.DT, t_end=K.T_END):
    phi0, _, g = K.setup("fig4_edge")
    held = phi0 > 0.5
    mask = nutrient_mask(nutrient)
    solve = QuasiStaticC(mask, g) if consumption == "printed" else FirstOrderC(mask, g)
    phi = phi0.copy()
    alpha = np.ones_like(phi)
    c = solve(phi)
    n = int(round(t_end / dt))
    rec = {"t": [0.0], "phi": [K.avg(phi)], "c": [K.avg(c)]}
    for s in range(1, n + 1):
        gx, gy, gz = K.grad_c(c)
        mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
        speed = np.where(mag > 1e-14, K.R * c / (K.K_M + c) / np.maximum(mag, 1e-14), 0.0)
        adv = K.upwind_dot(phi, (speed * gx, speed * gy, speed * gz))
        src = np.abs(adv) if growth == "abs" else -adv
        phi = phi + dt * (K.BETA * K.lap(phi) + K.K_A * alpha + src)
        phi = np.clip(phi, 0.0, 1.0)
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
    a = ap.parse_args(argv)
    paper = K.PAPER["fig4_edge"]
    try:
        results = json.loads(OUT.read_text())
    except (OSError, ValueError):
        results = {}
    for nut in a.nutrient:
        for sd in a.seed:
            t0 = time.time()
            rec, phi, c = run(sd, nut, a.consumption, a.growth)
            cmp = K.compare(rec, paper)
            key = f"seed={sd}/nutrient={nut}/consumption={a.consumption}/growth={a.growth}"
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

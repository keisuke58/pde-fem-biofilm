#!/usr/bin/env python3
"""Where does Klempt 2024 Fig. 7's factor of 5-10 actually live?

klempt2024_quantitative.py records the figure as not reproduced and names two
likely causes -- the one-node-thick initial disk, and the front speed of an
upwind FD scheme against the paper's Galerkin FEM -- without testing either.
This tests them, and two others, by changing one thing at a time from that
file's own baseline and reading phi(T*) off against the paper's curve.

Result (2026-09-30, fig7_high, phi at T* = 0.20; the paper is at 1.000):

    baseline, as klempt2024_quantitative runs it            0.105
    consumption 100x weaker (g = 1e6)                       0.108
    initial disk 3 nodes thick                              0.159
    initial disk 6 nodes thick                              0.238
    growth isotropic: drop the n_c projection               0.333
    isotropic and K_M = 0.01 instead of 1.0                 0.886

One recorded suspect is cleared directly and one is not tested; the gap turns
out to be two modelling choices. On the untested one: no grid-refinement study
was run here, so the upwind front speed is not ruled out by measurement -- only
crowded out, since items 3 and 4 below account for the factor on their own.

1. **Nutrient starvation is not it.** A hundredfold weaker consumption moves
   0.105 to 0.108. The colony is not short of nutrient; it is short of growth.

2. **The initial disk is not it either.** Six times the seed gives 0.238 --
   still a quarter of the paper -- and it scales with the seed rather than
   changing the rate, so no plausible initial condition closes a factor of ten.

3. **The n_c projection costs a factor of 3.2.** Eq. 34's growth carries
   n_gradphi . n_gradc. On the sides of the colony grad phi is radial while
   grad c is axial, so the product is ~0 and the colony grows upward as a
   column instead of spreading. Filling a cube needs the spreading.

4. **The Monod factor costs another 2.7.** Table 2 has K_M = 1.0 and the
   paper's own c runs from 1.0 down to a 0.49 plateau, so f = c/(K_M+c) never
   exceeds 0.5 and falls to 0.33. Setting K_M = 0.01 puts f near 1.

Together 0.105 -> 0.886 against the paper's 1.000, and the shape follows too.

**This is not a reproduction, and must not be read as one.** Both changes
contradict the paper as printed: the projection is in Eq. 34, and K_M = 1.0 is
in Table 2. Fitting two knobs until a curve matches is not evidence that the
knobs are right. What the runs establish is narrower and more useful -- the
discrepancy lives in these two terms and not in the discretisation, the seed
or the nutrient field, so either Fig. 7 was produced with something other than
the literal reading of Eq. 34 and Table 2, or c is normalised differently than
the figure's axis suggests. That is a question to put to the authors, not one
more thing to tune.

About 4 minutes per variant (21^3 nodes, 250 steps, a sparse solve each step).

    python JAXFEM/klempt2024_sensitivity.py
    python JAXFEM/klempt2024_sensitivity.py --only baseline isotropic
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import klempt2024_quantitative as K                       # noqa: E402

T_END = 0.25                      # the paper has the cube full at 0.20
PAPER_T = [0.01, 0.05, 0.10, 0.15, 0.20]
PAPER_PHI = [0.02, 0.12, 0.55, 0.90, 1.00]


def seed(nlayers):
    """The Fig. 7 disk (d = 5 um just above the bottom face), nlayers thick."""
    m = np.zeros_like(K.X)
    for k in range(1, nlayers + 1):
        m += (np.isclose(K.Z, k * K.H)
              & ((K.X - 10) ** 2 + (K.Y - 10) ** 2 <= 2.5 ** 2))
    return np.clip(m, 0.0, 1.0)


def run(phi, g, isotropic=False, k_m=None, dt=K.DT, t_end=T_END):
    """One trajectory. Everything not named here is klempt2024_quantitative's.

    A uniform c cannot be used as the "unlimited nutrient" control: growth is
    driven by n_c, so with no gradient there is no direction and the source is
    identically zero. Weak consumption is the usable stand-in, hence g = 1e6.
    """
    k_m = K.K_M if k_m is None else k_m
    mask = np.isclose(K.Z, 0.0)
    phi = phi.copy()
    alpha = np.ones_like(phi)
    c = K.solve_c(phi, g, mask, "first_order")
    ts, ps = [0.0], [K.avg(phi)]
    n = int(round(t_end / dt))
    for s in range(1, n + 1):
        f = K.R * c / (k_m + c)
        if isotropic:
            src = f * np.sqrt(sum(q ** 2 for q in K.grad_c(phi)))
        else:
            gx, gy, gz = K.grad_c(c)
            mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
            sp_ = np.where(mag > 1e-14, f / np.maximum(mag, 1e-14), 0.0)
            src = np.abs(K.upwind_dot(phi, (sp_ * gx, sp_ * gy, sp_ * gz)))
        phi = np.clip(phi + dt * (K.BETA * K.lap(phi) + K.K_A * alpha + src),
                      0.0, 1.0)
        alpha = alpha + dt * K.K_A * phi
        c = K.solve_c(phi, g, mask, "first_order")
        if s % 10 == 0:
            ts.append(s * dt)
            ps.append(K.avg(phi))
    return ts, ps


VARIANTS = {
    "baseline":        dict(nlayers=1, g=1e8, isotropic=False, k_m=None),
    "weak_consumption": dict(nlayers=1, g=1e6, isotropic=False, k_m=None),
    "disk3":           dict(nlayers=3, g=1e8, isotropic=False, k_m=None),
    "disk6":           dict(nlayers=6, g=1e8, isotropic=False, k_m=None),
    "isotropic":       dict(nlayers=1, g=1e8, isotropic=True,  k_m=None),
    "isotropic_km001": dict(nlayers=1, g=1e8, isotropic=True,  k_m=0.01),
}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="+", default=list(VARIANTS))
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("JAXFEM/klempt2024_results/sensitivity.json"))
    a = ap.parse_args(argv)

    print(f"grid N={K.N} h={K.H}  R={K.R} beta={K.BETA} K_M={K.K_M} dt={K.DT}")
    head = "  ".join(f"{t:5.2f}" for t in PAPER_T)
    print(f"{'':34s}{head}")
    print(f"{'PAPER (fig7_high)':34s}" + "  ".join(f"{v:5.3f}" for v in PAPER_PHI))

    out, t0 = {}, time.time()
    for name in a.only:
        if name not in VARIANTS:
            print(f"  unknown variant {name}")
            continue
        v = VARIANTS[name]
        ts, ps = run(seed(v["nlayers"]), v["g"], v["isotropic"], v["k_m"])
        at = np.interp(PAPER_T, ts, ps)
        print(f"{name:34s}" + "  ".join(f"{x:5.3f}" for x in at) +
              f"   (phi0 {ps[0]:.4f})", flush=True)
        out[name] = {**v, "t": PAPER_T, "phi": [round(float(x), 4) for x in at],
                     "phi0": round(float(ps[0]), 4)}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({"paper": dict(t=PAPER_T, phi=PAPER_PHI),
                                 "variants": out}, indent=1))
    print(f"wrote {a.out}  [{time.time()-t0:.0f} s]")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Local nutrient in the five-species point model (two-way step 1, Python).

As composition_local_nutrient.py (the partner's 8^3 model, 2 mm cube, seed of
32 elements, nutrient held at y = -1 mm and consumed in the seed, steady
d lap c = g phi c, Thiele number Lambda^2 = g L^2 / d scanned), but the point
model in every seed cell is the five-species model of the TMCMC calibration
(S. oralis, A. naeslundii, Veillonella, F. nucleatum, P. gingivalis; ecology_jax,
the code behind the ANSYS bridge), with c* = c*_0 c(x), c*_0 = 25 (the
calibration value, ecology_constants.py), alpha* = 0, eta_i = 1, Hill gate off.

Step 1 needs no new parameter: the MAP interaction matrices of the four
conditions are used as estimated. Step 2 (a composition-weighted front rate)
does not apply here: with eta_i = 1 and no growth-rate ratio in the
calibration, every species has the same rate.

Coupled scheme of the ANSYS runs (amount phi_cap = 0.9, rescaled every
coupling step, psi_i = 0.999 at the seed), with two choices that are
assumptions, not paper values:
  - seed shares 1/5 each (the TMCMC solver's default phi_init = 0.2);
  - time link s = 0.25: T* = 1 spans the calibration window (2500 steps of
    1e-4, the 21 days of the experiment), in 10 coupling steps of 0.1.

    python ansys_usermat/composition_local_nutrient_5sp.py --n 16
        -> assets/fig_composition_local_nutrient_5sp.png

MAP theta: tmcmc202601 data_5species/_runs/<condition>/theta_MAP.json
(theta_full, 26 Sep 2026); to be updated with the final calibration.

Result (5 Oct 2026, 8^3 grid, coupling step 0.1), range of the shares in the seed:
  Lambda = 1 (c 0.86-0.91): within 0.03 of the uniform c = 1 composition.
  Lambda = 3 (c 0.38-0.57): commensal static P. gingivalis 0.52-0.67 (0.73 at
    c = 1), A. naeslundii 0.03-0.21; dysbiotic static P. g. 0.59-0.74 (0.89).
  Lambda = 4 (c 0.23-0.44): commensal static P. g. 0.38-0.57; dysbiotic static
    0.46-0.64; commensal HOBIC A. n. 0.37-0.46 (0.56); dysbiotic HOBIC almost
    unchanged (Veillonella / P. g. 0.47-0.48 each).
The local nutrient gives the five-species composition a spatial variation of
up to about 0.3 in the shares, with no new parameter.
Caveat: at c = 1 these MAP values give P. gingivalis 0.73 in commensal static
(free point model from 0.2 each, 2500 steps of 1e-4: 0.81), not the
S. oralis-dominated commensal state. The coupled scheme is not the cause. The
end state depends strongly on the Hill gate, which is off here (decision
29 Sep): K_hill = 0.05 gives P. g. 0.71, K_hill = 0.15 gives 0.02 (A. naeslundii
0.64). Which species gains where should be read only once the calibration and
its settings (gate, initial shares, time) are final; the mechanism does not
depend on that.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import composition_local_nutrient as L  # noqa: E402
import figstyle  # noqa: E402
from ecology_constants import C_STAR  # noqa: E402

eco = L.T.eco
OUT = HERE.parent / "assets" / "fig_composition_local_nutrient_5sp.png"
SPECIES = ["S. oralis", "A. naeslundii", "Veillonella", "F. nucleatum", "P. gingivalis"]
S_5SP, PHI_CAP, NSTEP = 0.25, 0.9, 10
LAMBDAS = [1.0, 3.0, 4.0]

THETA_MAP = {   # theta_full, order a11 a12 a22 b1 b2 a33 a34 a44 b3 b4 a13 a14 a23 a24 a55 b5 a15 a25 a35 a45
    "commensal static": [1.536237, 1.388056, 1.919795, 2.150484, 2.231674, 1.029677, 1.041766, 1.623797,
                         2.754437, 0.593504, 1.84742, 1.427129, 2.694494, 1.071751, 2.621871, 1.88294,
                         0.797084, 2.015327, 1.365625, 2.788258],
    "commensal HOBIC": [1.505665, 1.754291, 1.928521, 0.908103, 2.473487, 0.874844, 1.910582, 0.918316,
                        0.997551, 2.034642, 2.259647, 0.482022, 2.265029, 0.200363, 0.886033, 1.19858,
                        0.492726, 2.623692, 2.087251, 1.66156],
    "dysbiotic static": [1.559555, 0.324685, 0.99858, 0.484752, 2.680076, 2.191778, 2.502476, 1.9431,
                         2.352477, 2.099383, 1.472655, 0.158679, 0.978861, 2.61594, 2.951558, 1.544449,
                         2.699247, 2.84245, 2.030109, 2.164549],
    "dysbiotic HOBIC": [0.520543, 1.820797, 0.744142, 2.079551, 2.415073, 0.023818, -0.44054, 2.80554,
                        4.988259, 0.891994, -0.439441, -0.253013, 1.653079, 1.420474, 0.009783, 0.190956,
                        1.564625, 1.247038, 17.335474, 4.578676],
}


def rescale(G, p):
    G = G.copy()
    G[:, :5] *= (p / G[:, :5].sum(1))[:, None]
    G[:, 5] = 1.0 - p
    return G


def shares(theta, cvals, dt=0.1):
    """Shares phi_i / sum phi_j at T* = 1 for every nutrient level in cvals."""
    dt_pm = S_5SP * dt
    n_sub = max(1, math.ceil(dt_pm / 1e-4 - 1e-9))
    th = jnp.asarray(theta, dtype=jnp.float64)
    mask = eco.active_mask(5)
    steps = jnp.arange(n_sub)
    run = jax.jit(jax.vmap(lambda g, c: eco._substep_scan(
        g, th, dt_pm / n_sub, steps, mask, c, jnp.float64(0.0), jnp.ones(5))[0]))
    cv = np.asarray(cvals, dtype=float)
    G = np.zeros((cv.size, 12)); G[:, :5] = 0.2; G[:, 6:11] = 0.999
    p = np.full(cv.size, PHI_CAP)
    cst = jnp.asarray(C_STAR * cv)
    for _ in range(int(round(1.0 / dt))):
        G = rescale(np.asarray(run(jnp.asarray(rescale(G, p)), cst)), p)
    return G[:, :5] / G[:, :5].sum(1, keepdims=True)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8, help="cells per edge (8 = the ANSYS mesh)")
    ap.add_argument("--dt", type=float, default=0.1, help="coupling step (default 0.1)")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    L.N, L.H = args.n, 2.0 / args.n
    figstyle.apply(size=11)
    S, D = L.grid_masks()
    res = {}
    for lam in LAMBDAS:
        cs = L.nutrient(S, D, lam)[S]
        print(f"Lambda = {lam}: c in the seed {cs.min():.2f} .. {cs.max():.2f}")
        # the uniform reference c = 1 first, then the seed's c
        for cond, th in THETA_MAP.items():
            sh = shares(th, np.concatenate([[1.0], cs]), args.dt)
            res[(lam, cond)] = (cs, sh[1:], sh[0])
            lo, hi = sh[1:].min(0), sh[1:].max(0)
            print(f"   {cond:17s} c=1: " + " ".join(f"{x:.3f}" for x in sh[0])
                  + " | range in the seed: " + " ".join(f"{a:.3f}-{b:.3f}" for a, b in zip(lo, hi)))
    fig, ax = plt.subplots(len(LAMBDAS), len(THETA_MAP), figsize=(14, 3.4 * len(LAMBDAS)),
                           sharex=True, sharey=True, squeeze=False)
    for r, lam in enumerate(LAMBDAS):
        for col, cond in enumerate(THETA_MAP):
            cs, sh, ref = res[(lam, cond)]
            o = np.argsort(cs)
            a = ax[r, col]
            for i in range(5):
                a.plot(cs[o], sh[o, i], ".", ms=3, color=f"C{i}", label=SPECIES[i] if r == col == 0 else None)
                a.axhline(ref[i], color=f"C{i}", lw=0.6, ls=":")
            if r == 0:
                a.set_title(cond, fontsize=11)
            if col == 0:
                a.set_ylabel(rf"$\Lambda={lam:g}$" + "\nshare $\\phi_i/\\sum_j\\phi_j$")
            if r == len(LAMBDAS) - 1:
                a.set_xlabel("local nutrient $c$ in the seed")
    fig.legend(loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.995))
    fig.suptitle(f"Local nutrient in the five-species point model (MAP of each condition), "
                 f"$c^*=25\\,c(x)$, {L.N}$^3$ grid; dotted: uniform $c=1$", y=0.955, fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


if __name__ == "__main__":
    main()

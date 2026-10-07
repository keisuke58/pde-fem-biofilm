#!/usr/bin/env python3
"""Two-way coupling, step 2, with four species (prototype in Python).

As two_way_step2.py, but the point model at every node is a four-species case
of Klempt et al. 2026 (Eq. 19-20, Table 2; klempt2026_cases.py, integrator of
klempt2026_reproduction.py, which reproduces the paper's Fig. 7-8). The front
growth rate of the field is the composition-weighted mean of species rates,

    r(x) = sum_i phi_i r_i / sum_j phi_j,   r_i = r (1/eta_i) / mean_j(1/eta_j),

with eta = (0.8, 1, 1.5, 2) from Eq. 20. That r_i scales like 1/eta_i is an
assumption (the same one that gives d = 1/3 for two species with eta = (1, 2));
for an equal mix r(x) = r. Points without biofilm keep r. Everything else as
in two_way_step2.py: Klempt 2024 test case 4.1 field, amount from the field,
phi_cap = 0.9, phi_min = 0.01, s = 0.15, coupling step 0.1. A new point starts
from the case's initial shares (case 1 equal, case 2 species 4 ten times the
others), psi_i = 0.999.

    python ansys_usermat/two_way_4sp.py [--dt 0.05] -> assets/fig_two_way_4sp.png

Result (5 Oct 2026, coupling step 0.1, 21^3), mean phi at T* = 1 and
biomass-weighted shares:
  case 1: one-way 0.704, two-way 0.744 (+6 %); shares 0.418/0.255/0.173/0.154
  case 2: one-way 0.704, two-way 0.744 (+6 %); shares 0.418/0.255/0.172/0.155
Largest pointwise change of phi 0.50 / 0.54. Both cases reach the same
composition (the four-species coexistence of the paper, species 1 largest)
within T* ~ 0.3, so r/r_bar is 1.08-1.14 everywhere in the biofilm and the
front is faster, as for case 3 of two species. With four species the
composition hardly varies in space here; a spatial pattern needs a local
input such as the nutrient (step 1). Convergence in the coupling step not
yet checked (two species: < 0.01 from 0.05 to 0.025).
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
import composition_transport_check as T  # noqa: E402
import figstyle  # noqa: E402
import klempt2026_reproduction as R  # noqa: E402
from klempt2026_cases import CASES as ALL  # noqa: E402

OUT = HERE.parent / "assets" / "fig_two_way_4sp.png"
CASES = {"4sp_case1": "case 1 (equal seed)", "4sp_case2": "case 2 (species 4 seeded 10x)"}
NS = 4


def stepper(case):
    """Advance a batch of states by s*dt in Klempt 2026 time (substeps of 1e-4)."""
    _, base = R.build(case)
    dt_pm = T.S * T.DT
    n_sub = max(1, math.ceil(dt_pm / 1e-4 - 1e-9))
    p = dict(base)
    p["dt_h"] = jnp.float64(dt_pm / n_sub)
    p["c"] = jnp.float64(case["c_star"])
    p["alpha"] = jnp.float64(case["alpha_star"])

    @jax.jit
    def run(G):
        def one(g):
            g, _ = jax.lax.scan(lambda g, _: (R.newton_step(g, p), None), g, jnp.arange(n_sub))
            return g
        return jax.vmap(one)(G)
    return lambda G: np.asarray(run(jnp.asarray(G)))


def seed_state(case, n):
    w = np.asarray(case["phi_init"], float)
    g = np.zeros((n, 12))
    g[:, :NS] = w / w.sum()
    g[:, 6:6 + NS] = R.PSI0
    return g


def rescale(g, p):
    g = g.copy()
    g[:, :NS] *= (p / g[:, :NS].sum(1))[:, None]
    g[:, 5] = 1.0 - p
    return g


def run(name, two_way):
    B, K = T.B, T.K
    s = T.S_MACRO
    case = ALL[name]
    adv = stepper(case)
    inv = 1.0 / np.asarray(case["eta"], float)
    fac_i = inv / inv.mean()                       # r_i / r
    phi, _, _ = K.setup("fig4_edge")
    n = phi.size
    G = np.zeros((n, 12)); started = np.zeros(n, bool)
    fac = np.ones(phi.shape)
    rec = {"t": [0.0], "phi": [phi.mean()], "share": [np.full(NS, np.nan)]}
    for k in range(1, int(round(1.0 / T.DT)) + 1):
        K.R = 100 * s * (fac if two_way else 1.0)
        K.BETA, K.K_A = 2 * s, 1e-3 * s
        _, phi, _ = B.run_setup(phi, B.nutrient_mask("edge"), 1e8, "ic", "first_order", "abs",
                                t_end=T.DT)
        p = np.minimum(phi.ravel(), T.PHI_CAP)
        new = (~started) & (p >= T.PHI_MIN)
        G[new] = seed_state(case, new.sum()); started |= new
        act = started & (p >= T.PHI_MIN)
        G[act] = rescale(adv(rescale(G[act], p[act])), p[act])
        sh = G[:, :NS] / np.maximum(G[:, :NS].sum(1, keepdims=True), 1e-300)
        fac = np.where(started, sh @ fac_i, 1.0).reshape(phi.shape)
        w = np.where(started, phi.ravel(), 0.0)
        rec["t"].append(k * T.DT); rec["phi"].append(phi.mean())
        rec["share"].append((w[:, None] * sh).sum(0) / max(w.sum(), 1e-300))
        print(f"  {name} {'two' if two_way else 'one'}-way T*={k * T.DT:.2f} mean phi {phi.mean():.4f}",
              flush=True)
    K.R = 100.0
    rec["share"] = np.array(rec["share"])
    return rec, phi, fac


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--dt", type=float, default=T.DT, help="coupling step (default 0.1)")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--save", help="write the results to this .npz")
    args = ap.parse_args(argv)
    T.DT = args.dt
    figstyle.apply(size=11)
    res = {(c, w): run(c, w) for c in CASES for w in (False, True)}
    if args.save:
        np.savez(args.save, **{f"{c}|{w}|{k}": np.asarray(v) for (c, w), (rec, phi, fac) in res.items()
                               for k, v in (("t", rec["t"]), ("phi_mean", rec["phi"]),
                                            ("share", rec["share"]), ("phi", phi), ("fac", fac))})
    for c in CASES:
        for w in (False, True):
            rec = res[(c, w)][0]
            sh = " ".join(f"{x:.3f}" for x in rec["share"][-1])
            print(f"{c} {'two' if w else 'one'}-way: mean phi at T*=1 {rec['phi'][-1]:.4f}, "
                  f"biomass-weighted shares {sh}")
        d = np.abs(res[(c, True)][1] - res[(c, False)][1]).max()
        print(f"   largest change of phi at T*=1: {d:.3f}")
    K = T.K
    kz = K.N // 2
    fig, ax = plt.subplots(2, 3, figsize=(13, 8.4))
    for r, c in enumerate(CASES):
        for col, w in enumerate((False, True)):
            _, phi, fac = res[(c, w)]
            a = ax[r, col]
            im = a.imshow(np.where(phi >= 0.05, fac, np.nan)[:, :, kz].T, origin="lower",
                          extent=(0, K.L, 0, K.L), cmap="viridis", vmin=0.8, vmax=1.2)
            a.contour(np.linspace(0, K.L, K.N), np.linspace(0, K.L, K.N), phi[:, :, kz].T, [0.5],
                      colors="black", linewidths=1)
            a.set_facecolor("#d9dde2"); a.grid(False)
            a.set_title("one-way" if not w else "two-way", fontsize=11)
            a.set_xticks([0, 10, 20]); a.set_yticks([0, 10, 20])
        ax[r, 0].set_ylabel(CASES[c] + "\n$y$ [$\\mu$m]")
        b = ax[r, 2]
        for w, ls in ((False, "-"), (True, "--")):
            rec = res[(c, w)][0]
            b.plot(rec["t"], rec["phi"], ls, color="k", label="mean $\\phi$" + (" (two-way)" if w else ""))
            for i in range(NS):
                b.plot(rec["t"], rec["share"][:, i], ls, color=f"C{i}",
                       label=f"$\\phi_{i + 1}/\\sum\\phi_j$" if not w else None)
        b.set_xlabel("$T^*$"); b.legend(fontsize=8, frameon=False, ncol=2)
    for a in ax[1, :2]:
        a.set_xlabel(r"$x$ [$\mu$m]")
    fig.colorbar(im, ax=ax[:, :2], label=r"$r/\bar r=\sum_i\phi_i r_i/(\bar r\sum_j\phi_j)$",
                 shrink=0.8, pad=0.02)
    fig.suptitle("Two-way coupling, step 2, four species (Klempt et al. 2026, Eq. 19-20): "
                 "$r_i\\propto 1/\\eta_i$ (assumption);\nKlempt 2024 test case 4.1 field, "
                 "mid-plane $z=10\\ \\mu$m, black: $\\phi=0.5$, grey: no biofilm", fontsize=11.5)
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


if __name__ == "__main__":
    main()

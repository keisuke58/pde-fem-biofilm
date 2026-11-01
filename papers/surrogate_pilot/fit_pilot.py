#!/usr/bin/env python3
"""Fit Gaussian processes to the pilot runs (run_pilot.py) and check them by leave-one-out.

Inputs x = (log10 alpha*_0, g l^2/d), scaled to [0, 1]. Two surrogates:
  living volume fraction of the seed at T* = 1, from all 32^3 runs (sets p and q), fitted on logit(L), L in (0, 1);
  seed von Mises stress at T* = 1, from the 64^3 runs (set f, the points of p), and for comparison from 32^3 (set p).
Kernel: constant x anisotropic RBF + white noise (scikit-learn). Writes summary.json and
assets/fig_surrogate_pilot.png: (a) living fraction over the parameter plane with the runs, (b) leave-one-out of the
living fraction, (c) leave-one-out of the 64^3 stress, (d) stress on 32^3 against 64^3 at the same points.

    python papers/surrogate_pilot/fit_pilot.py
"""
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from sklearn.gaussian_process import GaussianProcessRegressor  # noqa: E402
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

warnings.filterwarnings("ignore")
LO, HI = np.array([1.0, 0.5]), np.array([np.log10(300.0), 2.0])
EPS = 1e-3


def load(tags):
    X, L, V, names = [], [], [], []
    for t in tags:
        for p in sorted((HERE / "runs").glob(f"{t}[0-9][0-9][0-9]_*.json")):
            r = json.loads(p.read_text())["runs"]["separate_ab"]
            a0 = float(p.stem.split("_a")[1].split("_gc")[0]); gc = float(p.stem.split("_gc")[1])
            X.append([np.log10(a0), gc]); L.append(r["living_seed"][-1]); V.append(r["seed_vm_Pa"] * 1e4)
            names.append(p.stem[1:4])
    return (np.array(X) - LO) / (HI - LO), np.array(L), np.array(V), names


def gp():
    k = ConstantKernel(1.0, (1e-3, 1e3)) * RBF([0.3, 0.3], (0.02, 10.0)) + WhiteKernel(1e-6, (1e-10, 1e-1))
    return GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=5, random_state=0)


logit = lambda y: np.log(np.clip(y, EPS, 1 - EPS) / (1 - np.clip(y, EPS, 1 - EPS)))
expit = lambda z: 1 / (1 + np.exp(-z))


def loo(X, y, fwd=lambda v: v, back=lambda v: v):
    pred = np.zeros(len(X)); lo = np.zeros(len(X)); hi = np.zeros(len(X))
    for i in range(len(X)):
        m = np.arange(len(X)) != i
        g = gp().fit(X[m], fwd(y[m])); mu, s = g.predict(X[i:i + 1], return_std=True)
        pred[i], lo[i], hi[i] = back(mu[0]), back(mu[0] - 2 * s[0]), back(mu[0] + 2 * s[0])
    return pred, lo, hi


def main():
    figstyle.apply()
    Xl, L, _, _ = load("pq")                       # 32^3: living fraction
    Xp, _, V32, n32 = load("p")
    Xf, _, V64, n64 = load("f")                    # 64^3: stress
    pl, ll, hl = loo(Xl, L, logit, expit)
    pv, lv, hv = loo(Xf, V64)
    v32 = dict(zip(n32, V32)); pair = [(v32[k], v) for k, v in zip(n64, V64) if k in v32]
    summary = {
        "living_32": {"n": len(L), "loo_mean_abs": float(np.abs(pl - L).mean()), "loo_max_abs": float(np.abs(pl - L).max())},
        "vm_64": {"n": len(V64), "range_1e-4Pa": [float(V64.min()), float(V64.max())],
                  "loo_mean_abs": float(np.abs(pv - V64).mean()), "loo_max_abs": float(np.abs(pv - V64).max())},
        "vm_32_over_64": {"mean": float(np.mean([a / b for a, b in pair])), "min": float(min(a / b for a, b in pair)),
                          "max": float(max(a / b for a, b in pair))} if pair else None}
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))

    fig, axs = plt.subplots(1, 4, figsize=(8.0, 2.1), constrained_layout=True)
    g = gp().fit(Xl, logit(L))
    u = np.linspace(0, 1, 90); U, W = np.meshgrid(u, u)
    Z = expit(g.predict(np.c_[U.ravel(), W.ravel()])).reshape(U.shape)
    sc = lambda X: (10 ** (LO[0] + X[:, 0] * (HI[0] - LO[0])), LO[1] + X[:, 1] * (HI[1] - LO[1]))
    cs = axs[0].contourf(10 ** (LO[0] + U * (HI[0] - LO[0])), LO[1] + W * (HI[1] - LO[1]), Z, levels=np.linspace(0, 1, 11), cmap="Greens")
    nP = len(Xp); axs[0].scatter(*sc(Xl[:nP]), s=9, c="k", label="first 16"); axs[0].scatter(*sc(Xl[nP:]), s=9, c="#c92a2a", label="added")
    axs[0].set_xscale("log"); axs[0].set_xlabel(r"$\alpha^*_0$"); axs[0].set_ylabel(r"$g\ell^2/d$")
    fig.colorbar(cs, ax=axs[0], ticks=[0, 0.5, 1]).set_label(r"$\Sigma_i\bar\phi_i$ (GP)")
    for ax, y, p, lo, hi, lab in ((axs[1], L, pl, ll, hl, r"$\Sigma_i\bar\phi_i$ ($32^3$)"),
                                  (axs[2], V64, pv, lv, hv, r"$\sigma_\mathrm{vM}$ ($64^3$) [$10^{-4}$ Pa]")):
        a, b = y.min(), y.max(); pad = 0.06 * (b - a)
        ax.plot([a - pad, b + pad], [a - pad, b + pad], color="0.6", lw=0.8)
        ax.errorbar(y, p, yerr=[p - lo, hi - p], fmt="o", ms=3, lw=0.7, color="#3b5bdb")
        ax.set_xlabel("run: " + lab, fontsize=7.5); ax.set_ylabel("GP, left out", fontsize=7.5)
    if pair:
        a32, a64 = np.array(pair).T
        axs[3].scatter(a64, a32, s=10, color="#7048a8")
        m = [min(a32.min(), a64.min()), max(a32.max(), a64.max())]
        axs[3].plot(m, m, color="0.6", lw=0.8)
        axs[3].set_xlabel(r"$\sigma_\mathrm{vM}$ $64^3$ [$10^{-4}$ Pa]", fontsize=7.5); axs[3].set_ylabel(r"$\sigma_\mathrm{vM}$ $32^3$", fontsize=7.5)
    for ax, t in zip(axs, "abcd"):
        ax.text(0.04, 0.96, f"({t})", transform=ax.transAxes, va="top")
    out = ROOT / "assets" / "fig_surrogate_pilot.png"
    fig.savefig(out, dpi=400); fig.savefig(out.with_suffix(".pdf"))
    print("wrote", out)


if __name__ == "__main__":
    main()

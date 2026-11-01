#!/usr/bin/env python3
"""Fit a Gaussian process to the pilot runs (run_pilot.py) and check it by leave-one-out.

Inputs x = (log10 alpha*_0, g l^2/d), scaled to [0, 1]; outputs: seed living volume fraction and seed von Mises stress
at T* = 1. Kernel: constant x anisotropic RBF + white noise (scikit-learn). Writes summary.json and
assets/fig_surrogate_pilot.png (left: living fraction over the parameter plane from the GP with the runs; middle,
right: leave-one-out prediction against the run).

    python papers/surrogate_pilot/fit_pilot.py
"""
import json
import sys
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

LO, HI = np.array([1.0, 0.5]), np.array([np.log10(300.0), 2.0])


def load():
    X, Y = [], []
    for p in sorted((HERE / "runs").glob("p*.json")):
        r = json.loads(p.read_text())["runs"]["separate_ab"]
        a0 = float(p.stem.split("_a")[1].split("_gc")[0]); gc = float(p.stem.split("_gc")[1])
        X.append([np.log10(a0), gc]); Y.append([r["living_seed"][-1], r["seed_vm_Pa"] * 1e4])
    return (np.array(X) - LO) / (HI - LO), np.array(Y)


def gp():
    k = ConstantKernel(1.0, (1e-3, 1e3)) * RBF([0.3, 0.3], (0.02, 10.0)) + WhiteKernel(1e-6, (1e-10, 1e-2))
    return GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=5, random_state=0)


def main():
    figstyle.apply()
    X, Y = load()
    n = len(X)
    loo = np.zeros_like(Y); sd = np.zeros_like(Y)
    for j in range(2):
        for i in range(n):
            m = np.arange(n) != i
            g = gp().fit(X[m], Y[m, j])
            mu, s = g.predict(X[i:i + 1], return_std=True)
            loo[i, j], sd[i, j] = mu[0], s[0]
    err = np.abs(loo - Y)
    summary = {"n_runs": n, "loo_mean_abs_err": {"living": float(err[:, 0].mean()), "vm_1e-4Pa": float(err[:, 1].mean())},
               "loo_max_abs_err": {"living": float(err[:, 0].max()), "vm_1e-4Pa": float(err[:, 1].max())},
               "range": {"living": [float(Y[:, 0].min()), float(Y[:, 0].max())],
                         "vm_1e-4Pa": [float(Y[:, 1].min()), float(Y[:, 1].max())]}}
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))

    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.15), constrained_layout=True)
    g = gp().fit(X, Y[:, 0])
    u = np.linspace(0, 1, 80); U, V = np.meshgrid(u, u)
    Z = np.clip(g.predict(np.c_[U.ravel(), V.ravel()]).reshape(U.shape), 0, 1)
    a0 = 10 ** (LO[0] + U * (HI[0] - LO[0])); gc = LO[1] + V * (HI[1] - LO[1])
    cs = axs[0].contourf(a0, gc, Z, levels=np.linspace(0, 1, 11), cmap="Greens")
    axs[0].scatter(10 ** (LO[0] + X[:, 0] * (HI[0] - LO[0])), LO[1] + X[:, 1] * (HI[1] - LO[1]), s=10, c="k")
    axs[0].set_xscale("log"); axs[0].set_xlabel(r"$\alpha^*_0$"); axs[0].set_ylabel(r"$g\ell^2/d$")
    fig.colorbar(cs, ax=axs[0], label=r"$\Sigma_i\bar\phi_i$ (GP)", ticks=[0, 0.5, 1])
    for ax, j, lab in ((axs[1], 0, r"$\Sigma_i\bar\phi_i$"), (axs[2], 1, r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")):
        lo, hi = Y[:, j].min(), Y[:, j].max(); pad = 0.05 * (hi - lo + 1e-12)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="0.6", lw=0.8)
        ax.errorbar(Y[:, j], loo[:, j], yerr=2 * sd[:, j], fmt="o", ms=3, lw=0.8, color="#3b5bdb")
        ax.set_xlabel("run: " + lab); ax.set_ylabel("GP, left out: " + lab)
    for ax, t in zip(axs, "abc"):
        ax.text(0.04, 0.96, f"({t})", transform=ax.transAxes, va="top")
    out = ROOT / "assets" / "fig_surrogate_pilot.png"
    fig.savefig(out, dpi=400); fig.savefig(out.with_suffix(".pdf"))
    print("wrote", out)


if __name__ == "__main__":
    main()

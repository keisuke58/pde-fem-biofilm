#!/usr/bin/env python3
"""Uncertainty propagation with the surrogates of the pilot (Keio continuation, not thesis): assumed distributions of
alpha*_0 and g l^2/d are carried to the seed von Mises stress (Gaussian process on the 128^3 runs) and to the seed
living fraction (Gaussian process on logit(L) of the 32^3 runs). The distributions are illustrations, not data:
    alpha*_0 ~ log-normal, median 100, sigma of ln = ln 2 / 2 (95 % between 50 and 200), clipped to the design range
    g l^2/d  ~ uniform 0.7 .. 1.3
10 000 samples; the GP standard deviation is sampled as well (prediction plus GP noise), so the spread contains the
surrogate's own uncertainty.

    python papers/surrogate_pilot/propagate.py   -> assets/fig_surrogate_propagation.png/.pdf, propagation.json
"""
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
from fit_pilot import LO, HI, expit, gp, load, logit  # noqa: E402

warnings.filterwarnings("ignore")
N = 10_000


def main():
    figstyle.apply()
    rng = np.random.default_rng(0)
    a0 = np.clip(np.exp(np.log(100.0) + rng.normal(0.0, np.log(2.0) / 2, N)), 10.0, 300.0)
    gc = rng.uniform(0.7, 1.3, N)
    Xs = (np.c_[np.log10(a0), gc] - LO) / (HI - LO)

    Xl, L, _, _ = load("pq")                                  # living fraction, 32^3 (enough: 64^3 within 0.002)
    Xh, _, Vh, _ = load("gh")                                 # stress, 128^3 (16 points)
    gl = gp().fit(Xl, logit(L)); mu, sd = gl.predict(Xs, return_std=True)
    Ls = expit(mu + sd * rng.normal(size=N))
    gv = gp().fit(Xh, Vh); mu, sd = gv.predict(Xs, return_std=True)
    Vs = mu + sd * rng.normal(size=N)

    q = lambda v: [float(np.percentile(v, p)) for p in (5, 50, 95)]
    out = {"samples": N, "alpha0": "lognormal median 100, 95 % in 50..200", "g": "uniform 0.7..1.3",
           "living_5_50_95": q(Ls), "vm_1e-4Pa_5_50_95": q(Vs),
           "corr_living_vm": float(np.corrcoef(Ls, Vs)[0, 1])}
    (HERE / "propagation.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))

    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.0), constrained_layout=True)
    axs[0].hist(Ls, bins=50, color="#2b8a3e", alpha=0.8)
    axs[0].set_xlabel(r"$\Sigma_i\bar\phi_i$ in the seed"); axs[0].set_ylabel("samples")
    axs[1].hist(Vs, bins=50, color="#1f4e8c", alpha=0.8)
    axs[1].set_xlabel(r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")
    axs[2].plot(Ls[:2000], Vs[:2000], ".", ms=1.5, color="#7048a8", alpha=0.5)
    axs[2].set_xlabel(r"$\Sigma_i\bar\phi_i$"); axs[2].set_ylabel(r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")
    for ax, t in zip(axs, "abc"):
        ax.text(0.04, 0.96, f"({t})", transform=ax.transAxes, va="top")
    o = ROOT / "assets" / "fig_surrogate_propagation.png"
    fig.savefig(o, dpi=600); fig.savefig(o.with_suffix(".pdf"))
    print("wrote", o)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""fig_grid_examples.py -- grid convergence of all two-species examples (coupling paper, Section 3/4; thesis Ch. 4,
Table tab:ch5-gci): seed alpha - 1 (a) and seed von Mises stress (b) at T* = 1 against the grid (32^3, 64^3, 128^3),
outer step 0.0025, Python TVD, inherit rule. Dashed: Richardson extrapolation (Celik et al. 2008) where the three grids
are in the asymptotic range.

Data: the cloud runs results/2026-10-3d/*front1_dtc0.0025*inhliv*_{separate,mixed,mixed_ab}.json (no antibiotic, and the
Fritsch start fractions with antibiotic); the two antibiotic series with equal start amounts are from IKMHIWI03
(IKMHIWI03_TO_CLOUD.ja.md A-89; their JSONs stay on that machine), copied below.

    python papers/coupling_paper/fig_grid_examples.py -> papers/coupling_paper/fig_grid_examples.pdf
"""
import glob
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

RES = ROOT / "ansys_usermat/apdl/results/2026-10-3d"
N = (32, 64, 128)


def cloud(tag, key):
    out = []
    for n in N:
        pre = "" if n == 32 else f"nf{n}_"
        f = glob.glob(str(RES / f"coupled_living_consumption_{pre}front1_dtc0.0025_*inhliv*_tvd_wgmirror_fast_{tag}.json"))
        r = json.load(open(f[0]))["runs"][key]
        last = lambda q: r[q][-1] if isinstance(r[q], list) else r[q]  # noqa: E731
        out.append((last("seed_am1"), last("seed_vm_Pa")))
    return np.array(out)


SERIES = [  # label, (alpha-1, vM) per grid, colour, marker
    ("separate, antibiotic", np.array([[8.453e-4, 5.523e-4], [8.931e-4, 6.000e-4], [9.112e-4, 6.284e-4]]), "C0", "o"),
    ("separate", cloud("separate", "separate"), "C0", "s"),
    ("mixed, antibiotic", np.array([[7.42e-4, 2.00e-4], [8.04e-4, 2.82e-4], [8.32e-4, 3.30e-4]]), "C3", "o"),
    ("mixed", cloud("mixed", "mixed"), "C3", "s"),
    ("mixed (0.25, 0.30), antibiotic", cloud("mixed_ab", "mixed_ab"), "C2", "^"),
]


def gci3(f1, f2, f3, r=2.0):
    """Celik et al. 2008 with constant ratio r; finest grid first. None if not monotone in the asymptotic sense."""
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 / e21 <= 0:
        return None
    p = abs(np.log(abs(e32 / e21))) / np.log(r)
    if p < 0.5:
        return None
    return (r ** p * f1 - f2) / (r ** p - 1)


def main():
    figstyle.apply(size=8.5)
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.6), constrained_layout=True)
    h = 2.0 / np.array(N)
    for label, v, col, mk in SERIES:
        ls = "-" if "antibiotic" in label else "--"
        for k, ax in enumerate(axs):
            y = v[:, k] * 1e4
            ax.plot(h, y, color=col, marker=mk, ls=ls, lw=1, ms=4, label=label if k == 0 else None)
            ext = gci3(y[2], y[1], y[0])
            if ext is not None:
                ax.plot([h[2], 0], [y[2], ext], color=col, ls=":", lw=0.8)
                ax.plot([0], [ext], color=col, marker=mk, mfc="white", ms=4)
    for ax, yl, t in zip(axs, (r"seed $\alpha-1$ [$10^{-4}$]", r"seed von Mises [$10^{-4}$ Pa]"), "ab"):
        ax.set_xlabel(r"cell size $h$ [mm]")
        ax.set_ylabel(yl)
        ax.set_xlim(-0.003, 0.07)
        ax.text(-0.02, 1.02, f"({t})", transform=ax.transAxes, va="bottom", ha="left")
    axs[0].legend(fontsize=6.5, loc="lower left", frameon=False)
    out = HERE / "fig_grid_examples.pdf"
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()

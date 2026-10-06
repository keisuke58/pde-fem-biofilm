#!/usr/bin/env python3
"""Figure of keio_p1/results_2d.json.  python keio_p1/fig_nem_vs_galerkin.py -> assets/fig_nem_vs_galerkin.png"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ansys_usermat"))
import figstyle  # noqa: E402


def main():
    figstyle.apply(size=10)
    runs = [r for r in json.loads((HERE / "results_2d.json").read_text())["runs"] if not r["unstable"]]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    for beta, ls in ((0.02, "-"), (1e-4, "--")):
        for meth, col, lab in (("nem", "C0", "NEM, explicit"), ("galerkin", "C1", "Galerkin, implicit")):
            rr = sorted([r for r in runs if r["beta"] == beta and r["method"] == meth], key=lambda r: r["n"])
            h = np.array([2.0 / r["n"] for r in rr])
            b = r"$\beta=0.02$" if beta == 0.02 else r"$\beta=10^{-4}$"
            ax[0].loglog(h, [r["err_alpha"] for r in rr], ls, color=col, marker="o", ms=4, label=f"{lab}, {b}")
            ax[1].loglog(h, [r["step"] for r in rr], ls, color=col, marker="o", ms=4, label=f"{lab}, {b}")
            ax[2].loglog(h, [r["seconds"] for r in rr], ls, color=col, marker="o", ms=4, label=f"{lab}, {b}")
    ax[0].set(xlabel="element size $h$ [mm]", ylabel=r"relative error of seed $\alpha-1$")
    ax[1].set(xlabel="element size $h$ [mm]", ylabel=r"step of $\alpha-1$ at the seed surface / seed mean")
    ax[2].set(xlabel="element size $h$ [mm]", ylabel="run time [s] (Python, $T^*=1$)")
    for a in ax:
        a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    out = HERE.parent / "assets" / "fig_nem_vs_galerkin.png"
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    import matplotlib.ticker  # noqa: F401
    main()

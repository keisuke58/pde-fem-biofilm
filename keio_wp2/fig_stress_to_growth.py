#!/usr/bin/env python3
"""Figure of keio_wp2/results.json (stress_to_growth_2d.py).

    python keio_wp2/fig_stress_to_growth.py -> assets/fig_stress_to_growth.png
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ansys_usermat"))
import figstyle  # noqa: E402

K_ALPHA = 1e-3


def main():
    figstyle.apply(size=10)
    d = json.loads((HERE / "results.json").read_text())
    c = d["c_seed"]
    sweep = [r for r in d["runs"] if r["set"] == "sweep"]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    cmap = plt.get_cmap("viridis")
    for i, r in enumerate(sweep):
        h = np.array(r["hist"])
        pi = r["mu_star"] * c * K_ALPHA ** 2
        lab = r"$\mu^*=0$" if r["mu_star"] == 0 else (r"Table 2, $\Pi=%.0e$" % pi if r["mu_star"] > 1e9
                                                       else r"$\Pi=%.2g$" % pi)
        ax[0].plot(h[:, 0], h[:, 2], color=cmap(i / (len(sweep) - 1)), label=lab)
    ax[0].set(xlabel=r"time $T^*$", ylabel=r"mean $\phi$ in the seed")
    ax[0].legend(fontsize=7, frameon=False)
    pis = np.array([r["mu_star"] * c * K_ALPHA ** 2 for r in sweep[1:]])
    phis = np.array([r["hist"][-1][2] for r in sweep[1:]]) / sweep[0]["hist"][-1][2]
    ax[1].semilogx(pis, phis, "o-", color="C0")
    ax[1].axvline(1.0, color="0.5", ls=":", lw=1)
    ax[1].set(xlabel=r"$\Pi=\mu^*\,c\,(k_\alpha T^*)^2$", ylabel=r"seed $\phi$ at $T^*=1$ / without the term")
    for ms, mk in ((1e7, "-"), (d["mu_table2"], "--")):
        ref = [r for r in d["runs"] if r["set"] == "reference" and abs(r["mu_star"] - ms) < 1]
        if not ref:
            continue
        aref = ref[0]["hist"][-1][1]
        for j, sch in enumerate(["explicit", "semi", "iterated"]):
            rr = [r for r in d["runs"] if r["set"] == "scheme" and r["scheme"] == sch and abs(r["mu_star"] - ms) < 1]
            dts = np.array([r["dt"] for r in rr])
            err = np.array([abs(r["hist"][-1][1] - aref) / aref for r in rr])
            bad = np.array([r["unstable"] for r in rr])
            lab = f"{sch}, " + (r"$\mu^*=10^7$" if ms < 1e9 else "Table 2")
            ax[2].loglog(dts, err, mk, color=f"C{j}", marker="o", ms=4, label=lab)
            if bad.any():
                ax[2].loglog(dts[bad], err[bad], "x", color="k", ms=8)
    ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[2].set(xlabel=r"time step $\Delta t$", ylabel=r"relative error of seed $\alpha-1$")
    ax[2].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    out = HERE.parent / "assets" / "fig_stress_to_growth.png"
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__" and "--3d" not in sys.argv:
    main()


def fig3d():
    """seed phi at T* = 1 (relative to mu* = 0) against Pi, 2D and 3D (sphere, cube)."""
    figstyle.apply(size=10)
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    d2 = json.loads((HERE / "results.json").read_text())
    sw = [r for r in d2["runs"] if r["set"] == "sweep"]
    p0 = sw[0]["hist"][-1][2]
    ax.semilogx([r["mu_star"] * d2["c_seed"] * K_ALPHA ** 2 for r in sw[1:]],
                [r["hist"][-1][2] / p0 for r in sw[1:]], "o-", label="2D plane strain, c = 3.4 (seed)")
    d3 = json.loads((HERE / "results_3d_n16.json").read_text())
    for seed, mk in (("sphere", "s-"), ("cube", "^-")):
        rr = [r for r in d3["runs"] if r["seed"] == seed]
        crim = rr[0]["hist"][-1][6]
        q0 = rr[0]["hist"][-1][2]
        ax.semilogx([r["mu_star"] * crim * K_ALPHA ** 2 for r in rr[1:]], [r["hist"][-1][2] / q0 for r in rr[1:]],
                    mk, label=f"3D {seed}, c = {crim:.1f} (seed surface)")
    ax.axvline(1, color="0.5", ls=":", lw=1)
    ax.set(xlabel=r"$\Pi=\mu^*\,c\,(k_\alpha T^*)^2$", ylabel=r"seed $\phi$ at $T^*=1$ / without the term")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    out = HERE.parent / "assets" / "fig_stress_to_growth_3d.png"
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__" and "--3d" in sys.argv:
    fig3d()

#!/usr/bin/env python3
"""Figure of keio_wp2/results_tooth_bc.json (homeostatic_tooth.py --bc): the
bottom boundary condition, implant collar against tooth crown, beta = 0.02,
no feedback.

    python keio_wp2/fig_tooth_bc.py -> assets/fig_tooth_bc.png
"""
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

BOTTOMS = (("free", "bottom free", "C2"), ("uz", r"bottom $u_z=0$", "C0"), ("clamped", "bottom clamped", "C3"))


def main():
    figstyle.apply(size=10)
    d = json.loads((HERE / "results_tooth_bc.json").read_text())

    def pick(bulge, bottom):
        return next(r for r in d["runs"] if r["bulge"] == bulge and r["bottom"] == bottom
                    and r["P_h"] is None and r["nz"] == 64)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    for bottom, lab, c in BOTTOMS:
        for bulge, ls in ((1.0, "-"), (0.0, ":")):
            w = pick(bulge, bottom)["wall"]
            name = ("tooth, " if bulge else "implant, ") + lab
            ax[0].plot(-np.array(w["stt"]), w["z"], ls, color=c, label=name)
            ax[1].plot(np.array(w["tn"]), w["z"], ls, color=c, label=name)
    ax[0].set(xscale="log", xlabel=r"hoop compression on the surface $-\sigma_{\theta\theta}$ [Pa]",
              ylabel=r"$z$ [mm]", title="(a) hoop stress")
    ax[1].axvline(0, color="0.6", lw=0.8)
    ax[1].set(xlabel=r"normal stress on the surface $\sigma_{nn}$ [Pa] (tension $>0$)", ylabel=r"$z$ [mm]",
              title="(b) normal stress (pull-off)", xlim=(-6e-4, 3e-4))
    ax[1].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    ax[0].legend(fontsize=7, frameon=False, loc="upper right")
    x = np.arange(len(BOTTOMS))
    for k, (bulge, c, name) in enumerate(((0.0, "0.6", "implant collar"), (1.0, "C3", "tooth crown"))):
        v = [-min(pick(bulge, b)["wall"]["stt"]) for b, _, _ in BOTTOMS]
        ax[2].bar(x + (k - 0.5) * 0.36, v, 0.34, color=c, label=name)
        for xi, vi in zip(x, v):
            ax[2].text(xi + (k - 0.5) * 0.36, vi * 1.08, "%.1e" % vi, ha="center", fontsize=7)
    ax[2].set(yscale="log", xticks=x, xticklabels=["free", r"$u_z=0$", "clamped"], ylim=(5e-4, 3e-2),
              xlabel="bottom face", ylabel=r"largest $-\sigma_{\theta\theta}$ on the surface [Pa]",
              title="(c) peak hoop stress")
    ax[2].legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    out = HERE.parent / "assets" / "fig_tooth_bc.png"
    fig.savefig(out, dpi=200)
    print(out)


if __name__ == "__main__":
    main()

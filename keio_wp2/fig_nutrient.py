#!/usr/bin/env python3
"""Figure of keio_wp2/results_nutrient.json (homeostatic_nutrient.py), with
the runs without nutrient from results_tooth.json for reference.

    python keio_wp2/fig_nutrient.py -> assets/fig_nutrient_homeostatic.png
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

GEOM = {0.0: ("implant (nutrient from the top)", "C0"), 1.0: ("tooth (nutrient from the outer face)", "C3")}


def main():
    figstyle.apply(size=10)
    d = json.loads((HERE / "results_nutrient.json").read_text())
    d0 = json.loads((HERE / "results_tooth.json").read_text())

    def pick(bulge, P_h=None, local=True, w=0.5, g=1.0):
        return next(r for r in d["runs"] if r["bulge"] == bulge and r["P_h"] == P_h and r["local"] == local
                    and r["w"] == w and r["g"] == g)
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.9))

    for b, (name, c) in GEOM.items():
        h = np.array(pick(b)["hist"])
        ax[0].plot(h[:, 0], h[:, 2], "-", color=c, label=name.split(" (")[0] + r", $w=0.5$, $g=1$")
    for (b, w, g, ls, lab) in ((0.0, 0.0, 1.0, "--", r"implant, $w=0$"), (0.0, 0.5, 6.0, ":", r"implant, $g=6$")):
        h = np.array(pick(b, w=w, g=g)["hist"])
        ax[0].plot(h[:, 0], h[:, 2], ls, color="C0", label=lab)
    ax[0].set(xlabel=r"time $T^*$", ylabel=r"mean $\phi$ in the layer", ylim=(0, 1.05), title="(a) filling of the layer")
    ax[0].legend(fontsize=7, frameon=False, loc="lower right")

    for b, (name, c) in GEOM.items():
        for g, ls in ((1.0, "-"), (6.0, ":")):
            f = pick(b, g=g)["field"]
            ax[1].plot(np.array(f["c"])[0], f["z"], ls, color=c, label=name.split(" (")[0] + r", $g=%g$" % g)
    ax[1].set(xlabel=r"nutrient $c$ on the bonded surface, $T^*=1$", ylabel=r"$z$ [mm]", xlim=(0, 1.05),
              title="(b) nutrient")
    ax[1].legend(fontsize=7, frameon=False, loc="lower right")

    for b, (name, c) in GEOM.items():
        for P_h, ls, lab in ((None, "-", "no feedback"), (0.1, "--", r"$P_h=0.1$")):
            f = pick(b, P_h=P_h)["field"]
            ax[2].plot(-np.array(f["stt"]), f["z"], ls, color=c, label=name.split(" (")[0] + ", " + lab)
        r0 = next(r for r in d0["runs"] if r["beta"] == 0.02 and r["bulge"] == b and r["P_h"] is None)
        ax[2].plot(-np.array(r0["wall"]["stt"]), r0["wall"]["z"], ":", color=c,
                   label=name.split(" (")[0] + ", no nutrient (step 6)")
    ax[2].set(xscale="log", xlabel=r"hoop compression on the surface $-\sigma_{\theta\theta}$ [Pa]",
              ylabel=r"$z$ [mm]", title="(c) stress on the surface")
    ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[2].legend(fontsize=6.5, frameon=False, loc="center right")

    for b, (name, c) in GEOM.items():
        for P_h, ls, lab in ((None, "-", "no feedback"), (0.1, "--", r"$P_h=0.1$")):
            h = np.array(pick(b, P_h=P_h)["hist"])
            ax[3].plot(h[:, 0], h[:, 1], ls, color=c, label=name.split(" (")[0] + ", " + lab)
        r0 = next(r for r in d0["runs"] if r["beta"] == 0.02 and r["bulge"] == b and r["P_h"] is None)
        h0 = np.array(r0["hist"])
        ax[3].plot(h0[:, 0], h0[:, 1], ":", color=c, label=name.split(" (")[0] + ", no nutrient")
    ax[3].set(xlabel=r"time $T^*$", ylabel=r"mean $\alpha-1$ in the layer", title="(d) growth")
    ax[3].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax[3].legend(fontsize=6.5, frameon=False, loc="upper left")
    fig.tight_layout()
    out = HERE.parent / "assets" / "fig_nutrient_homeostatic.png"
    fig.savefig(out, dpi=200)
    print(out)


if __name__ == "__main__":
    main()

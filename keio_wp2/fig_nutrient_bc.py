#!/usr/bin/env python3
"""Figure of keio_wp2/results_nutrient_bc.json (homeostatic_nutrient.py --bc):
the bottom boundary condition with the nutrient and the front term, without
feedback and with P_h = 0.1.

    python keio_wp2/fig_nutrient_bc.py -> assets/fig_nutrient_bc.png
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

BOTTOMS = (("free", "bottom free", "C2"), ("uz", r"bottom $u_z=0$", "C0"), ("clamped", "bottom clamped", "C3"))


def main():
    figstyle.apply(size=10)
    d = json.loads((HERE / "results_nutrient_bc.json").read_text())

    def pick(bulge, bottom, P_h):
        return next(r for r in d["runs"] if r["bulge"] == bulge and r["bottom"] == bottom and r["P_h"] == P_h)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.0))
    x = np.arange(len(BOTTOMS))
    bars = ((0.0, None, "0.75", "implant, no feedback"), (1.0, None, "C3", "tooth, no feedback"),
            (0.0, 0.1, "0.45", r"implant, $P_h=0.1$"), (1.0, 0.1, "#8b1a1a", r"tooth, $P_h=0.1$"))
    for k, (b, P_h, c, lab) in enumerate(bars):
        v = [-min(pick(b, bot, P_h)["field"]["stt"]) for bot, _, _ in BOTTOMS]
        xs = x + (k - 1.5) * 0.2
        ax[0].bar(xs, v, 0.18, color=c, label=lab)
        for xi, vi in zip(xs, v):
            ax[0].text(xi, vi * 1.1, "%.1e" % vi, ha="center", fontsize=6, rotation=90)
    ax[0].set(yscale="log", xticks=x, xticklabels=["free", r"$u_z=0$", "clamped"], ylim=(1e-3, 1.5),
              xlabel="bottom face", ylabel=r"largest $-\sigma_{\theta\theta}$ on the surface [Pa]",
              title=r"(a) peak hoop stress, $T^*=1$")
    ax[0].legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
    for j, b in enumerate((1.0, 0.0)):
        for bot, lab, c in BOTTOMS:
            for P_h, ls in ((None, "-"), (0.1, "--")):
                f = pick(b, bot, P_h)["field"]
                ax[1 + j].plot(-np.array(f["stt"]), f["z"], ls, color=c,
                               label=lab + (", no feedback" if P_h is None else r", $P_h=0.1$"))
        ax[1 + j].set(xscale="log", xlabel=r"hoop compression on the surface $-\sigma_{\theta\theta}$ [Pa]",
                      ylabel=r"$z$ [mm]", title="(%s) %s" % ("bc"[j], "tooth crown" if b else "implant collar"))
        ax[1 + j].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax[1 + j].legend(fontsize=6.5, frameon=False, loc="center left", bbox_to_anchor=(0.3, 0.55))
    fig.tight_layout()
    out = HERE.parent / "assets" / "fig_nutrient_bc.png"
    fig.savefig(out, dpi=200)
    print(out)


if __name__ == "__main__":
    main()

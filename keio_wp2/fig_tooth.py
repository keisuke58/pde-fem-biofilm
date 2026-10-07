#!/usr/bin/env python3
"""Figure of keio_wp2/results_tooth.json (homeostatic_tooth.py): implant collar
against tooth crown, along the bonded surface, without feedback and with
P_h = 0.1.

    python keio_wp2/fig_tooth.py -> assets/fig_tooth_homeostatic.png
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

RI, HZ = 2.05, 2.0


def main():
    figstyle.apply(size=10)
    d = json.loads((HERE / "results_tooth.json").read_text())
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.9))
    z = np.linspace(0, HZ, 50)
    for b, c in ((0.0, "C0"), (1.0, "C3")):
        r = RI + b * np.sin(0.5 * np.pi * z / HZ)
        lab = "implant collar" if b == 0 else "tooth crown"
        ax[0].fill_betweenx(z, r, r + 0.25, color=c, alpha=0.35, lw=0, label=lab + " (biofilm)")
        ax[0].plot(r, z, color=c, lw=2)
    ax[0].set(xlabel=r"$r$ [mm]", ylabel=r"$z$ [mm]", title="(a) geometry (axisymmetric)", xlim=(1.9, 3.5))
    ax[0].set_aspect("equal")
    ax[0].legend(fontsize=7, frameon=False, loc="lower right")
    for j, beta in enumerate((0.02, 1e-3)):
        for b, c in ((0.0, "C0"), (1.0, "C3")):
            for P_h, local, ls in ((None, False, "-"), (0.1, True, "--")):
                r = next(x for x in d["runs"] if x["beta"] == beta and x["bulge"] == b
                         and x["P_h"] == P_h and x["local"] == local)
                w = r["wall"]
                lab = ("implant" if b == 0 else "tooth") + (", no feedback" if P_h is None else r", $P_h=0.1$ local")
                ax[1 + j].plot(-np.array(w["stt"]), w["z"], ls, color=c, label=lab)
                if j == 0:
                    ax[3].plot(w["alpha"], w["z"], ls, color=c, label=lab)
        ax[1 + j].set(xscale="log", xlabel=r"hoop compression on the surface $-\sigma_{\theta\theta}$ [Pa]",
                      ylabel=r"$z$ [mm]", title=r"(%s) $\beta = %g$ mm$^2/T^*$" % ("bc"[j], beta))
        ax[1 + j].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax[1 + j].legend(fontsize=7, frameon=False)
    ax[3].set(xlabel=r"$\alpha-1$ on the surface, $T^*=1$", ylabel=r"$z$ [mm]",
              title=r"(d) growth, $\beta = 0.02$ mm$^2/T^*$")
    ax[3].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    out = HERE.parent / "assets" / "fig_tooth_homeostatic.png"
    fig.savefig(out, dpi=200)
    print(out)


if __name__ == "__main__":
    main()

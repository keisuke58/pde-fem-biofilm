#!/usr/bin/env python3
"""Figure of keio_wp2/results_implant.json and results_implant_beta.json
(homeostatic_implant.py): the homeostatic growth law on the implant collar.

    python keio_wp2/fig_implant.py -> assets/fig_implant_homeostatic.png
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
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
import homeostatic_implant as H  # noqa: E402

LAB = {False: r"$p_h = P_h\,E\,k_\alpha T^*$", True: r"$p_h = P_h\,E(\phi^2+f)\,k_\alpha T^*$"}


def phi_profile(beta, dt=0.02, t_end=1.0):
    """phi across the layer at mid-height at T* (phi does not depend on p_h)."""
    M = H.Mesh()
    phi = M.seed.astype(float)
    nsub = int(np.ceil(dt / (0.8 * min(M.hr, M.hz) ** 2 / (4 * beta))))
    for _ in range(int(round(t_end / dt)) * nsub):
        phi = np.clip(phi + dt / nsub * (beta * M.lap(phi) + H.K_ALPHA), 0, 1)   # alpha ~ 1
    r = M.rc.reshape(M.nr, M.nz)[:, M.nz // 2]
    return r - H.RI, phi.reshape(M.nr, M.nz)[:, M.nz // 2]


def main():
    figstyle.apply(size=10)
    d = json.loads((HERE / "results_implant.json").read_text())
    db = json.loads((HERE / "results_implant_beta.json").read_text())
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.9))

    base = d["runs"][0]["hist"][-1][1]
    for c, local in (("C0", False), ("C3", True)):
        rs = [r for r in d["runs"] if r["P_h"] is not None and r["local"] == local]
        ax[0].plot([r["P_h"] for r in rs], [r["hist"][-1][1] / base for r in rs], "o-", color=c, label=LAB[local])
    ax[0].set(xscale="log", xlabel=r"$P_h$", ylabel=r"growth $\alpha-1$ / no feedback", ylim=(0, 1.05),
              title=r"(a) $\beta = 0.02$ mm$^2/T^*$")
    ax[0].legend(fontsize=7, frameon=False, loc="lower right")

    betas = sorted({r["beta"] for r in db["runs"]}, reverse=True)

    def pick(beta, P_h, local):
        return next(r["hist"][-1] for r in db["runs"]
                    if r["beta"] == beta and r["P_h"] == P_h and r["local"] == local)
    for c, local in (("C0", False), ("C3", True)):
        ax[1].plot(betas, [pick(b, 0.1, local)[1] / pick(b, None, False)[1] for b in betas], "o-",
                   color=c, label=LAB[local])
    ax[1].set(xscale="log", xlabel=r"$\beta$ [mm$^2/T^*$]", ylabel=r"growth $\alpha-1$ / no feedback",
              ylim=(0, 1.05), title=r"(b) $P_h = 0.1$")
    ax[1].legend(fontsize=7, frameon=False, loc="upper left")

    ax[2].plot(betas, [-pick(b, None, False)[5] for b in betas], "o-", color="k", label="no feedback")
    for c, local in (("C0", False), ("C3", True)):
        ax[2].plot(betas, [-pick(b, 0.1, local)[5] for b in betas], "o-", color=c, label=LAB[local])
    ax[2].set(xscale="log", yscale="log", xlabel=r"$\beta$ [mm$^2/T^*$]",
              ylabel=r"hoop compression on Ti $-\sigma_{\theta\theta}$ [Pa]", title=r"(c) stress on the titanium")
    ax[2].legend(fontsize=7, frameon=False)

    cmap = plt.get_cmap("viridis")
    for i, b in enumerate(betas):
        x, ph = phi_profile(b)
        ax[3].plot(x, ph, "o-", ms=3, color=cmap(i / (len(betas) - 1)), label=r"$\beta=%g$" % b)
    ax[3].set(xlabel=r"distance from the titanium [mm]", ylabel=r"$\phi$ at mid-height, $T^*=1$",
              title=r"(d) $\phi$ across the layer", ylim=(0, 1.05))
    ax[3].legend(fontsize=7, frameon=False)

    fig.tight_layout()
    out = HERE.parent / "assets" / "fig_implant_homeostatic.png"
    fig.savefig(out, dpi=200)
    print(out)


if __name__ == "__main__":
    main()

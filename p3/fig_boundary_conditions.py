"""fig_boundary_conditions.py -- boundary and initial conditions of the 3D seed problem of the unified model
(p3/uel/make_seed_inp.py, JAXFEM/fuel_front_3d_coupled.py): (a) the unit cube with the hemispherical seed,
(b) the mid-section y = 0.5 with the conditions on each face.

    python p3/fig_boundary_conditions.py      -> assets/fig_boundary_conditions.pdf, .png
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

SEED = "#1b6e5a"
NUT = "#3b6fb6"
SUB = "#8a8f8c"
INK = "#26302c"
R0, X0, Y0 = 0.25, 0.5, 0.5


def panel_a(ax):
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-58)
    for s in (ax.xaxis, ax.yaxis, ax.zaxis):
        s.pane.fill = False; s.pane.set_edgecolor("white")
    ax.grid(False)
    # substrate z = 0 and supply face x = 1
    ax.add_collection3d(Poly3DCollection([[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]], facecolor=SUB, alpha=0.35,
                                         edgecolor=INK, lw=0.6))
    ax.add_collection3d(Poly3DCollection([[(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)]], facecolor=NUT, alpha=0.30,
                                         edgecolor=INK, lw=0.6))
    # cube edges
    for a in (0, 1):
        for b in (0, 1):
            ax.plot([0, 1], [a, a], [b, b], color=INK, lw=0.6)
            ax.plot([a, a], [0, 1], [b, b], color=INK, lw=0.6)
            ax.plot([a, a], [b, b], [0, 1], color=INK, lw=0.6)
    # hemispherical seed
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 48), np.linspace(0, np.pi / 2, 16))
    ax.plot_surface(X0 + R0 * np.cos(u) * np.sin(v), Y0 + R0 * np.sin(u) * np.sin(v), R0 * np.cos(v),
                    color=SEED, alpha=0.75, lw=0, shade=True)
    ax.text(0.35, -0.05, -0.32, "substrate $z=0$: $\\mathbf{u}=\\mathbf{0}$", color=INK, fontsize=9, ha="center")
    ax.text(1.0, 0.3, 1.32, "supply $x=1$: $c=1$", color=NUT, fontsize=9, ha="center")
    ax.text(0.5, 0.5, 0.36, "seed", color=SEED, fontsize=9, ha="center")
    ax.set_xlabel("$x$", labelpad=-8); ax.set_ylabel("$y$", labelpad=-8); ax.set_zlabel("$z$", labelpad=-8)
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1]); ax.set_zticks([0, 1])
    ax.tick_params(pad=-3, labelsize=8)
    ax.set_title("(a) unit cube, hemispherical seed", loc="left", y=0.98)


def panel_b(ax):
    ax.set_aspect("equal"); ax.set_xlim(-0.33, 1.42); ax.set_ylim(-0.30, 1.22); ax.set_axis_off()
    # initial phase: smooth hemisphere phi = 0.9 (1 - tanh((r - R0)/w)) / 2
    x = np.linspace(0, 1, 300); z = np.linspace(0, 1, 300); Xg, Zg = np.meshgrid(x, z)
    phi = 0.45 * (1 - np.tanh((np.hypot(Xg - X0, Zg) - R0) / 0.03)) + 1e-4
    ax.contourf(Xg, Zg, phi, levels=np.linspace(0, 0.95, 11), cmap="Greens", alpha=0.9)
    ax.contour(Xg, Zg, phi, levels=[0.45], colors=SEED, linewidths=1.0)
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=INK, lw=0.9)
    # substrate hatching
    for xx in np.linspace(0, 1, 21):
        ax.plot([xx, xx - 0.04], [0, -0.05], color=SUB, lw=0.8)
    ax.plot([0, 1], [0, 0], color=INK, lw=2.0)
    # supply face
    ax.plot([1, 1], [0, 1], color=NUT, lw=3.0)
    for zz in np.linspace(0.12, 0.88, 5):
        ax.annotate("", xy=(1.0, zz), xytext=(1.13, zz), arrowprops=dict(arrowstyle="-|>", color=NUT, lw=0.9))
    T = dict(fontsize=8.6, color=INK, va="center")
    ax.text(0.5, -0.13, "$z=0$ (substrate): $\\mathbf{u}=\\mathbf{0}$;  no flux of $\\phi$, $c$", ha="center", **T)
    ax.text(1.16, 0.62, "$x=1$:\n$c=1$\n(supply)", ha="left", **dict(T, color=NUT))
    ax.text(1.16, 0.30, "traction-free,\nno flux of $\\phi$", ha="left", fontsize=8.0, color=INK, va="center")
    ax.text(-0.04, 0.5, "$x=0$, $y=0,1$, $z=1$:\ntraction-free;\nno flux of\n$\\phi$ and $c$", ha="right", **T)
    ax.text(0.5, 1.16, "initial: $\\phi=0.9$ in the seed ($r<0.25$, smooth edge 0.03), $10^{-4}$ outside;", ha="center",
            fontsize=8.2, color=INK)
    ax.text(0.5, 1.08, "$c=0.1$ (= $c_{\\rm stop}$) inside, composition from the calibrated Day-1 state", ha="center",
            fontsize=8.2, color=INK)
    ax.text(0.5, 0.18, "seed\n$\\phi\\approx0.9$", ha="center", fontsize=8.6, color="white")
    ax.annotate("", xy=(0.95, 0.40), xytext=(0.62, 0.40), arrowprops=dict(arrowstyle="<|-", color=NUT, lw=0.8, ls="--"))
    ax.text(0.79, 0.44, "nutrient diffuses in", ha="center", fontsize=7.8, color=NUT)
    ax.set_title("(b) mid-section $y=0.5$: boundary and initial conditions", loc="left")


if __name__ == "__main__":
    fig = plt.figure(figsize=(10.5, 4.6))
    a = fig.add_axes([0.0, 0.02, 0.40, 0.94], projection="3d")
    b = fig.add_axes([0.42, 0.02, 0.58, 0.94])
    panel_a(a); panel_b(b)
    out = ROOT / "assets" / "fig_boundary_conditions"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=170)
    print("wrote", out.with_suffix(".pdf"))

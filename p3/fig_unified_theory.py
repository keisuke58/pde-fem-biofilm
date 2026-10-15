"""fig_unified_theory.py -- four-panel figure of the unified two-level model for a paper:
(a) structure of the model, (b) the nutrient ceiling of the growth, (c) stress proportional to the
growth, (d) predicted seed stress against the stiffness.

    python p3/fig_unified_theory.py          -> assets/fig_unified_theory.pdf, .png

Data of (b)-(d): the user element (p3/uel/unified_uel.f90) on the 8^3 hemispherical seed, finite strain,
nu = 0.49 with the mean dilatation, Day 21, run through p3/uel/emulate_abaqus.py, which matches Abaqus 2024
to 4-5 digits (p3/uel/README.md). Values copied from the runs of 15 Oct 2026 (calibrate_growth.py).
Not from a paper: L0, Y, c_stop, D (the measured CH volume fixes L0 for a given Y).
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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

# (L0, seed alpha - 1, seed von Mises / E) per yield Y; 8^3, nu = 0.49, mean dilatation
RUNS = {
    1.0: [(1.926, 0.0346, 0.012986), (2.0, 0.0357, 0.013419), (3.597, 0.0572, 0.021770), (12.12, 0.1053, 0.045629),
          (13.76, 0.1067, 0.047813), (28.53, 0.0821, 0.055670), (40.31, 0.0516, 0.060875)],
    0.25: [(20.0, 0.0943, 0.034384), (50.0, 0.2057, 0.075620), (53.1, 0.2154, 0.079371)],
    0.1: [(50.0, 0.0995, 0.035916), (113.9, 0.2165, 0.077321), (120.0, 0.2271, 0.081122)],
}
TARGET = 1.80 ** (1 / 3) - 1           # measured CH volume, Day 3 -> 21: 1.80
NU = 0.49
COL = {1.0: "#1b6e5a", 0.25: "#c0662b", 0.1: "#5b5ea6"}


def box(ax, x, y, w, h, text, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02", fc=fc, ec="#33413b", lw=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.0, linespacing=1.35)


def arrow(ax, p, q, text="", off=(0, 0)):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=10, lw=0.9, color="#33413b",
                                 connectionstyle="arc3,rad=0.0"))
    if text:
        ax.text((p[0] + q[0]) / 2 + off[0], (p[1] + q[1]) / 2 + off[1], text, ha="center", va="center", fontsize=7.8,
                style="italic", color="#33413b")


def panel_a(ax):
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    box(ax, 0.02, 0.70, 0.40, 0.25, "composition (Klempt et al. 2026)\n$\\phi_i,\\ \\psi_i,\\ \\phi_0$ at each point\n"
        "living biomass $\\sum_i\\phi_i\\psi_i$", "#e3f0ea")
    box(ax, 0.58, 0.70, 0.40, 0.25, "phase (Klempt et al. 2024)\n$\\phi(\\mathbf{X},t)$, weight $\\phi^2$\n"
        "energy $\\beta|\\nabla\\phi|^2/2$", "#e8edf6")
    box(ax, 0.02, 0.36, 0.40, 0.22, "nutrient $c$\nenergy $c(\\ln c-1)$\n$\\dot c=D\\Delta c-Y\\dot\\xi$", "#eef3fa")
    box(ax, 0.58, 0.36, 0.40, 0.22, "growth, a reaction\n$\\dot\\xi=L_0\\,\\phi\\sum_i\\phi_i\\psi_i\\,(Y\\ln c-\\mu_b)$\n"
        "$\\alpha=1+\\xi$", "#f6efe4")
    box(ax, 0.20, 0.04, 0.60, 0.20, "mechanics: $\\mathbf{F}=\\mathbf{F}_e\\mathbf{F}_g$, $\\mathbf{F}_g=\\alpha\\mathbf{I}$\n"
        "$W=\\alpha^3(\\phi^2+f)\\,[\\mu(\\bar I_1-3)/2+K(J_e-1)^2/2]$", "#f3eef0")
    arrow(ax, (0.42, 0.82), (0.58, 0.82))
    ax.text(0.50, 0.99, "force on the phase $\\lambda h'(\\phi)(\\Psi-\\Psi_{\\rm ref})$", ha="center", va="center",
            fontsize=7.8, style="italic", color="#33413b")
    arrow(ax, (0.50, 0.965), (0.50, 0.83))
    arrow(ax, (0.30, 0.70), (0.66, 0.58), "catalyst", (-0.04, 0.0))
    arrow(ax, (0.42, 0.47), (0.58, 0.47), "affinity", (0, 0.035))
    arrow(ax, (0.58, 0.42), (0.42, 0.42), "consumption", (0, -0.035))
    arrow(ax, (0.78, 0.36), (0.62, 0.24), "$\\mathbf{F}_g$", (0.04, 0.0))
    arrow(ax, (0.80, 0.70), (0.80, 0.58), "$\\phi$", (0.03, 0.0))
    ax.text(0.5, -0.04, "one energy and one dissipation; residual = gradient of the discrete functional",
            ha="center", va="top", fontsize=7.8, color="#33413b")
    ax.set_title("(a) model", loc="left")


def panel_b(ax):
    for Y, rows in RUNS.items():
        r = sorted(rows)
        ax.plot([x[0] for x in r], [x[1] for x in r], "o-", ms=4, lw=1.2, color=COL[Y], label=f"$Y={Y:g}$")
    ax.axhline(TARGET, color="#555", lw=0.9, ls="--")
    ax.text(1.75, TARGET + 0.006, "measured CH volume, Day 3 to 21: $\\times1.80$", fontsize=8, color="#555")
    ax.set_xscale("log"); ax.set_xlabel("growth rate $L_0$"); ax.set_ylabel("seed $\\alpha-1$ at Day 21")
    ax.set_ylim(0, 0.25)
    ax.annotate("ceiling: fast growth\nstarves the seed", xy=(13.8, 0.105), xytext=(17, 0.012), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.8, color=COL[1.0]), color=COL[1.0])
    ax.legend(title="nutrient per growth", loc="upper left", bbox_to_anchor=(0.0, 0.84), frameon=False)
    ax.set_title("(b) nutrient limits the growth", loc="left")


def panel_c(ax):
    for Y, rows in RUNS.items():
        a = np.array([x[1] for x in rows]); s = np.array([x[2] for x in rows]); L = np.array([x[0] for x in rows])
        below = ~((Y == 1.0) & (L > 14))
        ax.plot(a[below], s[below], "o", ms=4.5, color=COL[Y], label=f"$Y={Y:g}$")
        if (~below).any():
            ax.plot(a[~below], s[~below], "o", ms=4.5, mfc="none", color=COL[Y], label="past the ceiling")
    pts = np.array([(x[1], x[2]) for Y, rows in RUNS.items() for x in rows if not (Y == 1.0 and x[0] > 14)])
    k = (pts[:, 0] @ pts[:, 1]) / (pts[:, 0] @ pts[:, 0])
    g = np.linspace(0, 0.24, 50)
    ax.plot(g, k * g, color="#222", lw=1.0)
    ax.text(0.13, k * 0.13 - 0.017, f"$\\sigma_{{vM}}\\approx{k:.2f}\\,E\\,(\\alpha-1)$", fontsize=9)
    ax.plot(g, g / (1 - NU) * 0.25, color="#999", lw=0.8, ls=":", label="fixed film $E(\\alpha-1)/(1-\\nu)$, $\\times0.25$")
    ax.set_xlabel("seed $\\alpha-1$"); ax.set_ylabel("seed von Mises stress $/E$")
    ax.set_xlim(0, 0.24); ax.set_ylim(0, 0.12)
    ax.legend(loc="upper left", frameon=False, fontsize=7.6)
    ax.set_title("(c) stress follows the growth", loc="left")
    return k


def panel_d(ax, k):
    E = np.logspace(0.5, 5, 100)
    s = k * TARGET * E
    ax.loglog(E, s, color="#1b6e5a", lw=1.6, label="seed von Mises, growth fitted to the volume")
    ax.loglog(E, TARGET / (1 - NU) * E, color="#999", lw=0.9, ls=":", label="laterally fixed film (upper bound)")
    ax.axhspan(1e5, 1e6, color="#c0662b", alpha=0.15, lw=0)
    ax.text(5, 2.2e5, "pressures reported to stop growth", fontsize=8, color="#8a4a1f")
    for Ev, lab, xt in ((10, "Klempt et al. 2024", (25, 0.45)), (600, "0.6 kPa", (1100, 3.5)),
                        (41e3, "Pattem et al. 2018", (4.5e3, 60))):
        ax.plot(Ev, k * TARGET * Ev, "o", color="#1b6e5a", ms=5)
        v = k * TARGET * Ev
        ax.annotate(f"{lab}\n{v:.2g} Pa" if v < 1e3 else f"{lab}\n{v / 1e3:.1f} kPa", xy=(Ev, v), xytext=xt,
                    fontsize=7.6, arrowprops=dict(arrowstyle="-", lw=0.6, color="#555"))
    ax.set_xlabel("biofilm stiffness $E$ [Pa]"); ax.set_ylabel("stress [Pa]")
    ax.set_xlim(3, 1e5); ax.set_ylim(0.3, 2e6)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.66), frameon=False, fontsize=7.8)
    ax.set_title("(d) predicted stress", loc="left")


if __name__ == "__main__":
    fig = plt.figure(figsize=(10.5, 8.2))
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    k = panel_c(fig.add_subplot(gs[1, 0]))
    panel_d(fig.add_subplot(gs[1, 1]), k)
    out = ROOT / "assets" / "fig_unified_theory"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=170)
    print("wrote", out.with_suffix(".pdf"), "k =", round(k, 3))

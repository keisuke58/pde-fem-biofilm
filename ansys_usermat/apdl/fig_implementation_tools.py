"""fig_implementation_tools.py -- how the coupled run with the calibrated point model is implemented and
which tools it uses: preparation (calibration, case file), the custom ANSYS executable (partner element,
material routine, point model in Fortran), and the post-processing.

    python ansys_usermat/apdl/fig_implementation_tools.py   -> assets/fig_implementation_tools.pdf, .png
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

INK = "#26302c"
USED = "#e6e8e7"      # used as it is
MINE = "#dce8f6"      # added in this work
EDGE_MINE = "#1e50a0"
TOOL = "#5d6b64"


def box(ax, x, y, w, h, title, body, tool, mine=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.015",
                                fc=MINE if mine else USED, ec=EDGE_MINE if mine else "#8a8f8c", lw=1.0))
    ax.text(x + w / 2, y + h - 0.035, title, ha="center", va="top", fontsize=9.2, fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h / 2 - 0.01, body, ha="center", va="center", fontsize=8.0, color=INK, linespacing=1.3)
    ax.text(x + w / 2, y + 0.022, tool, ha="center", va="bottom", fontsize=7.4, style="italic", color=TOOL)


def arrow(ax, p, q, text="", off=(0, 0.025), both=False):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>" if both else "-|>", mutation_scale=10, lw=1.0, color=INK))
    if text:
        ax.text((p[0] + q[0]) / 2 + off[0], (p[1] + q[1]) / 2 + off[1], text, ha="center", va="center", fontsize=7.8,
                color=INK)


def band(ax, y, h, label):
    ax.add_patch(FancyBboxPatch((0.005, y), 0.99, h, boxstyle="round,pad=0.004,rounding_size=0.012",
                                fc="none", ec="#b9c0bc", lw=0.8, ls="--"))
    ax.text(0.015, y + h - 0.012, label, ha="left", va="top", fontsize=8.4, fontweight="bold", color=TOOL)


if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()

    # 1. preparation
    band(ax, 0.75, 0.25, "1  preparation (once per condition)")
    box(ax, 0.05, 0.77, 0.25, 0.17, "input deck", "mesh, seed, nutrient supply,\ntime steps, rate $\\kappa$", "APDL")
    box(ax, 0.39, 0.77, 0.25, 0.17, "calibration", "TMCMC, five species,\nmeasured composition",
        "Python, JAX (GPU)")
    box(ax, 0.73, 0.77, 0.22, 0.17, "case file", "parameters $\\theta$ (MAP)\nand the Day-1 state",
        "Python script")
    arrow(ax, (0.64, 0.855), (0.73, 0.855))

    # 2. one ANSYS run
    band(ax, 0.27, 0.46, "2  one ANSYS run")
    ax.text(0.30, 0.718, "custom executable: ANSYS 2022 R2 with my Fortran linked in", ha="left", va="top",
            fontsize=8.0, style="italic", color=TOOL)
    box(ax, 0.05, 0.32, 0.25, 0.32, "partner element", "nonlocal field\n$\\phi$, $c$ at the nodes\n(Klempt 2024)",
        "Fortran, unchanged", mine=False)
    box(ax, 0.39, 0.32, 0.25, 0.32, "material routine", "per Gauss point:\n$\\Delta\\ln\\alpha=\\kappa\\,\\phi_{3D}"
        "\\int\\sum_i\\phi_i\\psi_i\\,dt$\n$\\mathbf{F}_g=\\alpha\\mathbf{I}$ $\\to$ $\\sigma$, tangent",
        "Fortran (USERMAT)")
    box(ax, 0.73, 0.32, 0.22, 0.32, "point model", "five species,\nimplicit, Newton,\none copy per\nGauss point",
        "Fortran, in the same exe")
    arrow(ax, (0.30, 0.52), (0.39, 0.52), "$\\phi_{3D}$")
    arrow(ax, (0.64, 0.52), (0.73, 0.52), "start, step", (0, 0.025))
    arrow(ax, (0.73, 0.45), (0.64, 0.45), "$\\int\\sum_i\\phi_i\\psi_i\\,dt$", (0, -0.035))
    ax.text(0.515, 0.292, "state of the point model kept as state variables of the Gauss point",
            ha="center", va="center", fontsize=7.6, color=TOOL)
    arrow(ax, (0.26, 0.77), (0.26, 0.64))
    arrow(ax, (0.84, 0.77), (0.84, 0.64))
    ax.text(0.855, 0.705, "path in an environment\nvariable", ha="left", va="center", fontsize=7.0, color=TOOL)

    # 3. post-processing
    band(ax, 0.0, 0.25, "3  results")
    box(ax, 0.05, 0.02, 0.25, 0.17, "post-processing", "stress, $\\alpha-1$, nutrient", "APDL macros")
    box(ax, 0.39, 0.02, 0.25, 0.17, "export", "seed means, composition", "Python $\\to$ JSON")
    box(ax, 0.73, 0.02, 0.22, 0.17, "figures, tables", "comparison with Eq. 36", "Python, matplotlib")
    arrow(ax, (0.30, 0.105), (0.39, 0.105)); arrow(ax, (0.64, 0.105), (0.73, 0.105))
    arrow(ax, (0.26, 0.32), (0.26, 0.19), "results file", (0.055, 0.0))

    ax.add_patch(FancyBboxPatch((0.70, -0.065), 0.012, 0.03, boxstyle="square,pad=0", fc=MINE, ec=EDGE_MINE, lw=0.8))
    ax.text(0.717, -0.05, "added in this work", fontsize=7.6, va="center", color=INK)
    ax.add_patch(FancyBboxPatch((0.86, -0.065), 0.012, 0.03, boxstyle="square,pad=0", fc=USED, ec="#8a8f8c", lw=0.8))
    ax.text(0.877, -0.05, "used as it is", fontsize=7.6, va="center", color=INK)
    ax.text(0.01, -0.05, "Build: Intel Fortran 2025 + Visual Studio, one link script.  "
            "Runs: PowerShell, one ANSYS run at a time, detached from the session.",
            fontsize=7.6, va="center", color=TOOL)
    ax.set_ylim(-0.08, 1.01)

    out = ROOT / "assets" / "fig_implementation_tools"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=170)
    print("wrote", out.with_suffix(".png"))

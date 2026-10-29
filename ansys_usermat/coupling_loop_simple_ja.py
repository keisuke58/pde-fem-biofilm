#!/usr/bin/env python3
"""coupling_loop_simple_ja.py -- one time step of the coupling as a simple loop, in Japanese, for readers who do not
work on biofilm models (Keio master's thesis and talks). Four boxes in a cycle (spreading and growth, nutrient and
antibiotic, the point model at every point, the living amount handed back) and a side box for the stress.

    python ansys_usermat/coupling_loop_simple_ja.py  -> assets/fig_coupling_loop_simple_ja.{pdf,png}
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.family": ["IPAGothic", "DejaVu Sans"], "mathtext.fontset": "stix",
                     "pdf.fonttype": 42, "font.size": 10})

C_FIELD, C_POINT, C_MECH = "#1f4e8c", "#c2571a", "#3d4a57"
FILL = {C_FIELD: "#e8eef7", C_POINT: "#fbeee5", C_MECH: "#eceff2"}


def box(ax, x, y, w, h, col, num, title, body, eq):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                fc=FILL[col], ec=col, lw=1.6))
    if num:
        ax.text(x + 0.16, y + h - 0.16, num, color="white", ha="center", va="center", fontsize=11, weight="bold",
                bbox=dict(boxstyle="circle,pad=0.25", fc=col, ec=col))
    ax.text(x + (0.34 if num else 0.08), y + h - 0.16, title, color=col, ha="left", va="center", fontsize=11.5,
            weight="bold")
    ax.text(x + 0.1, y + h - 0.36, body, ha="left", va="top", fontsize=9.5, linespacing=1.45)
    ax.text(x + w / 2, y + 0.07, eq, ha="center", va="bottom", fontsize=10, color="#333333")


def arrow(ax, p, q, col, rad=0.0, label=None, lpos=None):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=16, lw=1.8, color=col,
                                 connectionstyle=f"arc3,rad={rad}"))
    if label:
        ax.text(*lpos, label, ha="center", va="center", fontsize=9, color=col)


def main():
    fig, ax = plt.subplots(figsize=(6.3, 4.6))
    ax.set_xlim(0, 6.3); ax.set_ylim(0.6, 5.2); ax.axis("off")
    W, H = 2.55, 1.5
    # the four steps of one time step, clockwise
    box(ax, 0.15, 2.85, W, H, C_FIELD, "1", "膜が広がり、育つ",
        "バイオフィルムの量が\n栄養の濃い方へ広がり、\nその場で膨らむ（成長）",
        r"$\phi_{3D}$,  $\dot\alpha=k_\alpha\phi_{3D}$")
    box(ax, 3.55, 2.85, W, H, C_FIELD, "2", "栄養と抗菌薬が届く",
        "まわりから拡散してきて、\n生きている菌に\n食べられる",
        r"$c$,  $a$")
    box(ax, 3.55, 0.75, W, H, C_POINT, "3", "各点で菌が増減する",
        "計算の点ごとに「点モデル」。\n菌の種類ごとの量と、\n生きているか死んでいるか",
        r"$\phi_i$,  $\psi_i$  (点モデル)")
    box(ax, 0.15, 0.75, W, H, C_POINT, "4", "生きている量を返す",
        "生きている菌の量が、\n次のステップで\n栄養を食べる量になる",
        r"$\Sigma_i\,\bar\phi_i$  ($\bar\phi_i=\phi_i\psi_i$)")
    # cycle arrows
    arrow(ax, (2.75, 3.6), (3.5, 3.6), C_FIELD)
    arrow(ax, (4.83, 2.82), (4.83, 2.3), C_POINT, label="その場所の栄養・抗菌薬", lpos=(5.75, 2.56))
    arrow(ax, (3.5, 1.5), (2.75, 1.5), C_POINT)
    arrow(ax, (1.43, 2.3), (1.43, 2.82), C_FIELD, label="消費の量", lpos=(0.75, 2.56))
    ax.text(3.15, 2.56, "1 ステップ\n（くり返す）", ha="center", va="center", fontsize=10, color="#555555")
    # stress (one way, from the growth of box 1)
    ax.add_patch(FancyBboxPatch((0.05, 4.62), 2.85, 0.44, boxstyle="round,pad=0.02,rounding_size=0.06",
                                fc=FILL[C_MECH], ec=C_MECH, lw=1.4))
    ax.text(1.48, 4.84, r"応力：膨らみ $\alpha$ から計算（$\mathbf{F}_g=\alpha\mathbf{I}$）", ha="center", va="center",
            fontsize=9, color=C_MECH)
    arrow(ax, (1.43, 4.38), (1.43, 4.6), C_MECH)
    out = ROOT / "assets" / "fig_coupling_loop_simple_ja"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

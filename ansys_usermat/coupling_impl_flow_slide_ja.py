#!/usr/bin/env python3
"""coupling_impl_flow_slide_ja.py -- one outer step of the coupled model on one 16:9 slide: the time line (what is
evaluated when) and the six steps, each with its title and one equation; arrows say what they carry. The detailed
version with the ANSYS column is coupling_impl_flow_ja.py; its glossary goes to a page of its own here.

    python ansys_usermat/coupling_impl_flow_slide_ja.py
        -> assets/fig_coupling_impl_flow_slide_ja.{pdf,png}, assets/fig_coupling_glossary_ja.{pdf,png}
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from coupling_impl_flow_ja import GLOSS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.family": ["IPAGothic", "DejaVu Sans"], "mathtext.fontset": "stix", "pdf.fonttype": 42})

C = dict(field="#1f4e8c", nut="#2f7d6d", start="#8a5a9e", point="#c2571a", state="#3d4a57", mech="#555555")
F = dict(field="#e8eef7", nut="#e6f2ef", start="#f2ebf5", point="#fbeee5", state="#eef0f2", mech="#f0f0f0")

STEPS = [
    ("field", r"1. $\phi_{3D}$ と $\alpha$ を進める",
     r"$\dot\phi_{3D}=\beta\nabla^2\phi_{3D}+k_\alpha\alpha-r\,\frac{c}{k+c}\,\mathbf{n}_c\cdot\nabla\phi_{3D}$，"
     r"$\dot\alpha=k_\alpha\phi_{3D}$"),
    ("nut", r"2. 栄養 $c$ と抗菌薬 $a$",
     r"$\dot c=d\nabla^2c-g\,\phi_{3D}\,(\Sigma_i\bar\phi_i)\,c$　（$a$ も同じ式）"),
    ("start", r"3. 新しく届いた点を開始",
     r"$\phi_{3D}^{n+1}\geq\phi_{\min}$ の点：隣の組成で $\psi_i=0.999$ から"),
    ("point", r"4. 点モデル（各点で独立）",
     r"Klempt et al. 2026 の式のまま，$c^*=c^*_0\,c$，$\alpha^*=\alpha^*_0\,a$"),
    ("state", r"5. 生きている量",
     r"$\Sigma_i\bar\phi_i=\Sigma_i\,\phi_i\psi_i$"),
    ("mech", r"6. 応力（書き出す時刻だけ）",
     r"$\nabla\cdot\boldsymbol{\sigma}=\mathbf{0}$，$\boldsymbol{\sigma}=E_0(\phi_{3D}^2+f)\,\mathbb{C}:(\boldsymbol{\varepsilon}-(\alpha-1)\mathbf{I})$"),
]
ARROWS = [r"$\phi_{3D}^{n+1}$，$\alpha^{n+1}$",
          r"$c^{n+1}$，$a^{n+1}$",
          r"開始した点の $\mathbf{g}^n$",
          r"$\mathbf{g}^{n+1}$",
          r"$\phi_{3D}^{n+1}$，$\alpha^{n+1}$"]


def arrow(ax, p, q, col, ls="-", rad=0.0, lw=1.6, ms=13):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=col, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}"))


def timeline(ax, x, y, w, h):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.06", fc="white",
                                ec="#888888", lw=1.0))
    ax.text(x + 0.15, y + h - 0.12, r"外側の 1 ステップ $t_n\to t_{n+1}$（$\Delta t=0.0025\,T^*$）：何をいつ使うか",
            ha="left", va="top", fontsize=13, weight="bold", color="#333333")
    xa, xb = x + 3.6, x + w - 2.6
    yax = y + h - 0.62
    ax.plot([xa - 0.15, xb + 0.4], [yax, yax], color="#333333", lw=1.0)
    for xx, lab in ((xa, r"$t_n$"), (xb, r"$t_{n+1}$")):
        ax.plot([xx, xx], [y + 0.1, yax], color="#aaaaaa", lw=0.8, ls=":")
        ax.text(xx, yax + 0.04, lab, ha="center", va="bottom", fontsize=12)
    keys = ["field", "nut", "start", "point", "state"]
    labs = ["1. $\\phi_{3D}$，$\\alpha$", "2. $c$，$a$", "3. 開始", "4. 点モデル", "5. $\\Sigma_i\\bar\\phi_i$"]
    notes = ["$m$ 個の内側の刻み（RK4）", "後退オイラー 1 回", "$\\phi_{3D}^{n+1}$ で判定",
             "$n_{\\rm sub}$ 刻み，$c^*$ は補間", "次のステップの 2. へ（1 ステップ遅れ）"]
    ys = [yax - 0.4 - 0.38 * k for k in range(5)]
    for k, (key, lab, note, yy) in enumerate(zip(keys, labs, notes, ys)):
        col = C[key]
        ax.text(x + 0.25, yy, lab, ha="left", va="center", fontsize=12, color=col, weight="bold")
        if k == 0:
            m = 6
            xs = [xa + (xb - xa) * j / m for j in range(m + 1)]
            for x0, x1 in zip(xs[:-1], xs[1:]):
                arrow(ax, (x0, yy - 0.07), (x1, yy - 0.07), col, rad=-0.45, lw=1.2, ms=9)
        elif k == 1:
            arrow(ax, (xa, yy), (xb, yy), col, lw=1.8)
        elif k == 2:
            ax.plot([xb], [yy], marker="^", ms=9, color=col)
        elif k == 3:
            ax.plot([xa, xb], [yy, yy], color=col, lw=1.0)
            for j in range(21):
                xx = xa + (xb - xa) * j / 20
                ax.plot([xx, xx], [yy - 0.06, yy + 0.06], color=col, lw=0.8)
        else:
            ax.plot([xb], [yy], "o", ms=7, color=col)
        ax.text(xb + 0.25, yy, note, ha="left", va="center", fontsize=10.5, color=col)


def main():
    W, H = 13.33, 7.5
    fig, ax = plt.subplots(figsize=(W, H))
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    timeline(ax, 0.25, 4.7, W - 0.5, 2.7)
    # the six steps as strips
    x0, w, hb, gap = 1.75, 9.6, 0.48, 0.22
    y = 4.5
    tops = []
    for key, title, eq in STEPS:
        ax.add_patch(FancyBboxPatch((x0, y - hb), w, hb, boxstyle="round,pad=0.01,rounding_size=0.06",
                                    fc=F[key], ec=C[key], lw=1.3))
        ax.text(x0 + 0.15, y - hb / 2, title, ha="left", va="center", fontsize=13, weight="bold", color=C[key])
        ax.text(x0 + 3.55, y - hb / 2, eq, ha="left", va="center", fontsize=13.5)
        tops.append((y - hb, y))
        y -= hb + gap
    for k, ((b0, _), (_, t1)) in enumerate(zip(tops[:-1], tops[1:])):
        arrow(ax, (x0 + 0.5, b0), (x0 + 0.5, t1), "#444444", lw=1.5, ms=11)
        ax.text(x0 + 0.68, (b0 + t1) / 2, ARROWS[k], ha="left", va="center", fontsize=10.5, color="#222222")
    # 2 -> 4: c and a at both ends of the step (interpolated in the point model)
    (n0, n1), (p0, p1) = tops[1], tops[3]
    arrow(ax, (x0 - 0.02, (n0 + n1) / 2), (x0 - 0.02, (p0 + p1) / 2), C["nut"], rad=0.45, lw=1.5)
    ax.text(x0 - 0.42, (n0 + p1) / 2 - 0.12, "$c^n,c^{n+1}$\n$a^n,a^{n+1}$\n（補間）", ha="right", va="center",
            fontsize=10.5, color=C["nut"], linespacing=1.3)
    # 5 -> 2 of the next step
    (s0, s1), (q0, q1) = tops[4], tops[1]
    arrow(ax, (x0 + w + 0.02, (s0 + s1) / 2), (x0 + w + 0.02, (q0 + q1) / 2), C["point"], ls="--", rad=0.45, lw=1.5)
    ax.text(x0 + w + 0.6, (s1 + q0) / 2, "$(\\Sigma_i\\bar\\phi_i)^{n+1}$\n次のステップ\nの 2. へ", ha="left",
            va="center", fontsize=10.5, color=C["point"], linespacing=1.3)
    ax.text(x0, 0.12, "ANSYS も同じ順番（Gauss 点ごと）：場は USSFin（オリバーの要素），点モデルと $\\alpha$ は usermat。"
            "詳しい離散化は別図", ha="left", va="bottom", fontsize=10.5, color="#555555")
    out = ROOT / "assets" / "fig_coupling_impl_flow_slide_ja"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print("wrote", out.with_suffix(".pdf"))
    glossary()


def glossary():
    """the glossary of the detailed figure on a page of its own"""
    nl = [1 + m.count("\n") for _, m in GLOSS]
    h = 0.36 * sum(nl) + 0.12 * len(GLOSS) + 0.7
    fig, ax = plt.subplots(figsize=(13.33, h))
    ax.set_xlim(0, 13.33); ax.set_ylim(0, h); ax.axis("off")
    ax.text(0.2, h - 0.15, "用語（実装フロー図）", ha="left", va="top", fontsize=13, weight="bold", color="#333333")
    yy = h - 0.6
    for (term, mean), n in zip(GLOSS, nl):
        ax.text(0.3, yy, term, ha="left", va="top", fontsize=10.5, weight="bold", color="#333333")
        ax.text(3.6, yy, mean, ha="left", va="top", fontsize=10, color="#333333", linespacing=1.45)
        yy -= 0.36 * n + 0.12
    out = ROOT / "assets" / "fig_coupling_glossary_ja"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

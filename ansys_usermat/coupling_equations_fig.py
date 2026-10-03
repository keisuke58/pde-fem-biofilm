#!/usr/bin/env python3
"""The one-way coupling written out as equations: one increment at one
Gauss point. A companion to coupling_schematic.py (the overview); this one
shows each step the coupled ANSYS run takes, with the equation it uses and
where that equation comes from.

    python ansys_usermat/coupling_equations_fig.py        -> assets/fig_coupling_equations.png
    python ansys_usermat/coupling_equations_fig.py --ja   -> assets/fig_coupling_equations_ja.png

(1) The partner's element updates its nodal fields (Klempt et al. 2024,
Eq. 34/35). (2) At each Gauss point the material routine integrates Eq. 36
and computes the stress from F_g = alpha I. (3) For two species, the point
model of Klempt et al. 2026 is scaled to the field's amount, advanced by
s dt and scaled again; it only decides the composition. Nothing in (3) acts
back on (1) or (2). Values not taken from a paper are marked in the footer.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402

OUT = HERE.parent / "assets" / "fig_coupling_equations.png"
OUT_JA = HERE.parent / "assets" / "fig_coupling_equations_ja.png"
# Japanese: Times for Latin letters, a Mincho for kana and kanji (font fallback)
JA_SERIF = ["IPAexMincho", "Yu Mincho", "MS Mincho", "IPAMincho"]
INK, MUTED = "#1f2933", "#616e7c"
COL = {"field": ("#e8f1fb", "#1d4e89"), "mech": ("#eef1f4", "#3e4c59"),
       "micro": ("#fdf0e6", "#b45309")}
Y0, H, HEAD = 0.17, 0.70, 0.065
PX = {"field": (0.012, 0.300), "mech": (0.350, 0.300), "micro": (0.688, 0.300)}


ITALIC = ["italic"]


SPLIT = [False]          # set for Japanese: mathtext cannot fall back to a CJK font


def put(ax, x, y, s, ha="left", va="center", **kw):
    """ax.text, but a string that mixes text and $math$ is drawn piece by
    piece when SPLIT is set: mathtext renders the text outside $...$ in its
    own fonts, which have no kana or kanji."""
    if not SPLIT[0] or "$" not in s or s.isascii():
        return ax.text(x, y, s, ha=ha, va=va, **kw)
    if "\n" in s:
        rows = s.split("\n")
        step = 0.027 * kw.get("linespacing", 1.0) * kw.get("fontsize", 11) / 10.5
        y0 = y + step * (len(rows) - 1) / 2
        for i, row in enumerate(rows):
            put(ax, x, y0 - i * step, row, ha=ha, va=va, **kw)
        return None
    kw.pop("linespacing", None)
    pieces = [q for q in re.split(r"(\$[^$]*\$)", s) if q]
    fig = ax.figure
    rend = fig.canvas.get_renderer()
    axw = ax.get_window_extent(rend).width
    texts, widths = [], []
    for q in pieces:
        t = ax.text(0, y, q, ha="left", va="center_baseline", **kw)
        texts.append(t)
        widths.append(t.get_window_extent(rend).width / axw)
    total = sum(widths)
    cur = x - (total / 2 if ha == "center" else total if ha == "right" else 0)
    for t, w in zip(texts, widths):
        t.set_x(cur)
        cur += w
    return texts


def panel(ax, key, num, title, sub):
    x, w = PX[key]
    fill, acc = COL[key]
    ax.add_patch(FancyBboxPatch((x, Y0), w, H, boxstyle="round,pad=0,rounding_size=0.012",
                                fc=fill, ec=acc, lw=1.2, zorder=1))
    ax.add_patch(FancyBboxPatch((x, Y0 + H - HEAD), w, HEAD,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc=acc, ec=acc, lw=1.2, zorder=2))
    ax.add_patch(Rectangle((x, Y0 + H - HEAD), w, HEAD / 2, fc=acc, ec="none", zorder=2))
    ax.plot(x + 0.024, Y0 + H - HEAD / 2, "o", ms=21, mfc="white", mec="none", zorder=3)
    ax.text(x + 0.024, Y0 + H - HEAD / 2, num, ha="center", va="center", fontsize=12,
            weight="bold", color=acc, zorder=4)
    ax.text(x + 0.048, Y0 + H - HEAD / 2, title, ha="left", va="center", fontsize=13.5,
            weight="bold", color="white", zorder=4)
    put(ax, x + w / 2, Y0 + H - HEAD - 0.027, sub, ha="center", va="center",
            fontsize=10.5, style=ITALIC[0], color=MUTED, zorder=4)


def lines(ax, key, y, rows, dy=0.058):
    """rows: (kind, text). kind 'step' = numbered label, 'eq' = equation,
    'note' = small grey remark."""
    x, w = PX[key]
    acc = COL[key][1]
    for kind, txt in rows:
        if kind == "step":
            put(ax, x + 0.016, y, txt, ha="left", va="center", fontsize=11,
                    weight="bold", color=acc, zorder=4)
            y -= dy * 0.72
        elif kind == "eq":
            put(ax, x + 0.030, y, txt, ha="left", va="center", fontsize=12.5,
                    color=INK, zorder=4)
            y -= dy
        elif kind == "note":
            put(ax, x + 0.030, y + 0.012, txt, ha="left", va="center", fontsize=9.5,
                    color=MUTED, zorder=4)
            y -= dy * 0.62
        elif kind == "gap":
            y -= dy * 0.35
    return y


def arrow(ax, p, q, color=INK, ls="-", lw=1.6, rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=15, color=color,
                                 lw=lw, ls=ls, connectionstyle=f"arc3,rad={rad}", zorder=5))


EQ = {
    "c": r"$d\,\Delta c = g\,\phi,\qquad c = 1$",
    "phi": r"$\dot\phi = \beta\,\Delta\phi + k_\alpha\,\alpha"
           r" - r\,\dfrac{c}{k+c}\,\nabla\phi\cdot\dfrac{\nabla c}{|\nabla c|}$",
    "a36": r"$\dot\alpha = k_\alpha\,\phi$",
    "gp": r"$\phi_{3D} = \sum_I N_I\,\phi_I,\qquad 0 \leq \phi_{3D} \leq 1$",
    "aup": r"$\alpha_{n+1} = \alpha_n + k_\alpha\,\phi_{3D}\,\Delta t,\qquad \alpha(0) = 1$",
    "kin": r"$\mathbf{F}_g = \alpha_{n+1}\,\mathbf{I},\qquad"
           r"\mathbf{F}_e = \mathbf{F}\,\mathbf{F}_g^{-1},\qquad J_e = \det\mathbf{F}_e$",
    "E": r"$E(\phi_{3D}) = (\phi_{3D}^{\,2} + f)\,E_{bio},\qquad \nu = 0.49$",
    "sig": r"$\boldsymbol{\sigma} = \hat{\boldsymbol{\sigma}}(\mathbf{F}_e;\,E,\nu)$",
    "eq": r"$\int_\Omega \boldsymbol{\sigma} : \delta\boldsymbol{\varepsilon}\;dv = 0$",
    "hat": r"$\hat\phi = \mathrm{min}(\phi_{3D},\,\phi_{cap})$",
    "in": r"$\phi_i \leftarrow \phi_i\,\dfrac{\hat\phi}{\phi_1+\phi_2},\qquad \phi_0 = 1 - \hat\phi$",
    "pm1": r"$(\eta_{\phi,i}+\eta_i\psi_i^2)\,\dot\phi_i + \eta_i\phi_i\psi_i\,\dot\psi_i"
           r" = c^*\psi_i\sum_j A_{ij}\phi_j\psi_j - \gamma$",
    "pm2": r"$\eta_i(\phi_i\psi_i\,\dot\phi_i + \phi_i^2\,\dot\psi_i)"
           r" = c^*\phi_i\sum_j A_{ij}\phi_j\psi_j - \alpha^* b_i\psi_i$",
}

TEXT = {
    "en": {
        "title": r"One time increment $\Delta t$ at one Gauss point of the coupled ANSYS run",
        "p1": ("Element: growth field", "Klempt et al. 2024, partner's element (unchanged)"),
        "p2": ("Material routine: growth, stress", "added in this work"),
        "p3": ("Point model: composition", "Klempt et al. 2026, added in this work"),
        "field": [
            ("step", "a  nutrient, quasi-static (Eq. 35)"),
            ("eq", EQ["c"] + " on the source"),
            ("step", "b  biofilm field (Eq. 34)"),
            ("eq", EQ["phi"]),
            ("note", "diffusion, local source, front growth towards the nutrient"),
            ("eq", EQ["a36"] + "  (Eq. 36, inside the element)"),
            ("step", "c  value at the Gauss point"),
            ("eq", EQ["gp"]),
            ("note", r"$N_I$: shape functions, $\phi_I$: nodal values"),
        ],
        "mech": [
            ("step", "a  growth (Eq. 36, explicit)"),
            ("eq", EQ["aup"]),
            ("step", "b  kinematics"),
            ("eq", EQ["kin"]),
            ("note", r"growth is $\alpha - 1$ (stored as $\alpha - 1$ in the code)"),
            ("step", "c  stiffness and stress"),
            ("eq", EQ["E"]),
            ("eq", EQ["sig"] + "  (neo-Hookean)"),
            ("step", "d  equilibrium (ANSYS Newton iteration)"),
            ("eq", EQ["eq"]),
            ("note", "a free part grows without stress; constraint gives stress"),
            ("note", r"next increment: $\alpha_{n+1}$, $\phi_i$, $\psi_i$ carried as state variables"),
        ],
        "micro": [
            ("step", "a  amount the point model sees"),
            ("eq", EQ["hat"]),
            ("step", "b  scale in"),
            ("eq", EQ["in"]),
            ("step", r"c  advance by $s\,\Delta t$ ($n_{sub}$ implicit steps)"),
            ("eq", EQ["pm1"]),
            ("eq", EQ["pm2"]),
            ("note", r"$n_{sub} = \lceil s\,\Delta t / 10^{-4}\rceil$; barrier terms omitted"),
            ("step", "d  scale out and store"),
            ("eq", r"$\phi_1 + \phi_2 = \hat\phi$;  $\phi_i,\ \psi_i$ as state variables"),
            ("eq", r"output: $\phi_1/(\phi_1+\phi_2)$,  $\psi_i$"),
        ],
        "top": r"$\phi_{3D}$: the amount only",
        "none": "composition does not act back on (1) or (2): one-way coupling (two-way is the outlook)",
        "foot": r"Klempt et al. 2024 Table 2: $k_\alpha = 10^{-3}$ per $T^*$, $E_{bio} = 10$ Pa, $\nu = 0.49$.   "
                r"Klempt et al. 2026: $c^* = 100$, $\alpha^* = 10$, $A_{ij}$, $b_i$, $\eta_i$ per Table 1 case.   "
                "\nNot from a paper: "
                r"$s = 0.15$, $\phi_{cap} = 0.9$, $f = 10^{-3}$ (void stiffness), "
                r"and $\beta$, $r$, $k$, $d$, $g$ from the partner's example input.",
    },
    "ja": {
        "title": r"連成計算の1ステップ（時間増分 $\Delta t$、ある1つの積分点）",
        "p1": ("要素：成長場", "Klempt et al. 2024、共同研究者の要素（変更なし）"),
        "p2": ("材料ルーチン：成長と応力", "本研究で追加"),
        "p3": ("点モデル：菌種の組成", "Klempt et al. 2026、本研究で追加"),
        "field": [
            ("step", "a  栄養の拡散（Eq. 35、準静的）"),
            ("eq", EQ["c"] + "（栄養源）"),
            ("note", "栄養は成長よりずっと速く拡散するので、毎ステップ定常で解く"),
            ("step", "b  バイオフィルムの量（Eq. 34）"),
            ("eq", EQ["phi"]),
            ("note", "左から：拡散、その場での増加、栄養の方向への前線の成長"),
            ("eq", EQ["a36"] + "（Eq. 36、要素の内部）"),
            ("step", "c  積分点での値"),
            ("eq", EQ["gp"]),
            ("note", r"節点の値 $\phi_I$ を形状関数 $N_I$ で補間した、その点の量"),
        ],
        "mech": [
            ("step", "a  成長（Eq. 36、陽解法）"),
            ("eq", EQ["aup"]),
            ("note", r"バイオフィルムが多い点ほど $\alpha$ が速く増える"),
            ("step", "b  変形の分解"),
            ("eq", EQ["kin"]),
            ("note", r"全変形 = 弾性 × 成長。成長量は $\alpha-1$（コード内も $\alpha-1$ で保存）"),
            ("step", "c  剛性と応力"),
            ("eq", EQ["E"]),
            ("eq", EQ["sig"] + "（neo-Hooke）"),
            ("note", "応力は弾性部分 $\\mathbf{F}_e$ だけから生じる"),
            ("step", "d  つり合い（ANSYS の Newton 反復）"),
            ("eq", EQ["eq"]),
            ("note", "自由に膨らめれば応力はゼロ、周りに拘束されると応力が出る"),
        ],
        "micro": [
            ("step", "a  点モデルに渡す量"),
            ("eq", EQ["hat"]),
            ("note", r"満杯の領域では上限 $\phi_{cap}$ で頭打ち"),
            ("step", "b  量を場に合わせる"),
            ("eq", EQ["in"]),
            ("note", "組成（比）は保ったまま、合計だけを場の量にそろえる"),
            ("step", r"c  時間 $s\,\Delta t$ だけ進める（$n_{sub}$ 回の陰解法）"),
            ("eq", EQ["pm1"]),
            ("eq", EQ["pm2"]),
            ("note", r"$n_{sub} = \lceil s\,\Delta t / 10^{-4}\rceil$、境界項は省略"),
            ("step", "d  もう一度合わせて保存"),
            ("eq", r"$\phi_1 + \phi_2 = \hat\phi$、$\phi_i,\ \psi_i$ を状態変数に保存"),
            ("eq", r"出力：$\phi_1/(\phi_1+\phi_2)$、$\psi_i$"),
        ],
        "top": r"$\phi_{3D}$：量だけを渡す",
        "none": "組成は (1)、(2) に戻らない：片方向連成（双方向は今後の課題）",
        "foot": r"論文の値　Klempt et al. 2024 Table 2：$k_\alpha = 10^{-3}$ /$T^*$、$E_{bio} = 10$ Pa、$\nu = 0.49$。"
                r"　Klempt et al. 2026：$c^* = 100$、$\alpha^* = 10$、$A_{ij}$、$b_i$、$\eta_i$ は Table 1 の各ケース。"
                "\n論文にない値（仮定）　"
                r"$s = 0.15$、$\phi_{cap} = 0.9$、$f = 10^{-3}$（空隙の剛性）、"
                r"$\beta$、$r$、$k$、$d$、$g$ は共同研究者の入力例の値。",
    },
}


def use_japanese():
    """Put a Mincho after Times in the serif list; register it with
    matplotlib from fontconfig when matplotlib's own cache lacks it."""
    import subprocess
    from matplotlib import font_manager as fm
    known = {f.name for f in fm.fontManager.ttflist}
    for name in JA_SERIF:
        if name in known:
            break
        try:
            path = subprocess.run(["fc-match", "-f", "%{file}", name], capture_output=True,
                                  text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            continue
        if path and fm.FontProperties(fname=path).get_name() == name:
            fm.fontManager.addfont(path)
            break
    # fallback works across an explicit family list, not across the generic "serif"
    have = {f.name for f in fm.fontManager.ttflist}
    plt.rcParams["font.family"] = [f for f in figstyle.SERIF[:2] + JA_SERIF + figstyle.SERIF[2:]
                                   if f in have] or ["serif"]
    ITALIC[0] = "normal"
    SPLIT[0] = True


def main(lang="en"):
    figstyle.apply()
    T = TEXT[lang]
    if lang == "ja":
        use_japanese()
    fig = plt.figure(figsize=(16, 8.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    put(ax, 0.5, 0.968, T["title"], ha="center", va="center", fontsize=15, color=INK,
            weight="bold")
    panel(ax, "field", "1", *T["p1"])
    panel(ax, "mech", "2", *T["p2"])
    panel(ax, "micro", "3", *T["p3"])
    top = Y0 + H - HEAD - 0.075
    dy = 0.058 if lang == "en" else 0.052
    for key in ("field", "mech", "micro"):
        lines(ax, key, top, T[key], dy=dy)

    # flows
    xf, wf = PX["field"]; xm, wm = PX["mech"]; xc, wc = PX["micro"]
    ymid = Y0 + 0.30
    arrow(ax, (xf + wf, ymid), (xm, ymid))
    ax.text((xf + wf + xm) / 2, ymid + 0.025, r"$\phi_{3D}$", ha="center", fontsize=13, color=INK)
    ax.plot([xf + wf * 0.55, xf + wf * 0.55, xc + wc * 0.5], [Y0 + H, 0.902, 0.902],
            color=INK, lw=1.6, zorder=5)
    arrow(ax, (xc + wc * 0.5, 0.902), (xc + wc * 0.5, Y0 + H))
    put(ax, xm + wm / 2, 0.913 if lang == "en" else 0.921, T["top"], ha="center", fontsize=12, color=INK)
    # no feedback: a crossed, dashed line without arrowhead
    yb = Y0 - 0.045
    ax.plot([xc + wc * 0.5, xc + wc * 0.5, xf + wf * 0.5, xf + wf * 0.5], [Y0, yb, yb, Y0],
            color=MUTED, lw=1.3, ls="--", zorder=5)
    ax.plot(xm + wm * 0.5, yb, "o", ms=17, mfc="white", mec=MUTED, mew=1.3, zorder=6)
    ax.text(xm + wm * 0.5, yb, r"$\times$", ha="center", va="center", fontsize=15,
            color=MUTED, zorder=7)
    put(ax, xm + wm * 0.5, yb - 0.032, T["none"], ha="center", fontsize=11,
            style=ITALIC[0], color=MUTED)
    put(ax, 0.5, 0.035 if lang == "en" else 0.028, T["foot"], ha="center", va="center", fontsize=10.5, color=MUTED,
            linespacing=1.6)
    out = OUT_JA if lang == "ja" else OUT
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    main("ja" if "--ja" in sys.argv else "en")

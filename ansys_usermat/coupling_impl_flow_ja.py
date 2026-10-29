#!/usr/bin/env python3
"""coupling_impl_flow_ja.py -- the implementation of one outer step t_n -> t_{n+1} of the coupled model, precise:
state per point, discretisation, solver, order and lags; left the finite-volume code
(ansys_usermat/apdl/coupled_living_consumption.py), right the corresponding part of the ANSYS element.

    python ansys_usermat/coupling_impl_flow_ja.py -> assets/fig_coupling_impl_flow_ja.{pdf,png}
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.family": ["IPAGothic", "DejaVu Sans"], "mathtext.fontset": "stix", "pdf.fonttype": 42})

C = dict(state="#3d4a57", field="#1f4e8c", nut="#2f7d6d", point="#c2571a", start="#8a5a9e", mech="#555555")
F = dict(state="#eef0f2", field="#e8eef7", nut="#e6f2ef", point="#fbeee5", start="#f2ebf5", mech="#f0f0f0")

# (key, title, finite-volume text, ANSYS text, height)
STEPS = [
    ("state", r"状態（時刻 $t_n$、セルごと）",
     r"場：$\phi_{3D}^n,\ \alpha^n,\ c^n,\ a^n$　　点モデル：$\mathbf{g}^n=(\phi_1,\phi_2,\phi_0,\psi_1,\psi_2,\gamma)$、開始フラグ"
     "\n" r"生きている量 $(\Sigma_i\bar\phi_i)^n$，$\bar\phi_i=\phi_i\psi_i$（$\mathbf{g}^n$ から。始まっていないセルは $0$）",
     r"同じ量を Gauss 点の状態変数に持つ" "\n" r"（$\alpha-1$、$\mathbf{g}$ は usermat の ustatev）", 1.5),
    ("field", r"1. $\phi_{3D}$ と $\alpha$ を進める（RK4、内側の刻み $\Delta t/m$）",
     r"$\dot\phi_{3D}=\beta\nabla_h^2\phi_{3D}+k_\alpha\alpha-r\frac{c}{k+c}\,\mathbf{n}_c\cdot\nabla_h\phi_{3D}-\kappa_p\,[\phi_{3D}]_\pm$,"
     r"　$\dot\alpha=k_\alpha\phi_{3D}$" "\n"
     r"$\nabla_h^2$：7 点（壁は鏡映＝流れゼロ）。前線の項：MUSCL＋minmod の風上（TVD）" "\n"
     r"$\mathbf{n}_c=\nabla c/\sqrt{|\nabla c|^2+\varepsilon^2}$、$\varepsilon^2=10^{-6}$。$c$ は内側の各刻みで解き直す（前線が栄養に追従）" "\n"
     r"$m=\lceil\Delta t/\delta_{\max}\rceil$、$\delta_{\max}=\min\{0.15\,h^2/\beta,\ 0.4\,h/v,\ 2.2/(12\beta/h^2+2v/h)\}$、$v=r/2$（前線の最大速さ）",
     r"USSFin（オリバーの要素）：陽的に更新" "\n" r"$\nabla^2$、$\nabla$ は NEM" "\n"
     r"（近傍 30 点の重み付き最小二乗）" "\n" r"前線の項は風上側の MLS 勾配" "\n" r"$\alpha$：usermat で Gauss 点ごとに積分", 2.65),
    ("nut", r"2. 栄養 $c^{n+1}$ と抗菌薬 $a^{n+1}$（式 35、$\dot c$ 付き）",
     r"$\dot c=d\nabla^2c-g\,\phi_{3D}\,(\Sigma_i\bar\phi_i)\,c$ を後退オイラーで：$\dfrac{c^{n+1}-c^n}{d\,\Delta t}=\nabla_h^2c^{n+1}-\dfrac{g}{d}\,\phi_{3D}^{n+1}(\Sigma_i\bar\phi_i)^n c^{n+1}$、供給面で $c=1$" "\n"
     r"$a$ も同じ式（外の面で $a=1$、$T^*\geq1/3$）。AMG 前処理つき CG で解く。$(\Sigma_i\bar\phi_i)^n$ は 1 ステップ前の値" "\n"
     r"$d=10^8$ mm$^2/T^*$（--transient-c 1e8）。左辺 $\approx10^{-6}$；29 日までは左辺 0 の極限（差 $\leq0.1\,\%$）",
     r"同じ式 35 を要素の" "\n" r"2 つの栄養の場で、$\dot c$ 付き" "\n" r"（LIVING_CONS で消費に $\Sigma_i\bar\phi_i$ を掛ける）", 2.1),
    ("start", r"3. 新しく届いたセル（$\phi_{3D}^{n+1}\geq\phi_{\min}=10^{-3}$、未開始）",
     r"定着済みの隣の $\mathbf{g}$ を $\phi_{3D}$ の重みで平均し、量 $\phi_i$ を生きている量の比で配分、"
     r"$\psi_i=0.999$、$\gamma=0$" "\n"
     r"前線は 1 ステップに数セル進むので、埋まるまで層ごとに繰り返す（仮定）",
     r"ANSYS の計算は位置で決める規則" "\n" r"（受け継ぐ規則は Python と" "\n" r"Abaqus の要素）", 1.6),
    ("point", r"4. 点モデル（開始済みの各セル、互いに独立）",
     r"内側の刻み $k=1..n_{\rm sub}$、$w=k/n_{\rm sub}$：$c^*=c^*_0[(1-w)c^n+wc^{n+1}]$、$\alpha^*=\alpha^*_0[(1-w)a^n+wa^{n+1}]$" "\n"
     r"自分の時間で $s\Delta t/n_{\rm sub}=10^{-4}$ だけ進める：陰的オイラー、$\mathbf{Q}(\mathbf{g})=\mathbf{0}$ を Newton 法（解析的な接線）" "\n"
     r"Klempt et al. 2026 の式のまま（Fortran の ecology_native.f を ctypes で）→ $\mathbf{g}^{n+1}$",
     r"usermat が Gauss 点ごとに" "\n" r"同じ Fortran を呼ぶ" "\n" r"（eco_resjac、eco_solve）", 1.9),
    ("state", r"5. 生きている量の更新",
     r"$(\Sigma_i\bar\phi_i)^{n+1}=\Sigma_i\,\phi_i^{n+1}\psi_i^{n+1}$（次のステップの手順 2 で使う）",
     r"同じ（状態変数から）", 0.95),
    ("mech", r"6. 応力（書き出す時刻だけ。成長には戻さない）",
     r"$\phi_{3D}$ と $\alpha-1$ を $16^3$ の六面体要素ごとにセル平均 → $\mathbf{F}_g=\alpha\mathbf{I}$、ネオフック、"
     r"$E=E_0(\phi_{3D}^2+f)$、$\nu=0.49$" "\n" r"選択的低減積分で静的につり合いを解く（$32^3$、$64^3$ の要素でも種の応力は $0.2\,\%$ 以内）",
     r"usermat：$\mathbf{F}_e=\mathbf{F}\mathbf{F}_g^{-1}$、" "\n" r"接線は $\mathbf{F}$ の摂動（$10^{-7}$）", 1.6),
]

# what each arrow carries (between consecutive boxes)
ARROWS = [r"$\phi_{3D}^n,\ \alpha^n,\ c^n$",
          r"$\phi_{3D}^{n+1},\ \alpha^{n+1}$",
          r"$\phi_{3D}^{n+1}$（どのセルに届いたか）",
          r"新しく開始した点の $\mathbf{g}^n$",
          r"$\mathbf{g}^{n+1}$",
          r"$\phi_{3D}^{n+1},\ \alpha^{n+1}$（書き出す時刻）"]

# glossary under the flow (plain Japanese)
GLOSS = [
    ("ラプラシアン（7 点）", r"$\nabla_h^2u_i=(\sum_{j\in N(i)}u_j-6u_i)/h^2$（$N(i)$：6 つの隣）。2 次精度。拡散の項 $\beta\nabla^2\phi_{3D}$ の離散化"),
    ("鏡映（壁）", r"壁の外の仮想セルに $u_{\rm ghost}=u_{\rm in}$（壁の内側のセルの値） を置く。壁で $\partial u/\partial n=0$（流束ゼロ）を 2 次精度で満たす"),
    ("風上（upwind）", r"移流 $\mathbf{v}\cdot\nabla u$ のセル面の値を上流側から取る。1 次なら単調だが数値拡散 $\approx|\mathbf{v}|h/2$ を持つ"),
    ("MUSCL、minmod、TVD", r"MUSCL：セル内を勾配で線形に再構成し面の値を 2 次精度で作る。minmod$(a,b)$：同符号なら絶対値の小さい方、"
     "\n" r"異符号なら 0（極値で勾配を 0 にする）。これで全変動 $\sum_i|u_{i+1}-u_i|$ が増えない（TVD：新しい極値を作らない）"),
    ("前線の向きの正則化", r"$\mathbf{n}_c=\nabla c/\sqrt{|\nabla c|^2+\varepsilon^2}$、$\varepsilon^2=10^{-6}$。$|\nabla c|\gg\varepsilon$ で単位ベクトル、$|\nabla c|\to0$ で連続に 0"),
    ("RK4", r"4 段 4 次の陽的ルンゲ＝クッタ法。局所誤差 $O(\Delta t^5)$。拡散数 $\beta\Delta t/h^2$ とクーラン数 $|\mathbf{v}|\Delta t/h$ に安定の上限"),
    ("後退オイラー（式 35）", r"$\dot c$ を $(c^{n+1}-c^n)/\Delta t$ で置き、右辺を新しい時刻で評価する陰解法。どの $\Delta t$ でも安定で、1 ステップ 1 回の線形解。"
     "\n" r"拡散の時間 $\ell^2/d\approx10^{-8}T^*\ll\Delta t$ なので左辺は $10^{-6}$ 程度。左辺 0 の極限との差は $32^3$ で $\leq0.1\,\%$"),
    ("1 ステップ遅れ", r"手順 2 で $(\Sigma_i\bar\phi_i)^{n+1}$ の代わりに $(\Sigma_i\bar\phi_i)^n$ を使う（陽的な結合）。誤差は $O(\Delta t)$。2 回反復しても差は刻みの誤差以下"),
    ("疎行列、CG、AMG", r"$-\nabla_h^2+q$ の行列は対称正定値で各行の非零は 7 個。CG：対称正定値の連立方程式のクリロフ部分空間の反復法。"
     "\n" r"AMG：行列から粗い階層を作る前処理で、反復回数を格子の大きさにほぼよらなくする"),
    ("陰的オイラー法", r"$(\mathbf{g}^{k+1}-\mathbf{g}^k)/\delta t=\mathbf{f}(\mathbf{g}^{k+1})$。1 次精度、A 安定（刻みを大きくしても発散しない）"),
    ("Newton 法、ヤコビアン", r"$\mathbf{Q}(\mathbf{g})=\mathbf{0}$ を $\mathbf{g}\leftarrow\mathbf{g}-\mathbf{J}^{-1}\mathbf{Q}$、$\mathbf{J}=\partial\mathbf{Q}/\partial\mathbf{g}$ で解く。解の近くで 2 次収束"),
    ("ctypes", r"Python 標準の外部関数インターフェース。コンパイルした Fortran の共有ライブラリを"
     "\n" r"倍精度の配列のまま呼ぶ（ANSYS と同じルーチン）"),
    ("Gauss 点、状態変数", r"要素の積分に使う $2\times2\times2$ の Gauss–Legendre 点。状態変数：各点で増分をまたいで保存する内部変数（ustatev）"),
    ("NEM、MLS", r"近傍 30 個の Gauss 点の値に 2 次の Taylor 展開（未知数 9 個）を、Gauss 型の重み（幅 $\beta^*=0.2$ mm）の"
     "\n" r"最小二乗で当てはめ、その係数を $\nabla u$、$\nabla^2u$ とする（moving least squares）"),
    ("ネオフック", r"$W=C_{10}(\bar I_1-3)+(J-1)^2/D_1$：等容の部分と体積の部分に分けた超弾性。$C_{10}=\mu/2$"),
    (r"選択的低減積分（$\bar{\mathbf{B}}$ 法）", r"体積ひずみの項だけを要素で平均して積分し、$\nu\to0.5$ での体積ロッキング（剛性の過大評価）を避ける"),
    (r"$\mathbf{F}$ の摂動", r"$F_{ij}$ を相対 $10^{-7}$ だけずらし、応力の変化の前進差分で接線剛性 $\partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon}$ を作る（Sun et al. 2008）"),
]


def arrow(ax, p, q, col, ls="-", rad=0.0, lw=1.1):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, lw=lw, color=col, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}"))


def timeline(ax, x, y, w, h):
    """one outer step t_n -> t_{n+1}: when each quantity is evaluated"""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.05", fc="white",
                                ec="#777777", lw=0.9))
    ax.text(x + 0.12, y + h - 0.12, r"時間軸：外側の 1 ステップ $t_n\to t_{n+1}$ の中で、何をどの時刻で使うか", ha="left",
            va="top", fontsize=8.6, weight="bold", color="#333333")
    xa, xb = x + 4.7, x + w - 2.0
    yax = y + h - 0.95
    ax.plot([xa - 0.2, xb + 1.4], [yax, yax], color="#333333", lw=0.9)
    for xx, lab in ((xa, r"$t_n$"), (xb, r"$t_{n+1}$")):
        ax.plot([xx, xx], [y + 0.12, yax], color="#999999", lw=0.6, ls=":")
        ax.text(xx, yax + 0.06, lab, ha="center", va="bottom", fontsize=8)
    rows = [("field", r"手順 1：$\phi_{3D}$、$\alpha$（RK4）"),
            ("nut", r"手順 2：$c$、$a$（後退オイラー）"),
            ("point", r"手順 4：点モデル"),
            ("state", r"手順 5：$\Sigma_i\bar\phi_i$")]
    ys = [yax - 0.55 - 0.85 * k for k in range(4)]
    for (key, lab), yy in zip(rows, ys):
        ax.text(x + 0.2, yy, lab, ha="left", va="center", fontsize=7.4, color=C[key], weight="bold")
    # 1: m inner RK4 steps, c re-solved in each
    yy, col = ys[0], C["field"]
    m = 6
    xs = [xa + (xb - xa) * j / m for j in range(m + 1)]
    for x0, x1 in zip(xs[:-1], xs[1:]):
        arrow(ax, (x0, yy - 0.05), (x1, yy - 0.05), col, rad=-0.5, lw=0.9)
    ax.text((xa + xb) / 2, yy - 0.12, r"$m$ 個の内側の刻み $\delta=\Delta t/m$、各刻みの始めに $c$ を解き直す", ha="center", va="top",
            fontsize=6.5, color=col)
    # 2: one backward Euler step to t_{n+1}
    yy, col = ys[1], C["nut"]
    arrow(ax, (xa, yy), (xb, yy), col, lw=1.2)
    ax.plot([xb], [yy], "o", ms=3.5, color=col)
    ax.text(xa + 0.1, yy - 0.1, r"1 回の線形解（$\phi_{3D}^{n+1}$、$(\Sigma_i\bar\phi_i)^n$）→ $c^{n+1}$、$a^{n+1}$", ha="left",
            va="top", fontsize=6.5, color=col)
    # 4: n_sub point-model steps, c* and alpha* interpolated between c^n and c^{n+1}
    yy, col = ys[2], C["point"]
    n = 20
    for j in range(n + 1):
        xx = xa + (xb - xa) * j / n
        ax.plot([xx, xx], [yy - 0.06, yy + 0.06], color=col, lw=0.6)
    ax.plot([xa, xb], [yy, yy], color=col, lw=0.8)
    ax.plot([xa, xb], [yy + 0.14, yy + 0.3], color=C["nut"], lw=0.9)
    ax.text(xa - 0.08, yy + 0.14, r"$c^n$", ha="right", va="center", fontsize=6.6, color=C["nut"])
    ax.text(xb + 0.08, yy + 0.3, r"$c^{n+1}$", ha="left", va="center", fontsize=6.6, color=C["nut"])
    ax.text(xa + 0.1, yy - 0.1, r"$n_{\rm sub}=\lceil s\Delta t/10^{-4}\rceil$ 刻み。$c^*$、$\alpha^*$ は $c$、$a$ の線形補間（緑の線）",
            ha="left", va="top", fontsize=6.5, color=col)
    # 5: living fraction at t_{n+1}, used in step 2 of the next step
    yy, col = ys[3], C["state"]
    ax.plot([xb], [yy], "o", ms=3.5, color=col)
    arrow(ax, (xb + 0.08, yy + 0.04), (xb + 1.3, ys[1] - 0.04), C["point"], ls="--", rad=0.45, lw=0.9)
    ax.text(xb - 0.15, yy, r"$t_{n+1}$ の値を次のステップの手順 2 で使う（1 ステップ遅れ）", ha="right", va="center",
            fontsize=6.5, color=C["point"])


def main():
    W_L, W_R, X0, GAP = 12.4, 4.9, 1.0, 0.45
    GY, TH = 0.5, 4.6
    NL = [1 + m.count("\n") for _, m in GLOSS]
    GH = 0.42 * sum(NL) + 0.12 * len(GLOSS) + 0.7
    tot = TH + 0.35 + sum(s[4] for s in STEPS) + GY * (len(STEPS) - 1) + 0.9 + GH
    WT = X0 + W_L + GAP + W_R + 0.1
    sc = 8.3 / 17.95                                   # inches per unit, as before the left margin
    fig, ax = plt.subplots(figsize=(sc * WT, sc * tot))
    ax.set_xlim(0, WT); ax.set_ylim(0, tot); ax.axis("off")
    ax.text(X0 + W_L / 2, tot - 0.25, r"有限体積のコード（Python、外側の刻み $\Delta t=0.0025\,T^*$）", ha="center",
            va="center", fontsize=9.5, weight="bold")
    ax.text(X0 + W_L + GAP + W_R / 2, tot - 0.25, "ANSYS の要素で対応する所", ha="center", va="center", fontsize=9.5,
            weight="bold")
    yt = tot - 0.55
    timeline(ax, X0, yt - TH, W_L, TH)
    xr = X0 + W_L + GAP
    ax.add_patch(FancyBboxPatch((xr, yt - TH), W_R, TH, boxstyle="round,pad=0.01,rounding_size=0.05", fc="white",
                                ec="#777777", lw=0.9))
    ax.text(xr + 0.1, yt - 0.12, "量の置き場所", ha="left", va="top", fontsize=8.2, weight="bold", color="#333333")
    ax.text(xr + 0.1, yt - 0.5,
            r"Python：セルの中心" "\n" r"（$n^3$ 個、$h=2/n$ mm）。" "\n" r"場も点モデルもセルごと。" "\n"
            r"応力だけ $16^3$ の要素" "\n" r"（セル平均）。" "\n\n"
            r"ANSYS：要素の Gauss 点" "\n" r"（$2\times2\times2$ 個/要素）。" "\n" r"場は USSFin、点モデル" "\n"
            r"と $\alpha$ は usermat、" "\n" r"応力は同じ要素で。", ha="left", va="top", fontsize=6.9, linespacing=1.45, color="#333333")
    y = yt - TH - 0.35
    centers = []
    for key, title, left, right, h in STEPS:
        y0 = y - h
        for x, w in ((X0, W_L), (xr, W_R)):
            ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.01,rounding_size=0.05",
                                        fc=F[key], ec=C[key], lw=1.0))
        ax.text(X0 + 0.12, y - 0.12, title, ha="left", va="top", fontsize=8.6, weight="bold", color=C[key])
        ax.text(X0 + 0.12, y - 0.5, left, ha="left", va="top", fontsize=7.0, linespacing=1.5)
        ax.text(xr + 0.1, y - 0.12, right, ha="left", va="top", fontsize=6.6, linespacing=1.45, color="#333333")
        centers.append((y0, y))
        y = y0 - GY
    # what each arrow carries
    for k, ((b0, _), (_, t1)) in enumerate(zip(centers[:-1], centers[1:])):
        arrow(ax, (X0 + 0.6, b0 - 0.01), (X0 + 0.6, t1 + 0.01), "#444444", lw=1.2)
        ax.text(X0 + 0.78, (b0 + t1) / 2, ARROWS[k], ha="left", va="center", fontsize=7.0, color="#222222")
    # step 2 -> step 4: c and a at both ends of the step (the point model interpolates them)
    (n0, n1), (p0, p1) = centers[2], centers[4]
    arrow(ax, (X0 - 0.02, (n0 + n1) / 2), (X0 - 0.02, (p0 + p1) / 2), C["nut"], rad=0.35, lw=1.1)
    ax.text(0.28, ((n0 + n1) / 2 + (p0 + p1) / 2) / 2, r"$c^n,c^{n+1},a^n,a^{n+1}$（補間用）", rotation=90,
            ha="center", va="center", fontsize=6.8, color=C["nut"], bbox=dict(fc="white", ec="none", pad=0.3))
    # loop back to the next step (between the columns)
    xl = X0 + W_L + GAP / 2
    arrow(ax, (xl, centers[-2][0] + 0.1), (xl, centers[0][0] + 0.1), C["point"], ls="--", lw=1.0)
    ax.text(xl, (centers[-2][0] + centers[0][0]) / 2, r" $n\to n+1$：$(\Sigma_i\bar\phi_i)^{n+1}$ は次の手順 2 へ ",
            rotation=90, ha="center", va="center", fontsize=7.0, color=C["point"], bbox=dict(fc="white", ec="none", pad=0.5))
    # glossary
    yg = centers[-1][0] - 0.35
    ax.add_patch(FancyBboxPatch((X0, yg - GH + 0.2), W_L + GAP + W_R, GH - 0.2, boxstyle="round,pad=0.01,rounding_size=0.05",
                                fc="white", ec="#999999", lw=0.8))
    ax.text(X0 + 0.12, yg - 0.12, "用語", ha="left", va="top", fontsize=8.6, weight="bold", color="#333333")
    yy = yg - 0.5
    for (term, mean), nl in zip(GLOSS, NL):
        ax.text(X0 + 0.2, yy, term, ha="left", va="top", fontsize=7.0, weight="bold", color="#333333")
        ax.text(X0 + 3.9, yy, mean, ha="left", va="top", fontsize=6.8, color="#333333", linespacing=1.5)
        yy -= 0.42 * nl + 0.12
    out = ROOT / "assets" / "fig_coupling_impl_flow_ja"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

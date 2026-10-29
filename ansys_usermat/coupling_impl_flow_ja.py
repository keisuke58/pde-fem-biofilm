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
     "\n" r"生きている量 $L^n=\Sigma_i\,\phi_i\psi_i$（$\mathbf{g}^n$ から。始まっていないセルは $0$）",
     r"同じ量を Gauss 点の状態変数に持つ" "\n" r"（$\alpha-1$、$\mathbf{g}$ は usermat の ustatev）", 1.5),
    ("field", r"1. $\phi_{3D}$ と $\alpha$ を進める（RK4、内側の刻み $\Delta t/m$）",
     r"$\dot\phi_{3D}=\beta\nabla_h^2\phi_{3D}+k_\alpha\alpha-r\frac{c}{k+c}\,\mathbf{n}_c\cdot\nabla_h\phi_{3D}-\kappa_p\,[\phi_{3D}]_\pm$,"
     r"　$\dot\alpha=k_\alpha\phi_{3D}$" "\n"
     r"$\nabla_h^2$：7 点（壁は鏡映＝流れゼロ）。前線の項：MUSCL＋minmod の風上（TVD）" "\n"
     r"$\mathbf{n}_c=\nabla c/\sqrt{|\nabla c|^2+\varepsilon^2}$、$\varepsilon^2=10^{-6}$。$c$ は内側の各刻みで解き直す（前線が栄養に追従）",
     r"USSFin（オリバーの要素）：陽的に更新" "\n" r"$\nabla^2$、$\nabla$ は NEM" "\n"
     r"（近傍 30 点の重み付き最小二乗）" "\n" r"前線の項は風上側の MLS 勾配" "\n" r"$\alpha$：usermat で Gauss 点ごとに積分", 2.3),
    ("nut", r"2. 栄養 $c^{n+1}$ と抗菌薬 $a^{n+1}$（準定常）",
     r"$-\nabla_h^2c+\Lambda^2\,\phi_{3D}^{n+1}L^n\,c=0$、供給面で $c=1$（$\Lambda^2=g\ell^2/d=1$）" "\n"
     r"$a$ も同じ式（外の面で $a=1$、$T^*\geq1/3$）。疎行列を AMG 前処理つき CG で解く" "\n"
     r"$L^n$ は前のステップの値 → 消費は 1 ステップ遅れる（反復しても差は刻みの誤差以下）",
     r"同じ準定常の式を要素の" "\n" r"2 つの栄養の場で解く" "\n" r"（LIVING_CONS で消費に $L$ を掛ける）", 1.9),
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
     r"$L^{n+1}=\Sigma_i\,\phi_i^{n+1}\psi_i^{n+1}$（次のステップの手順 2 で使う）",
     r"同じ（状態変数から）", 0.95),
    ("mech", r"6. 応力（書き出す時刻だけ。成長には戻さない）",
     r"$\phi_{3D}$ と $\alpha-1$ を $16^3$ の六面体要素ごとにセル平均 → $\mathbf{F}_g=\alpha\mathbf{I}$、ネオフック、"
     r"$E=E_0(\phi_{3D}^2+f)$、$\nu=0.49$" "\n" r"選択的低減積分で静的につり合いを解く（$32^3$、$64^3$ の要素でも種の応力は $0.2\,\%$ 以内）",
     r"usermat：$\mathbf{F}_e=\mathbf{F}\mathbf{F}_g^{-1}$、" "\n" r"接線は $\mathbf{F}$ の摂動（$10^{-7}$）", 1.6),
]

# glossary under the flow (plain Japanese)
GLOSS = [
    ("セル", "立方体を細かく刻んだ小さな箱。値はここに置く"),
    ("ラプラシアン（7 点）", "まわりと比べてどれだけ高いか低いか。前後・左右・上下の 6 つと自分で計算。拡散に使う"),
    ("鏡映（壁）", "壁の外に鏡に映した自分がいるとみなす。壁から何も出入りしない条件になる"),
    ("風上（upwind）", "流れてくる側の値を使う計算。流れのある問題で安定する"),
    ("TVD、MUSCL、minmod", "風上の計算を 2 次の精度に上げつつ、縁でガタガタしないようブレーキをかける仕組み"),
    ("正則化 ε", "∇c ÷ |∇c| で 0 で割らないよう、分母に小さい数を足すこと"),
    ("RK4", "1 歩進むときに途中の傾きを 4 回測って平均する、精度の高い進め方"),
    ("準定常", "栄養は一瞬で行き渡るので、毎回「落ち着いた後の形」だけを計算する"),
    ("1 ステップ遅れ", "栄養の消費に、1 つ前の時刻の生きている量を使うこと"),
    ("疎行列、CG、AMG", "ほとんど 0 の大きな連立方程式を、少しずつ答えに近づけて解く（AMG は粗い格子で速める工夫）"),
    ("陰的オイラー法", "進んだ先の値で式を立てる進め方。手間はかかるが大きく進めても壊れない"),
    ("Newton 法、接線", "方程式の答えを接線をたどって探す方法。接線は「入力を少し変えたときの出力の変化」"),
    ("ctypes", "Python から Fortran のプログラムを呼ぶ仕組み"),
    ("Gauss 点、状態変数", "要素の中の計算点（1 要素に 8 個）と、そこで次の時刻へ持ち越す数値"),
    ("NEM、MLS", "近くの点の値に曲面を当てはめて、傾きや曲がりを出す方法（オリバーの要素）"),
    ("ネオフック", "ゴムのような材料のいちばん簡単なモデル"),
    ("選択的低減積分", "ほとんど縮まない材料で、要素が固くなりすぎないようにする工夫"),
    ("F の摂動", "変形を少しだけずらして答えの変化を測り、接線を数値で作ること"),
]


def main():
    W_L, W_R, X0, GAP = 12.4, 4.9, 0.1, 0.45
    GH = 0.42 * len(GLOSS) + 0.7
    tot = sum(s[4] for s in STEPS) + 0.32 * (len(STEPS) - 1) + 0.9 + GH
    fig, ax = plt.subplots(figsize=(8.3, 8.3 * tot / (X0 + W_L + GAP + W_R + 0.1)))
    ax.set_xlim(0, X0 + W_L + GAP + W_R + 0.1); ax.set_ylim(0, tot); ax.axis("off")
    ax.text(X0 + W_L / 2, tot - 0.25, r"有限体積のコード（Python、外側の刻み $\Delta t=0.0025\,T^*$）", ha="center",
            va="center", fontsize=9.5, weight="bold")
    ax.text(X0 + W_L + GAP + W_R / 2, tot - 0.25, "ANSYS の要素で対応する所", ha="center", va="center", fontsize=9.5,
            weight="bold")
    y = tot - 0.6
    centers = []
    for key, title, left, right, h in STEPS:
        y0 = y - h
        for x, w, txt, fs in ((X0, W_L, left, 7.6), (X0 + W_L + GAP, W_R, right, 7.2)):
            ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.01,rounding_size=0.05",
                                        fc=F[key], ec=C[key], lw=1.0))
        ax.text(X0 + 0.12, y - 0.12, title, ha="left", va="top", fontsize=8.6, weight="bold", color=C[key])
        ax.text(X0 + 0.12, y - 0.5, left, ha="left", va="top", fontsize=7.2, linespacing=1.5)
        ax.text(X0 + W_L + GAP + 0.1, y - 0.12, right, ha="left", va="top", fontsize=6.9, linespacing=1.45,
                color="#333333")
        centers.append((y0, y))
        y = y0 - 0.32
    for (b0, _), (_, t1) in zip(centers[:-1], centers[1:]):
        ax.add_patch(FancyArrowPatch((X0 + 0.6, b0 - 0.01), (X0 + 0.6, t1 + 0.01), arrowstyle="-|>",
                                     mutation_scale=10, lw=1.2, color="#444444"))
    # loop back to the next step (between the columns)
    xl = X0 + W_L + GAP / 2
    ax.add_patch(FancyArrowPatch((xl, centers[-2][0] + 0.1), (xl, centers[0][0] + 0.1), arrowstyle="-|>",
                                 mutation_scale=10, lw=1.0, color=C["point"], linestyle="--"))
    ax.text(xl, (centers[-2][0] + centers[0][0]) / 2, r" $n\to n+1$ ", rotation=90, ha="center", va="center",
            fontsize=7.5, color=C["point"], bbox=dict(fc="white", ec="none", pad=0.5))
    # glossary
    yg = centers[-1][0] - 0.35
    ax.add_patch(FancyBboxPatch((X0, yg - GH + 0.2), W_L + GAP + W_R, GH - 0.2, boxstyle="round,pad=0.01,rounding_size=0.05",
                                fc="white", ec="#999999", lw=0.8))
    ax.text(X0 + 0.12, yg - 0.12, "用語", ha="left", va="top", fontsize=8.6, weight="bold", color="#333333")
    for k, (term, mean) in enumerate(GLOSS):
        yy = yg - 0.5 - 0.42 * k
        ax.text(X0 + 0.2, yy, term, ha="left", va="top", fontsize=7.0, weight="bold", color="#333333")
        ax.text(X0 + 3.9, yy, mean, ha="left", va="top", fontsize=7.0, color="#333333")
    out = ROOT / "assets" / "fig_coupling_impl_flow_ja"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

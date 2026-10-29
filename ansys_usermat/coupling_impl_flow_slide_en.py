#!/usr/bin/env python3
"""coupling_impl_flow_slide_en.py -- English version of coupling_impl_flow_slide_ja.py: one outer step of the coupled model on one 16:9 slide: the time line (what is
evaluated when) and the six steps, each with its title and one equation; arrows say what they carry. The detailed
version with the ANSYS column is coupling_impl_flow_ja.py; its glossary goes to a page of its own here.

    python ansys_usermat/coupling_impl_flow_slide_en.py
        -> assets/fig_coupling_impl_flow_slide_en.{pdf,png}, assets/fig_coupling_glossary_en.{pdf,png}
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

import figstyle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
figstyle.apply()
plt.rcParams.update({"pdf.fonttype": 42})

C = dict(field="#1f4e8c", nut="#2f7d6d", start="#8a5a9e", point="#c2571a", state="#3d4a57", mech="#555555")
F = dict(field="#e8eef7", nut="#e6f2ef", start="#f2ebf5", point="#fbeee5", state="#eef0f2", mech="#f0f0f0")

STEPS = [
    ("field", r"1. Advance $\phi_{3D}$ and $\alpha$",
     r"$\dot\phi_{3D}=\beta\nabla^2\phi_{3D}+k_\alpha\alpha-r\,\frac{c}{k+c}\,\mathbf{n}_c\cdot\nabla\phi_{3D}$, "
     r"$\dot\alpha=k_\alpha\phi_{3D}$"),
    ("nut", r"2. Nutrient $c$, antibiotic $a$",
     r"$\dot c=d\nabla^2c-g\,\phi_{3D}\,(\Sigma_i\bar\phi_i)\,c$   (same equation for $a$)"),
    ("start", r"3. Start newly reached points",
     r"points with $\phi_{3D}^{n+1}\geq\phi_{\min}$: neighbours' composition, $\psi_i=0.999$"),
    ("point", r"4. Point model (each point)",
     r"Klempt et al. 2026 as published, $c^*=c^*_0\,c$, $\alpha^*=\alpha^*_0\,a$"),
    ("state", r"5. Living fraction",
     r"$\Sigma_i\bar\phi_i=\Sigma_i\,\phi_i\psi_i$"),
    ("mech", r"6. Stress (output times only)",
     r"$\nabla\cdot\boldsymbol{\sigma}=\mathbf{0}$, $\boldsymbol{\sigma}=E_0(\phi_{3D}^2+f)\,\mathbb{C}:(\boldsymbol{\varepsilon}-(\alpha-1)\mathbf{I})$"),
]
GLOSS = [
    ("Laplacian (7 points)", r"$\nabla_h^2u_i=(\sum_{j\in N(i)}u_j-6u_i)/h^2$ ($N(i)$: the 6 neighbours), second order; discretises $\beta\nabla^2\phi_{3D}$"),
    ("Mirror (walls)", r"a ghost cell outside the wall with $u_{\rm ghost}=u_{\rm in}$: zero flux $\partial u/\partial n=0$ at second order"),
    ("Upwind", r"the face value of the advection $\mathbf{v}\cdot\nabla u$ is taken from the upstream side; first order is monotone but adds diffusion $\approx|\mathbf{v}|h/2$"),
    ("MUSCL, minmod, TVD", r"MUSCL: a linear reconstruction in each cell gives second-order face values. minmod$(a,b)$: the smaller one if the signs agree,"
     "\n" r"else 0 (zero slope at extrema). The total variation $\sum_i|u_{i+1}-u_i|$ does not grow (TVD: no new extrema)"),
    ("Regularised front direction", r"$\mathbf{n}_c=\nabla c/\sqrt{|\nabla c|^2+\varepsilon^2}$, $\varepsilon^2=10^{-6}$: a unit vector for $|\nabla c|\gg\varepsilon$, continuously 0 for $|\nabla c|\to0$"),
    ("RK4", r"explicit four-stage Runge-Kutta, local error $O(\Delta t^5)$; stable below limits on $\beta\Delta t/h^2$ and the Courant number $|\mathbf{v}|\Delta t/h$"),
    ("Backward Euler (Eq. 35)", r"$\dot c$ replaced by $(c^{n+1}-c^n)/\Delta t$, right-hand side at the new time: stable for any $\Delta t$, one linear solve per step."
     "\n" r"Diffusion time $\ell^2/d\approx10^{-8}T^*\ll\Delta t$, so the time term is $\sim10^{-6}$; the limit without it differs by $\leq0.1\,\%$ ($32^3$)"),
    ("One-step lag", r"step 2 uses $(\Sigma_i\bar\phi_i)^n$ instead of $(\Sigma_i\bar\phi_i)^{n+1}$ (explicit coupling), error $O(\Delta t)$; two iterations change less than the step error"),
    ("Sparse matrix, CG, AMG", r"$-\nabla_h^2+q$ is symmetric positive definite with 7 non-zeros per row. CG: Krylov iteration for such systems."
     "\n" r"AMG: a preconditioner built from coarse levels of the matrix; the iteration count hardly grows with the grid"),
    ("Implicit Euler", r"$(\mathbf{g}^{k+1}-\mathbf{g}^k)/\delta t=\mathbf{f}(\mathbf{g}^{k+1})$: first order, A-stable"),
    ("Newton, Jacobian", r"$\mathbf{Q}(\mathbf{g})=\mathbf{0}$ solved by $\mathbf{g}\leftarrow\mathbf{g}-\mathbf{J}^{-1}\mathbf{Q}$, $\mathbf{J}=\partial\mathbf{Q}/\partial\mathbf{g}$; quadratic convergence near the solution"),
    ("ctypes", r"Python's foreign-function interface: calls the compiled Fortran library on double arrays (the same routine as in ANSYS)"),
    ("Gauss points, state variables", r"the $2\times2\times2$ Gauss-Legendre points of an element; state variables are kept per point between increments (ustatev)"),
    ("NEM, MLS", r"a second-order Taylor expansion (9 unknowns) fitted to the values of 30 neighbouring Gauss points by least squares with Gaussian"
     "\n" r"weights (width $\beta^*=0.2$ mm); its coefficients give $\nabla u$ and $\nabla^2u$ (moving least squares)"),
    ("Neo-Hookean", r"$W=C_{10}(\bar I_1-3)+(J-1)^2/D_1$: isochoric and volumetric parts, $C_{10}=\mu/2$ (Klempt et al. 2024, Eq. 20)"),
    (r"Selective reduced integration ($\bar{\mathbf{B}}$)", r"the volumetric strain is averaged over the element, avoiding volumetric locking for $\nu\to0.5$"),
    (r"Perturbation of $\mathbf{F}$", r"$F_{ij}$ shifted by $10^{-7}$ (relative); forward differences of the stress give the tangent $\partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon}$ (Sun et al. 2008)"),
]

ARROWS = [r"$\phi_{3D}^{n+1}$, $\alpha^{n+1}$",
          r"$c^{n+1}$, $a^{n+1}$",
          r"$\mathbf{g}^n$ of the started points",
          r"$\mathbf{g}^{n+1}$",
          r"$\phi_{3D}^{n+1}$, $\alpha^{n+1}$"]


def arrow(ax, p, q, col, ls="-", rad=0.0, lw=1.6, ms=13):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=col, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}"))


def timeline(ax, x, y, w, h):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.06", fc="white",
                                ec="#888888", lw=1.0))
    ax.text(x + 0.15, y + h - 0.12, r"One outer step $t_n\to t_{n+1}$ ($\Delta t=0.0025\,T^*$): what is used when",
            ha="left", va="top", fontsize=13, weight="bold", color="#333333")
    xa, xb = x + 3.6, x + w - 2.6
    yax = y + h - 0.62
    ax.plot([xa - 0.15, xb + 0.4], [yax, yax], color="#333333", lw=1.0)
    for xx, lab in ((xa, r"$t_n$"), (xb, r"$t_{n+1}$")):
        ax.plot([xx, xx], [y + 0.1, yax], color="#aaaaaa", lw=0.8, ls=":")
        ax.text(xx, yax + 0.04, lab, ha="center", va="bottom", fontsize=12)
    keys = ["field", "nut", "start", "point", "state"]
    labs = ["1. $\\phi_{3D}$, $\\alpha$", "2. $c$, $a$", "3. Start", "4. Point model", "5. $\\Sigma_i\\bar\\phi_i$"]
    notes = ["$m$ inner steps (RK4)", "one backward Euler step", "decided by $\\phi_{3D}^{n+1}$",
             "$n_{\\rm sub}$ steps, $c^*$ interpolated", "to 2. of next step (lag)"]
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
    ax.text(x0 - 0.42, (n0 + p1) / 2 - 0.12, "$c^n,c^{n+1}$\n$a^n,a^{n+1}$\n(interpolated)", ha="right", va="center",
            fontsize=10.5, color=C["nut"], linespacing=1.3)
    # 5 -> 2 of the next step
    (s0, s1), (q0, q1) = tops[4], tops[1]
    arrow(ax, (x0 + w + 0.02, (s0 + s1) / 2), (x0 + w + 0.02, (q0 + q1) / 2), C["point"], ls="--", rad=0.45, lw=1.5)
    ax.text(x0 + w + 0.6, (s1 + q0) / 2, "$(\\Sigma_i\\bar\\phi_i)^{n+1}$\nto 2. of the\nnext step", ha="left",
            va="center", fontsize=10.5, color=C["point"], linespacing=1.3)
    ax.text(x0, 0.12, "ANSYS: same order at the Gauss points (fields: USSFin; point model and $\\alpha$: usermat). Details: separate figure", ha="left", va="bottom", fontsize=10.5, color="#555555")
    out = ROOT / "assets" / "fig_coupling_impl_flow_slide_en"
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
    ax.text(0.2, h - 0.15, "Terms (implementation flow)", ha="left", va="top", fontsize=13, weight="bold", color="#333333")
    yy = h - 0.6
    for (term, mean), n in zip(GLOSS, nl):
        ax.text(0.3, yy, term, ha="left", va="top", fontsize=10.5, weight="bold", color="#333333")
        ax.text(3.4, yy, mean, ha="left", va="top", fontsize=10, color="#333333", linespacing=1.45)
        yy -= 0.36 * n + 0.12
    out = ROOT / "assets" / "fig_coupling_glossary_en"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

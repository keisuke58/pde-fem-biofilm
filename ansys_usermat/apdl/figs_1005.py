"""figs_1005.py -- figures for the 5 Oct meeting, from the ANSYS runs of
1-2 Oct with Klempt et al. 2024 Table 2 values (k_alpha = 1e-3 /T*, in both the
partner's field and the growth law), ONE species (second species field
started at 0 with rate 0), stored in results/2026-10-01_paper_values/
(git-ignored like every *.csv: local on IKMHIWI03; the same files are in
F:\biofilm_upf_wired as phi_trace_ds_pv_exact_1sp.csv and <name>_ds_kl_*.csv).
Stiffness also as in Table 2 / Eq. 20 (2 Oct): E = 10 Pa, nu = 0.49, weighted
with phi^2 (void floor 1e-3): --set YOUNG_BIO=1e-5 --set POISSON_BIO=0.49
--set YOUNG_VOID=-1e-3 (deck units MPa). Stresses are plotted in Pa.
Decks: make_wired_deck.py on the partner's stage-1 deck with
  --set K_LOCAL1=1e-3 --set K_LOCAL2=0 --set MY_BIOSTART2=0.0
  and --props 7=1e-3,28=1 (Eq. 36 in the material routine) or 28=3 (the
  partner's growth variable); --post both --post-elem 220.

    python ansys_usermat/apdl/figs_1005.py      -> assets/fig1005_*.png

1. exact solution: in a gradient-free region the partner's field (Klempt 2024
   Eq. 34/36) reduces to phi' = k alpha_K, alpha_K' = k phi, i.e.
   phi = sinh(k t), alpha_K = cosh(k t); the explicit update has its own
   discrete exact solution, which ANSYS follows closely (alpha_K to
   round-off; phi to ~1e-4 relative: element 1 is not perfectly
   gradient-free);
2. whole model, element by element: von Mises stress with Eq. 36 evaluated in
   the material routine vs with the partner's own growth variable;
3. the seeded element 220 over time: von Mises, mean stress, the largest
   neighbour von Mises, and alpha.
Result set 1 is dropped (the partner's routine returns zero stress on its
first calls while it sets up its operators).
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
RES = HERE / "results" / "2026-10-01_paper_values"
OUT = HERE.parents[1] / "assets"
K, DT = 1.0e-3, 0.1
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402
figstyle.apply(12)


def read_rows(p):
    with open(p, newline="") as f:
        return [[x.strip() for x in r] for r in csv.reader(f) if r]


def fig_exact():
    last = {}
    for r in read_rows(RES / "phi_trace_ds_pv_exact.csv"):
        if int(r[0]) == 1:
            last[int(r[3])] = (float(r[5]), float(r[6]))       # bio1, locbio1 at substep s
    s = np.array(sorted(last))
    t = (s - 1) * DT                                             # usermat sees the state after s-1 updates
    bio = np.array([last[i][0] for i in s])
    loc = np.array([last[i][1] for i in s])
    x, y, dx, dy = 0.0, 1.0, [], []
    for _ in s:
        dx.append(x); dy.append(y)
        x, y = x + DT * K * y, y + DT * K * x
    dx, dy = np.array(dx), np.array(dy)
    tt = np.linspace(0, t[-1], 200)

    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
    ax[0].plot(tt, np.sinh(K * tt), "k-", label=r"exact $\sinh(k_\alpha t)$")
    ax[0].plot(t, dx, "k:", lw=1.5, label="exact solution of the explicit update")
    ax[0].plot(t, bio, "o", color="C0", label="ANSYS")
    ax[0].set(xlabel=r"time $T^*$", ylabel=r"$\phi$", title=r"biofilm fraction $\phi$")
    ax[1].plot(tt, np.cosh(K * tt) - 1, "k-", label=r"exact $\cosh(k_\alpha t)-1$")
    ax[1].plot(t, dy - 1, "k:", lw=1.5, label="exact solution of the explicit update")
    ax[1].plot(t, loc - 1, "o", color="C1", label="ANSYS")
    ax[1].set(xlabel=r"time $T^*$", ylabel=r"$\alpha - 1$", title=r"growth $\alpha - 1$")
    for a in ax[:2]:
        a.legend(loc="upper left", fontsize=9)
        a.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    eb, el = np.abs(bio - dx), np.abs(loc - dy)
    ax[2].semilogy(t[1:], np.maximum(eb[1:], 1e-18), "o-", label=r"$\phi$")
    ax[2].semilogy(t[1:], np.maximum(el[1:], 1e-18), "s-", label=r"$\alpha - 1$")
    ax[2].set(xlabel=r"time $T^*$", ylabel="|ANSYS - discrete exact|",
              title="error vs. the update's exact solution")
    ax[2].legend(fontsize=10)
    fig.suptitle(rf"An element outside the seed, Klempt 2024 Eq. 34/36, $k_\alpha$ = {K:g} (Table 2), "
                 rf"$\Delta t$ = {DT:g}", y=1.03)
    fig.savefig(OUT / "fig1005_exact.png")
    plt.close(fig)
    return float(eb.max()), float(el.max())


def read_all(p):
    rows = read_rows(p)[1:]
    a = np.array([[float(x) for x in r] for r in rows])
    return a[:, 0].astype(int), a[:, 1] * 1e6, a[:, 2:5].mean(axis=1) * 1e6, a[:, 5]   # MPa -> Pa


def fig_whole():
    e1, q1, p1, a1 = read_all(RES / "all_stress_ds_pv_eq36.csv")
    e2, q2, p2, a2 = read_all(RES / "all_stress_ds_pv_partner.csv")
    assert (e1 == e2).all()
    cos = float(q1 @ q2 / (np.linalg.norm(q1) * np.linalg.norm(q2)))
    ratio = float(np.median(q2[q1 > 1e-3 * q1.max()] / q1[q1 > 1e-3 * q1.max()]))
    fig, ax0 = plt.subplots(figsize=(6.2, 4.6))
    ax = [ax0]
    m = (q1 > 0) & (q2 > 0)
    ax[0].loglog(q1[m], q2[m], ".", ms=3, alpha=0.5)
    lo, hi = q1[m].min(), q1[m].max()
    ax[0].plot([lo, hi], [lo, hi], "k--", lw=1, label="1 : 1")
    ax[0].plot([lo, hi], [ratio * lo, ratio * hi], "r-", lw=1, label=f"{ratio:.2f} : 1 (median)")
    ax[0].set(xlabel="von Mises [Pa], Eq. 36 in the material routine",
              ylabel="von Mises [Pa], partner's growth variable",
              title=f"{len(q1)} elements, cosine similarity {cos:.4f}")
    ax[0].legend(fontsize=10)
    fig.suptitle(r"Whole model at $T^*$ = 1.1, one species; Klempt 2024 Table 2: $k_\alpha$ = 1e-3, $E$ = 10 Pa, $\nu$ = 0.49, $E \propto \phi^2$", y=1.0, fontsize=10)
    fig.savefig(OUT / "fig1005_whole_model.png")
    plt.close(fig)
    return cos, ratio, float(q1.max()), float(q2.max())


def read_elem(p):
    a = np.array([[float(x) for x in r] for r in read_rows(p)[1:]])
    a = a[1:]                                                    # drop result set 1
    return a[:, 2], a[:, 9] * 1e6, a[:, 3:6].mean(axis=1) * 1e6, a[:, 11] * 1e6, a[:, 10]   # MPa -> Pa


def fig_elem():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    out = {}
    for name, f, ls in (("Eq. 36 in the material routine", "elem_stress_ds_pv_eq36.csv", "-"),
                        ("partner's growth variable", "elem_stress_ds_pv_partner.csv", "--")):
        t, q, p, nb, al = read_elem(RES / f)
        ax[0].plot(t, q, "C0" + ls, marker="o", ms=4, label=f"element 220, {name}")
        ax[0].plot(t, nb, "C3" + ls, marker="s", ms=4, label=f"largest neighbour, {name}")
        ax[1].plot(t, al, "C2" + ls, marker="o", ms=4, label=name)
        out[name] = (float(q[-1]), float(p[-1]), float(nb[-1]), float(al[-1]))
    ax[0].set(xlabel=r"time $T^*$", ylabel="von Mises [Pa]", title="seeded element 220 and its neighbours")
    ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax[0].legend(fontsize=8)
    ax[1].set(xlabel=r"time $T^*$", ylabel=r"$\alpha - 1$", title=r"growth $\alpha - 1$ of element 220")
    ax[1].legend(fontsize=9)
    ax[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.suptitle(r"One species; Klempt 2024 Table 2: $k_\alpha$ = 1e-3, $E$ = 10 Pa, $\nu$ = 0.49, $E \propto \phi^2$; result set 1 dropped", y=1.02, fontsize=11)
    fig.savefig(OUT / "fig1005_element220.png")
    plt.close(fig)
    return out


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    print("exact: max |ANSYS - discrete| phi %.2e, alpha_K %.2e" % fig_exact())
    print("whole: cosine %.4f, median ratio partner/ours %.3f, max SEQV ours %.4g, partner %.4g" % fig_whole())
    for k, v in fig_elem().items():
        print(f"elem 220, {k}: SEQV {v[0]:.4g}, p {v[1]:.4g}, nbr max {v[2]:.4g}, alpha {v[3]:.4g}")

#!/usr/bin/env python3
"""Klempt et al. 2024: Eq. 36 against Table 1 step 3a, the two growth laws
for alpha printed in the paper.

  Eq. 36:          d alpha / dt = k_alpha phi
                   -> growth wherever there is biofilm, in proportion to time;
  Table 1 step 3a: (alpha_{n+1} - alpha_n) / ((1 + alpha_{n+1}) dt)
                     - (k_alpha / eta_alpha) (phi_{n+1} - phi_n) / dt = 0
                   -> growth only where phi changes; integrated,
                   (1 + a) = (1 + a_0) exp(k (phi - phi_0)).

Table 1 says it solves Eq. 36, so one of the two is a misprint. Here k_alpha /
eta_alpha is read as k_alpha (parameters divided by eta lose their bar, so
k_alpha already contains the division), and alpha in Table 1 as the growth
alpha - 1 of F_g = alpha I. Field: test case 4.1 as fitted in
klempt2024_timescale.py (Table 2, edge-line source, time scale 1.5), 50 steps
of 0.02 up to T* = 1. The ANSYS seed is the limit phi = 1 at all times.

    python JAXFEM/klempt2024_alpha_law.py -> assets/fig_klempt2024_alpha_law.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402

OUT = ROOT / "assets" / "fig_klempt2024_alpha_law.png"
CACHE = ROOT / "figdata" / "klempt2024_alpha_law_snapshots.npz"   # the field, so a redraw needs no rerun
K_ALPHA, S_MACRO, N_STEP = 1e-3, 1.5, 50


def snapshots():
    if CACHE.exists() and "--recompute" not in sys.argv:
        return np.load(CACHE)["snaps"]
    K.R, K.BETA, K.K_A = 100 * S_MACRO, 2 * S_MACRO, 1e-3 * S_MACRO
    phi0, _, _ = K.setup("fig4_edge")
    snaps = [phi0.copy()]
    for _ in range(N_STEP):
        _, phi, _ = B.run_setup(snaps[-1], B.nutrient_mask("edge"), 1e8, "ic",
                                "first_order", "abs", t_end=1.0 / N_STEP)
        snaps.append(phi)
    snaps = np.array(snaps)
    CACHE.parent.mkdir(exist_ok=True)
    np.savez_compressed(CACHE, snaps=snaps.astype(np.float32))
    return snaps


def laws(snaps):
    dt = 1.0 / N_STEP
    a36 = K_ALPHA * np.cumsum(0.5 * (snaps[1:] + snaps[:-1]) * dt, axis=0)    # trapezoid
    at1 = np.expm1(K_ALPHA * (snaps[1:] - snaps[0]))                             # exact
    return np.concatenate([np.zeros_like(snaps[:1]), a36]), np.concatenate([np.zeros_like(snaps[:1]), at1])


def main():
    figstyle.apply(size=12)
    S = snapshots()
    a36, at1 = laws(S)
    t = np.linspace(0, 1, N_STEP + 1)
    w = S[-1] >= 0.5                                   # biofilm at T* = 1
    seed = S[0] >= 0.5                                 # biofilm from the start
    print(f"nodes with phi >= 0.5 at T* = 1: {w.sum()}, of them from the start: {(w & seed).sum()}")
    for name, a in (("Eq. 36", a36), ("Table 1", at1)):
        print(f"{name:8s} mean alpha-1 at T* = 1 over phi>=0.5: {a[-1][w].mean():.3e};"
              f" initial biofilm: {a[-1][w & seed].mean():.3e}; new biofilm: {a[-1][w & ~seed].mean():.3e}")
    print(f"ANSYS seed (phi = 1 const.) at T* = 1: Eq. 36 {K_ALPHA:.1e}, Table 1 0")

    kz = K.N // 2
    with plt.rc_context({"font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 9, "xtick.labelsize": 7.5,
                         "ytick.labelsize": 7.5, "legend.fontsize": 7.5}):   # drawn at the thesis text width
        fig = plt.figure(figsize=(6.3, 2.9))
        ax = [fig.add_axes([0.075, 0.33, 0.19, 0.53]), fig.add_axes([0.30, 0.33, 0.19, 0.53]),
              fig.add_axes([0.665, 0.33, 0.315, 0.53])]
        cax = fig.add_axes([0.505, 0.33, 0.012, 0.53])
        vmax = max(a36[-1].max(), at1[-1].max())
        for a, (name, f) in zip(ax[:2], (("(a) Eq. 36", a36), ("(b) Table 1, step 3a", at1))):
            im = a.imshow(f[-1][:, :, kz].T, origin="lower", extent=(0, K.L, 0, K.L), cmap="viridis",
                          vmin=0, vmax=vmax)
            a.contour(np.linspace(0, K.L, K.N), np.linspace(0, K.L, K.N), S[-1][:, :, kz].T, [0.5],
                      colors="white", linewidths=0.8)
            a.set_title(name)
            a.set_xlabel(r"$x$ [$\mu$m]"); a.grid(False)
            a.set_xticks([0, 10, 20]); a.set_yticks([0, 10, 20])
        ax[0].set_ylabel(r"$y$ [$\mu$m]"); ax[1].set_yticklabels([])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(r"$\alpha-1$ at $T^*=1$")
        cb.formatter.set_powerlimits((0, 0)); cb.update_ticks()
        b = ax[2]
        for name, f, ls in (("Eq. 36", a36, "-"), ("Table 1", at1, "--")):
            b.plot(t, [f[k][w & seed].mean() for k in range(len(t))], ls, color="C0", lw=1.1,
                   label=f"{name}, initial biofilm")
            b.plot(t, [f[k][w & ~seed].mean() for k in range(len(t))], ls, color="C1", lw=1.1,
                   label=f"{name}, new biofilm")
        b.plot(t, K_ALPHA * t, ":", color="0.3", lw=1.1, label=r"ANSYS seed, Eq. 36")
        b.set_xlabel(r"$T^*$"); b.set_ylabel(r"mean $\alpha-1$")
        b.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        b.set_title("(c) mean growth")
        b.legend(frameon=False, handlelength=1.8, ncol=2, loc="upper center", bbox_to_anchor=(0.42, -0.27),
                 columnspacing=0.8)
        fig.savefig(OUT, dpi=250); fig.savefig(OUT.with_suffix(".pdf"))
    print("wrote", OUT)


if __name__ == "__main__":
    main()

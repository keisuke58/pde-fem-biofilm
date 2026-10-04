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
K_ALPHA, S_MACRO, N_STEP = 1e-3, 1.5, 50


def snapshots():
    K.R, K.BETA, K.K_A = 100 * S_MACRO, 2 * S_MACRO, 1e-3 * S_MACRO
    phi0, _, _ = K.setup("fig4_edge")
    snaps = [phi0.copy()]
    for _ in range(N_STEP):
        _, phi, _ = B.run_setup(snaps[-1], B.nutrient_mask("edge"), 1e8, "ic",
                                "first_order", "abs", t_end=1.0 / N_STEP)
        snaps.append(phi)
    return np.array(snaps)


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
    fig = plt.figure(figsize=(14, 4.6))
    ax = [fig.add_axes([0.05, 0.14, 0.22, 0.68]), fig.add_axes([0.30, 0.14, 0.22, 0.68]),
          fig.add_axes([0.66, 0.14, 0.32, 0.68])]
    cax = fig.add_axes([0.535, 0.14, 0.01, 0.68])
    vmax = max(a36[-1].max(), at1[-1].max())
    for a, (name, f) in zip(ax[:2], (("Eq. 36", a36), ("Table 1, step 3a", at1))):
        im = a.imshow(f[-1][:, :, kz].T, origin="lower", extent=(0, K.L, 0, K.L), cmap="viridis",
                      vmin=0, vmax=vmax)
        a.contour(np.linspace(0, K.L, K.N), np.linspace(0, K.L, K.N), S[-1][:, :, kz].T, [0.5],
                  colors="white", linewidths=1)
        a.set_title(rf"{name}: $\alpha-1$ at $T^*=1$")
        a.set_xlabel(r"$x$ [$\mu$m]"); a.grid(False)
        a.set_xticks([0, 10, 20]); a.set_yticks([0, 10, 20])
    ax[0].set_ylabel(r"$y$ [$\mu$m]"); ax[1].set_yticklabels([])
    fig.colorbar(im, cax=cax, label=r"$\alpha-1$ [-]")
    b = ax[2]
    for name, f, ls in (("Eq. 36", a36, "-"), ("Table 1", at1, "--")):
        b.plot(t, [f[k][w & seed].mean() for k in range(len(t))], ls, color="C0",
               label=f"{name}, biofilm from the start")
        b.plot(t, [f[k][w & ~seed].mean() for k in range(len(t))], ls, color="C1",
               label=f"{name}, new biofilm")
    b.plot(t, K_ALPHA * t, ":", color="0.3", label=r"ANSYS seed, Eq. 36 ($\phi=1$)")
    b.set_xlabel(r"$T^*$"); b.set_ylabel(r"mean $\alpha-1$ [-]")
    b.legend(fontsize=9, frameon=False)
    fig.suptitle(r"Klempt et al. 2024: the two $\alpha$ laws on the test case 4.1 field "
                 r"($k_\alpha = 10^{-3}$, mid-plane $z = 10\ \mu$m, white: $\phi = 0.5$)", fontsize=12, y=0.97)
    fig.savefig(OUT, dpi=200)
    print("wrote", OUT)


if __name__ == "__main__":
    main()

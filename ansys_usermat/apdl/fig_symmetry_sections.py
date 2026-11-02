#!/usr/bin/env python3
"""Mirror symmetry of the separate example with antibiotic (64^3, T* = 1), for the question of 2 Nov 2026 (meeting
with Assoc. Prof. Muramatsu: the 3D views looked unsymmetric). Two sections of the same run:
  (a) the plane y = const through the middle of the seed: mirror symmetric about x = 0 and about z = 0;
  (b) the plane z = 0: symmetric about x = 0, not about y = 0 (the nutrient is held at y = -1 mm only).
Fields: share of species A phi_1/(phi_1 + phi_2), living volume fraction sum_i phibar_i, nutrient c, von Mises stress
(16^3 stress mesh). The largest mirror differences are printed in the title of each row.

    python ansys_usermat/apdl/fig_symmetry_sections.py   -> assets/fig_symmetry_sections.png/.pdf
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
import wp3_reference as W  # noqa: E402

NPZ = ROOT / "figdata" / ("coupled_living_separate_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1"
                          "_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast_separate_ab.npz")
T = "1.0"


def main():
    figstyle.apply()
    z = np.load(NPZ)
    seed = W.seed_mask(64, W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json"))
    ys = np.flatnonzero(seed.any(axis=(0, 2)))
    jy = int(round(ys.mean()))                      # cell layer through the middle of the seed
    phi = z[f"phi_{T}"]
    F = [("share", np.where(phi > 0.5, z[f"shareA_{T}"], np.nan), "RdBu", (0, 1), r"$\phi_1/(\phi_1+\phi_2)$"),
         ("living", np.where(phi > 0.5, z[f"L_{T}"], np.nan), "Greens", (0, 1), r"$\Sigma_i\bar\phi_i$"),
         ("c", z[f"c_{T}"], "Blues", None, r"$c$"),
         ("vm", 1e4 * np.nan_to_num(z[f"vm16_{T}"]), "magma_r", None, r"$\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")]

    def mirr(u, ax):
        u = np.nan_to_num(u)
        return np.abs(u - np.flip(u, ax)).max() / (np.abs(u).max() + 1e-30)

    dxz = max(max(mirr(u, 0), mirr(u, 2)) for _, u, *_ in F)
    dy = min(mirr(u, 1) for _, u, *_ in F[2:])
    fig, axs = plt.subplots(2, 4, figsize=(6.3, 3.6), constrained_layout=True)
    ext = (-1, 1, -1, 1)
    for j, (_, u, cm, lim, lab) in enumerate(F):
        n = u.shape[0]
        jj = int(jy * n / 64)
        A = u[:, jj, :].T                           # rows z, columns x
        B = u[:, :, n // 2 - 1: n // 2 + 1].mean(2).T if n % 2 == 0 else u[:, :, n // 2].T   # rows y, columns x
        vmin, vmax = lim if lim else (np.nanmin(u), np.nanmax(u))
        for i, (S, yl) in enumerate(((A, "$z$ [mm]"), (B, "$y$ [mm]"))):
            ax = axs[i, j]
            im = ax.imshow(S, origin="lower", extent=ext, cmap=cm, vmin=vmin, vmax=vmax, interpolation="nearest")
            ax.axvline(0, color="k", ls="--", lw=0.6)
            if i == 0:
                ax.axhline(0, color="k", ls="--", lw=0.6)
            ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
            ax.set_xlabel("$x$ [mm]", labelpad=1)
            if j == 0:
                ax.set_ylabel(yl, labelpad=1)
            else:
                ax.set_yticklabels([])
            ax.set_aspect("equal")
        cb = fig.colorbar(im, ax=axs[:, j], orientation="horizontal", location="bottom", shrink=0.9, pad=0.02)
        cb.set_label(lab)
    yv = -1 + (jy + 0.5) * 2 / 64
    axs[0, 0].text(-0.62, 0.5, f"(a) $y = {yv:.2f}$ mm\n(seed middle)", transform=axs[0, 0].transAxes, rotation=90,
                   va="center", ha="center", fontsize=7.5)
    axs[1, 0].text(-0.62, 0.5, "(b) $z = 0$", transform=axs[1, 0].transAxes, rotation=90, va="center", ha="center",
                   fontsize=7.5)
    fig.text(0.5, 1.0, f"dashed: mirror planes. Largest mirror difference about $x = 0$ and $z = 0$: {dxz:.0e} (relative)",
             ha="center", va="bottom", fontsize=7.5)
    o = ROOT / "assets" / "fig_symmetry_sections.png"
    fig.savefig(o, dpi=600, bbox_inches="tight"); fig.savefig(o.with_suffix(".pdf"), bbox_inches="tight")
    print(f"wrote {o}; largest mirror difference about x = 0 and z = 0: {dxz:.1e}; about y = 0 (smallest of c, vM): {dy:.2f}")


if __name__ == "__main__":
    main()

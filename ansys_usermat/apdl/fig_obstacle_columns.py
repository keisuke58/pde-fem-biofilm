#!/usr/bin/env python3
"""fig_obstacle_columns.py -- growth against four rigid bars (Klempt et al. 2024, Sect. 4.3, Fig. 9) with the agreed
coupling, 16^3, main constants (coupled_living_consumption.py --obstacle columns --nf 16 --snap-every 0.25 ...).
Rows: separate and mixed example without antibiotic. Columns: phi_3D at T* = 0.25 and 1, nutrient c at T* = 1, von
Mises stress at T* = 1 (16^3 stress mesh), on the mid-plane z = 0 (the bars along z cross it, black; the bars along x
lie beside the plane). Lengths on the deck's 2 mm cube (1 mm = 10 um on the paper's 20 um cube).

    python ansys_usermat/apdl/fig_obstacle_columns.py -> assets/fig_obstacle_columns.pdf/.png
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

NF = int(sys.argv[sys.argv.index("--nf") + 1]) if "--nf" in sys.argv else 16
TAG = ("_nf16" if NF == 16 else "") + "_front1_phimin0.001_v2_fortran_ansysbc_a100_gab1_gc1_tab0.333333_snaps_obstcolumns"


def mid(u):
    k = u.shape[2] // 2
    return 0.5 * (u[:, :, k - 1] + u[:, :, k])


def main():
    figstyle.apply(size=8)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    fig, axs = plt.subplots(2, 4, figsize=(6.3, 3.5), constrained_layout=True)
    for r, (key, lab) in enumerate((("separate", "separate"), ("mixed", "mixed"))):
        f = ROOT / "figdata" / f"coupled_living_{key}{TAG}.npz"
        if NF != 16:                             # the 32^3 runs were made one example at a time (--only)
            f = ROOT / "figdata" / f"coupled_living_{key}{TAG}_{key}.npz"
        d = np.load(f)
        ob = mid(d["obst_1.0"].astype(float)) > 0.5
        vmax_vm = np.nanmax(d["vm16_1.0"])
        for c, (fld, T, cm, lim, name) in enumerate((("phi", 0.25, "viridis", (0, 1), r"$\phi_{3D}$, $T^*=0.25$"),
                                                     ("phi", 1.0, "viridis", (0, 1), r"$\phi_{3D}$, $T^*=1$"),
                                                     ("c", 1.0, "Blues", (0, 1), r"$c$, $T^*=1$"),
                                                     ("vm16", 1.0, "magma", (0, vmax_vm), r"$\sigma_\mathrm{vM}$ [Pa], $T^*=1$"))):
            a = axs[r, c]
            u = mid(d[f"{fld}_{T}"])
            k_ = ob.shape[0] // u.shape[0]       # the stress lives on the 16^3 mesh
            ob_ = ob if k_ <= 1 else ob.reshape(u.shape[0], k_, u.shape[1], k_).mean((1, 3)) > 0.5
            e = np.linspace(-1, 1, u.shape[0] + 1)
            u = np.where(ob_, np.nan, u)
            im = a.pcolormesh(e, e, np.ma.masked_invalid(u).T, cmap=cm, vmin=lim[0], vmax=lim[1], shading="flat",
                              rasterized=True)
            a.pcolormesh(e, e, np.ma.masked_where(~ob_, ob_).T, cmap=ListedColormap(["k"]), shading="flat")
            a.set_aspect("equal"); a.set_xticks([-1, 0, 1]); a.set_yticks([-1, 0, 1]); a.grid(False)
            a.tick_params(labelsize=6)
            cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=5.5)
            if r == 0:
                a.set_title(f"({'abcd'[c]}) {name}", fontsize=7)
            if c == 0:
                a.set_ylabel(f"{lab}\n$y$ [mm]", fontsize=7)
            if r == 1:
                a.set_xlabel("$x$ [mm]", fontsize=7)
    out = ROOT / "assets" / ("fig_obstacle_columns" + ("" if NF == 16 else f"_nf{NF}"))
    fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"), dpi=300)
    print("wrote", out.with_suffix(".png"))


if __name__ == "__main__":
    main()

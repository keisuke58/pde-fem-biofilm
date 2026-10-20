"""fig_coupled_examples.py -- the two examples of the agreed coupling (species A in the seed, B in four corner blocks;
A and B mixed in the seed), with and without the antibiotic from T* = 0.5, on the mid-plane z = 0 at T* = 1.
Python (coupled_living_consumption.py --v2 --phi-min 0.1 --front 0.1 --fortran --ansys-bc), i.e. the boundary
conditions of the ANSYS decks and the Fortran point model of the ANSYS runs; synthetic data (Fritsch et al. 2025,
case I). Columns: share of A where the point model runs, living volume fraction sum phibar_i, nutrient c,
antibiotic a; the contour phi_3D = 0.5.

    python ansys_usermat/apdl/fig_coupled_examples.py   -> assets/fig_coupled_examples.pdf/.png
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

FD = ROOT / "figdata"
TAG = "front0.1_phimin0.1_v2_fortran_ansysbc"
RUNS = [("separate", "A in the seed, B in the corners"), ("separate_ab", "the same, antibiotic from $T^*=0.5$"),
        ("mixed", "A and B mixed in the seed"), ("mixed_ab", "the same, antibiotic from $T^*=0.5$")]


def load(key):
    for name in (f"coupled_living_{key}_{TAG}_{key}.npz", f"coupled_living_{key}_{TAG}.npz"):
        p = FD / name
        if p.exists():
            return np.load(p)
    raise SystemExit(f"missing {key}")


def mid(u):
    k = u.shape[2] // 2
    return 0.5 * (u[:, :, k - 1] + u[:, :, k])


def main():
    figstyle.apply(size=8)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = [("shareA", r"share of A", "coolwarm_r", (0, 1)), ("L", r"$\Sigma_i\bar\phi_i$", "Greens", (0, 1)),
            ("c", r"nutrient $c$", "Blues", (0, 1)), ("a", r"antibiotic $a$", "Reds", (0, 1))]
    fig, ax = plt.subplots(len(RUNS), len(cols), figsize=(6.3, 6.6), squeeze=False)
    for i, (key, lab) in enumerate(RUNS):
        d = load(key)
        phi = d["phi_1.0"]
        n = phi.shape[0]
        e = np.linspace(-1, 1, n + 1)
        xc = 0.5 * (e[1:] + e[:-1])
        for j, (f, t, cm, lim) in enumerate(cols):
            a = ax[i, j]
            u = mid(d[f"{f}_1.0"])
            im = a.pcolormesh(e, e, np.ma.masked_invalid(u).T, cmap=cm, vmin=lim[0], vmax=lim[1], shading="flat",
                              rasterized=True)
            a.contour(xc, xc, mid(phi).T, levels=[0.5], colors="k", linewidths=0.6)
            a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([]); a.grid(False)
            if i == 0:
                a.set_title(t, fontsize=8)
            if j == 0:
                a.set_ylabel(lab, fontsize=7)
            if i == len(RUNS) - 1:
                cb = fig.colorbar(im, ax=a, orientation="horizontal", fraction=0.06, pad=0.04)
                cb.ax.tick_params(labelsize=6)
    fig.tight_layout(h_pad=0.4, w_pad=0.3)
    out = ROOT / "assets" / "fig_coupled_examples"
    fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"), dpi=300)
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

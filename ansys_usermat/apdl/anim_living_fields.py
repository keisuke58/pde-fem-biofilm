#!/usr/bin/env python3
"""anim_living_fields.py -- animation of an ANSYS coupling run (40 times, fields_<run>_dense.json) on the mid-plane
z = 0 for the colloquium: biofilm phi_3D, living volume fraction, share of species A, nutrient c, antibiotic a,
von Mises stress. Writes an animated GIF (no ffmpeg needed) and a PDF of four frames.

    python ansys_usermat/apdl/anim_living_fields.py u1_16_mixeq_ab [u1_16_sep_ab ...]
      -> assets/anim_<run>.gif, assets/frames_<run>.pdf
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

PANELS = [("phi", r"biofilm $\phi_{3D}$", "Greys", (0, 1)),
          ("L", r"living $\sum_i\bar\phi_i$", "Greens", (0, 1)),
          ("shareA", r"share of A", "coolwarm_r", (0, 1)),
          ("nut", r"nutrient $c$", "Blues", (0, 1)),
          ("ab", r"antibiotic $a$", "Reds", (0, 1)),
          ("seqv", r"$\sigma_{vM}$ [mPa]", "magma", None)]


def midplane(d):
    x, y, z = (np.round(np.array(d[k]), 6) for k in ("cx", "cy", "cz"))
    zs = np.unique(np.abs(z))
    zmid = zs.min()
    ax = np.unique(x)
    n = len(ax)
    idx = {v: i for i, v in enumerate(ax)}
    lo = [k for k in range(len(z)) if np.isclose(z[k], -zmid)]
    hi = [k for k in range(len(z)) if np.isclose(z[k], zmid)]

    def grid(u):
        g = np.zeros((n, n)); cnt = np.zeros((n, n))
        for k in lo + hi:
            i, j = idx[x[k]], idx[y[k]]
            if np.isfinite(u[k]):
                g[i, j] += u[k]; cnt[i, j] += 1
        return np.where(cnt > 0, g / np.maximum(cnt, 1), np.nan)
    e = np.concatenate([[ax[0] - (ax[1] - ax[0]) / 2], (ax[1:] + ax[:-1]) / 2, [ax[-1] + (ax[1] - ax[0]) / 2]])
    return grid, e


def main(runs):
    figstyle.apply(size=8)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    for run in runs:
        d = json.loads((HERE / "results/2026-10-living" / f"fields_{run}_dense.json").read_text())
        grid, e = midplane(d)
        t = d["t"]
        data = {}
        for key, *_ in PANELS:
            arr = np.array(d[key], float)
            if key == "seqv":
                arr = arr * 1e9               # MPa -> mPa
            if key in ("L", "shareA"):
                arr = np.where(np.array(d["phi"], float) >= 0.1, arr, np.nan)
            data[key] = [grid(arr[k]) for k in range(len(t))]
        vmax_s = np.nanmax([np.nanmax(g) for g in data["seqv"]])
        fig, axs = plt.subplots(2, 3, figsize=(6.3, 4.4))
        ims = []
        for a, (key, title, cmap, lim) in zip(axs.ravel(), PANELS):
            lo, hi = lim if lim else (0, vmax_s)
            im = a.pcolormesh(e, e, data[key][0].T, cmap=cmap, vmin=lo, vmax=hi, shading="flat")
            a.set_title(title, fontsize=8); a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([]); a.grid(False)
            cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=6)
            ims.append((im, key))
        sup = fig.suptitle("", fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, 0.95))

        def draw(k):
            for im, key in ims:
                im.set_array(data[key][k].T.ravel())
            sup.set_text(f"$T^*$ = {t[k]:.3f}" + ("  (antibiotic from 0.5)" if t[k] > 0.5 else ""))
            return [im for im, _ in ims] + [sup]
        anim = FuncAnimation(fig, draw, frames=len(t), blit=False)
        out = ROOT / "assets" / f"anim_{run}.gif"
        anim.save(out, writer=PillowWriter(fps=5), dpi=110)
        # four frames as a static figure
        ks = [int(np.argmin(np.abs(np.array(t) - s))) for s in (0.25, 0.5, 0.75, 1.0)]
        fig2, axs2 = plt.subplots(len(PANELS), 4, figsize=(6.3, 8.6))
        for r, (key, title, cmap, lim) in enumerate(PANELS):
            lo, hi = lim if lim else (0, vmax_s)
            for c, k in enumerate(ks):
                a = axs2[r, c]
                im = a.pcolormesh(e, e, data[key][k].T, cmap=cmap, vmin=lo, vmax=hi, shading="flat")
                a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([]); a.grid(False)
                if r == 0:
                    a.set_title(f"$T^*={t[k]:.2f}$", fontsize=8)
                if c == 0:
                    a.set_ylabel(title, fontsize=7)
            fig2.colorbar(im, ax=axs2[r, :], fraction=0.02, pad=0.01).ax.tick_params(labelsize=6)
        fig2.savefig(ROOT / "assets" / f"frames_{run}.pdf")
        plt.close("all")
        print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1:] or ["u1_16_mixeq_ab", "u1_16_sep_ab"])

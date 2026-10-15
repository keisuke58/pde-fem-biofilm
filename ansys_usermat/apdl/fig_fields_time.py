"""fig_fields_time.py -- fields on the mid-plane of the 2 mm cube over time, laid out like Klempt's SIIRI poster
(rows: time, columns: variables): phase phi of the field, nutrient c, growth alpha - 1 and von Mises / E for
Eq. 36, and von Mises / E for growth from the living biomass. Calibrated five species (CH, MAP of the measured
data), 16^3, element values of the layer next to z = 0 (A-15: fields_<run>.json; phi = SVAR 41, the phase the
material routine receives; Gauss-point means).

    python ansys_usermat/apdl/fig_fields_time.py   -> assets/fig_fields_time.pdf, .png
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
sys.path.insert(0, str(ROOT / "ansys_usermat/apdl"))
import figstyle  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from fig_cal5_fields import cube  # noqa: E402

R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
E = 10.0  # Pa


def load(fn):
    d = json.load(open(R / fn))
    cx, cy, cz = (np.array(d[k]) for k in ("cx", "cy", "cz"))
    zs = np.unique(np.round(cz, 6)); z0 = zs[zs > 0].min()
    xs = np.unique(np.round(cx, 6)); m = len(xs)
    sel = np.where(np.isclose(cz, z0))[0]
    ix, iy = np.searchsorted(xs, np.round(cx[sel], 6)), np.searchsorted(xs, np.round(cy[sel], 6))
    out = {}
    for k, scale in (("phi", 1.0), ("nut", 1.0), ("alpha", 1e3), ("seqv", 1e6 / E * 1e5)):
        a = np.array(d[k]); g = np.full((len(d["t"]), m, m), np.nan)
        g[:, iy, ix] = a[:, sel] * scale
        out[k] = g
    h = xs[1] - xs[0]
    return d["t"], out, xs[0] - h / 2, xs[-1] + h / 2


def seed_mask():
    """initial seed elements (component BIOFILM1) on the plane next to z = 0"""
    d = json.load(open(R / "w16_CH_g1_s015_raw.json")); a = d["all_stress"]
    c = np.array([a["cx"], a["cy"], a["cz"]]).T; s = np.isin(a["elem"], d["seed_BIOFILM1"])
    xs = np.unique(np.round(c[:, 0], 6)); z0 = np.unique(np.round(c[:, 2], 6)); z0 = z0[z0 > 0].min()
    m = np.zeros((len(xs), len(xs)), bool)
    for i in np.where(s & np.isclose(c[:, 2], z0))[0]:
        m[np.searchsorted(xs, round(c[i, 1], 6)), np.searchsorted(xs, round(c[i, 0], 6))] = True
    return m


if __name__ == "__main__":
    seed = seed_mask()
    t, a, lo, hi = load("fields_w16_CH_g1_s015_raw.json")
    _, b, _, _ = load("fields_w16_CH_cal_k0814_raw.json")
    cols = ((a["phi"], "phase $\\phi$", (0, 1)),
            (a["nut"], "nutrient $c$", None),
            (a["alpha"], "$\\alpha-1$ ($10^{-3}$)", None),
            (a["seqv"], "von Mises$/E$ ($10^{-5}$)", "vm"))
    vmax_vm = np.nanmax(a["seqv"])
    fig = plt.figure(figsize=(10.5, 6.2))
    TS = [0, 1, 3]                                   # T* = 0.25, 0.5, 1 (Day 6, 11, 21)
    cax3 = fig.add_axes([0.0, 0.3, 0.22, 0.45], projection="3d"); cube(cax3)
    for v, lab in (((1.3, -1, -1), "$x$"), ((1, 1.3, -1), "$y$"), ((-1, -1, 1.3), "$z$")):
        cax3.text(*v, lab, fontsize=8)
    fig.text(0.11, 0.24, "each panel: plane $z=0$ seen from above\n(white dashed: initial seed elements;\ngreen: nutrient supply, $y=-1$ mm)",
             ha="center", va="top", fontsize=9)
    x0, y0, W, H, dx, dy = 0.3, 0.70, 0.16, 0.25, 0.175, 0.28
    for j, (data, lab, rng) in enumerate(cols):
        if rng == "vm":
            vmin, vmax = 0.0, vmax_vm
        elif rng is None:
            vmin, vmax = np.nanmin(data), np.nanmax(data)
        else:
            vmin, vmax = rng
        for r, tr in enumerate(TS):
            ax = fig.add_axes([x0 + j * dx, y0 - r * dy, W, H])
            im = ax.imshow(data[tr], origin="lower", extent=(lo, hi, lo, hi), cmap="viridis", vmin=vmin, vmax=vmax)
            ax.contour(np.linspace(lo, hi, seed.shape[0] + 1)[:-1] + (hi - lo) / seed.shape[0] / 2,
                       np.linspace(lo, hi, seed.shape[0] + 1)[:-1] + (hi - lo) / seed.shape[0] / 2,
                       seed.astype(float), levels=[0.5], colors="w", linewidths=0.8, linestyles="--")
            ax.plot([lo, hi], [lo, lo], color="#2e8b57", lw=3.0, solid_capstyle="butt", clip_on=False)
            if j == 0 and r == len(TS) - 1:
                ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1]); ax.tick_params(labelsize=6.5, length=2, pad=1)
                ax.set_xlabel("$x$ [mm]", fontsize=7, labelpad=0); ax.set_ylabel("$y$ [mm]", fontsize=7, labelpad=0)
            else:
                ax.set_xticks([]); ax.set_yticks([])
            if r == 0:
                ax.set_title(lab, fontsize=10.5)
            if j == 0:
                ax.text(-0.3, 0.5, f"Day {1 + 20 * t[tr]:g}\n$T^*={t[tr]:g}$", transform=ax.transAxes, ha="right",
                        va="center", fontsize=10.5)
        cax = fig.add_axes([x0 + j * dx + 0.008, y0 - (len(TS) - 1) * dy - 0.075, W - 0.016, 0.02])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.ax.tick_params(labelsize=8.5)
        cb.ax.locator_params(nbins=4)
    out = ROOT / "assets" / "fig_fields_time"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=170)
    print("wrote", out.with_suffix(".png"))

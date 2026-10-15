"""fig_cal5_fields.py -- fields of the calibrated five-species runs (CH, MAP of the measured data, 16^3, T* = 1)
on the mid-plane z = 0 of the 2 mm cube, laid out like Klempt's SIIRI poster: a 3D sketch of the cube with the
seed and the nutrient supply face, then one row per growth law (Eq. 36 / living biomass) and one column per field
(nutrient c, alpha - 1, von Mises stress; the exported phi columns are
point-model states, not the field, so they are not shown). Element values, the layer of elements next to z = 0.

    python ansys_usermat/apdl/fig_cal5_fields.py   -> assets/fig_cal5_fields.pdf, .png
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
RUNS = (("w16_CH_g1_s015_raw.json", "Eq. 36\n($\\phi$ of the field)"),
        ("w16_CH_cal_k0814_raw.json", "living\nbiomass"))


def plane(fn):
    d = json.load(open(R / fn))
    a, n = d["all_stress"], d["nut_field"]
    ca = np.array([a["cx"], a["cy"], a["cz"]]).T
    cn = np.array([n["cx"], n["cy"], n["cz"]]).T
    zs = np.unique(np.round(ca[:, 2], 6)); z0 = zs[np.argmin(np.abs(zs - 0.0) + (zs < 0))]   # first layer above 0
    xs = np.unique(np.round(ca[:, 0], 6)); m = len(xs)
    out = {k: np.full((m, m), np.nan) for k in ("phi", "c", "alpha", "vm")}
    ix = lambda v: np.searchsorted(xs, np.round(v, 6))
    for i in np.where(np.isclose(ca[:, 2], z0))[0]:
        out["alpha"][ix(ca[i, 1]), ix(ca[i, 0])] = a["alpha"][i] * 1e3
        out["vm"][ix(ca[i, 1]), ix(ca[i, 0])] = a["seqv"][i] * 1e6 / 10.0 * 1e5   # MPa -> Pa, / E = 10 Pa, in 1e-5
    for i in np.where(np.isclose(cn[:, 2], z0))[0]:
        out["phi"][ix(cn[i, 1]), ix(cn[i, 0])] = n["phi1"][i] + n["phi2"][i]
        out["c"][ix(cn[i, 1]), ix(cn[i, 0])] = n["nut1"][i]
    return out, xs[0] - (xs[1] - xs[0]) / 2, xs[-1] + (xs[1] - xs[0]) / 2


def cube(ax):
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=20, azim=-60); ax.set_axis_off()
    for a in (-1, 1):
        for b in (-1, 1):
            ax.plot([-1, 1], [a, a], [b, b], color="#555", lw=0.6)
            ax.plot([a, a], [-1, 1], [b, b], color="#555", lw=0.6)
            ax.plot([a, a], [b, b], [-1, 1], color="#555", lw=0.6)
    s = 0.5
    f = lambda pts: Poly3DCollection([pts], alpha=0.85)
    for pts in ([(-s, -s, -s), (s, -s, -s), (s, s, -s), (-s, s, -s)], [(-s, -s, s), (s, -s, s), (s, s, s), (-s, s, s)],
                [(-s, -s, -s), (s, -s, -s), (s, -s, s), (-s, -s, s)], [(-s, s, -s), (s, s, -s), (s, s, s), (-s, s, s)],
                [(-s, -s, -s), (-s, s, -s), (-s, s, s), (-s, -s, s)], [(s, -s, -s), (s, s, -s), (s, s, s), (s, -s, s)]):
        p = f(pts); p.set_facecolor("#c0392b"); p.set_edgecolor("#7b241c"); ax.add_collection3d(p)
    g = Poly3DCollection([[(-1, -1, -1), (1, -1, -1), (1, -1, 1), (-1, -1, 1)]], alpha=0.45)
    g.set_facecolor("#2e8b57"); ax.add_collection3d(g)
    z = Poly3DCollection([[(-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0)]], alpha=0.15)
    z.set_facecolor("#1e50a0"); z.set_edgecolor("#1e50a0"); ax.add_collection3d(z)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.text2D(0.5, -0.02, "seed (red, sketch), nutrient supply $y=-1$ mm (green),\nplane $z=0$ shown (blue)",
              transform=ax.transAxes, ha="center", va="top", fontsize=7.6)
    ax.set_title("2 mm cube, $16^3$", fontsize=9)


if __name__ == "__main__":
    data = [plane(fn)[0] for fn, _ in RUNS]
    _, lo, hi = plane(RUNS[0][0])
    cols = (("c", "nutrient $c$", "viridis"),
            ("alpha", "$\\alpha-1$ ($10^{-3}$)", "viridis"), ("vm", "von Mises$/E$ ($10^{-5}$)", "viridis"))
    fig = plt.figure(figsize=(10.5, 4.8))
    cube(fig.add_axes([0.0, 0.18, 0.24, 0.7], projection="3d"))
    W, H, x0, y0, dx = 0.2, 0.37, 0.33, 0.5, 0.225
    for j, (k, lab, cm) in enumerate(cols):
        vmin = min(np.nanmin(d[k]) for d in data); vmax = max(np.nanmax(d[k]) for d in data)
        for r, d in enumerate(data):
            ax = fig.add_axes([x0 + j * dx, y0 - r * (H + 0.06), W, H])
            im = ax.imshow(d[k], origin="lower", extent=(lo, hi, lo, hi), cmap=cm, vmin=vmin, vmax=vmax)
            ax.set_xticks([]); ax.set_yticks([])
            if r == 0:
                ax.set_title(lab, fontsize=9)
            if j == 0:
                ax.text(-0.12, 0.5, RUNS[r][1], transform=ax.transAxes, ha="right", va="center", fontsize=8.6)
        cax = fig.add_axes([x0 + j * dx + 0.01, 0.02, W - 0.02, 0.025])
        fig.colorbar(im, cax=cax, orientation="horizontal").ax.tick_params(labelsize=7)
    fig.text(x0 - 0.005, 0.965, "mid-plane $z=0$, $T^*=1$ (Day 21), commensal flow condition, calibrated five species",
             fontsize=8.6, ha="left")
    out = ROOT / "assets" / "fig_cal5_fields"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=170)
    print("wrote", out.with_suffix(".png"))

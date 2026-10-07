#!/usr/bin/env python3
"""figs_beta.py -- the diffusion coefficient beta of the phi field (Klempt et al.
2024 Eq. 34) in the partner's element: seed stresses and the composition behind
the front, from the run JSON files of 5 Oct (no ANSYS needed).

    python ansys_usermat/apdl/figs_beta.py
        -> assets/fig_beta_seed_stress.png, assets/fig_beta_share_front.png

beta = 1e-4 is the element's example input; beta = 0.02 mm^2/T* is Table 2's
beta = 2 converted to the 2 mm cube (KLEMPT2024_REPRODUCTION.md sec. 8).

Figure 1, one species, Klempt 2024 Table 2 stiffness (E = 10 Pa, nu = 0.49):
mid-plane sections of alpha - 1 and of the mean stress p on 16^3 for both
beta; p by layer around the seed (layer k: centroid k * 0.25 mm outside the
seed) and alpha - 1 on a line through the seed, on 8^3 and 16^3.
Figure 2, two species, Klempt et al. 2026 case 6, point model in every
element: mid-plane sections of the share phi_1/(phi_1+phi_2) where the point
model ran (grey: never ran) and of the nutrient c. The consumption values
(4, 6) are example inputs of the element, not paper values; the local nutrient
enters the point model as c* = c*_0 min(c/c_ref, 1), an assumption of this work.
The two-species decks keep the element's example stiffness, so Figure 2 shows
no stresses.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import colors  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402

RES = HERE / "results" / "2026-10-05_ansys"
ASSETS = HERE.parents[1] / "assets"


def load(run):
    rec = json.loads((RES / f"{run}.json").read_text())
    a = rec["all_stress"]
    cen = np.stack([a["cx"], a["cy"], a["cz"]], 1)
    elem = np.array(a["elem"], int)
    seed = np.isin(elem, rec["seed_BIOFILM1"])
    p = (np.array(a["sx"]) + np.array(a["sy"]) + np.array(a["sz"])) / 3 * 1e6   # MPa -> Pa
    out = {"elem": elem, "cen": cen, "seed": seed, "p": p, "alpha": np.array(a["alpha"])}
    if "nut_field" in rec:
        nf = rec["nut_field"]
        pos = {e: i for i, e in enumerate(elem)}
        order = np.array([pos[e] for e in nf["elem"]])
        for k in ("nut1", "phi1", "phi2"):
            v = np.full(len(elem), np.nan)
            v[order] = nf[k]
            out[k] = v
    return out


def section(r, key, z=None):
    """values on the element layer nearest to z = 0 (positive side), as an (x, y) image."""
    cz = np.unique(np.round(r["cen"][:, 2], 6))
    z = cz[cz > 0].min() if z is None else z
    m = np.isclose(r["cen"][:, 2], z)
    xs, ys = np.unique(np.round(r["cen"][m, 0], 6)), np.unique(np.round(r["cen"][m, 1], 6))
    img = np.full((len(ys), len(xs)), np.nan)
    smask = np.zeros_like(img, bool)
    v = r[key] if isinstance(key, str) else key
    for c, val, s in zip(r["cen"][m], v[m], r["seed"][m]):
        i, j = np.searchsorted(ys, round(c[1], 6)), np.searchsorted(xs, round(c[0], 6))
        img[i, j] = val
        smask[i, j] = s
    h = xs[1] - xs[0]
    ext = (xs[0] - h / 2, xs[-1] + h / 2, ys[0] - h / 2, ys[-1] + h / 2)
    return img, smask, ext, z


def outline(ax, smask, ext):
    ny, nx = smask.shape
    hx, hy = (ext[1] - ext[0]) / nx, (ext[3] - ext[2]) / ny
    for i in range(ny):
        for j in range(nx):
            if not smask[i, j]:
                continue
            x0, y0 = ext[0] + j * hx, ext[2] + i * hy
            for di, dj, seg in ((0, -1, ([x0, x0], [y0, y0 + hy])), (0, 1, ([x0 + hx, x0 + hx], [y0, y0 + hy])),
                                (-1, 0, ([x0, x0 + hx], [y0, y0])), (1, 0, ([x0, x0 + hx], [y0 + hy, y0 + hy]))):
                ii, jj = i + di, j + dj
                if not (0 <= ii < ny and 0 <= jj < nx and smask[ii, jj]):
                    ax.plot(*seg, color="k", lw=0.8)


def layer_means(r):
    cen, seed = r["cen"], r["seed"]
    h = np.min(np.diff(np.unique(np.round(cen[:, 0], 9))))
    d = np.full(len(cen), np.inf)
    for b in cen[seed]:
        d = np.minimum(d, np.max(np.maximum(np.abs(cen - b) - h / 2, 0), axis=1))
    lay = np.ceil(d / 0.25 - 1e-9).astype(int)
    lay[seed] = 0
    ks = [k for k in range(0, 4) if (lay == k).any()]
    return ks, [r["p"][lay == k].mean() for k in ks]


def line_x(r, key="alpha"):
    cen = r["cen"]
    cy = np.unique(np.round(cen[:, 1], 6)); y0 = cy[cy > 0].min()
    cz = np.unique(np.round(cen[:, 2], 6)); z0 = cz[cz > 0].min()
    m = np.isclose(cen[:, 1], y0) & np.isclose(cen[:, 2], z0)
    o = np.argsort(cen[m, 0])
    return cen[m, 0][o], r[key][m][o]


def fig_stress():
    runs = {("1e-4", 8): "ds8_ref", ("1e-4", 16): "ds16_pv_eq36",
            ("0.02", 8): "ds8_beta002_dt2", ("0.02", 16): "ds16_beta002_dt2"}
    R = {k: load(v) for k, v in runs.items()}
    fig, ax = plt.subplots(2, 3, figsize=(14, 8.2), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    lab = {"1e-4": r"$\beta = 10^{-4}$ (example input)", "0.02": r"$\beta = 0.02$ mm$^2$/T$^*$ (Klempt 2024, converted)"}
    for row, b in enumerate(("1e-4", "0.02")):
        r = R[(b, 16)]
        img, sm, ext, z = section(r, "alpha")
        im = ax[row, 0].imshow(img, origin="lower", extent=ext, cmap="viridis")
        outline(ax[row, 0], sm, ext)
        cb = fig.colorbar(im, ax=ax[row, 0], shrink=0.85, label=r"$\alpha - 1$")
        cb.formatter.set_powerlimits((-2, 2)); cb.update_ticks()
        img, sm, ext, _ = section(r, "p")
        a = np.nanmax(np.abs(img))
        im = ax[row, 1].imshow(img, origin="lower", extent=ext, cmap="RdBu_r",
                               norm=colors.TwoSlopeNorm(vcenter=0, vmin=-a, vmax=a))
        outline(ax[row, 1], sm, ext)
        cb = fig.colorbar(im, ax=ax[row, 1], shrink=0.85, label=r"mean stress $p$ [Pa]")
        cb.formatter.set_powerlimits((-2, 2)); cb.update_ticks()
        for c in (0, 1):
            ax[row, c].set(xlabel="x [mm]", ylabel="y [mm]")
            ax[row, c].grid(False)
        ax[row, 0].set_title(f"{lab[b]}\n" + r"$\alpha - 1$, 16$^3$, section $z$ = " + f"{z:.3f} mm", fontsize=10)
        ax[row, 1].set_title(r"mean stress $p$ (blue: compression)", fontsize=10)
    col = {"1e-4": "#1f77b4", "0.02": "#d62728"}
    for b, ls in (("1e-4", "--"), ("0.02", "-")):
        for n, mk, fill, ms in ((16, "s", True, 6), (8, "o", False, 9)):     # 8^3 hollow, on top
            kw = dict(color=col[b], marker=mk, ms=ms, mfc=col[b] if fill else "white", mew=1.4,
                      label=f"{lab[b].split(' (')[0]}, {n}$^3$")
            ks, ps = layer_means(R[(b, n)])
            ax[0, 2].plot(ks, ps, ls, **kw)
            x, al = line_x(R[(b, n)])
            ax[1, 2].plot(x, al, ls, **{**kw, "ms": ms * 0.6})
    ax[0, 2].set_yscale("symlog", linthresh=1e-6)
    ax[0, 2].axhline(0, color="k", lw=0.6)
    ax[0, 2].set(xticks=[0, 1, 2, 3], xticklabels=["seed", "layer 1", "layer 2", "layer 3"],
                 ylabel=r"mean $p$ [Pa]", title="mean stress by layer (0.25 mm each)")
    ax[1, 2].set(xlabel="x [mm] (line through the seed)", ylabel=r"$\alpha - 1$", title=r"$\alpha - 1$ across the seed")
    for a_ in (ax[0, 2], ax[1, 2]):
        a_.legend(fontsize=8.5)
    fig.suptitle(r"Seed stress for two values of $\beta$ (one species, $E$ = 10 Pa, $\nu$ = 0.49, $T^*$ = 1.1)", fontsize=12)
    fig.tight_layout()
    out = ASSETS / "fig_beta_seed_stress.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


def fig_share():
    panels = (("ds16_c6_nut_g4", "share", r"$\beta = 10^{-4}$, consumption 4"),
              ("ds16_c6_all_b002", "share", r"$\beta = 0.02$, no local nutrient"),
              ("ds16_c6_nut_g6_b002", "share", r"$\beta = 0.02$, consumption 6"),
              ("ds16_c6_nut_g6_b002", "nut1", r"$\beta = 0.02$, consumption 6: nutrient $c$"))
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.6))
    for k, (run, key, title) in enumerate(panels):
        r = load(run)
        if key == "share":
            s = r["phi1"] + r["phi2"]
            ran = (s > 1e-6) & (s < 0.95)
            v = np.where(ran, r["phi1"] / np.where(ran, s, 1), np.nan)
            img, sm, ext, z = section(r, v)
            cmap = matplotlib.colormaps["coolwarm"].copy(); cmap.set_bad("0.85")
            im = ax[k].imshow(img, origin="lower", extent=ext, cmap=cmap, vmin=0, vmax=0.5)
            fig.colorbar(im, ax=ax[k], shrink=0.8, label=r"$\phi_1/(\phi_1+\phi_2)$")
        else:
            img, sm, ext, z = section(r, "nut1")
            im = ax[k].imshow(img, origin="lower", extent=ext, cmap="Greens", vmin=0, vmax=1)
            fig.colorbar(im, ax=ax[k], shrink=0.8, label=r"nutrient $c$")
        outline(ax[k], sm, ext)
        ax[k].set(xlabel="x [mm]", ylabel="y [mm]", title=title)
        ax[k].title.set_fontsize(10.5)
        ax[k].grid(False)
    fig.suptitle("Composition in the point model, Klempt et al. 2026 case 6, 16$^3$, $T^*$ = 1, section $z$ = "
                 f"{z:.3f} mm; nutrient held at y = −1 mm; grey: point model never ran; consumption: example input",
                 fontsize=11)
    fig.tight_layout()
    out = ASSETS / "fig_beta_share_front.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    figstyle.apply()
    fig_stress()
    fig_share()

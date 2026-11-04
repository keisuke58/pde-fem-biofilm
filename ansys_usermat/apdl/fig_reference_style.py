"""fig_reference_style.py -- four figures of the kind the references use and this work did not yet have
(4 Nov 2026), drawn at the printed width (16 cm), Times via figstyle, vector PDF and 400 dpi PNG.

(a) fig_profiles_y:      profiles along the nutrient direction y through the seed centre at several times
                         (phi_3D, c, sum_i phibar_i, alpha - 1), separate and mixed starts with antibiotic, 64^3
                         (as Feng et al. 2021, Fig. 13, and Klempt et al. 2024, Fig. 6).
(b) fig_hydrostatic:     hydrostatic stress sigma_m = tr(sigma)/3 at T* = 1 on the mid-planes z = 0 and x = 0, red tension,
                         blue compression (as Klempt et al. 2024, Fig. 8), separate and mixed, stress mesh 16^3.
(c) fig_cost_breakdown:  computing time by part, uniform grid against the two- and three-level grids, same machine
                         (as von Zabiensky et al. 2026, Figs. 10 and 11).
(e) fig_spacetime_y:      space-time diagram of phi_3D and sum_i phibar_i along the same line (front position over time).
(d) fig_limiter_grid:    seed von Mises stress of the mixed start with antibiotic against the grid for five
                         reconstructions of the front term, with the extrapolated values.

    python ansys_usermat/apdl/fig_reference_style.py   -> assets/fig_<name>.pdf, .png
"""
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(8)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

OUT = ROOT / "assets"
FD = ROOT / "figdata"
R3 = ROOT / "ansys_usermat/apdl/results/2026-10-3d"
W = 16 / 2.54
BASE = "front1_phimin0.001_v2_fortran_ansysbc_inhliv_{mix}a100_gab1_gc1_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast"
CASES = (("separate_ab", "", "separate start"), ("mixed_ab", "mixeq_", "mixed start"))
T_AB = 1 / 3


def save(fig, name):
    for ext, kw in (("pdf", {}), ("png", {"dpi": 400})):
        fig.savefig(OUT / f"fig_{name}.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)
    print("wrote", OUT / f"fig_{name}.pdf")


def snaps(key, mix):
    f = FD / f"coupled_living_{key}_nf64_{BASE.format(mix=mix)}_{key}.npz"
    d = np.load(f)
    times = sorted({float(k.split("_")[-1]) for k in d.files if k.startswith("phi_")})
    return d, times


def fig_profiles():
    want = (0.05, 0.1, 0.25, 0.5, 1.0)
    rows = (("phi", r"$\phi_{3D}$"), ("c", r"$c$"), ("L", r"$\Sigma_i\,\bar\phi_i$"), ("alpha", r"$\alpha-1$ [$10^{-4}$]"))
    fig, ax = plt.subplots(2, 4, figsize=(W, 0.52 * W), sharex=True)
    cmap = plt.get_cmap("viridis")
    for r, (key, mix, title) in enumerate(CASES):
        d, times = snaps(key, mix)
        ts = [min(times, key=lambda t: abs(t - w)) for w in want]
        n = d["phi_1.0"].shape[0]
        y = -1 + (np.arange(n) + 0.5) * 2 / n
        mid = slice(n // 2 - 1, n // 2 + 1)
        for k, t in enumerate(ts):
            col = cmap(k / (len(ts) - 1) * 0.9)
            for j, (v, lab) in enumerate(rows):
                name = f"{v}_{t:g}" if f"{v}_{t:g}" in d.files else f"{v}_{t}"
                if name not in d.files:
                    continue
                u = d[name][mid, :, mid].mean((0, 2))
                if v == "alpha":
                    u = (u - 1) * 1e4
                ax[r, j].plot(y, u, color=col, lw=1.0, label=fr"$T^*={t:g}$")
        for j, (v, lab) in enumerate(rows):
            a = ax[r, j]
            a.axvspan(-0.25, 0.25, color="0.9", lw=0, zorder=0)
            a.set_title(f"({'abcdefgh'[4 * r + j]}) {lab}", fontsize=8, loc="left")
            a.tick_params(labelsize=7)
            if r == 1:
                a.set_xlabel(r"$y$ [mm]")
        ax[r, 0].set_ylabel(title)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.legend(h, l, loc="upper center", ncol=len(l), fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.0))
    save(fig, "profiles_y")


def fig_hydro():
    fig, ax = plt.subplots(2, 2, figsize=(0.78 * W, 0.74 * W))
    P = []
    for key, mix, _ in CASES:
        d, _ = snaps(key, mix)
        P.append(d["p16_1.0"] * 1e4)
    lim = max(np.abs(p).max() for p in P)
    norm = TwoSlopeNorm(0, -lim, lim)
    n = P[0].shape[0]
    e = np.linspace(-1, 1, n + 1)
    for r, (p, (_, _, title)) in enumerate(zip(P, CASES)):
        mz = 0.5 * (p[:, :, n // 2 - 1] + p[:, :, n // 2])      # plane z = 0: (x, y)
        mx = 0.5 * (p[n // 2 - 1] + p[n // 2])                  # plane x = 0: (y, z)
        for c_, (m, xl, yl, pl) in enumerate(((mz.T, "x", "y", "z = 0"), (mx, "z", "y", "x = 0"))):
            a = ax[r, c_]
            im = a.pcolormesh(e, e, m, cmap="RdBu_r", norm=norm, shading="flat")
            a.set_aspect("equal")
            a.set_xlabel(f"${xl}$ [mm]")
            a.set_ylabel(f"${yl}$ [mm]")
            a.set_title(f"({'abcd'[2 * r + c_]}) {title}, plane ${pl}$", fontsize=8, loc="left")
            a.tick_params(labelsize=7)
            a.grid(False)
    fig.subplots_adjust(wspace=0.45, hspace=0.4)
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.03)
    cb.set_label(r"$\sigma_m=\mathrm{tr}\,\sigma/3$ [$10^{-4}$ Pa] (red: tension)")
    save(fig, "hydrostatic")


def fig_cost():
    # wall times [s] from the logs of the runs (one core each, same cloud machine; README of results/2026-11-amr)
    runs = [("uniform $128^3$", 2230, 1290, 3732 - 2230 - 1290, "3.587 / 8.495"),
            ("two levels, base $64^3$", 576, 798, 1567 - 576 - 798, "3.644 / 8.531"),
            ("uniform $64^3$", None, None, 280, "2.815 / 8.040"),
            ("two levels, base $32^3$", 61, 76, 165 - 61 - 76, "2.813 / 8.040"),
            ("three levels, base $16^3$", 74, 83, 215 - 74 - 83, "2.813 / 8.039")]
    fig, axs = plt.subplots(1, 2, figsize=(W, 0.3 * W), gridspec_kw=dict(width_ratios=(1.15, 1)))
    cols = ("#3b6ea5", "#d08c2b", "#9aa5b1")
    for a, grp, xmax, ttl in ((axs[0], runs[:2], 5200, "(a) $128^3$ cells at the edge"),
                              (axs[1], runs[2:], 520, "(b) $64^3$ cells at the edge")):
        y = np.arange(len(grp))[::-1].astype(float)
        for yi, (lab, nu, pt, rest, val) in zip(y, grp):
            if nu is None:
                a.barh(yi, rest, color="white", edgecolor="0.3", hatch="///", height=0.6)
            else:
                left = 0
                for v, c in zip((nu, pt, rest), cols):
                    a.barh(yi, v, left=left, color=c, height=0.6)
                    left += v
            tot = (nu or 0) + (pt or 0) + rest
            a.text(tot + 0.015 * xmax, yi + 0.12, f"{tot:.0f} s", va="center", fontsize=7)
            a.text(tot + 0.015 * xmax, yi - 0.2, val, va="center", fontsize=6, color="0.35")
        a.set_yticks(y)
        a.set_yticklabels([r[0] for r in grp], fontsize=7.5)
        a.set_xlabel("wall time [s], one core")
        a.set_xlim(0, xmax)
        a.set_ylim(-0.6, len(grp) - 0.4)
        a.set_title(ttl, fontsize=8, loc="left")
        a.spines[["top", "right"]].set_visible(False)
        a.grid(axis="y", visible=False)
    for c, l in zip(cols, ("nutrient", "point model", "fields and other")):
        axs[0].barh(np.nan, 0, color=c, label=l)
    axs[0].barh(np.nan, 0, color="white", edgecolor="0.3", hatch="///", label="total only")
    h, l = axs[0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.legend(h, l, loc="upper center", ncol=4, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.0))
    save(fig, "cost_breakdown")


def _gci(f3, f2, f1, r21=1.5, r32=2.0):
    from scipy.optimize import brentq
    e32, e21 = f3 - f2, f2 - f1
    s = np.sign(e32 / e21)
    p = brentq(lambda p: p - abs(np.log(abs(e32 / e21)) + np.log((r21 ** p - s) / (r32 ** p - s))) / np.log(r21), 0.01, 20)
    return p, (r21 ** p * f1 - f2) / (r21 ** p - 1)


def fig_limiters():
    lim = {"minmod": ("minmod", "o", "#444444"), "vanleer": ("van Leer", "v", "#2a9d8f"), "mc": ("MC", "s", "#3b6ea5"),
           "superbee": ("superbee", "D", "#c0392b"), "weno5": ("WENO5", "^", "#8e44ad")}
    data = {k: {} for k in lim}
    for f in glob.glob(str(R3 / "*front1_dtc0.0025*mixeq*tvd_*mixed_ab.json")):
        if any(t in f for t in ("k26c4", "obst", "transc", "stiff", "eps2", "b1edge", "bcorner", "seedx")):
            continue
        m = re.search(r"tvd_(superbee|mc|vanleer|weno5)_", f)
        k = m.group(1) if m else "minmod"
        nf = re.search(r"nf(\d+)", f)
        nf = int(nf.group(1)) if nf else 32
        r = json.load(open(f))["runs"]["mixed_ab"]
        v = r["seed_vm_Pa"]
        data[k][nf] = (v[-1] if isinstance(v, list) else v) * 1e4
    # minmod: the rows of the grid table (tab:gci of the paper); its files carry other variants with the same name parts
    data["minmod"] = {32: 2.00, 64: 2.82, 128: 3.30, 192: 3.52, 256: 3.63}
    fig, ax = plt.subplots(figsize=(0.62 * W, 0.42 * W))
    for k, (lab, mk, c) in lim.items():
        g = sorted(data[k])
        if not g:
            continue
        ax.plot([1 / n for n in g], [data[k][n] for n in g], marker=mk, color=c, lw=1.0, ms=3.5, label=lab)
        if all(n in data[k] for n in (64, 128, 192)):
            p, ext = _gci(data[k][64], data[k][128], data[k][192])
            ax.plot([0], [ext], marker=mk, mfc="white", color=c, ms=4.5, clip_on=False)
            ax.annotate(f"{ext:.2f} ($p$ = {p:.2f})", (0, ext), xytext=(4, 0), textcoords="offset points",
                        fontsize=6.5, color=c, va="center")
    ax.set_xlim(0, 1 / 30)
    ax.set_xticks([0, 1 / 256, 1 / 128, 1 / 64, 1 / 32])
    ax.set_xticklabels(["0", "1/256", "1/128", "1/64", "1/32"])
    ax.set_xlabel(r"$h/(2\,\mathrm{mm})=1/N$")
    ax.set_ylabel(r"seed von Mises [$10^{-4}$ Pa]")
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "limiter_grid")
    return data


def fig_spacetime():
    """space-time diagram along the line x = z = 0: phi_3D and sum_i phibar_i over y and T*"""
    fig, ax = plt.subplots(2, 2, figsize=(W, 0.5 * W), sharex=True, sharey=True)
    for r, (key, mix, title) in enumerate(CASES):
        d, times = snaps(key, mix)
        n = d["phi_1.0"].shape[0]
        mid = slice(n // 2 - 1, n // 2 + 1)
        e = np.linspace(-1, 1, n + 1)
        ts = np.array(times)
        te = np.concatenate([[ts[0] - 0.5 * (ts[1] - ts[0])], 0.5 * (ts[1:] + ts[:-1]), [ts[-1] + 0.5 * (ts[-1] - ts[-2])]])
        for c_, (v, lab, cm) in enumerate((("phi", r"$\phi_{3D}$", "viridis"), ("L", r"$\Sigma_i\,\bar\phi_i$", "magma"))):
            M = np.array([d[f"{v}_{t:g}"][mid, :, mid].mean((0, 2)) if f"{v}_{t:g}" in d.files else
                          d[f"{v}_{t}"][mid, :, mid].mean((0, 2)) for t in ts])
            a = ax[r, c_]
            im = a.pcolormesh(e, te, M, cmap=cm, vmin=0, vmax=1, shading="flat", rasterized=True)
            a.axhline(T_AB, color="w", lw=0.8, ls="--")
            for yv in (-0.25, 0.25):
                a.axvline(yv, color="w", lw=0.5, ls=":")
            a.set_title(f"({'abcd'[2 * r + c_]}) {lab}, {title}", fontsize=8, loc="left")
            a.grid(False)
            a.tick_params(labelsize=7)
            if r == 1:
                a.set_xlabel(r"$y$ [mm]")
            if c_ == 0:
                a.set_ylabel(r"$T^*$")
            if r == 0:
                fig.colorbar(im, ax=ax[:, c_], shrink=0.9, pad=0.02)
    save(fig, "spacetime_y")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    fig_profiles()
    fig_hydro()
    fig_cost()
    print(fig_limiters())
    fig_spacetime()

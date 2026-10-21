#!/usr/bin/env python3
"""fig_contact_zone.py -- what happens where the A and B regions grow and touch (Oliver, 20 Oct: open question).

Separate example (A in the seed, B in four corner blocks), Python with the boundary conditions of the ANSYS decks
and the Fortran point model, T* = 2, fields every 0.25 (coupled_living_consumption.py --tend 2 --snap-every 0.25).
Two start rules for points the biofilm reaches:
  neighbours : the phi-weighted mean state of the colonised neighbours (composition carried with the biofilm),
  position   : species B within 0.5625 mm of the vertical edges, A elsewhere (the rule of the ANSYS runs).
(a) share of A on the mid-plane z = 0 at T* = 0.5, 1, 1.5, 2, grey where phi_3D < 0.1, contour phi_3D = 0.5;
(a') for the run with antibiotic also the living volume fraction (green), which falls in the B regions while
     the share of A there stays small: A does not take over the space of the dead B;
(b) share of A along the diagonal of the mid-plane from the seed centre to a corner; the contact position
    (share 0.5) and the width of the mixed zone (share between 0.1 and 0.9) are printed and saved.

    python ansys_usermat/apdl/fig_contact_zone.py   -> assets/fig_contact_zone.pdf/.png, results JSON
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

FD = ROOT / "figdata"
TAG = "_front0.1_phimin0.1_v2_fortran_ansysbc{sp}_tend2_snaps_transc1_{key}"   # Eq. 35 with its time derivative, d = 1 mm^2/T*
RUNS = [("neighbours", "separate", ""), ("position", "separate", "_startpos"), ("neighbours, antibiotic", "separate_ab", "")]
TIMES = (0.5, 1.0, 1.5, 2.0)


def load(key, sp):
    name = f"coupled_living_{key}" + TAG.format(sp=sp, key=key) + ".npz"
    p = FD / name
    return np.load(p) if p.exists() else None


def mid(u):
    k = u.shape[2] // 2
    return 0.5 * (u[:, :, k - 1] + u[:, :, k])


def diag(u):
    n = u.shape[0]
    i = np.arange(n // 2, n)
    return u[i, i]


def contact(sh, s):
    ok = np.isfinite(sh)
    if ok.sum() < 3:
        return np.nan, np.nan
    sh, s = sh[ok], s[ok]
    cross = np.flatnonzero((sh[:-1] - 0.5) * (sh[1:] - 0.5) <= 0)
    pos = np.nan
    if len(cross):
        j = cross[0]
        pos = s[j] + (0.5 - sh[j]) * (s[j + 1] - s[j]) / (sh[j + 1] - sh[j])
    mixed = s[(sh > 0.1) & (sh < 0.9)]
    width = (mixed.max() - mixed.min() + (s[1] - s[0])) if len(mixed) else 0.0
    return pos, width


def main():
    figstyle.apply(size=8)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    have = [(lab, key, sp, load(key, sp)) for lab, key, sp in RUNS]
    have = [h for h in have if h[3] is not None]
    rows = [(lab, key, sp, d, "shareA") for lab, key, sp, d in have]
    rows += [(lab + r": $\Sigma_i\bar\phi_i$", key, sp, d, "L") for lab, key, sp, d in have if key.endswith("_ab")]
    if not have:
        raise SystemExit("no runs yet")
    fig = plt.figure(figsize=(6.3, 1.7 * len(rows) + 2.0))
    gs = fig.add_gridspec(len(rows) + 1, len(TIMES), height_ratios=[1] * len(rows) + [1.1])
    out = {}
    for r, (lab, key, sp, d, fld) in enumerate(rows):
        n = d["phi_1.0"].shape[0]
        e = np.linspace(-1, 1, n + 1); xc = 0.5 * (e[1:] + e[:-1])
        out[lab] = {}
        for c, T in enumerate(TIMES):
            a = fig.add_subplot(gs[r, c])
            ph, sh = d[f"phi_{T}"], d[f"{fld}_{T}"]
            m = np.where(mid(ph) >= 0.1, mid(sh), np.nan)
            a.set_facecolor("0.85")
            im_ = a.pcolormesh(e, e, np.ma.masked_invalid(m).T, cmap="coolwarm_r" if fld == "shareA" else "Greens",
                               vmin=0, vmax=1, shading="flat", rasterized=True)
            if fld == "shareA":
                im = im_
            a.contour(xc, xc, mid(ph).T, levels=[0.5], colors="k", linewidths=0.5)
            a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([]); a.grid(False)
            if r == 0:
                a.set_title(f"$T^*={T:g}$", fontsize=8)
            if c == 0:
                a.set_ylabel(lab, fontsize=7)
            s = np.sqrt(2) * xc[n // 2:]
            if fld == "shareA":
                pos, width = contact(diag(m), s)
                out[lab][T] = dict(contact_mm=None if np.isnan(pos) else float(pos), mixed_width_mm=float(width))
            else:
                out[lab][T] = dict(diag_living=[None if not np.isfinite(v) else float(v) for v in diag(m)])
    cb = fig.colorbar(im, ax=fig.axes, fraction=0.02, pad=0.01)
    cb.set_label("share of A", fontsize=7); cb.ax.tick_params(labelsize=6)
    ax = fig.add_subplot(gs[len(rows), :])
    styles = ["-", "--", ":"]
    for (lab, key, sp, d), ls in zip(have, styles):
        n = d["phi_1.0"].shape[0]
        e = np.linspace(-1, 1, n + 1); xc = 0.5 * (e[1:] + e[:-1]); s = np.sqrt(2) * xc[n // 2:]
        for T, col in zip(TIMES, ["0.75", "0.5", "0.3", "0.0"]):
            m = np.where(mid(d[f"phi_{T}"]) >= 0.1, mid(d[f"shareA_{T}"]), np.nan)
            ax.plot(s, diag(m), ls, color=col, lw=0.9, label=f"{lab}, $T^*={T:g}$" if T in (1.0, 2.0) else None)
    ax.set_xlabel("distance from the seed centre along the diagonal [mm]")
    ax.set_ylabel("share of A"); ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=6, ncol=2, frameon=False)
    fig.savefig(ROOT / "assets" / "fig_contact_zone.pdf"); fig.savefig(ROOT / "assets" / "fig_contact_zone.png", dpi=300)
    (HERE / "results/2026-10-3d/contact_zone.json").write_text(json.dumps(out, indent=1) + "\n")
    for lab, v in out.items():
        if "Sigma" not in lab:
            print(lab, {T: (x["contact_mm"] and round(x["contact_mm"], 3), round(x["mixed_width_mm"], 3)) for T, x in v.items()})


if __name__ == "__main__":
    main()

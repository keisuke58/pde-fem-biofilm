"""Figures of the measured biofilm volume of Heine et al. (aggregated, per
condition and day; data/heine_volume/*.csv) and of the growth fit used for the
time link (RESEARCH_IDEAS.ja.md, section 10).

    python ansys_usermat/apdl/volume_figs.py
writes assets/fig_vol_{volume,live,species,fit,stress}.pdf
"""
import csv, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle
figstyle.apply()
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = ROOT / "data" / "heine_volume"
OUT = ROOT / "assets"
COND = [("Commensal", "Static", "CS"), ("Commensal", "HOBIC", "CH"),
        ("Dysbiotic", "Static", "DS"), ("Dysbiotic", "HOBIC", "DH")]
COL = {"CS": "#1f77b4", "CH": "#2ca02c", "DS": "#d4a017", "DH": "#e8620c"}


def table(name):
    out = {}
    for r in csv.DictReader(open(D / name)):
        key = next(c for a, b, c in COND if a == r["condition"] and b == r["cultivation"])
        out.setdefault(key, []).append(r)
    return out


def series(rows, val, sd=None):
    rows = sorted(rows, key=lambda r: int(r["day"]))
    d = np.array([int(r["day"]) for r in rows])
    v = np.array([float(r[val]) for r in rows])
    s = np.array([float(r[sd]) for r in rows]) if sd else None
    n = np.array([int(r["n_biofilms"]) for r in rows]) if "n_biofilms" in rows[0] else None
    return d, v, s, n


def fig_volume():
    T = table("volume_summary.csv")
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
    for i, grp in enumerate((("CS", "CH"), ("DS", "DH"))):
        for c in grp:
            d, v, s, n = series(T[c], "mean_biofilm", "sd_biofilm")
            ax[i].errorbar(d, v, yerr=s, marker="o", ms=4, capsize=3, color=COL[c], label=c)
        ax[i].set(xlabel="day", ylabel="biofilm volume ($10^6\\,\\mu$m$^3$ per field)",
                  title=["commensal", "dysbiotic"][i], ylim=(0, 1.1), xticks=[1, 3, 6, 10, 15, 21])
        ax[i].legend(frameon=False)
    fig.tight_layout(); fig.savefig(OUT / "fig_vol_volume.pdf")


def fig_live():
    T = table("live_fraction_summary.csv")
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for _, _, c in COND:
        d, v, s, _ = series(T[c], "mean_biofilm", "sd_biofilm")
        ax.errorbar(d, v, yerr=s, marker="o", ms=4, capsize=3, color=COL[c], label=c)
    ax.set(xlabel="day", ylabel="live fraction (LIVE/DEAD)", ylim=(0, 1.05), xticks=[1, 3, 6, 10, 15, 21])
    ax.legend(frameon=False, ncol=4, loc="lower center", fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_vol_live.pdf")


def fig_species():
    T = table("live_volume_species.csv")
    sp = ["So", "An", "Vei", "Fn", "Pg"]
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.0), sharey=True)
    for ax, (_, _, c) in zip(axs, COND):
        rows = T[c]
        days = sorted({int(r["day"]) for r in rows})
        bottom = np.zeros(len(days))
        for s in sp:
            v = np.array([float(next(r["median"] for r in rows if int(r["day"]) == d and r["species"] == s))
                          for d in days])
            ax.bar(range(len(days)), v, bottom=bottom, label=s)
            bottom += v
        ax.set(xticks=range(len(days)), xticklabels=days, xlabel="day", title=c)
    axs[0].set_ylabel("live volume per species ($10^6\\,\\mu$m$^3$)")
    axs[-1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_vol_species.pdf")


def growth_fit(c):
    d, v, s, n = series(table("volume_summary.csv")[c], "mean_biofilm", "sd_biofilm")
    se = s / np.sqrt(n)
    i0 = int(np.argmin(v)); m = d >= d[i0]
    dd, vv, ss = d[m], v[m], se[m]
    best = None
    for V0 in np.linspace(0.5 * vv[0], 1.5 * vv[0], 201):
        for r in np.linspace(0, 0.2, 2001):
            q = (((V0 * (1 + r * (dd - dd[0])) ** 3 - vv) / ss) ** 2).sum()
            if best is None or q < best[0]:
                best = (q, V0, r)
    return dd[0], best[1], best[2], best[0], len(dd) - 2


def fig_fit():
    T = table("volume_summary.csv")
    fig, ax = plt.subplots(figsize=(5.2, 3.3))
    for c in ("CH", "DS"):
        d, v, s, _ = series(T[c], "mean_biofilm", "sd_biofilm")
        ax.errorbar(d, v, yerr=s, fmt="o", ms=4, capsize=3, color=COL[c], label=f"{c} measured")
        d0, V0, r, q, dof = growth_fit(c)
        t = np.linspace(d0, 21, 100)
        ax.plot(t, V0 * (1 + r * (t - d0)) ** 3, color=COL[c],
                label=f"{c} fit: $\\alpha-1$ = {r:.4f}/day")
    t = np.linspace(3, 21, 100)
    ax.plot(t, 0.24 * (1 + 6.75e-4 / 21 * (t - 3)) ** 3, "k--", lw=1,
            label="element, $k_\\alpha=10^{-3}$, $T^*=1$ = 21 days")
    ax.set(xlabel="day", ylabel="biofilm volume ($10^6\\,\\mu$m$^3$)", ylim=(0, 1.0), xticks=[1, 3, 6, 10, 15, 21])
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig_vol_fit.pdf")


def fig_stress():
    eps = np.linspace(0, 0.22, 100)
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    for E, lab in ((41e3, "$E$ = 41 kPa"), (600.0, "$E$ = 0.6 kPa"), (10.0, "$E$ = 10 Pa (Klempt 2024)")):
        ax.plot(eps, E * eps / (1 - 0.49), label=lab)
    ax.axvspan(0.158, 0.196, color="0.85", label="measured growth strain")
    ax.set(yscale="log", xlabel="growth strain $\\varepsilon_g=\\ln J/3$", ylabel="film stress bound (Pa)", ylim=(1e-1, 1e5))
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_vol_stress.pdf")


if __name__ == "__main__":
    fig_volume(); fig_live(); fig_species(); fig_fit(); fig_stress()
    for c in ("CH", "DS"):
        print(c, growth_fit(c))

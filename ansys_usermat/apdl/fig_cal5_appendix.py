"""fig_cal5_appendix.py -- calibrated five-species point model (commensal flow condition, CH, MAP of the
measured data) in the partner element, for the appendix of the 19 Oct deck:
(a) shares phi_i/sum(phi) at a seed Gauss point with the coupling of Eq. 36, (b) living fractions phi_i psi_i,
(c) seed means of alpha - 1 and von Mises stress for Eq. 36 and for growth from the living biomass, 8^3 and 16^3.

    python ansys_usermat/apdl/fig_cal5_appendix.py   -> assets/fig_cal5_appendix.pdf, .png
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

R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
NAMES = ["So", "An", "Vei", "Fn", "Pg"]


def trace(fn):
    """Point-model state at the first traced Gauss point, the converged row of each substep."""
    pm = json.load(open(R / fn)); C = {c: i for i, c in enumerate(pm["columns"])}
    rows = pm["last_rows"] if "seed_trace" not in pm else pm["seed_trace"]
    r0 = rows[0]
    rows = [r for r in rows if r[C["elem"]] == r0[C["elem"]] and r[C["ip"]] == r0[C["ip"]]]
    d = {}
    for r in rows:
        d[(int(r[C["ldstep"]]), int(r[C["isubst"]]))] = r
    rows = [d[k] for k in sorted(d)]
    t = np.arange(1, len(rows) + 1) / len(rows)
    phi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(1, 6)] for r in rows])
    psi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(7, 12)] for r in rows])
    return t, phi, psi


def seed_means(fn):
    r = json.load(open(R / fn)); a = r["all_stress"]
    s = np.isin(np.array(a["elem"]), r["seed_BIOFILM1"])
    return np.mean(np.array(a["alpha"])[s]), np.mean(np.array(a["seqv"])[s]) * 1e6   # MPa -> Pa


if __name__ == "__main__":
    t, phi, psi = trace("w8_CH_g1_s015_raw_pm.json")
    tl, phil, psil = trace("w8_CH_cal_k0814_raw_pm.json")
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.0), gridspec_kw=dict(width_ratios=[1, 1, 1.05]))
    share = phi / phi.sum(axis=1, keepdims=True)
    for i in range(5):
        l, = ax[0].plot(t, share[:, i], lw=1.4, label=NAMES[i])
        ax[1].plot(tl, phil[:, i] * psil[:, i], lw=1.4, color=l.get_color())
    ax[0].set(xlabel="$T^*$", ylabel="$\\phi_i/\\sum_j\\phi_j$", title="(a) shares, Eq. 36 coupling", ylim=(0, 1))
    ax[0].legend(fontsize=7.5, frameon=False, ncol=2, loc="upper right")
    ax[1].plot(tl, (phil * psil).sum(axis=1), color="#222", lw=1.6, label="living $\\sum_i\\phi_i\\psi_i$")
    ax[1].set(xlabel="$T^*$", ylabel="$\\phi_i\\psi_i$", title="(b) living fractions, point model")
    ax[1].legend(fontsize=7.5, frameon=False, loc="upper right")

    vals = {}
    for n in (8, 16):
        vals[n] = (seed_means(f"w{n}_CH_g1_s015_raw.json"), seed_means(f"w{n}_CH_cal_k0814_raw.json"))
    x = np.arange(2); w = 0.36
    for k, (lab, col) in enumerate((("Eq. 36", "#4c72b0"), ("living biomass", "#dd8452"))):
        v = [vals[n][k][1] * 1e4 for n in (8, 16)]
        b = ax[2].bar(x + (k - 0.5) * w, v, w * 0.92, color=col, label=lab)
        for xi, vi, n in zip(x, v, (8, 16)):
            if k == 1:
                d = vals[n][1][1] / vals[n][0][1] - 1
                ax[2].text(xi + 0.5 * w, vi + 0.08, f"{100 * d:+.1f} %", ha="center", fontsize=8)
    ax[2].set_xticks(x, ["$8^3$", "$16^3$"])
    ax[2].set(ylabel="seed von Mises ($10^{-4}$ Pa)", title="(c) seed stress, $T^*=1$", ylim=(0, 7.2))
    ax[2].legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.tight_layout()
    out = ROOT / "assets" / "fig_cal5_appendix"
    fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"), dpi=170)
    for n in (8, 16):
        (a0, v0), (a1, v1) = vals[n]
        print(n, f"alpha-1 {a0:.5e} {a1:.5e}  vM {v0:.5e} {v1:.5e}  {100 * (v1 / v0 - 1):+.2f} %")
    print("final shares (Eq. 36):", np.round(share[-1], 3), " living sum start/end:",
          round((phil * psil).sum(1)[0], 3), round((phil * psil).sum(1)[-1], 3))

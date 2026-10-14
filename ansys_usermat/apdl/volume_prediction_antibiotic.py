"""Predicted volume of the commensal flow biofilm with an antibiotic, for growth from phi
(scheme a) and from the living biomass (scheme c), at growth rates fitted to the measured
volume (x1.80 from Day 3 to Day 21, data/heine_volume/volume_summary.csv).

Uses the point-model traces of the element runs (first seed Gauss point, results/2026-10-cal5).
The point model does not depend on alpha (one-way coupling), so int phi dt and
int sum phi_i psi_i dt from these traces hold for any growth rate; only the rate is rescaled:
  alpha(T*) = 1 + k * int_0^T* g dT*,  V/V0 = alpha^3,  (a) g = phi_used, (c) g = sum phi_i psi_i.
The rate is fixed by V(1)/V(0) = 1.80 without antibiotic.  The mechanics at these large growths is
not computed here (the element stops at T* ~ 0.35 for kappa 2.8, see A-8)."""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle
figstyle.apply()
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
TARGET = 1.80


def trace(fn):
    pm = json.load(open(R / fn)); C = {k: i for i, k in enumerate(pm["columns"])}
    rows = pm["seed_trace"]; rows = [r for r in rows if r[C["ip"]] == rows[0][C["ip"]] and r[C["elem"]] == rows[0][C["elem"]]]
    d = {int(r[C["isubst"]]): r for r in rows}; k = sorted(d); rows = [d[i] for i in k]
    t = np.array([0.0] + list(np.array(k) / max(k)))
    phi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(1, 6)] for r in rows])
    psi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(7, 12)] for r in rows])
    pu = np.array([float(r[C["phi_used"]]) for r in rows])
    live = np.concatenate([[(phi * psi).sum(1)[0]], (phi * psi).sum(1)])
    pu = np.concatenate([[pu[0]], pu])
    return t, pu, live


def cumint(t, g):
    return np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(t))])


def rate(I1):
    return (TARGET ** (1 / 3) - 1) / I1


runs = [("w8_CH_cal_k0882_pm.json", "w8_CH_g1_s015_pm.json", "no antibiotic", "k", "-"),
        ("w8_CH_cal_k0882_ab05_pm.json", "w8_CH_g1_s015_ab05_pm.json", "all $b_i=1$ from $T^*=0.5$", "#4c72b0", "--"),
        ("w8_CH_cal_k0882_ab05_ana_pm.json", "w8_CH_g1_s015_ab05_ana_pm.json", "anaerobes from $T^*=0.5$", "#c44e52", "--"),
        ("w8_CH_cal_k0882_ab_pm.json", "w8_CH_g1_s015_ab_pm.json", "all $b_i=1$ from $T^*=0$", "#4c72b0", ":")]

t0, _, lv0 = trace(runs[0][0]); kc = rate(cumint(t0, lv0)[-1])
ta0, pu0, _ = trace(runs[0][1]); ka = rate(cumint(ta0, pu0)[-1])
print(f"rates fitted to x{TARGET} (per T*, at this Gauss point; not the element's prop(7)): (a) k = {ka:.3f}, (c) kappa = {kc:.3f}")
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.7), sharey=True)
for fc, fa, lab, c, ls in runs:
    tc, _, lv = trace(fc); ta, pu, _ = trace(fa)
    Va = (1 + ka * cumint(ta, pu)) ** 3; Vc = (1 + kc * cumint(tc, lv)) ** 3
    ax[0].plot(ta, Va, color=c, ls=ls, label=lab); ax[1].plot(tc, Vc, color=c, ls=ls, label=lab)
    print(f"{lab:28s} V(1)/V(0): (a) {Va[-1]:.3f}  (c) {Vc[-1]:.3f}")
for a, ttl in zip(ax, ["(a) growth from $\\phi$", "(c) growth from the living biomass"]):
    a.axvline(0.5, color="0.5", lw=0.7)
    a.set(xlabel="$T^*$ (1 = Day 21)", title=ttl)
ax[0].set_ylabel("volume $V/V_0=\\alpha^3$")
ax[1].legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(ROOT / "assets/fig_volume_prediction_antibiotic.pdf")

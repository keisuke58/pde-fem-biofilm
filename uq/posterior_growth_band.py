"""posterior_growth_band.py -- step 1 of the stress-interval idea (KEIO_PLAN.ja.md section 5): how much the growth
of a seed varies over the final posterior of chapter 3, when growth follows the living amount (scheme (c),
d alpha/dt / alpha = kappa phi_3D sum_i phi_i psi_i).

Data: data/posterior_final/<C>_particles.npz (EML export of 19 Oct 2026): 200 particles per condition, point-model
trajectories on 0-21 days. Seed packed (phi_3D = 1), T* = t / 21 (s = 0.2375, as the ANSYS runs of scheme (c)).
Then ln alpha(T*) = kappa * I(T*), I(T*) = int_0^T* sum phi psi dT'. For the small growth of the runs (alpha - 1 ~ 1e-3)
alpha - 1 = kappa I, and the seed stress is proportional to alpha - 1 (one seed, linear range). So the spread of I over
the particles is the spread of the growth and of the seed stress, relative to the MAP, without a finite element run.
With Eq. 36 (the thesis scheme) the growth does not depend on the composition and has no spread.

    python3 uq/posterior_growth_band.py   -> uq/posterior_growth_band.json, assets/fig_posterior_growth_band.png
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(8.5)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CONDS = [("CS", "#1e50a0"), ("CH", "#2a9d8f"), ("DS", "#e09f3e"), ("DH", "#c0392b")]
DATA = ROOT / "data" / "posterior_final"


def cumtrapz(y, x):
    return np.concatenate([np.zeros(y.shape[:-1] + (1,)), np.cumsum(0.5 * (y[..., 1:] + y[..., :-1]) * np.diff(x), -1)], -1)


out = {}
fig, ax = plt.subplots(1, 3, figsize=(6.3, 2.4), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
for k, (c, col) in enumerate(CONDS):
    d = np.load(DATA / f"{c}_particles.npz")
    t = d["t_days"]; T = t / 21.0
    L = d["sum_phipsi"]; Lm = d["map_sum_phipsi"]
    I = cumtrapz(L, T); Im = cumtrapz(Lm, T)
    rel = I[:, -1] / Im[-1]                       # growth (and seed stress) at T* = 1 relative to the MAP
    q = np.percentile(rel, [5, 50, 95])
    qL = np.percentile(L, [5, 50, 95], axis=0)
    out[c] = dict(I_map=float(Im[-1]), rel_5_50_95=[round(float(v), 3) for v in q],
                  spread_90_percent=round(float(q[2] - q[0]), 3),
                  live_day21_5_50_95=[round(float(v), 3) for v in np.percentile(L[:, -1], [5, 50, 95])])
    ax[0].fill_between(t, qL[0], qL[2], color=col, alpha=0.25, lw=0)
    ax[0].plot(t, Lm, color=col, lw=1.1, label=c)
    ax[1].fill_between(T, np.percentile(I, 5, 0), np.percentile(I, 95, 0), color=col, alpha=0.25, lw=0)
    ax[1].plot(T, Im, color=col, lw=1.1)
    parts = ax[2].violinplot([rel], positions=[k], widths=0.7, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(col); b.set_alpha(0.45)
    ax[2].plot([k - 0.25, k + 0.25], [q[1]] * 2, color=col, lw=1.2)
ax[0].set(xlabel="day (0 = initial state)", ylabel=r"$\Sigma\phi_i\psi_i$", title="(a) living amount", ylim=(0, 1.02))
ax[0].legend(frameon=False, ncol=2, fontsize=7.5, loc="lower left")
ax[1].set(xlabel="$T^*$", ylabel=r"$\int_0^{T^*}\Sigma\phi_i\psi_i\,dT'$", title="(b) its integral")
ax[2].axhline(1.0, color="0.6", lw=0.8)
ax[2].set(xticks=range(4), xticklabels=[c for c, _ in CONDS], ylabel="growth / MAP at $T^*=1$",
          title="(c) seed growth, stress")
fig.tight_layout(w_pad=1.0)
fig.savefig(ROOT / "assets" / "fig_posterior_growth_band.png", dpi=250)
(Path(__file__).with_name("posterior_growth_band.json")).write_text(json.dumps(out, indent=1))
for c in out:
    print(c, out[c])

"""fig_limiters_general.py -- textbook figure of the slope limiters used for the front term (minmod, van Leer, MC, superbee).

(a) Sweby diagram: the limiter psi(r) against the ratio of successive gradients r, with the second-order TVD region
    (between minmod and superbee, Sweby 1984).
(b) 1D linear advection u_t + u_x = 0 on a periodic interval, 100 cells, MUSCL reconstruction with each limiter, upwind
    flux, SSP-RK3, CFL 0.4: a square pulse and a smooth bump after one period, against the exact profile and first-order
    upwind. Minmod smears the edge most, superbee keeps it sharpest (and steepens the smooth bump).

    python ansys_usermat/apdl/fig_limiters_general.py -> assets/fig_limiters_general.pdf/.png
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402
figstyle.apply()

LIM = {
    "minmod": lambda r: np.maximum(0, np.minimum(1, r)),
    "van Leer": lambda r: (r + np.abs(r)) / (1 + np.abs(r)),
    "MC": lambda r: np.maximum(0, np.minimum(np.minimum(2 * r, 0.5 * (1 + r)), 2)),
    "superbee": lambda r: np.maximum(0, np.maximum(np.minimum(2 * r, 1), np.minimum(r, 2))),
}
COL = {"minmod": "#1c7ed6", "van Leer": "#2b8a3e", "MC": "#e8590c", "superbee": "#c2255c", "upwind": "0.55"}


def rhs(u, lim, h):
    du_m = u - np.roll(u, 1); du_p = np.roll(u, -1) - u
    if lim is None:
        ul = u
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(np.abs(du_p) > 1e-14, du_m / du_p, 0.0)
        ul = u + 0.5 * LIM[lim](r) * du_p          # left state at the face i+1/2 (speed +1)
    f = ul
    return -(f - np.roll(f, 1)) / h


def advect(u0, lim, n_per=1, cfl=0.4):
    n = u0.size; h = 1.0 / n; dt = cfl * h; steps = int(round(n_per / dt)); u = u0.copy()
    for _ in range(steps):                          # SSP-RK3
        u1 = u + dt * rhs(u, lim, h)
        u2 = 0.75 * u + 0.25 * (u1 + dt * rhs(u1, lim, h))
        u = u / 3 + 2 / 3 * (u2 + dt * rhs(u2, lim, h))
    return u


n = 100
x = (np.arange(n) + 0.5) / n
u0 = np.where((x > 0.1) & (x < 0.4), 1.0, 0.0) + np.where((x > 0.55) & (x < 0.9), np.sin(np.pi * (x - 0.55) / 0.35) ** 2, 0.0)

fig, ax = plt.subplots(1, 2, figsize=(6.3, 2.5), gridspec_kw=dict(width_ratios=[1, 1.6]))
r = np.linspace(0, 3, 400)
ax[0].fill_between(r, LIM["minmod"](r), LIM["superbee"](r), color="0.9", lw=0)
for k, f in LIM.items():
    ax[0].plot(r, f(r), color=COL[k], lw=1.3, label=k)
ax[0].set_xlim(0, 3); ax[0].set_ylim(0, 2.2)
ax[0].set_xlabel(r"$r$ (ratio of successive gradients)", fontsize=8); ax[0].set_ylabel(r"$\psi(r)$", fontsize=8)
ax[0].text(2.5, 1.04, "second-order\nTVD region", fontsize=7, ha="center", va="bottom", color="0.35")
ax[0].legend(frameon=False, fontsize=7, loc="upper left"); ax[0].tick_params(labelsize=7)
ax[0].set_title("(a) Sweby diagram", fontsize=8.5, loc="left")

ax[1].plot(x, u0, color="k", lw=1.0, ls=":", label="exact")
ax[1].plot(x, advect(u0, None), color=COL["upwind"], lw=1.0, label="first-order upwind")
for k in LIM:
    ax[1].plot(x, advect(u0, k), color=COL[k], lw=1.2, label=k)
ax[1].set_xlim(0, 1); ax[1].set_ylim(-0.42, 1.15)
ax[1].set_yticks(np.arange(0, 1.01, 0.2)); ax[1].set_xlabel(r"$x$", fontsize=8); ax[1].set_ylabel(r"$u$", fontsize=8); ax[1].tick_params(labelsize=7)
ax[1].legend(frameon=False, fontsize=6.5, loc="lower center", ncol=3)
ax[1].set_title("(b) advection after one period (100 cells, CFL 0.4)", fontsize=8.5, loc="left")
fig.tight_layout()
A = HERE.parents[1] / "assets"
fig.savefig(A / "fig_limiters_general.pdf"); fig.savefig(A / "fig_limiters_general.png", dpi=400)
print("wrote", A / "fig_limiters_general.pdf")

"""fig_fritsch_case1.py -- the point model alone for Fritsch et al. 2025, case I (two species, true parameters), the
species data of the coupled examples.

A = [[1, 0.1], [0.1, 1]], b = (1, 2), eta = (1, 2), start phi = (0.25, 0.30), c* = 100, 1000 steps of 1e-4, so
T* = 1 is the whole run (time link s = 0.1). Three antibiotic settings, nutrient held at c* = 100 (c = 1):
(a) no antibiotic; (b) the coupled examples: alpha* = 100 from T* = 1/3 (Klempt 2026, case 4); (c) Fritsch's own case I,
alpha* = 10 from the start. Solid: volume fraction phi_i; dashed: living volume fraction phibar_i = phi_i psi_i.

    python ansys_usermat/apdl/fig_fritsch_case1.py -> assets/fig_fritsch_case1.pdf/.png
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(ROOT / "JAXFEM"))
import figstyle  # noqa: E402
import klempt2026_reproduction as KR  # noqa: E402
from fritsch2025_cases import CASES  # noqa: E402

figstyle.apply()
base = dict(CASES["fritsch_2sp"])
KR.DT = base["dt"]
STEPS = base["steps"]
SETTINGS = [("(a) no antibiotic", 0.0),
            (r"(b) $\alpha^*=100$ from $T^*=1/3$", lambda k: 100.0 if k > STEPS / 3 else 0.0),
            (r"(c) $\alpha^*=10$ from $T^*=0$", 10.0)]
COL = ("#1c6fb0", "#c0392b")
LAB = ("A", "B")

fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.35), sharey=True, constrained_layout=True)
for ax, (title, al) in zip(axs, SETTINGS):
    case = dict(base, alpha_star=al)
    r = KR.simulate(case)
    t = np.arange(r["phi"].shape[0]) / STEPS                 # T*
    phibar = r["phi"] * r["psi"]
    for i in range(2):
        ax.plot(t, r["phi"][:, i], color=COL[i], lw=1.3, label=rf"$\phi_{i+1}$ ({LAB[i]})")
        ax.plot(t, phibar[:, i], color=COL[i], lw=1.1, ls="--", label=rf"$\bar\phi_{i+1}$")
    ax.plot(t, phibar.sum(1), color="0.25", lw=1.0, ls=":", label=r"$\Sigma_i\bar\phi_i$")
    if callable(al):
        ax.axvline(1 / 3, color="0.6", lw=0.7, ls=(0, (1, 2)))
    ax.set_title(title, fontsize=8, loc="left")
    ax.set_xlabel(r"$T^*$", fontsize=8); ax.set_xlim(0, 1); ax.set_ylim(0, 1.0); ax.tick_params(labelsize=7)
    print(title, "end phi", np.round(r["phi"][-1], 3), "psi", np.round(r["psi"][-1], 3), "sum phibar", round(float(phibar[-1].sum()), 3))
axs[0].set_ylabel("volume fraction", fontsize=8)
axs[0].legend(frameon=False, fontsize=6.3, loc="center right")
A = ROOT / "assets"
fig.savefig(A / "fig_fritsch_case1.pdf"); fig.savefig(A / "fig_fritsch_case1.png", dpi=400)
print("wrote", A / "fig_fritsch_case1.png")

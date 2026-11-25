"""fig_point_model_survival.py -- the point model alone (Fritsch et al. 2025, case I, true parameters), to interpret the
coupled examples (Oliver, 25 Nov 2026). (a) Living volume fraction sum_i phibar_i under the antibiotic alpha* = 100 from
T* = 1/3 for different start compositions and nutrient levels c (c* = 100 c, held): only a point dominated by A with
c close to 1 keeps living biomass. (b) Share phi_1/(phi_1+phi_2) without antibiotic from Fritsch's start (0.25, 0.30):
A overtakes B because eta = (1, 2); with equal or swapped eta, B wins.

    python ansys_usermat/apdl/fig_point_model_survival.py -> assets/fig_point_model_survival.pdf/.png
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
base = CASES["fritsch_2sp"]
KR.DT = base["dt"]
N = base["steps"]
ab = lambda k: 100.0 if k > N / 3 else 0.0   # noqa: E731
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.5), constrained_layout=True)
runs = [("A only, $c = 1$", [0.55, 1e-6], 1.0, "#1c3faa", "-"), ("A only, $c = 0.85$", [0.55, 1e-6], 0.85, "#4c6ef5", "-"),
        ("A only, $c = 0.7$", [0.55, 1e-6], 0.7, "#91a7ff", "-"), ("B only, $c = 1$", [1e-6, 0.55], 1.0, "#c92a2a", "--"),
        ("A, B equal, $c = 1$", [0.275, 0.275], 1.0, "#7048a8", ":")]
for lab, pi, c, col, ls in runs:
    r = KR.simulate(dict(base, c_star=100.0 * c, alpha_star=ab, phi_init=pi))
    t = np.arange(r["phi"].shape[0]) / N
    lv = (r["phi"] * r["psi"]).sum(1)
    axs[0].plot(t, lv, color=col, ls=ls, lw=1.2, label=lab)
    print(f"(a) {lab:22s} living at T*=1: {lv[-1]:.3f}")
axs[0].axvline(1 / 3, color="0.5", ls=(0, (1, 2)), lw=0.8)
axs[0].set_ylabel(r"$\Sigma_i\bar\phi_i$", fontsize=8); axs[0].set_ylim(0, 1.02)
axs[0].set_title(r"(a) antibiotic $\alpha^* = 100$ from $T^* = 1/3$", fontsize=8, loc="left")
axs[0].legend(frameon=False, fontsize=6.3, loc="lower left")
for lab, eta, col, ls in (("$\\eta = (1, 2)$, as in the examples", [1.0, 2.0], "0.15", "-"), ("$\\eta = (1, 1)$", [1.0, 1.0], "#e8590c", "--"),
                          ("$\\eta = (2, 1)$", [2.0, 1.0], "#2f9e44", ":")):
    r = KR.simulate(dict(base, alpha_star=0.0, eta=eta))
    t = np.arange(r["phi"].shape[0]) / N
    sh = r["phi"][:, 0] / r["phi"].sum(1)
    axs[1].plot(t, sh, color=col, ls=ls, lw=1.2, label=lab)
    print(f"(b) eta {eta}: share A at T*=1/3, 1: {sh[N // 3]:.3f} {sh[-1]:.3f}")
axs[1].set_ylabel(r"$\phi_1/(\phi_1+\phi_2)$", fontsize=8); axs[1].set_ylim(0, 1.02)
axs[1].set_title(r"(b) no antibiotic, $c = 1$, start $(0.25, 0.30)$", fontsize=8, loc="left")
axs[1].legend(frameon=False, fontsize=6.3, loc="upper left")
for ax in axs:
    ax.set_xlim(0, 1); ax.set_xlabel(r"$T^*$", fontsize=8); ax.tick_params(labelsize=7)
A = ROOT / "assets"
fig.savefig(A / "fig_point_model_survival.pdf"); fig.savefig(A / "fig_point_model_survival.png", dpi=400)
print("wrote", A / "fig_point_model_survival.png")

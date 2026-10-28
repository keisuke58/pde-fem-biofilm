"""make_dieckow_loo.py -- thesis Section 3.6: leave-one-patient-out prediction of the in-vivo guild compositions
(10 patients of Dieckow et al.; gLV, eHOMD + pFBA sign prior, the best model of Table 3.6, mean RMSE 0.049).

Data: figures/data/dieckow_loo_fits.json, read back from the vector graphics of the internship figure
dieckow_loo_fits_thesis.pdf (bar heights and markers; every week sums to 1 and the printed per-patient RMSE is
reproduced to three digits). Writes
    figures/dieckow_loo_patients.pdf   observed W1-W3 and predicted W2, W3 per patient
    figures/dieckow_loo_scatter.pdf    predicted against observed, all patients and guilds

    python3 thesis/figures/make_dieckow_loo.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(9)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

D = json.load(open(HERE / "data" / "dieckow_loo_fits.json"))
SHORT = ["Actino.", "Bacilli", "Bactero.", "β-Prot.", "Clostr.", "Corio.", "Fuso.", "γ-Prot.", "Negat.", "Other"]
COL = ["#97683a", "#22b055", "#eb2424", "#26b0e0", "#99c457", "#b59e7d", "#fac20a", "#1470b5", "#703b99", "#7f7f7f"]
pats = list(D["patients"])


def get(p):
    o = np.array(D["patients"][p]["observed_W1_W3"]); o /= o.sum(1, keepdims=True)
    q = np.array(D["patients"][p]["loo_pred_W2_W3"]); q /= q.sum(1, keepdims=True)
    return o, q


# ---- figure 1: per patient
fig, axs = plt.subplots(2, 5, figsize=(6.3, 3.9), sharey=True)
x = np.arange(10)
for ax, p in zip(axs.flat, pats):
    o, q = get(p)
    for w, a in enumerate((0.35, 0.6, 0.9)):
        ax.bar(x + (w - 1) * 0.27, o[w], width=0.27, color=COL, alpha=a, linewidth=0)
    ax.plot(x, q[0], "x", color="k", ms=3.5, mew=0.9)
    ax.plot(x + 0.27, q[1], "+", color="k", ms=4.5, mew=0.9)
    rm = np.sqrt(np.mean((q - o[1:]) ** 2))
    ax.set_title(f"patient {p}", fontsize=9, pad=2)
    ax.text(0.97, 0.95, f"RMSE {rm:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
    ax.set_xticks(x, SHORT, rotation=90, fontsize=7)
    ax.set_ylim(0, 1); ax.set_yticks([0, 0.5, 1])
    ax.tick_params(axis="x", length=0, pad=1.5)
for ax in axs[:, 0]:
    ax.set_ylabel("guild fraction")
hand = [Patch(color="0.55", alpha=a, label=l) for a, l in ((0.35, "week 1 (start)"), (0.6, "week 2"), (0.9, "week 3"))]
hand += [Line2D([], [], ls="", marker="x", color="k", ms=4, label="predicted week 2"),
         Line2D([], [], ls="", marker="+", color="k", ms=5, label="predicted week 3")]
fig.legend(handles=hand, ncol=5, fontsize=7.5, frameon=False, loc="lower center", handletextpad=0.4, columnspacing=1.0)
fig.tight_layout(rect=(0, 0.05, 1, 1), w_pad=0.3, h_pad=0.6)
fig.savefig(HERE / "dieckow_loo_patients.pdf")
fig.savefig(HERE / "dieckow_loo_patients.png", dpi=170)

# ---- figure 2: predicted vs observed
fig, ax = plt.subplots(figsize=(3.4, 3.2))
allo, allq = [], []
for p in pats:
    o, q = get(p)
    for w, mk in ((0, "o"), (1, "s")):
        ax.scatter(o[w + 1], q[w], c=COL, marker=mk, s=14, edgecolors="none", alpha=0.85)
    allo.append(o[1:]); allq.append(q)
allo = np.concatenate(allo).ravel(); allq = np.concatenate(allq).ravel()
ax.plot([0, 1], [0, 1], color="0.4", lw=0.8, ls="--")
r2 = 1 - np.sum((allq - allo) ** 2) / np.sum((allo - allo.mean()) ** 2)
ax.text(0.04, 0.96, f"RMSE {np.sqrt(np.mean((allq - allo) ** 2)):.3f}\n$R^2$ {r2:.2f}", transform=ax.transAxes,
        va="top", fontsize=8)
ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="observed fraction", ylabel="predicted fraction (left-out patient)",
       aspect="equal")
hand = [Line2D([], [], ls="", marker="o", color="0.4", ms=4, label="week 2"),
        Line2D([], [], ls="", marker="s", color="0.4", ms=4, label="week 3")]
hand += [Line2D([], [], ls="", marker="o", color=c, ms=4, label=s) for c, s in zip(COL, SHORT)]
ax.legend(handles=hand, fontsize=6.5, frameon=False, loc="lower right", ncol=2, handletextpad=0.2, columnspacing=0.6)
fig.tight_layout()
fig.savefig(HERE / "dieckow_loo_scatter.pdf")
fig.savefig(HERE / "dieckow_loo_scatter.png", dpi=170)
print("pooled RMSE", round(float(np.sqrt(np.mean((allq - allo) ** 2))), 4), "R2", round(float(r2), 3),
      "mean per-patient RMSE", round(float(np.mean([D["patients"][p]["rmse_printed"] for p in pats])), 4))

"""fig_p3_posterior.py -- figure of uq/p3_posterior.py (CH, 30 random posterior particles and the MAP, 2x6x8):
growth, band stress and cervical ratio without and with the homeostatic-pressure law.

    python3 uq/fig_p3_posterior.py -> assets/fig_p3_posterior.png
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

r = json.load(open(ROOT / "uq" / "p3_posterior_CH_2x6x8_r30.json"))
keys = [("core_alpha_minus_1", r"(a) core $\alpha-1$"), ("band_vM_over_E", "(b) band von Mises$/E$"),
        ("ratio", "(c) band / core")]
fig, ax = plt.subplots(1, 3, figsize=(6.3, 2.3))
for a, (k, ttl) in zip(ax, keys):
    for j, (pr, col, lab) in enumerate(((0.0, "#c0392b", "no law"), (0.05, "#1e50a0", "$p_h/E=0.05$"))):
        v = np.array([x[k] for x in r if x["p_ref"] == pr and x["particle"] != "MAP"])
        m = [x[k] for x in r if x["p_ref"] == pr and x["particle"] == "MAP"][0]
        parts = a.violinplot([v], positions=[j], widths=0.7, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(col); b.set_alpha(0.45)
        a.plot(np.full(v.size, j) + np.random.default_rng(1).uniform(-0.12, 0.12, v.size), v, ".", color=col, ms=2.5)
        a.plot(j, m, "*", color="k", ms=6, label="MAP" if j == 0 else None)
    a.set_xticks([0, 1], ["no law", "$p_h/E=0.05$"])
    a.set_title(ttl, fontsize=9)
ax[2].axhline(1.0, color="0.6", lw=0.8)
ax[0].legend(frameon=False, loc="upper right", fontsize=7.5)
fig.tight_layout(w_pad=1.0)
fig.savefig(ROOT / "assets" / "fig_p3_posterior.png", dpi=250)
print("wrote assets/fig_p3_posterior.png")

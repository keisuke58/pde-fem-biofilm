"""fig_threshold.py -- the cervical concentration ratio against p_h/E (a) and against p_h/p_0 (b) for four cases
that change p_0, the band pressure without the law (threshold_study.py). In (b) the curves collapse and cross 1 at
p_h/p_0 = 1.0-1.4. Mesh 4x12x16 (filled) and 2x6x8 (open).

    python3 p3/uel/fig_threshold.py   -> assets/fig_threshold.pdf, .png
Data: p3/uel/results/threshold_summary.json.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

rows = json.load(open(ROOT / "p3/uel/results/threshold_summary.json"))
CASES = [("main", "main case", "#1e50a0", "s"), ("t05", "half the growth time", "#2a9d8f", "D"),
         ("b05", "half the crown bulge", "#e09f3e", "^"), ("bf", "bottom face clamped", "#c0392b", "o")]

fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
for case, lab, col, mk in CASES:
    for nr, filled in ((4, True), (2, False)):
        rs = sorted([r for r in rows if r["case"] == case and r["nr"] == nr and r["m"]], key=lambda r: r["m"])
        if not rs:
            continue
        kw = dict(color=col, marker=mk, ms=5, lw=1.4 if filled else 0.9, ls="-" if filled else ":",
                  mfc=col if filled else "none", label=lab if filled else None)
        ax[0].plot([r["p_ref"] for r in rs], [r["ratio"] for r in rs], **kw)
        ax[1].plot([r["m"] for r in rs], [r["ratio"] for r in rs], **kw)
for a in ax:
    a.axhline(1.0, color="0.6", lw=0.8)
    a.set_ylabel("band / core von Mises")
ax[0].set(xlabel="$p_h/E$", title="(a) against $p_h/E$")
ax[1].set(xlabel="$p_h/p_0$  ($p_0$: band pressure without the law)", title="(b) against $p_h/p_0$")
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, fontsize=7.5, frameon=False, loc="lower center", ncol=4)
fig.tight_layout(rect=(0, 0.08, 1, 1))
out = ROOT / "assets" / "fig_threshold"
fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"), dpi=170)
print("wrote", out.with_suffix(".png"))

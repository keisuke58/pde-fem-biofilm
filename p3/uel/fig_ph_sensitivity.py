"""fig_ph_sensitivity.py -- sensitivity of the cervical result to the homeostatic pressure p_h/E, which the
literature does not fix (LITERATURE_P3.md). Concentration ratio and band stress against p_h/E = 0.02-0.2, with the
law-off value as a horizontal line, for the outer face free and fixed (nutrient from the crown, bulge 1, bottom
u_z = 0, explicit gate, mesh 4x12x16; Gauss-point means).

    python3 p3/uel/fig_ph_sensitivity.py   -> assets/fig_ph_sensitivity.pdf, .png
Data: p3/uel/results/implant_cervical_summary.json (cervical_study.py).
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

rows = json.load(open(ROOT / "p3/uel/results/implant_cervical_summary.json"))
NR = 4
PREFS = [0.02, 0.05, 0.1, 0.2]


def pick(outer, pref):
    for r in rows:
        if (r["nr"] == NR and r["nut"] == "outer" and r["bulge"] == 1.0 and r["bottom"] == "uz" and r["outer"] == outer
                and abs(r["p_ref"] - pref) < 1e-9 and r["gate"] == "explicit"):
            return r
    return None


fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
table = []
for outer, col, ls, mk in [("free", "#1e50a0", "-", "s"), ("fixed", "#c0392b", "--", "o")]:
    pts = [(p, pick(outer, p)) for p in PREFS]
    pts = [(p, r) for p, r in pts if r]
    off = pick(outer, 0.0)
    lab = f"outer face {outer}"
    ax[0].plot([p for p, _ in pts], [r["ratio"] for _, r in pts], color=col, ls=ls, marker=mk, lw=1.4, ms=5,
               label=lab + ", law on")
    ax[1].plot([p for p, _ in pts], [r["band_vM_over_E"] for _, r in pts], color=col, ls=ls, marker=mk, lw=1.4, ms=5)
    if off:
        ax[0].axhline(off["ratio"], color=col, ls=":", lw=1.0, label=lab + ", no law")
        ax[1].axhline(off["band_vM_over_E"], color=col, ls=":", lw=1.0)
    table.append((outer, off and (round(off["ratio"], 2), round(off["band_vM_over_E"], 4)),
                  [(p, round(r["ratio"], 2), round(r["band_vM_over_E"], 4)) for p, r in pts]))
for a in ax:
    a.set_xscale("log")
    a.set_xticks(PREFS, [str(p) for p in PREFS])
    a.minorticks_off()
    a.set_xlabel("$p_h/E$")
ax[0].axhline(1.0, color="0.6", lw=0.8)
ax[0].set(ylabel="band / core von Mises", title="(a) cervical concentration")
ax[1].set(ylabel="band mean von Mises$/E$", title="(b) stress in the cervical band", yscale="log")
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, fontsize=7.5, frameon=False, loc="lower center", ncol=4)
fig.tight_layout(rect=(0, 0.07, 1, 1))
out = ROOT / "assets" / "fig_ph_sensitivity"
fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"), dpi=170)
print("wrote", out.with_suffix(".png"))
for t in table:
    print(*t)

"""fig_cervical.py -- the cervical-line result of the tooth-like sector for the P3 paper: concentration ratio
(band mean vM/E over core mean vM/E) and band stress with and without the homeostatic-pressure law, over the mesh
levels, for the bottom face u_z = 0 and the outer face free or fixed (nutrient from the crown, bulge 1, seed
0.125 mm, E = 1, p_h/E = 0.05; Gauss-point means; explicit gate).

    python3 p3/uel/fig_cervical.py        -> assets/fig_cervical.pdf, .png
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


def pick(nr, outer, pref):
    for r in rows:
        if (r["nr"] == nr and r["nut"] == "outer" and r["bulge"] == 1.0 and r["bottom"] == "uz" and r["outer"] == outer
                and abs(r["p_ref"] - pref) < 1e-9 and r["gate"] == "explicit"):
            return r
    return None


levels = [(2, "2×6×8"), (4, "4×12×16"), (6, "6×18×24")]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
style = {("free", 0.0): ("#c0392b", "-", "o"), ("free", 0.05): ("#1e50a0", "-", "s"),
         ("fixed", 0.0): ("#c0392b", "--", "o"), ("fixed", 0.05): ("#1e50a0", "--", "s")}
for (outer, pref), (col, ls, mk) in style.items():
    pts = [(i, pick(nr, outer, pref)) for i, (nr, _) in enumerate(levels)]
    pts = [(i, r) for i, r in pts if r]
    lab = f"outer face {'free' if outer == 'free' else 'fixed'}, " + ("no law" if pref == 0 else "$p_h/E=0.05$")
    ax[0].plot([i for i, _ in pts], [r["ratio"] for _, r in pts], color=col, ls=ls, marker=mk, label=lab, lw=1.4, ms=5)
    ax[1].plot([i for i, _ in pts], [r["band_vM_over_E"] for _, r in pts], color=col, ls=ls, marker=mk, lw=1.4, ms=5)
for a in ax:
    a.set_xticks(range(3), [l for _, l in levels])
    a.set_xlabel("mesh (r × θ × z elements)")
ax[0].axhline(1.0, color="0.6", lw=0.8)
ax[0].set(ylabel="band / core von Mises", title="(a) cervical concentration")
ax[1].set(ylabel="band mean von Mises$/E$", title="(b) stress in the cervical band", yscale="log")
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, fontsize=7.5, frameon=False, loc="lower center", ncol=4)
fig.tight_layout(rect=(0, 0.07, 1, 1))
out = ROOT / "assets" / "fig_cervical"
fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"), dpi=170)
print("wrote", out.with_suffix(".png"))
for (outer, pref) in style:
    print(outer, pref, [(nr, round(pick(nr, outer, pref)["ratio"], 2), round(pick(nr, outer, pref)["band_vM_over_E"], 4))
                       for nr, _ in levels if pick(nr, outer, pref)])

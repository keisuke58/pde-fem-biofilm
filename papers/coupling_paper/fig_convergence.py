"""fig_convergence.py -- Figure 2 of the coupling paper: seed alpha - 1 (a) and seed von Mises stress (b) at T* = 1 against
the mesh (h = 2 mm / n), one line per scheme and code, from the GCI data of ansys_usermat/apdl/gci_a52.py
(results/2026-10-living/gci_a52.json; A-52, A-57 to A-64). Separate example with antibiotic. Grey band: the finest
grids of the two second-order schemes (Python TVD 128^3, Abaqus SUPG 48^3).

    python papers/coupling_paper/fig_convergence.py [path/to/gci_a52.json] -> papers/coupling_paper/fig_convergence.pdf
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply()
src = Path(sys.argv[1]) if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else ROOT / "ansys_usermat/apdl/results/2026-10-living/gci_a52.json"
res = json.load(open(src, encoding="utf-8"))
STYLE = {  # label, colour, marker, linestyle
    "Python (r h/4)": ("Python, constant $r h/4$", "C0", "o", "--"),
    "Abaqus": ("Abaqus, constant $r h/4$", "C1", "s", "--"),
    "ANSYS": ("ANSYS, constant $r h/4$", "C2", "^", "--"),
    "Python (TVD)": ("Python, TVD", "C0", "o", "-"),
    "Python (LW + sensor)": ("Python, LW + sensor", "C0", "D", ":"),
    "Abaqus (SUPG)": ("Abaqus, SUPG", "C1", "s", "-"),
    "ANSYS (living10)": ("ANSYS, LW + sensor", "C2", "^", ":"),
    "ANSYS (upwind)": ("ANSYS, upwind MLS gradients", "C2", "v", "-"),
}
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.9), constrained_layout=True)
for ax, q, lab, scale in ((axs[0], "seed_am1", r"seed mean of $\alpha-1$ [$10^{-4}$]", 1e4),
                          (axs[1], "seed_vm_Pa", r"seed mean of $\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]", 1e4)):
    for code, r in res.items():
        vals = r.get("values", {})
        ns = sorted(int(n) for n in vals)
        ns = [n for n in ns if n >= 16]
        if not ns or code not in STYLE:
            continue
        name, col, mk, ls = STYLE[code]
        ax.plot([2.0 / n for n in ns], [vals[str(n)][q] * scale for n in ns], marker=mk, color=col, ls=ls, ms=4, lw=1.0, label=name)
    # the converged value: the band of the finest grids of the two second-order schemes (Python TVD, Abaqus SUPG)
    fin = [r["values"][str(max(int(n) for n in r["values"]))][q] * scale for c_, r in res.items()
           if c_ in ("Python (TVD)", "Abaqus (SUPG)") and r.get("values")]
    if fin:
        ax.axhspan(min(fin), max(fin), color="0.85", lw=0, zorder=0)
    ax.set_xlabel("$h$ [mm]")
    ax.set_ylabel(lab)
    ax.set_xscale("log")
    ax.set_xticks([2 / 16, 2 / 32, 2 / 64, 2 / 128])
    ax.set_xticklabels(["$2/16$", "$2/32$", "$2/64$", "$2/128$"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.invert_xaxis()
    ax.grid(True, lw=0.3, alpha=0.5)
axs[0].text(-0.2, 1.02, "(a)", transform=axs[0].transAxes)
axs[1].text(-0.2, 1.02, "(b)", transform=axs[1].transAxes)
fig.legend(*axs[1].get_legend_handles_labels(), fontsize=7.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.0), frameon=False)
out = HERE / "fig_convergence.pdf"
fig.savefig(out, bbox_inches="tight")
fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=200)
print("wrote", out)

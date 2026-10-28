#!/usr/bin/env python3
"""fig_mms.py -- Figure 4 of the coupling paper: manufactured solution (abaqus_living/mms.py, A-50) of the phi and c
equations with the front term, L2 error of phi (a) and of c (b) at T* = 0.05 against the cell size h = 2 mm / n
(n = 8, 16, 32) for the Python finite-volume code and the Abaqus Galerkin element; dotted: slope 2.

Data: Python from results/2026-10-3d/coupled_living_consumption_{,nf16_,nf8_}front1_dtc0.005_phimin0.001_v2_tend0.05
_a100_gab1_gc1_tab0.333333_mms_artdiff_wgmirror_mixed.json (key runs.mixed.mms, `--mms --tend 0.05 --front 1
--artdiff --wall-grad mirror`, outer step 0.005); Abaqus from the table of report A-50 in
ansys_usermat/apdl/IKMHIWI03_TO_CLOUD.ja.md (`make_living_inp.py --mms`, PROPS(19) = 1, increment 2.5e-4; the Abaqus
job files stay on IKMHIWI03). The observed orders are printed.

    python papers/coupling_paper/fig_mms.py -> papers/coupling_paper/fig_mms.pdf
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

RES = ROOT / "ansys_usermat/apdl/results/2026-10-3d"
N = (8, 16, 32)
H = 2.0 / np.array(N)
ABAQUS = {"phi": [3.26e-3, 1.02e-3, 2.90e-4], "c": [2.68e-3, 6.69e-4, 1.62e-4]}   # A-50, L2 errors at T* = 0.05


def python_errors():
    phi, c = [], []
    for n in N:
        pre = "" if n == 32 else f"nf{n}_"
        f = RES / f"coupled_living_consumption_{pre}front1_dtc0.005_phimin0.001_v2_tend0.05_a100_gab1_gc1_tab0.333333_mms_artdiff_wgmirror_mixed.json"
        m = json.load(open(f))["runs"]["mixed"]["mms"]
        assert abs(m["h"] - 2.0 / n) < 1e-12 and abs(m["T"] - 0.05) < 1e-12, (f, m)
        phi.append(m["l2_phi"]); c.append(m["l2_c"])
    return {"phi": phi, "c": c}


def orders(e):
    e = np.asarray(e)
    return np.log(e[:-1] / e[1:]) / np.log(2.0)


def main():
    figstyle.apply(size=8.5)
    py = python_errors()
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.6), constrained_layout=True)
    for ax, q, lab, t in zip(axs, ("phi", "c"), (r"$L_2$ error of $\phi_{3D}$", r"$L_2$ error of $c$"), "ab"):
        for code, e, col, mk in (("Python, finite volumes", py[q], "C0", "o"), ("Abaqus, Galerkin", ABAQUS[q], "C1", "s")):
            ax.loglog(H, e, color=col, marker=mk, ls="-", lw=1, ms=4, label=code)
            print(f"{code:24s} {q}: errors {', '.join(f'{v:.2e}' for v in e)}; observed orders {np.round(orders(e), 2)}")
        ref = 0.6 * min(py[q][-1], ABAQUS[q][-1]) * (H / H[-1]) ** 2
        ax.loglog(H, ref, color="0.4", ls=":", lw=0.9, label="slope 2")
        ax.set_xlabel(r"cell size $h$ [mm]")
        ax.set_ylabel(lab)
        ax.set_xticks(H)
        ax.set_xticklabels([f"$2/{n}$" for n in N])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlim(H[-1] / 1.4, H[0] * 1.4)
        ax.grid(True, which="both", lw=0.3, alpha=0.4)
        ax.text(-0.18, 1.02, f"({t})", transform=ax.transAxes, va="bottom", ha="left")
    axs[0].legend(fontsize=7, loc="lower right", frameon=False)
    out = HERE / "fig_mms.pdf"
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()

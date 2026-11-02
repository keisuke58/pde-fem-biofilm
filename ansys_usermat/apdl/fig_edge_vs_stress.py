#!/usr/bin/env python3
"""Does the seed stress follow the sharpness of the biofilm edge? (Keio continuation / paper discussion, not thesis.)

For every saved 64^3 run of the cube examples (figdata/coupled_living_*nf64*.npz, T* = 1), two measures of the edge of
the seed are compared with the seed von Mises stress of the same run:
  contrast of phi_3D: mean phi_3D in the seed minus mean phi_3D in the shell of one element (0.25 mm) around it;
  contrast of growth: the same for alpha - 1;
  stiffness-weighted growth gradient: mean of (phi_3D^2 + f) |grad alpha| over the seed and the shell (the stiffness
  E = E_0 (phi_3D^2 + f) of the runs), the measure of incompatible growth used in panel (a).
The stress is the seed mean on the 16^3 stress mesh (seed_vm_Pa of the run's JSON, recomputed from vm16 here).

    python ansys_usermat/apdl/fig_edge_vs_stress.py   -> assets/fig_edge_vs_stress.png/.pdf, results/2026-10-3d/edge_vs_stress.json
"""
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import binary_dilation  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
import wp3_reference as W  # noqa: E402

CELLS8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
SEED = W.seed_mask(64, CELLS8)
SHELL = binary_dilation(SEED, np.ones((3, 3, 3), bool), iterations=8) & ~SEED   # one element = 8 cells at 64^3
SEED16 = W.seed_mask(16, CELLS8)


def group(name):
    law = "living" if "growlivE" in name else "published" if "growliv_" not in name else None
    if law is None:
        return None
    if "front0.1" in name:
        g = "slow front (r/10)"
    elif "csine" in name or "cconst" in name:
        g = "nutrient sine / held 0.5"
    elif re.search(r"_a(?!100_)[0-9.]+_gab1_gc", name):
        g = "pilot points (alpha*_0, g l^2/d)"
    else:
        g = "cube examples"
    if "abpulse" in name or "absine" in name:
        g = "antibiotic schedules"
    if any(s in name for s in ("eps2", "dtc", "snaps")):
        return None                                # step and regularisation checks, not separate cases
    return law, g


def main():
    figstyle.apply()
    rows = []
    for f in sorted((ROOT / "figdata").glob("coupled_living_*_nf64_*.npz")):
        gr = group(f.name)
        if gr is None:
            continue
        z = np.load(f)
        if "phi_1.0" not in z.files or "vm16_1.0" not in z.files or "alpha_1.0" not in z.files:
            continue
        phi, am1, vm = z["phi_1.0"], z["alpha_1.0"] - 1.0, z["vm16_1.0"]
        g = np.sqrt(sum(d ** 2 for d in np.gradient(am1, 2.0 / 64)))
        band = SEED | SHELL
        rows.append(dict(name=f.stem.replace("coupled_living_", ""), law=gr[0], group=gr[1],
                         mixed="mixed" in f.name, antibiotic="_ab_nf" in f.name,
                         dphi=float(phi[SEED].mean() - phi[SHELL].mean()),
                         egrad_1e4=float(1e4 * ((phi ** 2 + 1e-3) * g)[band].mean()),
                         dgrowth_1e4=float(1e4 * (am1[SEED].mean() - am1[SHELL].mean())),
                         growth_1e4=float(1e4 * am1[SEED].mean()),
                         vm_1e4Pa=float(1e4 * np.nan_to_num(vm)[SEED16].mean())))
    out = HERE / "results/2026-10-3d/edge_vs_stress.json"
    for law in ("published", "living"):
        r = [x for x in rows if x["law"] == law]
        if len(r) > 2:
            for k in ("egrad_1e4", "dphi", "dgrowth_1e4", "growth_1e4"):
                print(law, k, "corr with vM:", round(float(np.corrcoef([x[k] for x in r], [x["vm_1e4Pa"] for x in r])[0, 1]), 3), len(r))
    out.write_text(json.dumps(rows, indent=1) + "\n")

    marks = {"cube examples": "o", "slow front (r/10)": "s", "nutrient sine / held 0.5": "D",
             "pilot points (alpha*_0, g l^2/d)": ".", "antibiotic schedules": "^"}
    cols = {"published": "#1f4e8c", "living": "#c92a2a"}
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.3), constrained_layout=True)
    for x in rows:
        kw = dict(marker=marks[x["group"]], ms=3 if x["group"].startswith("pilot") else 4.5, ls="",
                  mfc=cols[x["law"]] if not x["mixed"] else "white", mec=cols[x["law"]], mew=0.8)
        axs[0].plot(x["egrad_1e4"], x["vm_1e4Pa"], **kw)
        axs[1].plot(x["dgrowth_1e4"], x["vm_1e4Pa"], **kw)
        axs[2].plot(x["growth_1e4"], x["vm_1e4Pa"], **kw)
    axs[0].set_xlabel(r"$(\phi_{3D}^2+f)\,|\nabla\alpha|$ [$10^{-4}$/mm]"); axs[0].set_ylabel(r"seed $\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")
    axs[1].set_xlabel(r"$\alpha-1$: seed $-$ shell [$10^{-4}$]")
    axs[2].set_xlabel(r"$\alpha-1$ in the seed [$10^{-4}$]")
    for ax, t in zip(axs, "abc"):
        ax.text(0.04, 0.96, f"({t})", transform=ax.transAxes, va="top")
    h = [plt.Line2D([], [], marker=m, ls="", color="k", mfc="k" if g != "mixed" else "white", ms=4)
         for g, m in marks.items()]
    h += [plt.Line2D([], [], marker="o", ls="", color=c, ms=4) for c in cols.values()]
    h += [plt.Line2D([], [], marker="o", ls="", color="k", mfc="white", ms=4)]
    fig.legend(h, list(marks) + ["published growth law", "living-biomass growth law", "mixed start"],
               loc="outside lower center", ncol=4, fontsize=6.5, frameon=False)
    o = ROOT / "assets" / "fig_edge_vs_stress.png"
    fig.savefig(o, dpi=600); fig.savefig(o.with_suffix(".pdf"))
    print("wrote", o, len(rows), "runs")


if __name__ == "__main__":
    main()

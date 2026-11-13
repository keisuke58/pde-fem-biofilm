"""fig_mech_K_sweep.py -- K sweep of the mechanical term -K phi_3D (I:C_e - 3) of Klempt 2024 Eq. 30 (fifa, part 12).

Seed means at T* = 1 against K = 0 on the same grid, 64^3 / 128^3 / 192^3, superbee, --transient-c 1e8, outer step 0.025:
(a) change of alpha - 1, (b) change of the von Mises stress; separate and mixed (equal shares), both with antibiotic.
Runs without a JSON diverged (separate: K = 3e5 and 1e6; mixed: 1e6) and are marked on the axis.

    python ansys_usermat/apdl/fig_mech_K_sweep.py -> assets/fig_mech_K_sweep.pdf/.png
"""
import glob
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402
figstyle.apply()

R = HERE / "results" / "2026-10-3d"
EX = {"separate": ("separate_ab", "_inhliv_a100", "#1c7ed6"), "mixed": ("mixed_ab", "_inhliv_mixeq_a100", "#e8590c")}
MK = {64: "o", 128: "s", 192: "^"}


def load(ex):
    key, mid, _ = EX[ex]
    out = {}
    pat = str(R / f"coupled_living_consumption_nf*_front1_phimin0.001_v2_fortran_ansysbc{mid}_gab1_gc1_tab0.333333*"
                  f"transc1e+08_tvd_superbee_wgmirror_fast_{key}.json")
    for f in glob.glob(pat):
        m = re.search(r"nf(\d+)_front1_phimin0\.001_v2_fortran_ansysbc" + re.escape(mid) +
                      r"_gab1_gc1_tab0\.333333(?:_mech(\d+))?_transc", f)
        if not m:
            continue
        run = json.load(open(f))["runs"][key]
        out[(int(m.group(1)), float(m.group(2) or 0))] = (run["seed_am1"], run["seed_vm_Pa"])
    return out


fig, ax = plt.subplots(1, 2, figsize=(6.3, 2.6), sharex=True)
for ex, (_, _, col) in EX.items():
    d = load(ex)
    for nf in (64, 128, 192):
        ks = sorted(k for (n, k) in d if n == nf and k > 0)
        if (nf, 0.0) not in d or not ks:
            continue
        a0, v0 = d[(nf, 0.0)]
        for j, a in enumerate(ax):
            y = [100 * (d[(nf, k)][j] / (a0, v0)[j] - 1) for k in ks]
            a.plot(ks, y, color=col, marker=MK[nf], ms=3.5, lw=1.0, mfc="w" if nf != 192 else col,
                   ls={64: ":", 128: "--", 192: "-"}[nf], label=f"{ex}, ${nf}^3$")
    print(ex, {k: (f"{100*(v[0]/d[(k[0],0.0)][0]-1):+.2f}", f"{100*(v[1]/d[(k[0],0.0)][1]-1):+.2f}") for k, v in sorted(d.items()) if k[1] > 0})
for a, lab, t in zip(ax, (r"change of $\alpha-1$ [%]", "change of von Mises [%]"), ("(a)", "(b)")):
    a.set_xscale("log"); a.axhline(0, color="0.6", lw=0.6); a.axhspan(-1, 1, color="0.9", lw=0)
    a.set_xlabel(r"$K$ [$1/T^*$]", fontsize=8); a.set_ylabel(lab, fontsize=8); a.tick_params(labelsize=7)
    a.set_title(t, fontsize=8.5, loc="left"); a.set_xlim(6e2, 1.6e6)
    a.axvline(3e5, color=EX["separate"][2], lw=0.6, ls=(0, (1, 2))); a.axvline(1e6, color="k", lw=0.6, ls=(0, (1, 2)))
ax[1].text(3e5, 3.5, "separate\ndiverges", fontsize=6, ha="right", va="bottom", color=EX["separate"][2])
ax[1].text(1e6, 3.5, "both\ndiverge", fontsize=6, ha="right", va="bottom")
ax[0].legend(frameon=False, fontsize=6.3, loc="lower left")
fig.tight_layout()
A = HERE.parents[1] / "assets"
fig.savefig(A / "fig_mech_K_sweep.pdf"); fig.savefig(A / "fig_mech_K_sweep.png", dpi=400)
print("wrote", A / "fig_mech_K_sweep.png")

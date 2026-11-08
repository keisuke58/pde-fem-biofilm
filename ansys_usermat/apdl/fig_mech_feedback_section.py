"""fig_mech_feedback_section.py -- where the mechanical term of Klempt 2024 Eq. 30 acts (Keio continuation, not the thesis).

--mech-feedback K adds -K phi (I:C_e - 3) = -2 K phi tr(eps_e) to d phi_3D/dt: where the biofilm is compressed
(tr eps_e < 0, pressure p > 0) phi_3D grows faster, where it is stretched it grows slower. The figure compares a run with
K and the run without, section z = 0 at T*: (a) mean stress -p of the run without K (16^3 stress mesh, tension positive; the snapshots store p16 with tension positive),
(b) the change of phi_3D, (c) the change of alpha - 1, (d) the von Mises stress with and without K along x at y = 0.

    python ansys_usermat/apdl/fig_mech_feedback_section.py [--nf 32] [--K 100000] [--t 1.0]
needs the snapshots of the two runs (--snap-every 0.5 --snap-stress), see the header of the K sweep in CLOUD_TO_FIFA part 12.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402
figstyle.apply()
ap = argparse.ArgumentParser()
ap.add_argument("--nf", type=int, default=32); ap.add_argument("--obst", action="store_true"); ap.add_argument("--K", default="100000"); ap.add_argument("--t", default="1.0"); ap.add_argument("--exp", action="store_true"); ap.add_argument("--m1", default=None); ap.add_argument("--clip", action="store_true")
A = ap.parse_args()
nf = "" if A.nf == 32 else f"_nf{A.nf}"
base = f"coupled_living_separate_ab{nf}_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_snapsvm{{m}}{{c}}_transc1e+08_tvd_superbee_wgmirror_fast{{o}}_separate_ab.npz"
O = "_obstcolumns" if A.obst else ""
C = "_clip" if A.clip else ""
Z0 = np.load(ROOT / "figdata" / base.format(m="", o=O, c=C))
Z1 = np.load(ROOT / "figdata" / base.format(m=(A.m1 if A.m1 else f"_mech{A.K}" + ("exp" if A.exp else "")), o=O, c=C))
t = A.t
n = Z0[f"phi_{t}"].shape[0]; k = n // 2
nm = Z0[f"p16_{t}"].shape[0]; km = nm // 2
sec = lambda a, kk: 0.5 * (a[:, :, kk - 1] + a[:, :, kk])          # z = 0 lies between the two middle layers
p0 = 1e4 * sec(Z0[f"p16_{t}"], km)
dphi = sec(Z1[f"phi_{t}"] - Z0[f"phi_{t}"], k)
dal = 1e4 * sec(Z1[f"alpha_{t}"] - Z0[f"alpha_{t}"], k)
ext = (-1, 1, -1, 1)
fig, ax = plt.subplots(1, 4, figsize=(6.3, 1.9))
for a_, f, cm, lab in ((ax[0], p0, "RdBu_r", r"(a) mean stress, no $K$ [$10^{-4}$ Pa]"),
                       (ax[1], dphi, "PuOr", r"(b) $\Delta\phi_{3D}$"),
                       (ax[2], dal, "PuOr", r"(c) $\Delta(\alpha-1)$ [$10^{-4}$]")):
    v = float(np.nanmax(np.abs(f))) or 1.0
    im = a_.imshow(f.T, origin="lower", extent=ext, cmap=cm, vmin=-v, vmax=v)
    a_.set_title(lab, fontsize=7.5, loc="left"); a_.set_xlabel("$x$ [mm]", fontsize=7); a_.tick_params(labelsize=6)
    cb = fig.colorbar(im, ax=a_, fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=6)
ax[0].set_ylabel("$y$ [mm]", fontsize=7)
xs = np.linspace(-1, 1, nm, endpoint=False) + 1 / nm
for Z, lab, ls in ((Z0, "without", "-"), (Z1, f"K = {float(A.K):g}", "--")):
    ax[3].plot(xs, 1e4 * sec(Z[f"vm16_{t}"], km)[:, km], ls, lw=1.1, label=lab)
ax[3].set_title(r"(d) $\sigma_\mathrm{vM}$ at $y=0$ [$10^{-4}$ Pa]", fontsize=7.5, loc="left")
ax[3].set_xlabel("$x$ [mm]", fontsize=7); ax[3].tick_params(labelsize=6); ax[3].legend(frameon=False, fontsize=6)
fig.tight_layout()
out = ROOT / "assets" / f"fig_mech_feedback_section{nf}{'_bars' if A.obst else ''}_K{A.K}{'exp' if A.exp else ''}{'_clip' if A.clip else ''}"
fig.savefig(out.with_suffix(".png"), dpi=300); fig.savefig(out.with_suffix(".pdf"))
print("wrote", out, "dphi range", float(dphi.min()), float(dphi.max()), "p range", float(p0.min()), float(p0.max()))

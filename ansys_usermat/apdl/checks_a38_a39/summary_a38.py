"""T* = 0.25 ... 1: Python 8^3 variants vs ANSYS u3a_8_sep_ab_dt0005 / u6_8_sep_ab_dt0005 (separate_ab)."""
import glob, json
import numpy as np

R3 = "ansys_usermat/apdl/results/2026-10-3d/"
RL = "ansys_usermat/apdl/results/2026-10-living/"
base = "coupled_living_consumption_nf8_front1{dtc}_phimin0.001_v2_ansysbc_startpos_a100_gab1_gc1_tab0.333333_snaps_artdiff{wg}_separate_ab"
for run in ("u3a_8_sep_ab_dt0005", "u3a_8_sep_ab_dt00025", "u3a_8_sep_ab_dt000125", "u6t_8_sep_ab_dt0005", "u6s_8_sep_ab_dt0005"):
    fa = json.load(open(RL + f"fields_{run}.json"))
    seed = np.isin(np.array(fa["elem"]), fa["seed"])
    m = [np.mean(p) for p in fa["phi"]]
    k = -1
    print(f"| ANSYS {run} | " + " | ".join(f"{x:.3f}" for x in m) + f" | {np.array(fa['L'][k])[seed].mean():.3f} | "
          f"{np.array(fa['nut'][k])[seed].mean():.3f} | {np.array(fa['alpha'][k])[seed].mean():.2e} |")
for dtc, wg in (("", ""), ("", "_wgmirror"), ("", "_wgzero"), ("_dtc0.005", "_wgmirror"), ("_dtc0.005", "_wgzero")):
    b = base.format(dtc=dtc, wg=wg)
    r = json.load(open(R3 + b + ".json"))["runs"]["separate_ab"]
    z = np.load(glob.glob("figdata/coupled_living_separate_ab" + b[len("coupled_living_consumption"):] + ".npz")[0])
    m = [float(z[f"phi_{T}"].mean()) for T in (0.25, 0.5, 0.75, 1.0)]
    print(f"| Python{dtc or ' dtc 0.025'} {wg[3:] or 'onesided'} | " + " | ".join(f"{x:.3f}" for x in m)
          + f" | {r['living_seed'][-1]:.3f} | {r['c_seed'][-1]:.3f} | {r['seed_am1']:.2e} |")

"""a14_summary.py -- request A-14: scheme (c) against (a) in the Abaqus cube
(make_cube_inp.py, five calibrated species, native point model) for each mesh,
from the last increment of F:\\abaqus_work\\comp_a14_cube<n>_CH_{a,c}\\*.dat as
compare_ansys.py reads it (centroidal S, SDV84 = alpha - 1, the seed from
<job>.seed.txt). Next to it the same measures of the ANSYS runs
(results/2026-10-cal5/w<n>_CH_g1_s015_raw.json and w<n>_CH_cal_k0814_raw.json).

    python abaqus_composition/a14_summary.py [8 16 24]

Writes abaqus_composition/results_a14/a14_cube_CH.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from compare_ansys import last_tables  # noqa: E402

WORK = Path(r"F:\abaqus_work")
R = HERE.parent / "ansys_usermat" / "apdl" / "results" / "2026-10-cal5"


def abaqus(n, k):
    job = f"a14_cube{n}_CH_{k}"
    dat = WORK / f"comp_{job}" / f"{job}.dat"
    t = last_tables(dat)
    S = next(x for x in t if x["head"][:1] == ["S11"])["rows"]
    V = next(x for x in t if "SDV84" in x["head"])
    iv = V["head"].index("SDV84")
    seed = [int(x) for x in (dat.parent / f"{job}.seed.txt").read_text().split()]
    s = np.array([S[e] for e in seed]) * 1e6                     # MPa -> Pa
    sx, sy, sz, sxy, sxz, syz = (s[:, i] for i in range(6))
    vm = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3 * (sxy ** 2 + sxz ** 2 + syz ** 2))
    a = np.array([V["rows"][e][iv] for e in seed])
    return dict(job=job, seed_elements=len(seed), alpha_minus_1=float(a.mean()), vM=float(vm.mean()))


def ansys(n, run):
    d = json.loads((R / f"w{n}_{run}.json").read_text())
    s = set(int(e) for e in d["seed_BIOFILM1"]); a = d["all_stress"]
    m = np.array([int(e) in s for e in a["elem"]])
    return dict(run=f"w{n}_{run}", seed_elements=int(m.sum()), alpha_minus_1=float(np.array(a["alpha"], float)[m].mean()),
                vM=float(1e6 * np.array(a["seqv"], float)[m].mean()))


def main(meshes):
    rows = []
    for n in meshes:
        r = dict(n=n, abaqus_a=abaqus(n, "a"), abaqus_c=abaqus(n, "c"),
                 ansys_a=ansys(n, "CH_g1_s015_raw"), ansys_c=ansys(n, "CH_cal_k0814_raw"))
        for p in ("abaqus", "ansys"):
            r[f"{p}_diff_alpha"] = r[f"{p}_c"]["alpha_minus_1"] / r[f"{p}_a"]["alpha_minus_1"] - 1
            r[f"{p}_diff_vM"] = r[f"{p}_c"]["vM"] / r[f"{p}_a"]["vM"] - 1
        rows.append(r)
        print(f"{n:2d}: Abaqus a {r['abaqus_a']['alpha_minus_1']:.5e} c {r['abaqus_c']['alpha_minus_1']:.5e} "
              f"({100 * r['abaqus_diff_alpha']:+.2f} %), vM a {r['abaqus_a']['vM']:.4e} c {r['abaqus_c']['vM']:.4e} "
              f"({100 * r['abaqus_diff_vM']:+.2f} %) | ANSYS alpha {100 * r['ansys_diff_alpha']:+.3f} %, "
              f"vM {100 * r['ansys_diff_vM']:+.2f} %")
    out = dict(request="A-14", condition="CH", theta="Tmcmc202601 6e18e0b", kappa=8.135097e-3,
               measures="seed means (Abaqus: <job>.seed.txt, centroid of the last increment; ANSYS: seed_BIOFILM1); "
                        "(c)/(a) - 1 within each program", meshes=rows)
    dst = HERE / "results_a14" / "a14_cube_CH.json"
    dst.parent.mkdir(exist_ok=True)
    dst.write_text(json.dumps(out, indent=1))
    print("wrote", dst)


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]] or [8, 16, 24])

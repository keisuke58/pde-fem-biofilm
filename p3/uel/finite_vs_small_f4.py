"""finite_vs_small_f4.py -- finite against small strain on the tooth-like collar with the synthetic four-species case of
Fritsch et al. 2025 (anim_cervical_3d.synthetic_eco; no experimental data), where growth is large (alpha - 1 up to
0.4). Same deck as anim_cervical_3d.py (4x12x16, bulge 1, bottom u_z = 0, outer face free, nutrient from the outer
face, seed 0.125 mm), without the gate and with p_h/E = 0.05, each with finite=True and finite=False of
make_implant_uel_inp.run. Measures of cervical_study.measures (band, core, ratio) plus the seed mean alpha - 1.

    python p3/uel/finite_vs_small_f4.py   -> p3/uel/results/finite_vs_small_f4.json and a table on stdout
"""
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "finite_vs_small_f4.json"

WORKER = r'''
import json, sys
sys.path.insert(0, sys.argv[1])
import make_implant_uel_inp as mi
from cervical_study import measures
p_ref, finite = float(sys.argv[2]), sys.argv[3] == "1"
nodes, elems, layer, ijk = mi.sector(4, 12, 16, 1.0)
r = mi.run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=finite, seed_layers=2, p_ref=p_ref, outer="free",
           bottom="uz", nut="outer")
r.update(nr=4, nt=12, nz=16, bulge=1.0, bottom="uz", outer="free", nut="outer", p_ref=p_ref)
s = measures(r)
s.update(seed_alpha_minus_1=r["seed_alpha_minus_1"], finite=finite, wall_time_s=r["wall_time_s"])
print("RESULT " + json.dumps(s))
'''

def synthetic_eco(path):
    """the same file as anim_cervical_3d.synthetic_eco (Fritsch et al. 2025, Case II M3, species 5 absent)"""
    sys.path.insert(0, str(HERE.parents[1] / "JAXFEM"))
    from fritsch2025_cases import CASES
    c = CASES["fritsch_4sp_M3"]
    A, b = c["A"], c["b"]
    th = [0.0] * 20
    th[0], th[1], th[2], th[3], th[4] = A[0][0], A[0][1], A[1][1], b[0], b[1]
    th[5], th[6], th[7], th[8], th[9] = A[2][2], A[2][3], A[3][3], b[2], b[3]
    th[10], th[11], th[12], th[13] = A[0][2], A[0][3], A[1][2], A[1][3]
    th[14], th[15] = 1.0, 1.0
    lines = ["5", f"{c['c_star']} 0.0", " ".join(["1.0"] * 5), " ".join(f"{v:.17g}" for v in th),
             "phi_init " + " ".join(f"{v:.17g}" for v in c["phi_init"] + [0.0]), "newton 60 1e-26"]
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    scr = Path(os.environ.get("UQ_SCRATCH", "/tmp/fvs_f4")); scr.mkdir(parents=True, exist_ok=True)
    eco = scr / "uel_eco_f4.txt"; synthetic_eco(eco)
    env = dict(os.environ, BIOFILM_ECO_CASE=str(eco))
    jobs = [(p, f) for p in (0.0, 0.05) for f in (1, 0)]
    procs = [subprocess.Popen([sys.executable, "-c", WORKER, str(HERE), str(p), str(f)], env=env, cwd=scr,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True) for p, f in jobs]
    rows = []
    for pr in procs:
        out = pr.communicate()[0]
        rows.append(json.loads([l for l in out.splitlines() if l.startswith("RESULT ")][-1][7:]))
    OUT.write_text(json.dumps(rows, indent=1))
    keys = ["seed_alpha_minus_1", "core_vM_over_E", "band_vM_over_E", "ratio", "band_max_vM_over_E", "seed_max_p_over_E"]
    print("p_h/E  strain  " + "  ".join(keys))
    for r in rows:
        print(f"{r['p_ref']:<5}  {'finite' if r['finite'] else 'small '}  " + "  ".join(f"{r[k]:.4f}" for k in keys))

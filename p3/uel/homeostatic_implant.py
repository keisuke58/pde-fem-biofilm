"""homeostatic_implant.py -- the homeostatic-pressure growth law on the implant sector (make_implant_uel_inp.py):
runs the matrix mesh x bulge x outer-face condition x p_h through the emulation and tabulates it.

    BIOFILM_ECO_CASE=<dir>/uel_eco.txt python homeostatic_implant.py [--levels 1 2] [--p-ref 0 0.05] [--outer free fixed]
                                                                    [--bulge 0 1] [--ninc 60] [--gate-implicit] [--table-only]

Levels f = 1, 2, 3 are nr x nt x nz = f (2 x 6 x 8) with the seed fixed at 0.125 mm (--seed-layers f). A case whose
results/implant_<tag>.json exists is not rerun (--table-only: no runs at all). p_ref is p_h / E (PROPS(19), form 0:
p_h = p_ref, the "E" form of keio_wp2). The table gives, per case, the seed-core means (implant_convergence.core:
seed elements with 0.25 <= z <= 1.75, 15 <= theta <= 75 degrees) of alpha - 1, vM/E and p/E, the top-row corner
value (largest vM/E in the row of elements at the top face, a stiff-soft corner that the mesh study found singular),
the largest p/E over the seed and the share of seed elements whose end pressure is at or above p_h (growth stopped).
Writes results/implant_homeostatic_summary.json (the table rows) and prints the table. --gate-implicit runs the
cases with the implicit gate (tag `_gi`; column "gate": explicit or implicit).
"""
import argparse
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from implant_convergence import CORE_TH, CORE_Z  # noqa: E402

OUT = HERE / "results"


def tag_of(f, bulge, p_ref, outer, ninc=None, gate_implicit=False):
    t = f"{2 * f}x{6 * f}x{8 * f}" + ("_b1" if bulge else "") + (f"_sl{f}" if f > 1 else "")
    t += (f"_ph{p_ref:g}" + ("_gi" if gate_implicit else "") if p_ref > 0 else "") + ("_of" if outer == "fixed" else "") + (f"_n{ninc}" if ninc else "")
    return t


def summarise(r):
    el = r["elements"]
    th = np.degrees(np.arctan2(el["cy"], el["cx"])); z = np.array(el["cz"]); lay = np.array(el["layer"])
    xi = np.array(el["alpha_minus_1"]); vm = np.array(el["vM_over_E"]); p = np.array(el["p_over_E"])
    seed = lay < r["seed_layers"]
    core = seed & (z >= CORE_Z[0]) & (z <= CORE_Z[1]) & (th >= CORE_TH[0]) & (th <= CORE_TH[1])
    top = z > z.max() - 1e-9
    it = int(np.where(top)[0][vm[top].argmax()])
    p_ref = r.get("p_ref", 0.0)
    return dict(nr=r["nr"], bulge=r["bulge"], outer=r.get("outer", "free"), p_ref=p_ref, gate="implicit" if r.get("gate_implicit") else "explicit",
                increments=r["increments_nominal"] if "increments_nominal" in r else 60,
                core_alpha_minus_1=float(xi[core].mean()), core_vM_over_E=float(vm[core].mean()), core_p_over_E=float(p[core].mean()),
                seed_alpha_minus_1=float(xi[seed].mean()), seed_vM_over_E=float(vm[seed].mean()),
                corner_vM_over_E=float(vm[it]), corner_alpha_minus_1=float(xi[it]), corner_theta=float(th[it]),
                max_vM_over_E=r["max_vM_over_E"], seed_max_p_over_E=float(p[seed].max()),
                stopped_share=float((p[seed] >= p_ref).mean()) if p_ref > 0 else 0.0,
                cutbacks=r["incrementation"]["cutbacks"], wall_time_s=r.get("wall_time_s"))


def table(rows):
    hdr = ("| mesh | bulge | outer face | p_h/E | gate | increments | core alpha-1 | core vM/E | core p/E | seed alpha-1 | corner vM/E "
           "| corner alpha-1 | max p/E (seed) | stopped share | s |")
    out = [hdr, "|" + "---|" * 16]
    for s in rows:
        out.append(f"| {s['nr']}x{3 * s['nr']}x{4 * s['nr']} | {s['bulge']:g} | {s['outer']} | {s['p_ref'] if s['p_ref'] else 'off'} "
                   f"| {s['gate'] if s['p_ref'] else '-'} | {s['increments']} "
                   f"| {s['core_alpha_minus_1']:.4f} | {s['core_vM_over_E']:.4f} | {s['core_p_over_E']:.4f} | {s['seed_alpha_minus_1']:.4f} "
                   f"| {s['corner_vM_over_E']:.3f} | {s['corner_alpha_minus_1']:.4f} | {s['seed_max_p_over_E']:.3f} | {s['stopped_share']:.2f} "
                   f"| {s['wall_time_s'] or 0:.0f} |")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--p-ref", type=float, nargs="+", default=[0.0, 0.05])
    ap.add_argument("--outer", nargs="+", default=["free", "fixed"])
    ap.add_argument("--bulge", type=float, nargs="+", default=[0.0, 1.0])
    ap.add_argument("--ninc", type=int, default=None)
    ap.add_argument("--gate-implicit", action="store_true", help="the implicit gate (make_implant_uel_inp.py --gate-implicit; tag _gi)")
    ap.add_argument("--table-only", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    if not a.table_only:
        import make_implant_uel_inp as mi
        for f in a.levels:
            for bulge in a.bulge:
                for outer in a.outer:
                    for p_ref in a.p_ref:
                        tag = tag_of(f, bulge, p_ref, outer, a.ninc, a.gate_implicit and p_ref > 0)
                        path = OUT / f"implant_{tag}.json"
                        if path.exists():
                            continue
                        nodes, elems, layer, ijk = mi.sector(2 * f, 6 * f, 8 * f, bulge)
                        r = mi.run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, seed_layers=f, p_ref=p_ref, outer=outer, ninc=a.ninc,
                                   gate_implicit=a.gate_implicit and p_ref > 0)
                        r.update(nr=2 * f, nt=6 * f, nz=8 * f, bulge=bulge, nu=0.49, kinematics="finite strain, mean dilatation",
                                 elements_total=len(elems))
                        path.write_text(json.dumps(r))
                        s = summarise(r)
                        print(f"{tag}: core alpha-1 {s['core_alpha_minus_1']:.4f}, vM/E {s['core_vM_over_E']:.4f}, p/E {s['core_p_over_E']:.4f}, "
                              f"corner vM/E {s['corner_vM_over_E']:.3f}, {r['wall_time_s']:.0f} s", flush=True)
    rows = []
    for fpath in sorted(glob.glob(str(OUT / "implant_*.json"))):
        m = re.match(r"implant_(\d+)x(\d+)x(\d+)(_b1)?(_sl\d+)?(_ph[\d.]+)?(_gi)?(_of)?(_n\d+)?\.json$", Path(fpath).name)
        if not m:
            continue
        r = json.load(open(fpath))
        if "p_over_E" not in r["elements"] or r.get("seed_layers", 1) != int(m.group(1)) // 2:
            continue
        rows.append(summarise(r))
    rows.sort(key=lambda s: (s["nr"], s["bulge"], s["outer"] != "free", s["p_ref"], s["gate"] != "explicit", s["increments"]))
    (OUT / "implant_homeostatic_summary.json").write_text(json.dumps(rows, indent=0))
    print(table(rows))


if __name__ == "__main__":
    main()

"""threshold_study.py -- does the cervical concentration vanish at p_h = p_0 in other cases too? p_0 is the mean
pressure in the cervical band without the law (the pressure growth alone builds there). In the main case
(cervical_study.py: bulge 1, bottom u_z = 0, outer face free, 4x12x16) p_0/E = 0.214 and the ratio crosses 1 at
p_h/E = 0.2, i.e. p_h/p_0 = 0.93 (README, "Sensitivity to p_h/E"). Here the same test in cases that change p_0:

  b05   bulge 0.5 (crown slope 20 degrees instead of 38)
  bf    bottom face u = 0 (clamped) instead of u_z = 0
  t05   half the growth time (the GROW step stopped at T/2): a younger biofilm with less growth

For each case: the law-off run gives p_0 (band mean p/E), then p_h = m p_0 for m = 0.5, 1, 2. Nutrient from the
outer face, outer face free, explicit gate, seed 0.125 mm, 4x12x16 (or --nr 2 for 2x6x8).

    BIOFILM_ECO_CASE=<dir>/uel_eco.txt python threshold_study.py [--cases b05 bf t05 main] [--nr 4] [--mult 0.5 1 2]
Results: results/threshold_<case>_<mesh>_<m>.json, summary results/threshold_summary.json.
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_implant_uel_inp as mi  # noqa: E402
from cervical_study import measures  # noqa: E402
from driver import props_base  # noqa: E402
import numpy as np  # noqa: E402

OUT = HERE / "results"
CASES = {"b05": dict(bulge=0.5, bottom="uz", s=1.0), "bf": dict(bulge=1.0, bottom="fixed", s=1.0),
         "t05": dict(bulge=1.0, bottom="uz", s=1.0, half=True),
         "main": dict(bulge=1.0, bottom="uz", s=1.0)}


def set_time_link(s):
    def element_props(nu, p_ref=0.0, p_form=0, gate_implicit=False):
        P = props_base(s=s, D=0.2, L0=2.0, lam=0.0, psiref=0.0, nu=nu)
        if p_ref > 0:
            P = np.r_[P, 1e30, 0.0, p_ref, float(p_form + (10 if gate_implicit else 0))]
        return P
    mi.element_props = element_props


def one(case, nr, p_ref, label):
    c = CASES[case]
    path = OUT / f"threshold_{case}_{nr}x{3 * nr}x{4 * nr}_{label}.json"
    if path.exists():
        return json.loads(path.read_text())
    set_time_link(c["s"])
    f = nr // 2
    nodes, elems, layer, ijk = mi.sector(nr, 3 * nr, 4 * nr, c["bulge"])
    r = mi.run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, seed_layers=f, p_ref=p_ref, outer="free",
               bottom=c["bottom"], nut="outer", **({"nstop": mi.NINC // 2} if c.get("half") else {}))
    r.update(nr=nr, nt=3 * nr, nz=4 * nr, bulge=c["bulge"], nu=0.49, time_link_s=c["s"], case=case, label=label)
    s = measures(r)
    s.update(case=case, label=label, time_link_s=c["s"])
    path.write_text(json.dumps(dict(summary=s, run={k: v for k, v in r.items() if k != "elements"})))
    print(f"{case} {label}: p_h/E {p_ref:.4f}, ratio {s['ratio']:.2f}, band vM/E {s['band_vM_over_E']:.4f}, "
          f"band p/E {s['band_p_over_E']:.4f}, core alpha-1 {s['core_alpha_minus_1']:.4f}, cutbacks {s['cutbacks']}, "
          f"{r['wall_time_s']:.0f} s", flush=True)
    return dict(summary=s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="+", default=["b05", "bf", "t05"])
    ap.add_argument("--nr", type=int, default=4)
    ap.add_argument("--mult", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    a = ap.parse_args()
    rows = []
    for case in a.cases:
        off = one(case, a.nr, 0.0, "off")["summary"]
        p0 = off["band_p_over_E"]
        rows.append(dict(off, p0=p0, m=None))
        for m in a.mult:
            s = one(case, a.nr, m * p0, f"m{m:g}")["summary"]
            rows.append(dict(s, p0=p0, m=m))
    old = []
    sp = OUT / "threshold_summary.json"
    if sp.exists():
        old = [r for r in json.loads(sp.read_text()) if not any(r["case"] == n["case"] and r["label"] == n["label"]
                                                               and r["nr"] == n["nr"] for n in rows)]
    sp.write_text(json.dumps(old + rows, indent=0))
    print("\n| case | mesh | p_0/E | m = p_h/p_0 | ratio | band vM/E |\n|---|---|---|---|---|---|")
    for r in old + rows:
        print(f"| {r['case']} | {r['mesh']} | {r['p0']:.3f} | {r['m'] if r['m'] else 'off'} | {r['ratio']:.2f} "
              f"| {r['band_vM_over_E']:.4f} |")


if __name__ == "__main__":
    main()

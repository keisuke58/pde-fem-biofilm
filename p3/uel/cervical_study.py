"""cervical_study.py -- the cervical line of the tooth-like sector (make_implant_uel_inp.py --bulge 1) with a
gingiva-side condition on the bottom face, with and without the homeostatic-pressure growth law (KEIO_PLAN.ja.md
section 1, P3 (e); the keio_wp2 tooth of stages 6 and 7).

    BIOFILM_ECO_CASE=<dir>/uel_eco.txt python cervical_study.py [--levels 1 2] [--bottom uz fixed] [--outer free fixed]
                                                                 [--p-ref 0 0.05 0.02] [--nut top outer] [--gate-implicit] [--table-only]

Where the cervical line is. In the WP2 tooth (keio_wp2/homeostatic_tooth.py) the stress concentration sits at z = 0
of the first element ring on the bonded face: the crown's surface r_in(z) = r_i + b sin(pi z / 2h) is steepest there
(slope b pi / 2h, 38 degrees for b = 1 mm) and meets the bottom face, whose condition stands for the gingiva or the
bone (u_z = 0, u = 0, or free). The sector here has the same surface, so the cervical line is the bottom edge of the
bonded inner face, z = 0, r = r_i.

Measures (E = 1 in the decks, so vM/E is the von Mises stress itself; p = -tr(sigma)/3, Cauchy, positive in
compression, p/E likewise; per-element values are Gauss-point means):
  - cervical band: seed elements (the inner 0.125 mm) whose centroid has 0 <= z <= 0.25 mm and 15 <= theta <= 75
    degrees (one element row at 2x6x8, two at 4x12x16, three at 6x18x24: a fixed region like the core, away from the
    sector cuts);
  - core: the seed core of implant_convergence.py, seed elements with 0.25 <= z <= 1.75 and 15 <= theta <= 75;
  - ratio = mean vM/E of the band / mean vM/E of the core, the cervical concentration;
  - band max = the largest element vM/E in the band, and the largest p/E over the seed.
The bottom edge is a corner of the elastic problem (a rigid oblique face meets a constrained bottom face), so its
largest value need not converge with the mesh; the ratio uses band means, and the band max is reported beside it.

Runs: bulge 1, seed fixed at 0.125 mm (--seed-layers f), bottom u_z = 0 ("uz") or u = 0 ("fixed"), outer face free
or u = 0, law off or p_h/E = 0.05 and 0.02 (form 0), levels f = 1, 2 (3 for the main case), and the nutrient on the
top face (the implant deck; the bottom row then hardly grows within T) or on the outer face (saliva on a crown, the
WP2 tooth of stage 7, `--nut outer`, tag `_no`), where the growth is uniform along z and the cervical line can be
loaded as in keio_wp2. A case whose
results/implant_<tag>.json exists is not rerun; bottom-free runs of homeostatic_implant.py are tabulated too when
present. Writes results/implant_cervical_summary.json and prints the table, the spread of the ratio over the four
boundary-condition combinations and the change between levels. --gate-implicit runs the law-on cases with the
implicit gate (tag `_gi`, column "gate"); the spread lines use the explicit runs.
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
BAND_Z = (0.0, 0.25)            # the cervical band in z (mm); theta as the core
BOTTOM_TAG = {"free": "", "uz": "_bz", "fixed": "_bf"}


def tag_of(f, bottom, outer, p_ref, bulge=1.0, nut="top", gate_implicit=False):
    t = f"{2 * f}x{6 * f}x{8 * f}" + ("_b1" if bulge else "") + (f"_sl{f}" if f > 1 else "")
    t += f"_ph{p_ref:g}" + ("_gi" if gate_implicit else "") if p_ref > 0 else ""
    return t + ("_of" if outer == "fixed" else "") + BOTTOM_TAG[bottom] + ("_no" if nut == "outer" else "")


def masks(r):
    """band and core element masks of a result dict (make_implant_uel_inp.run)"""
    el = r["elements"]
    th = np.degrees(np.arctan2(el["cy"], el["cx"])); z = np.array(el["cz"])
    seed = np.array(el["layer"]) < r.get("seed_layers", 1)
    sector = (th >= CORE_TH[0]) & (th <= CORE_TH[1])
    band = seed & sector & (z >= BAND_Z[0]) & (z <= BAND_Z[1])
    core = seed & sector & (z >= CORE_Z[0]) & (z <= CORE_Z[1])
    return band, core, seed


def measures(r):
    el = r["elements"]
    xi = np.array(el["alpha_minus_1"]); vm = np.array(el["vM_over_E"]); p = np.array(el["p_over_E"])
    th = np.degrees(np.arctan2(el["cy"], el["cx"])); z = np.array(el["cz"]); rr = np.hypot(el["cx"], el["cy"])
    band, core, seed = masks(r)
    ib = int(np.where(band)[0][vm[band].argmax()]) if band.any() else -1
    cm = float(vm[core].mean()) if core.any() else float("nan")
    bm = float(vm[band].mean()) if band.any() else float("nan")
    return dict(nr=r["nr"], mesh=f"{r['nr']}x{r['nt']}x{r['nz']}", bulge=r.get("bulge", 0.0), bottom=r.get("bottom", "free"),
                outer=r.get("outer", "free"), nut=r.get("nut", "top"), p_ref=r.get("p_ref", 0.0),
                gate="implicit" if r.get("gate_implicit") else "explicit",
                band_elements=int(band.sum()), core_elements=int(core.sum()),
                core_alpha_minus_1=float(xi[core].mean()) if core.any() else float("nan"), core_vM_over_E=cm,
                core_p_over_E=float(p[core].mean()) if core.any() else float("nan"),
                band_alpha_minus_1=float(xi[band].mean()) if band.any() else float("nan"), band_vM_over_E=bm,
                band_p_over_E=float(p[band].mean()) if band.any() else float("nan"),
                ratio=bm / cm if core.any() and cm > 0 else float("nan"),
                band_max_vM_over_E=float(vm[ib]) if ib >= 0 else float("nan"),
                band_max_at_r_theta_z=[float(rr[ib]), float(th[ib]), float(z[ib])] if ib >= 0 else None,
                seed_max_p_over_E=float(p[seed].max()), max_vM_over_E=r["max_vM_over_E"],
                cutbacks=r.get("incrementation", {}).get("cutbacks"), wall_time_s=r.get("wall_time_s"))


def table(rows):
    out = ["| nutrient | mesh | bottom | outer | p_h/E | gate | core alpha-1 | core vM/E | band mean vM/E | ratio | band max vM/E | max p/E (seed) |",
           "|" + "---|" * 12]
    for s in rows:
        out.append(f"| {s['nut']} | {s['mesh']} | {s['bottom']} | {s['outer']} | {s['p_ref'] if s['p_ref'] else 'off'} "
                   f"| {s['gate'] if s['p_ref'] else '-'} | {s['core_alpha_minus_1']:.4f} "
                   f"| {s['core_vM_over_E']:.4f} | {s['band_vM_over_E']:.4f} | {s['ratio']:.2f} | {s['band_max_vM_over_E']:.3f} "
                   f"| {s['seed_max_p_over_E']:.3f} |")
    return "\n".join(out)


def spread(rows, bottoms=("uz", "fixed"), outers=("free", "fixed")):
    """per level and p_h: the ratio over the boundary-condition combinations, its spread (max - min) / mean, and the
    same for the band mean vM/E; then the change of the ratio from one level to the next per case"""
    out, by = [], {}
    for s in rows:
        if s["gate"] == "implicit":
            continue                                        # the spread lines compare the explicit runs
        by[(s["nut"], s["nr"], s["p_ref"], s["bottom"], s["outer"])] = s
    nuts = sorted({s["nut"] for s in rows}, reverse=True)
    levels = sorted({s["nr"] for s in rows}); prefs = sorted({s["p_ref"] for s in rows})
    for nut in nuts:
        for nr in levels:
            for pr in prefs:
                sel = [by[k] for k in ((nut, nr, pr, b, o) for b in bottoms for o in outers) if k in by]
                if len(sel) < 2:
                    continue
                rat = np.array([s["ratio"] for s in sel]); bm = np.array([s["band_vM_over_E"] for s in sel])
                out.append(f"  nutrient {nut}, {sel[0]['mesh']}, p_h/E {pr if pr else 'off'}: ratio "
                           + ", ".join(f"{s['bottom']}/{s['outer']} {s['ratio']:.2f}" for s in sel)
                           + f"; spread {(rat.max() - rat.min()) / rat.mean():.0%} (band mean vM/E spread {(bm.max() - bm.min()) / bm.mean():.0%})")
        for pr in prefs:
            for b in bottoms:
                for o in outers:
                    seq = [by[(nut, nr, pr, b, o)] for nr in levels if (nut, nr, pr, b, o) in by]
                    if len(seq) >= 2:
                        out.append(f"  nutrient {nut}, {b}/{o}, p_h/E {pr if pr else 'off'}: ratio " + " -> ".join(f"{s['ratio']:.2f}" for s in seq)
                                   + " (" + ", ".join(f"{(s2['ratio'] - s1['ratio']) / s1['ratio']:+.0%}" for s1, s2 in zip(seq, seq[1:])) + ")")
    return "\n".join(out)


def load_rows(bulge=1.0):
    rows = []
    for fpath in sorted(glob.glob(str(OUT / "implant_*.json"))):
        m = re.match(r"implant_(\d+)x(\d+)x(\d+)(_b1)?(_sl\d+)?(_ph[\d.]+)?(_gi)?(_of)?(_b[zf])?(_no)?\.json$", Path(fpath).name)
        if not m or bool(m.group(4)) != bool(bulge):
            continue
        r = json.load(open(fpath))
        if "p_over_E" not in r["elements"] or r.get("seed_layers", 1) != int(m.group(1)) // 2:
            continue
        rows.append(measures(r))
    order = {"free": 0, "uz": 1, "fixed": 2}
    rows.sort(key=lambda s: (s["nut"] != "top", s["nr"], order[s["bottom"]], s["outer"] != "free", -s["p_ref"] if s["p_ref"] else 0,
                             s["gate"] != "explicit"))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--bottom", nargs="+", default=["uz", "fixed"])
    ap.add_argument("--outer", nargs="+", default=["free", "fixed"])
    ap.add_argument("--p-ref", type=float, nargs="+", default=[0.0, 0.05, 0.02])
    ap.add_argument("--nut", nargs="+", default=["top", "outer"], help="face with c = 1: top (implant deck) and/or outer (crown)")
    ap.add_argument("--gate-implicit", action="store_true", help="the implicit gate for the law-on cases (tag _gi)")
    ap.add_argument("--table-only", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    if not a.table_only:
        import make_implant_uel_inp as mi
        for f in a.levels:
            for nut in a.nut:
                for bottom in a.bottom:
                    for outer in a.outer:
                        for p_ref in a.p_ref:
                            tag = tag_of(f, bottom, outer, p_ref, nut=nut, gate_implicit=a.gate_implicit and p_ref > 0)
                            path = OUT / f"implant_{tag}.json"
                            if path.exists():
                                continue
                            nodes, elems, layer, ijk = mi.sector(2 * f, 6 * f, 8 * f, 1.0)
                            r = mi.run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, seed_layers=f, p_ref=p_ref, outer=outer,
                                       bottom=bottom, nut=nut, gate_implicit=a.gate_implicit and p_ref > 0)
                            r.update(nr=2 * f, nt=6 * f, nz=8 * f, bulge=1.0, nu=0.49, kinematics="finite strain, mean dilatation",
                                     elements_total=len(elems))
                            path.write_text(json.dumps(r))
                            s = measures(r)
                            print(f"{tag}: core vM/E {s['core_vM_over_E']:.4f}, band {s['band_vM_over_E']:.4f}, ratio {s['ratio']:.2f}, "
                                  f"band max {s['band_max_vM_over_E']:.3f}, {r['wall_time_s']:.0f} s", flush=True)
    rows = load_rows()
    (OUT / "implant_cervical_summary.json").write_text(json.dumps(rows, indent=0))
    print(table(rows))
    print("\nspread over the boundary conditions and change between levels:")
    print(spread(rows))


if __name__ == "__main__":
    main()

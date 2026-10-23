"""gci_a52.py -- A-52: grid convergence index (Roache 1994; Celik et al. 2008, J. Fluids Eng. 130, 078001) of the
separate example with antibiotic (front 1, alpha*_0 = 100, antibiotic from T* = 1/3) at T* = 1, for each implementation
on its own meshes, then the comparison of the implementations within their numerical uncertainty (ASME V&V 20 style).

  Python (finite volumes, beta + r h/4 as ANSYS, wall gradient mirror, increment 0.00125): 8 .. 96
  Python (TVD front term, no added diffusion; reference without stabilisation):               16, 32, 64, 128
  Abaqus (Galerkin hex8 element, beta + r h/4, increment 0.00125):                            8, 16, 24, 32, 48 (quarter)
  Python (LW + sensor, explicit as the element, A-59): 16 .. 64; Abaqus SUPG (quarter): 16 .. 48; ANSYS living10: 16, 24
  ANSYS (Oliver's element, living8 build, beta + r h/4):                                      8, 16, 24 (24: increment 0.0005)

The stabilisation r h/4 vanishes with h, so all series approach the same continuum solution (first order expected
with r h/4, second order for the TVD series). Quantities: seed means of alpha - 1 and von Mises (Pa).

    python gci_a52.py [--out results/2026-10-living/gci_a52]
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PY = HERE / "results/2026-10-3d"
AB = HERE.parent.parent / "abaqus_living"
AN = HERE / "results/2026-10-living"
OUT = Path(sys.argv[sys.argv.index("--out") + 1]).resolve() if "--out" in sys.argv else AN / "gci_a52"
Q = ("seed_am1", "seed_vm_Pa")


def py_val(pattern):
    fs = sorted(glob.glob(str(PY / pattern)))
    if not fs:
        return None
    r = json.load(open(fs[-1]))["runs"]["separate_ab"]
    return {"seed_am1": r["seed_am1"], "seed_vm_Pa": r["seed_vm_Pa"], "file": Path(fs[-1]).name}


def ab_val(name):
    p = AB / f"{name}_summary.json"
    if not p.exists():
        return None
    s = json.load(open(p))[-1]
    if abs(s["t"] - 1.0) > 1e-3:
        return None
    return {"seed_am1": s["seed_am1"], "seed_vm_Pa": s["seed_vm_Pa"], "file": p.name}


def an_val(run):
    p = AN / f"fields_{run}.json"
    if not p.exists():
        return None
    f = json.load(open(p))
    if abs(f["t"][-1] - 1.0) > 1e-3:
        return None
    s = np.isin(np.array(f["elem"]), f["seed"])
    return {"seed_am1": float(np.array(f["alpha"][-1])[s].mean()),
            "seed_vm_Pa": float(np.array(f["seqv"][-1])[s].mean()) * 1e6, "file": p.name}


SERIES = {
    "Python (r h/4)": {n: py_val(f"*nf{n}_front1_dtc0.00125_*ansysbc_startpos_a100_gab1_gc1_tab0.333333_artdiff_wgmirror_fast_separate_ab.json")
                       for n in (8, 16, 24, 48, 64, 96)}
                      | {32: py_val("coupled_living_consumption_front1_dtc0.00125_*ansysbc_startpos_a100_gab1_gc1_tab0.333333_artdiff_wgmirror_fast_separate_ab.json")},
    "Python (TVD)": {16: py_val("*nf16_front1_phimin*ansysbc_startpos_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_separate_ab.json"),
                     32: py_val("coupled_living_consumption_front1_phimin*ansysbc_startpos_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_separate_ab.json"),
                     64: py_val("*nf64_front1_phimin*ansysbc_startpos_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_separate_ab.json"),
                     128: py_val("*nf128_front1_phimin*ansysbc_startpos_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_separate_ab.json")},
    "Abaqus": {8: ab_val("lv8_sep_ab"), 16: ab_val("lv16_sep_ab_dt000125"), 24: ab_val("lv24_sep_ab"), 32: ab_val("lv32_sep_ab") or ab_val("lv32q_sep_ab"), 48: ab_val("lv48q_sep_ab")},
    "ANSYS": {8: an_val("u8c_8_sep_ab_dt0005"), 16: an_val("u8_16_sep_ab"), 24: an_val("u8_24_sep_ab_dt0005")},
    # A-57..A-59: the second-order leaning schemes
    "Python (LW + sensor)": {n: py_val(f"*nf{n}_front1_dtc0.00125_*ansysbc_startpos_a100_gab1_gc1_tab0.333333_lw_sensor2_wgmirror_fast_separate_ab.json")
                             for n in (16, 24, 48, 64)}
                            | {32: py_val("coupled_living_consumption_front1_dtc0.00125_*ansysbc_startpos_a100_gab1_gc1_tab0.333333_lw_sensor2_wgmirror_fast_separate_ab.json")},
    "Abaqus (SUPG)": {n: ab_val(f"lv{n}qs2_sep_ab") for n in (16, 24, 32, 48)},
    "ANSYS (living10)": {16: an_val("lw16_16_sep_ab"), 24: an_val("lw24_24_sep_ab")},
    "ANSYS (upwind)": {8: an_val("up8_8_sep_ab"), 16: an_val("up16_16_sep_ab")},
}


def gci3(n1, n2, n3, f1, f2, f3):
    """Celik et al. 2008: grids 1 (finest) .. 3, h = 2/n; returns p, f_ext, relative GCI of grid 1, oscillatory flag"""
    r21, r32 = n1 / n2, n2 / n3
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0:
        return None
    s = np.sign(e32 / e21)
    p = 1.0
    for _ in range(200):
        q = np.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(np.log(abs(e32 / e21)) + q) / np.log(r21)
        if abs(p_new - p) < 1e-10:
            break
        p = p_new
    f_ext = (r21 ** p * f1 - f2) / (r21 ** p - 1)
    gci = 1.25 * abs((f1 - f2) / f1) / (r21 ** p - 1)
    return {"p": p, "f_ext": f_ext, "gci": gci, "oscillatory": bool(s < 0), "grids": [n1, n2, n3]}


def gci2(n1, n2, f1, f2, p):
    """two grids, assumed order p, safety factor 3 (Roache)"""
    r = n1 / n2
    return {"p": p, "f_ext": (r ** p * f1 - f2) / (r ** p - 1), "gci": 3.0 * abs((f1 - f2) / f1) / (r ** p - 1),
            "oscillatory": None, "grids": [n1, n2], "assumed_p": True}


res = {}
for code, ser in SERIES.items():
    have = sorted(n for n, v in ser.items() if v)
    res[code] = {"values": {n: ser[n] for n in have}}
    for q in Q:
        out = []
        if len(have) >= 3:
            for k in range(len(have) - 2):
                n3, n2, n1 = have[k:k + 3]
                g = gci3(n1, n2, n3, ser[n1][q], ser[n2][q], ser[n3][q])
                if g:
                    out.append(g)
        if len(have) >= 2:   # always also the finest pair with p = 1 (r h/4 is first order), Fs = 3
            n2, n1 = have[-2:]
            out.append(gci2(n1, n2, ser[n1][q], ser[n2][q], 1.0))
        res[code][q] = out

lines = ["# A-52 GCI（自動）", "",
         "T\\* = 1、分けた例（抗菌薬あり）、種の平均。h = 2/n。GCI は最も細かい格子の値の相対的な誤差の幅（Celik 2008、Fs = 1.25、二つの格子だけのときは Fs = 3、p = 1 と仮定）。", ""]
for q in Q:
    lines += [f"## {q}", "", "| 実装 | 格子 | 値（細かい順） | p | 外挿値 | GCI | 振動 |", "|---|---|---|---|---|---|---|"]
    for code, r in res.items():
        vals = r["values"]
        for g in r[q]:
            v = ", ".join(f"{vals[n][q]:.4e}" for n in g["grids"])
            lines.append(f"| {code} | {'/'.join(map(str, g['grids']))} | {v} | {g['p']:.2f}{' (仮定)' if g.get('assumed_p') else ''} | "
                         f"{g['f_ext']:.4e} | {100 * g['gci']:.1f} % | {'' if g['oscillatory'] is None else ('はい' if g['oscillatory'] else 'いいえ')} |")
    lines.append("")
    # comparison of the finest values within their bands: |E| against u = sqrt(U_A^2 + U_B^2), U = GCI * |f1|
    best = {}
    for code, r in res.items():
        # the finest triple if it is in the asymptotic range (monotone, 0.5 <= p <= 3), else the finest pair, Fs = 3
        tri = [g for g in r[q] if len(g["grids"]) == 3]
        g = tri[-1] if tri and not tri[-1]["oscillatory"] and 0.5 <= tri[-1]["p"] <= 3.0 else (r[q][-1] if r[q] else None)
        if g:
            n1 = g["grids"][0]; f1 = r["values"][n1][q]
            best[code] = (n1, f1, g["gci"] * abs(f1), g["f_ext"])
    lines += [f"比べ方：最も細かい格子の値 f1 ± U（U = GCI·|f1|）。|E| ≤ u = √(U_A² + U_B²) なら数値の誤差の幅の中で一致。", "",
              "| A | B | f1 A | f1 B | E/f1 B | u/f1 B | |E|/u | 外挿値の差 |", "|---|---|---|---|---|---|---|---|"]
    ks = list(best)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            a, b = best[ks[i]], best[ks[j]]
            E = a[1] - b[1]; u = np.hypot(a[2], b[2])
            lines.append(f"| {ks[i]} ({a[0]}) | {ks[j]} ({b[0]}) | {a[1]:.4e} | {b[1]:.4e} | {100 * E / b[1]:+.1f} % | "
                         f"{100 * u / abs(b[1]):.1f} % | {abs(E) / u:.2f} | {100 * (a[3] - b[3]) / b[3]:+.1f} % |")
    lines.append("")
missing = [f"{c} {n}" for c, s in SERIES.items() for n, v in s.items() if not v]
if missing:
    lines += ["まだない値：" + "、".join(missing), ""]
OUT.with_suffix(".md").write_text("\n".join(lines), encoding="utf-8")
OUT.with_suffix(".json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
print("\n".join(lines))
print(f"FILE:{OUT.with_suffix('.md').relative_to(HERE.parent.parent).as_posix()}")
print(f"FILE:{OUT.with_suffix('.json').relative_to(HERE.parent.parent).as_posix()}")

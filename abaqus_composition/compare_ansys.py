"""compare_ansys.py -- an Abaqus cube run (make_cube_inp.py) next to the ANSYS run
it was built from, with the measures of summarize_runs_json.py.

    python abaqus_composition/compare_ansys.py JOB.dat RUN.json

Reads the last increment's centroidal S11..S23 and SDV84 (alpha - 1), TEMP (phi)
from the .dat (*EL PRINT, POSITION=CENTROIDAL), takes the seed from
JOB.seed.txt next to the input (make_cube_inp.py), and prints for both
programs: seed mean von Mises and mean stress p, p in layers 1-3 with the share
in tension, alpha - 1 in the seed interior and on its surface and the step
between them, and the mean phi of the seed and of layer 1. Pa.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ansys_usermat" / "apdl"))
from summarize_runs_json import layers, surface  # noqa: E402

NUM = re.compile(r"^\s+(\d+)\s+((?:[-+]?\d+\.\d*(?:E[-+]\d+)?\s*)+)$")


def last_tables(dat):
    lines = Path(dat).read_text(errors="replace").splitlines()
    start = max(i for i, l in enumerate(lines) if re.search(r"INCREMENT\s+\d+ SUMMARY", l))
    tables, cur = [], None
    for l in lines[start:]:
        if "ELEMENT  FOOT-" in l:
            cur = {"head": l.split()[2:], "rows": {}}
            tables.append(cur)
            continue
        m = NUM.match(l)
        if m and cur is not None:
            cur["rows"][int(m.group(1))] = [float(x) for x in m.group(2).split()]
    return tables


def measures(elem, cen, seed_ids, s6, alpha, phi):
    seed = np.isin(elem, seed_ids)
    lay, h = layers(cen, seed)
    sx, sy, sz, sxy, sxz, syz = (s6[:, i] for i in range(6))
    p = (sx + sy + sz) / 3
    vm = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3 * (sxy ** 2 + sxz ** 2 + syz ** 2))
    surf = surface(cen, seed, h)
    inner = seed & ~surf
    out = {"seed vM": vm[seed].mean(), "seed p": p[seed].mean()}
    for k in (1, 2, 3):
        out[f"layer {k} p"] = p[lay == k].mean()
        out[f"layer {k} tension"] = float(np.mean(p[lay == k] > 0))
    ai, as_ = alpha[inner].mean(), alpha[surf].mean()
    out.update({"alpha interior": ai, "alpha surface": as_, "alpha step %": (ai - as_) / ai * 100})
    if phi is not None:
        out.update({"phi seed": phi[seed].mean(), "phi layer 1": phi[lay == 1].mean()})
    return out


def main(dat, js):
    t = last_tables(dat)
    S = next(x for x in t if x["head"][:1] == ["S11"])["rows"]
    V = next(x for x in t if "SDV84" in x["head"])
    iv, it = V["head"].index("SDV84"), V["head"].index("TEMP") if "TEMP" in V["head"] else None
    el = np.array(sorted(S))
    n = round(len(el) ** (1 / 3))
    h = 2.0 / n
    ijk = [((e - 1) % n, (e - 1) // n % n, (e - 1) // n // n) for e in el]
    cen = np.array([[-1 + h * (i + 0.5), -1 + h * (j + 0.5), -1 + h * (k + 0.5)] for i, j, k in ijk])
    seed_ids = [int(x) for x in Path(dat).with_suffix(".seed.txt").read_text().split()] \
        if Path(dat).with_suffix(".seed.txt").exists() else \
        [int(x) for x in Path(dat).parent.joinpath(Path(dat).stem + ".seed.txt").read_text().split()]
    s6 = np.array([S[e] for e in el]) * 1e6
    alpha = np.array([V["rows"][e][iv] for e in el])
    phi = np.array([V["rows"][e][it] for e in el]) if it is not None else None
    ab = measures(el, cen, seed_ids, s6, alpha, phi)

    rec = json.loads(Path(js).read_text())
    a = rec["all_stress"]
    ael = np.array(a["elem"], int)
    acen = np.stack([a["cx"], a["cy"], a["cz"]], 1)
    as6 = np.stack([a["sx"], a["sy"], a["sz"], a.get("sxy", np.zeros(len(ael))),
                    a.get("sxz", np.zeros(len(ael))), a.get("syz", np.zeros(len(ael)))], 1) * 1e6
    an = measures(ael, acen, rec["seed_BIOFILM1"], as6, np.array(a["alpha"]), None)
    an["seed vM"] = float(np.mean(np.array(a["seqv"])[np.isin(ael, rec["seed_BIOFILM1"])]) * 1e6)
    print(f"{'measure':20s} {'Abaqus':>12s} {'ANSYS':>12s} {'ratio':>8s}")
    for k in ab:
        x, y = ab[k], an.get(k)
        r = f"{x / y:8.4f}" if y not in (None, 0) else ""
        print(f"{k:20s} {x:12.4e} " + (f"{y:12.4e} {r}" if y is not None else ""))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

"""maze_compare_a67.py -- the maze (Klempt et al. 2024, Sect. 4.4, Fig. 11) in the Abaqus element (make_maze_inp.py,
SUPG, 4 um cells, one layer) against the Python finite-volume code (maze_coupled_2d.py, same 4 um cells): number of
colonised cells (phi >= phi_min) and the closest phi >= 0.5 to the maze centre (um) over time, without and with the
antibiotic.

    python maze_compare_a67.py F:\\abaqus_work\\maze60b\\living_gp.txt F:\\abaqus_work\\maze60b_ab\\living_gp.txt
        -> ansys_usermat/apdl/results/2026-10-maze/maze_abaqus_vs_python_a67.json and a table on stdout
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "ansys_usermat" / "apdl" / "results" / "2026-10-maze"
PHI_MIN = 1e-3


def blocks(path):
    """(t, array) per output step of living_gp.txt (element gp x y z phi c a alpha-1 living vM p shareA started)"""
    t, cur = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("# step"):
                if cur:
                    yield t, np.array(cur)
                t, cur = float(line.split()[-1]), []
            elif line[0] != "#":
                cur.append(line.split())
    if cur:
        yield t, np.array(cur)


def measures(path):
    out = {}
    for t, a in blocks(path):
        a = a.astype(float)
        el, x, y, phi = a[:, 0].astype(int), a[:, 2], a[:, 3], a[:, 5]
        if not out:
            measures.centre = ((x.min() + x.max()) / 2, (y.min() + y.max()) / 2)
        x0, y0 = measures.centre
        phie = np.bincount(el, phi) / np.maximum(np.bincount(el), 1)
        r = np.hypot(x - x0, y - y0) * 1e3
        on = a[:, 13] > 0.5                                    # point model started
        out[round(t, 3)] = dict(colonised=int((phie[np.bincount(el) > 0] >= PHI_MIN).sum()),
                                reach_um=float(r[phi >= 0.5].min()) if (phi >= 0.5).any() else None,
                                living_mean=float(a[on, 9].mean()) if on.any() else 0.0)
    return out


res = {}
for p, ab in zip(sys.argv[1:3], (False, True)):
    A = measures(p)
    P = json.load(open(RES / f"maze_h4_tend60{'_ab' if ab else ''}_k24fig.json"))["hist"]
    tp = np.array(P["t"][1:])
    rows = []
    for t in sorted(A):
        if t < 1 or abs(t / 6 - round(t / 6)) > 1e-6:          # Abaqus output every 6 T*
            continue
        k = int(np.argmin(abs(tp - t)))
        rows.append(dict(t=t, abaqus=A[t], python=dict(colonised=P["colonised"][k + 1], reach_um=P["reach_um"][k], living_mean=P["living_mean"][k])))
    res["antibiotic" if ab else "no_antibiotic"] = rows
    print(f"\n{'with' if ab else 'without'} antibiotic: T*, colonised cells (Abaqus / Python), closest phi >= 0.5 to the "
          f"centre in um (Abaqus / Python), mean living fraction of the colonised points (Abaqus / Python)")
    for r in rows:
        ra, rp = r["abaqus"]["reach_um"], r["python"]["reach_um"]
        print(f"  {r['t']:5.1f}  {r['abaqus']['colonised']:6d} / {r['python']['colonised']:6d}   "
              f"{'-' if ra is None else f'{ra:6.1f}'} / {'-' if rp is None else f'{rp:6.1f}'}   "
              f"{r['abaqus']['living_mean']:.3f} / {r['python']['living_mean']:.3f}")
(RES / "maze_abaqus_vs_python_a67.json").write_text(json.dumps(res, indent=1))
print("wrote", RES / "maze_abaqus_vs_python_a67.json")

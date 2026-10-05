"""compare_klempt.py -- domain means of phi and c over time from an Abaqus run of
make_klempt_inp.py, next to the finite-difference reproduction
(JAXFEM/klempt2024_quantitative.py, growth="printed", zero-order consumption).

    python abaqus_composition/compare_klempt.py JOB.dat [--case fig7_high] [--ref ref.json]

The mean of the centroid values of a trilinear field over equal hexahedra is
its volume mean, the same quantity as the reproduction's trapezoid average.
--ref: a json of the reproduction (t, phi, c); without it the reproduction is
run here (about 40 s for T* = 0.2).
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def series(dat):
    t, phi, c, cur, rows, tab_t = [], [], [], None, {}, None
    for l in Path(dat).read_text(errors="replace").splitlines():
        m = re.search(r"STEP TIME COMPLETED\s+([0-9.Ee+-]+)", l)
        if m:
            cur = float(m.group(1).rstrip(","))
            continue
        if "ELEMENT  FOOT-" in l:
            if rows and tab_t is not None:
                v = np.array(list(rows.values()))
                t.append(tab_t); phi.append(v[:, 0].mean()); c.append(v[:, 1].mean())
            rows, tab_t = {}, cur                    # the increment summary comes before its table
            continue
        m = re.match(r"^\s+(\d+)\s+((?:[-+]?\d+\.\d*(?:E[-+]\d+)?\s*){2})$", l)
        if m:
            rows[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    if rows and tab_t is not None:
        v = np.array(list(rows.values()))
        t.append(tab_t); phi.append(v[:, 0].mean()); c.append(v[:, 1].mean())
    return np.array(t), np.array(phi), np.array(c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dat")
    ap.add_argument("--case", default="fig7_high")
    ap.add_argument("--ref", default=None)
    a = ap.parse_args()
    t, phi, c = series(a.dat)
    if a.ref:
        ref = json.loads(Path(a.ref).read_text())
    else:
        sys.path.insert(0, str(ROOT / "JAXFEM"))
        import klempt2024_quantitative as k
        ref = k.run(a.case, "zero_order", growth="printed", t_end=float(t[-1]))
    rt = np.array(ref["t"])
    print(f"{'T*':>6s} {'phi Abaqus':>11s} {'phi FD':>9s} {'ratio':>7s} {'c Abaqus':>9s} {'c FD':>8s}")
    for ti, p, cc in zip(t, phi, c):
        i = int(np.argmin(np.abs(rt - ti)))
        print(f"{ti:6.3f} {p:11.4f} {ref['phi'][i]:9.4f} {p / ref['phi'][i]:7.3f} {cc:9.4f} {ref['c'][i]:8.4f}")


if __name__ == "__main__":
    main()

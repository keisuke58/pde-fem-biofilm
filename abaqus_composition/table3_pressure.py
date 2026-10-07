"""table3_pressure.py -- Klempt 2024 Table 3 (hydrostatic stress of test case
4.1) from an Abaqus run of make_klempt_inp.py --case fig4_corner --s-every.

    python abaqus_composition/table3_pressure.py JOB.dat --E 10

Per printed increment: max(alpha - 1) and the range of p = tr(sigma)/3
(tension > 0) over all elements and over the biofilm (centroid phi > 0.5),
in the stress unit of the run and divided by E. Table 3's legend is
+3e-4 ... -6.5e-4 MPa; JAXFEM/klempt2024_pressure_scale.py (small strain,
edge source) gave p/E -3.4e-5 ... +1.5e-5 at T* = 0.05.
The pressure of BIOFILM_STRESS_CORE carries the spherical term of
DEVIATOR_SCALING_FINDING.md, about (1 - 2 nu)/(1 + nu) = 1.3 % of the
volumetric stress at nu = 0.49.
"""
import argparse
import re
from pathlib import Path

import numpy as np

NUM = re.compile(r"^\s+(\d+)\s+((?:[-+]?\d+\.\d*(?:E[-+]\d+)?\s*)+)$")


def tables(dat):
    out, cur_t, cur = {}, None, None
    for l in Path(dat).read_text(errors="replace").splitlines():
        m = re.search(r"STEP TIME COMPLETED\s+([0-9.Ee+-]+)", l)
        if m:
            cur_t = float(m.group(1).rstrip(","))
            continue
        if "ELEMENT  FOOT-" in l:
            cur = {"head": l.split()[2:], "rows": {}}
            out.setdefault(round(cur_t, 6), []).append(cur)
            continue
        m = NUM.match(l)
        if m and cur is not None:
            cur["rows"][int(m.group(1))] = [float(x) for x in m.group(2).split()]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dat")
    ap.add_argument("--E", type=float, default=10.0)
    a = ap.parse_args()
    print(f"{'T*':>5s} {'max a-1':>9s} {'p min':>10s} {'p max':>10s} {'p/E min':>10s} {'p/E max':>10s}"
          f" {'biofilm p min':>13s} {'max':>10s}")
    for t, tabs in sorted(tables(a.dat).items()):
        S = next((x for x in tabs if "S11" in x["head"]), None)
        V = next((x for x in tabs if "TEMP" in x["head"]), None)
        if S is None or V is None:
            continue
        el = sorted(S["rows"])
        s = np.array([S["rows"][e] for e in el])
        h = S["head"]
        p = (s[:, h.index("S11")] + s[:, h.index("S22")] + s[:, h.index("S33")]) / 3
        a1 = s[:, h.index("SDV84")]
        phi = np.array([V["rows"][e][V["head"].index("TEMP")] for e in el])
        b = phi > 0.5
        print(f"{t:5.2f} {a1.max():9.2e} {p.min():10.2e} {p.max():10.2e} {p.min() / a.E:10.2e} {p.max() / a.E:10.2e}"
              f" {p[b].min():13.2e} {p[b].max():10.2e}")


if __name__ == "__main__":
    main()

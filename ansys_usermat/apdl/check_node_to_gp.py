#!/usr/bin/env python3
"""check_node_to_gp.py -- NOTE (A-20, 20 Oct): not applicable to the partner element. In the neighboured element
method the fields phi, c, a live at the integration points (and surface points); there are no nodal values and no
interpolation, the material routine reads its own point. Kept for elements with nodal fields.

 are the nutrient c and the antibiotic a that the material routine receives at a Gauss point
the element's nodal values interpolated with its shape functions? (A-19, Oliver's question of 20 Oct.)

Input (CSV, from IKMHIWI03, any result set of a run of the coupling, e.g. lv8_sep at T* = 0.5):
  nodes.csv : node, x, y, z, c, a            nodal values of the element's two nutrient unknowns
  elems.csv : elem, n1, ..., n8              connectivity in the ANSYS SOLID185 order (I J K L M N O P)
  gauss.csv : elem, ip, c_gp, a_gp           the values the material routine read at each integration point
                                             (SVAR 51, 52 with ERESX,NO, or from the trace), all 8 per element
For every element the 8 trilinear interpolations at the natural coordinates (+-1/sqrt 3)^3 are compared with the 8
Gauss-point values. The comparison sorts both sets of 8, so it does not depend on how the integration points are
numbered. A second comparison against the element mean shows whether the routine gets the point values or only an
element average.

    python ansys_usermat/apdl/check_node_to_gp.py DIR
Pass: max |c_gp - interpolation| at round-off (~1e-12 relative); an element-mean hand-over would instead give
|c_gp - mean| ~ 0 and |c_gp - interpolation| ~ the in-element variation.
"""
import csv
import sys
from pathlib import Path

import numpy as np

# SOLID185 corner nodes I..P in natural coordinates (r, s, t)
CORNERS = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
G = 1.0 / np.sqrt(3.0)
GPS = np.array([[r, s, t] for t in (-G, G) for s in (-G, G) for r in (-G, G)])


def shape(rst):
    return 0.125 * np.prod(1.0 + CORNERS * rst[None, :], axis=1)


N_AT_GP = np.array([shape(p) for p in GPS])        # 8 x 8


def read(path):
    with open(path) as f:
        rows = [r for r in csv.reader(f) if r and not r[0].strip().startswith("#")]
    if rows and not rows[0][0].strip().lstrip("-").replace(".", "").isdigit():
        rows = rows[1:]
    return rows


def main():
    d = Path(sys.argv[1])
    nodes = {int(float(r[0])): (float(r[4]), float(r[5])) for r in read(d / "nodes.csv")}
    elems = {int(float(r[0])): [int(float(x)) for x in r[1:9]] for r in read(d / "elems.csv")}
    gp = {}
    for r in read(d / "gauss.csv"):
        gp.setdefault(int(float(r[0])), []).append((float(r[2]), float(r[3])))
    worst = {"c": 0.0, "a": 0.0}
    worst_mean = {"c": 0.0, "a": 0.0}
    spread = {"c": 0.0, "a": 0.0}
    n_el = 0
    for e, con in elems.items():
        if e not in gp or len(gp[e]) != 8:
            continue
        for k, name in enumerate(("c", "a")):
            nv = np.array([nodes[n][k] for n in con])
            interp = np.sort(N_AT_GP @ nv)
            got = np.sort(np.array([v[k] for v in gp[e]]))
            worst[name] = max(worst[name], float(np.abs(got - interp).max()))
            worst_mean[name] = max(worst_mean[name], float(np.abs(got - nv.mean()).max()))
            spread[name] = max(spread[name], float(interp.max() - interp.min()))
        n_el += 1
    print(f"{n_el} elements with 8 Gauss-point values")
    for name in ("c", "a"):
        print(f"  {name}: max |GP - shape-function interpolation| = {worst[name]:.3e};  "
              f"max |GP - element mean of the nodes| = {worst_mean[name]:.3e};  largest in-element variation "
              f"{spread[name]:.3e}")
    ok = all(worst[n] <= 1e-10 * max(1.0, spread[n]) + 1e-12 for n in worst)
    print("PASS: the Gauss-point values are the nodal values interpolated with the shape functions" if ok else
          "CHECK: the Gauss-point values are not the shape-function interpolation (see the element-mean column)")


if __name__ == "__main__":
    main()

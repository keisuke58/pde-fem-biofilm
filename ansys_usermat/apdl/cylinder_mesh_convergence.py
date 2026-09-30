#!/usr/bin/env python3
"""cylinder_mesh_convergence.py --- compare t_growth_cylinder_shell*.dat
results across meshes at the same physical locations.

Each result file (growth_cylinder_result*.txt: PRNSOL,U + NLIST + growth-
layer PRESOL,S,COMP) is reduced to quantities that exist independently of the
mesh, then compared on a common grid:

  * outer surface (r=R_OUT), mid-length: radial displacement u_r(theta),
    interpolated onto theta = 0..60 deg;
  * growth-layer nodal-averaged stress at the outer surface and at the bonded
    interface (r=R_MID), mid-length: hoop sigma_tt(theta) and von Mises;
  * growth-layer von Mises over all element-nodes: mean and max.

Usage:
    python cylinder_mesh_convergence.py LABEL=path [LABEL=path ...]

List files coarse to fine; relative differences are reported against the
finest one. The result files are large and live in the F:-based working
directory, not in the repo.
"""
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_cylinder_bulge import parse as parse_u_and_coords  # noqa: E402

R_MID, R_OUT, LEN, ARC = 4.3, 4.4, 3.0, 60.0
THETA = np.linspace(0.0, ARC, 61)
_NUM = re.compile(r"-?\d+\.\d+(?:E[+-]\d+)?")


def parse_growth_layer_stress(text):
    """node -> mean (sx, sy, sz, sxy, syz, sxz) over the growth-layer elements
    that list it, plus every element-node von Mises value."""
    acc = defaultdict(lambda: np.zeros(6))
    cnt = defaultdict(int)
    vms = []
    in_block = False
    for line in text.splitlines():
        if "ELEMENT NODAL STRESS LISTING" in line:
            in_block = True
            continue
        if not in_block:
            continue
        if "PRINT" in line and "ELEMENT SOLUTION" in line:
            continue
        parts = line.split(None, 1)
        if len(parts) == 2 and parts[0].isdigit():
            nums = _NUM.findall(parts[1])
            if len(nums) == 6:
                s = np.array([float(v) for v in nums])
                nid = int(parts[0])
                acc[nid] += s
                cnt[nid] += 1
                vms.append(von_mises(s))
    return {n: acc[n] / cnt[n] for n in acc}, np.array(vms)


def von_mises(s):
    sx, sy, sz, sxy, syz, sxz = s
    return math.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)
                     + 3.0 * (sxy ** 2 + syz ** 2 + sxz ** 2))


def hoop(s, x, y):
    sx, sy, _sz, sxy, _syz, _sxz = s
    r = math.hypot(x, y)
    c, sn = x / r, y / r
    return sx * sn * sn + sy * c * c - 2.0 * sxy * sn * c


def ring(coord, r0, z0=LEN / 2, r_tol=1e-3, z_tol=1e-3):
    out = []
    for nid, (x, y, z) in coord.items():
        if abs(math.hypot(x, y) - r0) < r_tol and abs(z - z0) < z_tol:
            out.append((math.degrees(math.atan2(y, x)), nid))
    out.sort()
    return out


def on_grid(pairs):
    th = np.array([p[0] for p in pairs])
    v = np.array([p[1] for p in pairs])
    return np.interp(THETA, th, v)


def reduce(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    disp, coord = parse_u_and_coords(text)
    stress, vms = parse_growth_layer_stress(text)

    outer = ring(coord, R_OUT)
    inter = ring(coord, R_MID)
    ur = []
    for th, n in outer:
        x, y, _ = coord[n]
        ux, uy, _ = disp[n]
        ur.append((th, (ux * x + uy * y) / math.hypot(x, y)))

    def stress_ring(nodes):
        tt, vm = [], []
        for th, n in nodes:
            if n in stress:
                x, y, _ = coord[n]
                tt.append((th, hoop(stress[n], x, y)))
                vm.append((th, von_mises(stress[n])))
        return on_grid(tt), on_grid(vm), len(tt)

    tt_out, vm_out, n_out = stress_ring(outer)
    tt_int, vm_int, n_int = stress_ring(inter)
    return {
        "nodes": len(coord), "ring_nodes": (len(outer), n_out, n_int),
        "u_r": on_grid(ur), "tt_out": tt_out, "vm_out": vm_out,
        "tt_int": tt_int, "vm_int": vm_int,
        "vm_mean": float(vms.mean()), "vm_max": float(vms.max()),
    }


def rel_l2(a, ref):
    return float(np.linalg.norm(a - ref) / np.linalg.norm(ref))


def main():
    runs = []
    for arg in sys.argv[1:]:
        label, path = arg.split("=", 1)
        print(f"reading {label}: {path}", flush=True)
        runs.append((label, reduce(path)))
    if not runs:
        print(__doc__)
        sys.exit(1)

    ref = runs[-1][1]
    print(f"\n{'mesh':8s} {'nodes':>8s} {'max u_r':>12s} {'u_r@30':>12s} "
          f"{'tt_out@30':>12s} {'tt_int@30':>12s} {'vm mean':>12s} {'vm max':>12s}")
    i30 = int(np.argmin(abs(THETA - 30.0)))
    for label, r in runs:
        print(f"{label:8s} {r['nodes']:8d} {r['u_r'].max():12.5e} {r['u_r'][i30]:12.5e} "
              f"{r['tt_out'][i30]:12.5e} {r['tt_int'][i30]:12.5e} "
              f"{r['vm_mean']:12.5e} {r['vm_max']:12.5e}")
    print(f"\nrelative L2 difference of the theta-profiles vs {runs[-1][0]}:")
    print(f"{'mesh':8s} {'u_r':>9s} {'tt_out':>9s} {'vm_out':>9s} {'tt_int':>9s} {'vm_int':>9s}")
    for label, r in runs[:-1]:
        print(f"{label:8s} " + " ".join(f"{rel_l2(r[k], ref[k]):9.2e}"
              for k in ("u_r", "tt_out", "vm_out", "tt_int", "vm_int")))
    print("\nring node counts (outer, outer with stress, interface with stress):")
    for label, r in runs:
        print(f"  {label}: {r['ring_nodes']}")


if __name__ == "__main__":
    main()

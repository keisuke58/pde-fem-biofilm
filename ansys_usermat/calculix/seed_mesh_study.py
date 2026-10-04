#!/usr/bin/env python3
"""Seed stress against mesh and element type, in CalculiX.

The same mechanical problem as ansys_usermat/apdl/mesh_study_seed.py (Check 3
of the 5 Oct slides), solved by an independent finite element code:
  - the partner's 2 mm cube; the seed = the 32 elements of BIOFILM1 of the
    8^3 deck, as a region in space;
  - growth alpha - 1 = 1.1e-3 in the seed as an isotropic eigenstrain
    (*EXPANSION 1.1e-3, temperature 0 -> 1), none elsewhere;
  - stiffness E (phi^2 + f): 10.01 Pa in the seed, 0.01 Pa elsewhere, nu = 0.49;
  - only rigid-body motion constrained, at three corner nodes.
Element types: C3D8 (full 2x2x2 integration, no B-bar), C3D8I (incompatible
modes), C3D20R (quadratic, reduced). Reported as in mesh_study_seed.py: seed
von Mises (mean, max of the element averages), seed mean stress, largest von
Mises outside the seed. Small strain, linear elastic.

    python ansys_usermat/calculix/seed_mesh_study.py --ccx /path/to/ccx_2.21 \
        --n 8 16 32 --types C3D8 C3D8I C3D20R
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "apdl"))
sys.path.insert(0, str(HERE.parent))
import mesh_study_seed as S  # noqa: E402  (seed boxes from the partner's deck)

E, NU, FLOOR, GROWTH = S.E, S.NU, S.FLOOR, S.GROWTH
# local node order of C3D8 / C3D20 (CalculiX = Abaqus), on the doubled grid
HEX8 = [(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0), (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2)]
MID = [(1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0), (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
       (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1)]


def write_inp(path: Path, n: int, etype: str, boxes: np.ndarray) -> np.ndarray:
    quad = etype.startswith("C3D20")
    h = 2.0 / n
    m = 2 * n + 1                                          # doubled grid, spacing h/2
    lid = lambda i, j, k: (i * m + j) * m + k + 1
    loc = HEX8 + (MID if quad else [])
    conn, used = [], set()
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = [lid(2 * i + a, 2 * j + b, 2 * k + d) for a, b, d in loc]
                conn.append(c); used.update(c)
    ax = -1 + h * (np.arange(n) + 0.5)
    cx, cy, cz = np.meshgrid(ax, ax, ax, indexing="ij")
    cen = np.stack([cx.ravel(), cy.ravel(), cz.ravel()], 1)
    inseed = np.zeros(n ** 3, bool)
    for b in boxes:
        inseed |= np.all(np.abs(cen - b) < 0.125 + 1e-9, axis=1)
    corner = lambda p: lid(*[0 if c < 0 else 2 * n for c in p])
    with open(path, "w") as f:
        f.write("*NODE,NSET=NALL\n")
        for i in range(m):
            for j in range(m):
                for k in range(m):
                    if lid(i, j, k) in used:
                        f.write(f"{lid(i, j, k)},{-1 + i * h / 2:.10g},{-1 + j * h / 2:.10g},{-1 + k * h / 2:.10g}\n")
        for name, sel in (("SEED", inseed), ("VOID", ~inseed)):
            f.write(f"*ELEMENT,TYPE={etype},ELSET={name}\n")
            for e in np.flatnonzero(sel):
                c = conn[e]
                f.write(f"{e + 1}," + ",".join(map(str, c[:15])) + ("," if len(c) > 15 else "") + "\n")
                if len(c) > 15:
                    f.write(",".join(map(str, c[15:])) + "\n")
        for name, Ee, a in (("SEED", E * (1 + FLOOR), GROWTH), ("VOID", E * FLOOR, 0.0)):
            f.write(f"*MATERIAL,NAME=M{name}\n*ELASTIC\n{Ee:.10g},{NU}\n*EXPANSION,ZERO=0.\n{a:.10g}\n"
                    f"*SOLID SECTION,ELSET={name},MATERIAL=M{name}\n")
        c0, c1, c2 = (corner(p) for p in S.CORNERS)
        f.write("*INITIAL CONDITIONS,TYPE=TEMPERATURE\nNALL,0.\n*STEP\n*STATIC\n*BOUNDARY\n"
                f"{c0},1,3,0.\n{c1},2,3,0.\n{c2},3,3,0.\n*TEMPERATURE\nNALL,1.\n"
                "*EL PRINT,ELSET=SEED\nS\n*EL PRINT,ELSET=VOID\nS\n*END STEP\n")
    return inseed


def read_stress(dat: Path, ne: int) -> np.ndarray:
    """Element averages of the integration-point stresses (Cauchy, Voigt order of CalculiX)."""
    acc = np.zeros((ne, 6)); cnt = np.zeros(ne)
    for line in open(dat):
        p = line.split()
        if len(p) == 8 and p[0].isdigit() and p[1].isdigit():
            e = int(p[0]) - 1
            acc[e] += [float(x) for x in p[2:]]; cnt[e] += 1
    return acc / cnt[:, None]


def run(ccx: str, n: int, etype: str, boxes, work: Path):
    job = work / f"seed_{etype}_{n}"
    inseed = write_inp(job.with_suffix(".inp"), n, etype, boxes)
    t0 = time.time()
    subprocess.run([ccx, job.name], cwd=work, check=True, stdout=subprocess.DEVNULL,
                   env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    sig = read_stress(job.with_suffix(".dat"), n ** 3)
    p = sig[:, :3].mean(1)
    s = sig[:, :3] - p[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return inseed, vm, p, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccx", required=True)
    ap.add_argument("--n", type=int, nargs="+", default=[8, 16, 32])
    ap.add_argument("--types", nargs="+", default=["C3D8", "C3D8I", "C3D20R"])
    ap.add_argument("--work", default=None, help="directory for the jobs (default: a temporary one)")
    args = ap.parse_args()
    boxes = S.seed_boxes()
    work = Path(args.work or tempfile.mkdtemp(prefix="ccx_seed_"))
    work.mkdir(parents=True, exist_ok=True)
    print(f"seed = {len(boxes)} elements of the 8^3 deck; growth {GROWTH:g}, E = {E:g} Pa, nu = {NU}")
    print(f"{'type':>7} {'mesh':>6} {'seed vM mean':>13} {'seed vM max':>12} {'seed p mean':>12}"
          f" {'outside vM max':>15}   [Pa]")
    for etype in args.types:
        for n in args.n:
            ins, vm, p, dt = run(args.ccx, n, etype, boxes, work)
            print(f"{etype:>7} {n:>4}^3 {vm[ins].mean():13.3e} {vm[ins].max():12.3e} {p[ins].mean():12.3e}"
                  f" {vm[~ins].max():15.3e}   ({dt:.0f} s)", flush=True)


if __name__ == "__main__":
    main()

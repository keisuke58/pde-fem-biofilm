#!/usr/bin/env python3
"""Stress of the step-2 fields (ROADMAP_TWO_WAY.md) in CalculiX.

Step 2 lets the composition change how fast the biofilm spreads. This script
asks whether that change reaches the stress. It takes the fields that
ansys_usermat/two_way_step2.py --save writes (phi and alpha - 1 = k_alpha int
phi dt at T* = 1, Klempt 2024 test case 4.1, 21^3 nodes on the 20 um cube) and
solves the mechanics once per case and d:
  - 20^3 C3D8I elements on the grid of the field (1 um);
  - growth alpha - 1 as an isotropic eigenstrain: nodal temperature = alpha - 1,
    expansion coefficient 1;
  - stiffness E (phi^2 + f), E = 10 Pa, nu = 0.49, f = 1e-3 (Klempt 2024 Eq. 20
    weighting and Table 2), phi of the element = mean of its 8 nodes, binned to
    steps of 0.01;
  - only rigid-body motion constrained (three corner nodes). The paper gives no
    mechanical boundary condition for case 4.1; this is an assumption.
Small strain, linear elastic. Reported: von Mises and mean stress, averaged over
the biofilm (elements with phi >= 0.5) and as the 99th percentile.

    python ansys_usermat/calculix/two_way_stress.py --ccx ccx_2.21 --npz s.npz
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np

E, NU, FLOOR = 10.0, 0.49, 1e-3
HEX8 = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def write_inp(path: Path, phi: np.ndarray, am1: np.ndarray, L: float = 20.0) -> np.ndarray:
    N = phi.shape[0]; n = N - 1; h = L / n
    lid = lambda i, j, k: (i * N + j) * N + k + 1
    conn = np.array([[lid(i + a, j + b, k + c) for a, b, c in HEX8]
                     for i in range(n) for j in range(n) for k in range(n)])
    pe = phi.ravel()[conn - 1].mean(1)
    lev = np.round(np.clip(pe, 0, 1) * 100).astype(int)
    with open(path, "w") as f:
        f.write("*NODE,NSET=NALL\n")
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    f.write(f"{lid(i, j, k)},{i * h:g},{j * h:g},{k * h:g}\n")
        for l in np.unique(lev):
            f.write(f"*ELEMENT,TYPE=C3D8I,ELSET=P{l}\n")
            for e in np.flatnonzero(lev == l):
                f.write(f"{e + 1}," + ",".join(map(str, conn[e])) + "\n")
            Ee = E * ((l / 100) ** 2 + FLOOR)
            f.write(f"*MATERIAL,NAME=M{l}\n*ELASTIC\n{Ee:.10g},{NU}\n*EXPANSION,ZERO=0.\n1.\n"
                    f"*SOLID SECTION,ELSET=P{l},MATERIAL=M{l}\n")
        f.write("*ELSET,ELSET=EALL\n" + "\n".join(f"P{l}" for l in np.unique(lev)) + "\n")
        f.write("*INITIAL CONDITIONS,TYPE=TEMPERATURE\nNALL,0.\n*STEP\n*STATIC\n*BOUNDARY\n"
                f"{lid(0, 0, 0)},1,3,0.\n{lid(n, 0, 0)},2,3,0.\n{lid(0, n, 0)},3,3,0.\n*TEMPERATURE\n")
        a = am1.ravel()
        for i in range(N ** 3):
            f.write(f"{i + 1},{a[i]:.10g}\n")
        f.write("*EL PRINT,ELSET=EALL\nS\n*END STEP\n")
    return pe


def stresses(dat: Path, ne: int):
    acc = np.zeros((ne, 6)); cnt = np.zeros(ne)
    for line in open(dat):
        p = line.split()
        if len(p) == 8 and p[0].isdigit() and p[1].isdigit():
            e = int(p[0]) - 1
            acc[e] += [float(x) for x in p[2:]]; cnt[e] += 1
    sig = acc / cnt[:, None]
    pm = sig[:, :3].mean(1)
    s = sig[:, :3] - pm[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm, pm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccx", required=True)
    ap.add_argument("--npz", required=True, help="output of two_way_step2.py --save")
    ap.add_argument("--work", default=None)
    args = ap.parse_args()
    z = np.load(args.npz)
    work = Path(args.work or tempfile.mkdtemp(prefix="ccx_tw_"))
    work.mkdir(parents=True, exist_ok=True)
    runs = sorted({tuple(k.split("|")[:2]) for k in z.files if k.endswith("|alpha_m1")})
    print(f"{'case':>10} {'d':>6} {'alpha-1 max':>11} {'biofilm vol':>11} {'vM mean':>10} {'vM p99':>10}"
          f" {'p mean':>11} {'p p1':>11}   [Pa], biofilm = phi >= 0.5")
    for c, d in runs:
        phi, am1 = z[f"{c}|{d}|phi"], z[f"{c}|{d}|alpha_m1"]
        job = work / f"tw_{c}_d{round(float(d) * 1000):04d}"
        pe = write_inp(job.with_suffix(".inp"), phi, am1)
        subprocess.run([args.ccx, job.name], cwd=work, check=True, stdout=subprocess.DEVNULL)
        vm, pm = stresses(job.with_suffix(".dat"), pe.size)
        b = pe >= 0.5
        print(f"{c:>10} {float(d):6.3f} {am1.max():11.3e} {b.mean():11.3f} {vm[b].mean():10.3e}"
              f" {np.percentile(vm[b], 99):10.3e} {pm[b].mean():11.3e} {np.percentile(pm[b], 1):11.3e}",
              flush=True)


if __name__ == "__main__":
    main()

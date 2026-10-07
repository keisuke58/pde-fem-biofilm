#!/usr/bin/env python3
"""Run the one-species coupling against Klempt 2024's own PDE, without ANSYS.

The PDE (``klempt2024_quantitative.py``) evolves phi, c and -- by Eq. 36 --
alpha on a 3D grid. This script records, at a few grid points, the phi
history the PDE produces and the alpha it integrates from it, then feeds the
same phi history through the Fortran routine meant for the partner's element
(``growth_from_phi.f`` -> ``BIOFILM_GROWTH_VISCO_V01``, via
``crosscheck/one_species_driver.f``) and checks that the two agree.

Two independent implementations, one shared variable: this is the cross-check
ONE_SPECIES_COUPLING.md said the one-species route makes possible. It is also
the answer an ANSYS run of the same history has to reproduce.

The stress printed is the material-point stress at a fully constrained point
(F = I): what the routine returns for that alpha, not an FE equilibrium.

    python JAXFEM/one_species_pde_crosscheck.py [--case fig7_high] [--t-end 1.0]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "JAXFEM"))
sys.path.insert(0, str(_ROOT / "ansys_usermat"))
sys.path.insert(0, str(_ROOT / "ansys_usermat" / "apdl"))

import klempt2024_quantitative as kq        # noqa: E402
import closed_form_reference as cf          # noqa: E402

AU = _ROOT / "ansys_usermat"
E_BIO, NU_BIO, E_VOID, NU_VOID = 1000.0, 0.30, 1.0, 0.30


def pde_histories(case, points, dt=kq.DT, t_end=kq.T_END, variant="printed"):
    """Mirror kq.run's loop step for step, recording phi and alpha at points."""
    phi, mask, g = kq.setup(case)
    alpha = np.ones_like(phi)
    c = kq.solve_c(phi, g, mask, variant)
    n = int(round(t_end / dt))
    hist = {p: {"phi": [], "alpha_K": []} for p in points}
    for _ in range(n):
        gx, gy, gz = kq.grad_c(c)
        mag = np.sqrt(gx**2 + gy**2 + gz**2)
        speed = np.where(mag > 1e-14,
                         kq.R * c / (kq.K_M + c) / np.maximum(mag, 1e-14), 0.0)
        v = (speed * gx, speed * gy, speed * gz)
        src = -kq.upwind_dot(phi, v)
        phi = phi + dt * (kq.BETA * kq.lap(phi) + kq.K_A * alpha + src)
        phi = np.clip(phi, 0.0, 1.0)
        alpha = alpha + dt * kq.K_A * phi
        c = kq.solve_c(phi, g, mask, variant)
        for p in points:
            hist[p]["phi"].append(float(phi[p]))
            hist[p]["alpha_K"].append(float(alpha[p]))
    return hist


def build_driver(tmp: Path) -> Path:
    fc = shutil.which("gfortran")
    cc = shutil.which("cc") or shutil.which("gcc")
    if fc is None or cc is None:
        raise RuntimeError("gfortran and a C compiler are needed")
    o = {n: tmp / f"{n}.o" for n in ("hook", "core", "shim")}
    subprocess.run([fc, "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(AU / "coupling" / "usermat_py_hook.f"), "-o",
                    str(o["hook"])], check=True, cwd=tmp)
    subprocess.run([fc, "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(AU / "usermat_biofilm.f"), "-o", str(o["core"])],
                   check=True, cwd=tmp)
    subprocess.run([cc, "-c", "-fPIC", str(AU / "coupling" / "biofilm_py_eval.c"),
                    "-o", str(o["shim"])], check=True)
    exe = tmp / "one_species"
    subprocess.run([fc, "-ffixed-line-length-132", "-I", str(tmp),
                    str(AU / "crosscheck" / "one_species_driver.f"),
                    str(AU / "growth_from_phi.f"), str(AU / "biofilm_material_v01.f"),
                    str(o["core"]), str(o["hook"]), str(o["shim"]), "-o", str(exe)],
                   check=True)
    return exe


def run_fortran(exe, phi, k_alpha, dt):
    stdin = ("1 0 0 0 1 0 0 0 1\n"
             f"{E_BIO:.17e} {E_VOID:.17e} {NU_BIO:.17e} {NU_VOID:.17e} "
             f"1.0 {k_alpha:.17e} {dt:.17e} {len(phi)}\n"
             + "\n".join(f"{p:.17e}" for p in phi) + "\n")
    r = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True,
                       check=True)
    t = r.stdout.split()
    n = len(phi)
    return np.array([float(x) for x in t[:n]]), np.array([float(x) for x in t[n:n + 6]])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="fig7_high", choices=list(kq.PAPER))
    ap.add_argument("--t-end", type=float, default=kq.T_END)
    a = ap.parse_args(argv)

    m = kq.N // 2
    # colony seed (z = H, centre), a point above it, one at the edge, one far
    pts = [(m, m, 1), (m, m, 4), (m + 3, m, 1), (kq.N - 2, kq.N - 2, kq.N - 2)]
    labels = ["seed centre", "above seed", "seed edge", "far corner"]
    print(f"case {a.case}, N={kq.N}, dt={kq.DT}, t_end={a.t_end}, "
          f"k_alpha={kq.K_A}")
    hist = pde_histories(a.case, pts, t_end=a.t_end)

    exe = build_driver(Path(tempfile.mkdtemp()))
    D1 = 2.0 / (E_BIO / (3.0 * (1.0 - 2.0 * NU_BIO)))
    C10 = 0.5 * E_BIO / (2.0 * (1.0 + NU_BIO))
    worst = 0.0
    print(f"\n{'point':<12} {'phi end':>9} {'alpha PDE':>12} {'alpha Fortran':>14}"
          f" {'|diff|':>9} {'sigma11 (Pa)':>13} {'closed form':>12}")
    for lab, p in zip(labels, pts):
        h = hist[p]
        a_f, sig = run_fortran(exe, h["phi"], kq.K_A, kq.DT)
        a_pde = np.array(h["alpha_K"]) - 1.0          # Klempt -> this repo
        d = float(np.max(np.abs(a_f - a_pde)))
        worst = max(worst, d)
        cfs = (cf.constrained_stress(a_f[-1], D1) + cf.spurious_term(a_f[-1], C10))[0]
        print(f"{lab:<12} {h['phi'][-1]:9.4f} {a_pde[-1]:12.6e} {a_f[-1]:14.6e}"
              f" {d:9.1e} {sig[0]:13.6e} {cfs:12.6e}")
    print(f"\nlargest alpha difference over every step and point: {worst:.2e}")
    ok = worst < 1e-12
    print("PASS: the Fortran routine reproduces the PDE's Eq. 36" if ok
          else "FAIL: the two implementations of Eq. 36 disagree")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

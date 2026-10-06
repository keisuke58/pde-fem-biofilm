"""ecology_native.f (the point model in Fortran, no material server) against
ecology_jax.ecology_substeps, the server's computation: the same states,
cases, sub-steps and local nutrient, compared to round-off.

Needs gfortran and jax; skipped without them.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CP = ROOT / "ansys_usermat" / "coupling"
sys.path[:0] = [str(CP), str(ROOT / "ansys_usermat"), str(ROOT)]

DRIVER = """
      program drv
      use biofilm_py_bridge
      implicit none
      integer :: nact, nsub, hc, st
      double precision :: cs, als, eta(5), th(20), g(12), dt, crel
      double precision :: gn(12), pint
      logical :: ok, hcl
      do
        read(*, *, iostat=st) nact, cs, als, eta, hc, th, g, dt, nsub,
     &                        crel
        if (st .ne. 0) exit
        hcl = hc .ne. 0
        call eco_native_set(nact, cs, als, eta, th, hcl)
        call biofilm_ecology_hook_c(g, th, dt, nsub, crel, gn, pint, ok)
        write(*, '(13ES25.16E3, I3)') gn, pint, merge(1, 0, ok)
      end do
      end program drv
"""

pytest.importorskip("jax")


@pytest.fixture(scope="module")
def driver(tmp_path_factory):
    gf = shutil.which("gfortran")
    if gf is None:
        pytest.skip("gfortran not available")
    d = tmp_path_factory.mktemp("econ")
    (d / "drv.f").write_text(DRIVER)
    exe = d / ("drv.exe" if os.name == "nt" else "drv")
    subprocess.run([gf, "-O2", "-ffixed-line-length-132", str(CP / "ecology_native.f"), str(d / "drv.f"),
                    "-o", str(exe)], cwd=d, check=True)
    return exe


def cases():
    import material_server as ms
    import ecology_jax as ej
    out = []
    rng = np.random.default_rng(1)
    for name in ("2sp_case6", "2sp_case3", "4sp_case1", None):
        if name:
            ms.set_case(name)
            th, hp, n = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"], ms.ECOLOGY_ACTIVE
        else:
            ms.set_case(None)
            ms.set_active_species(5)          # set_case(None) keeps the last n
            th = [float(x) for x in np.asarray(__import__("jax_hamilton_0d_5species_demo").THETA_DEMO)]
            hp, n = None, 5
        for k in range(4):
            g = np.array(ej.default_initial_state(), dtype=float)
            if k:
                phi = rng.uniform(0.02, 0.3, 5)
                g[0:5] = phi
                g[5] = 1.0 - phi.sum()
                g[6:11] = rng.uniform(0.3, 0.999, 5)
                g[11] = rng.uniform(-1, 1)
            # (0.002, 1): with THETA_DEMO and c* = 25 one step of 0.01 from the
            # default state is beyond what 6 Newton iterations converge (both
            # versions end at the clip bounds, in round-off-dependent places)
            dt, nsub = [(0.002, 1), (0.025, 5), (0.1, 20), (0.01, 3)][k]
            crel = [-1.0, 1.0, 0.35, 0.0][k]
            out.append((name, n, hp, th, g, dt, nsub, crel))
    return out


def reference(name, n, hp, th, g, dt, nsub, crel):
    import ecology_jax as ej
    h = hp
    if crel >= 0:
        base = hp if hp is not None else {"c": ej.C_STAR, "alpha": ej.ALPHA_STAR, "eta": [1.0] * 5}
        h = dict(base, c=float(base["c"]) * crel)
    gn, pint = ej.ecology_substeps(g, th, dt, nsub, n, h)
    return np.asarray(gn, dtype=float), float(pint)


def test_native_equals_server_point_model(driver):
    import ecology_jax as ej
    rows, refs = [], []
    for name, n, hp, th, g, dt, nsub, crel in cases():
        cs = hp["c"] if hp else ej.C_STAR
        als = hp["alpha"] if hp else ej.ALPHA_STAR
        eta = hp["eta"] if hp else [1.0] * 5
        vals = [n, cs, als, *eta, 1 if hp else 0, *th, *g, dt, nsub, crel]
        rows.append(" ".join(repr(float(v)) if not isinstance(v, int) else str(v) for v in vals))
        refs.append(reference(name, n, hp, th, g, dt, nsub, crel))
    out = subprocess.run([str(driver)], input="\n".join(rows) + "\n", capture_output=True, text=True,
                         check=True).stdout.split("\n")
    worst = 0.0
    for line, (gn, pint) in zip(out, refs):
        v = [float(x) for x in line.split()]
        assert int(v[13]) == 1
        got = np.array(v[:12])
        err = np.max(np.abs(got - gn) / np.maximum(np.abs(gn), 1e-3))
        err = max(err, abs(v[12] - pint) / max(abs(pint), 1e-12))
        worst = max(worst, err)
    assert worst < 1e-10, worst
    print("largest relative difference", worst)

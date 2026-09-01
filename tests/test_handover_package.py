"""The package handed to the partner group must build on its own and answer
identically to the sources it was extracted from.

Two things this guards. First, that the extraction really is self-contained --
the constitutive core is lifted out of usermat_biofilm.f, whose `usermat`
routine `use`s the Python-coupling bridge, so it would be easy for the package
to quietly acquire a dependency on the coupling layer or the C shim. Here it is
compiled with nothing else present.

Second, that it has not drifted. The package is generated rather than kept as a
copy, but a generator can still fall behind its input; comparing outputs at
zero tolerance is what makes "extracted, not rewritten" a checked claim rather
than a comment.

Requires gfortran; skipped where absent.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_MAKE = _ROOT / "handover" / "make_handover.py"
_AU = _ROOT / "ansys_usermat"
_FC = shutil.which("gfortran")
_CC = shutil.which("cc") or shutil.which("gcc")

pytestmark = pytest.mark.skipif(
    _FC is None or _CC is None or not _MAKE.exists(),
    reason="gfortran/cc or the handover generator is unavailable")

F_TEST = np.array([[1.06, 0.02, 0.0], [0.0, 0.97, 0.01], [0.0, 0.0, 0.98]])
I3 = np.eye(3)

# (biofilm, growth, eta, dt, c01r, mtype) -- spans both material paths, the
# elastic and viscous branches, and a cut-back.
CASES = [
    (1.0, 0.20, 5.0, 1.0e-4, 0.15, 1.0),
    (0.0, 0.20, 5.0, 1.0e-4, 0.15, 1.0),
    (0.5, 0.00, 0.0, 1.0e-2, 0.00, 0.0),
    (1.0, 0.35, 0.05, 1.0e-4, 0.30, 1.0),
]


def _stdin(F, Fv, E, EL, nu, nuL, biofilm, growth, eta, dt, c01r, mtype):
    return (" ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
            " ".join(f"{Fv[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
            f"{E:.17e} {EL:.17e} {nu:.17e} {nuL:.17e} {biofilm:.17e} "
            f"{growth:.17e} {eta:.17e} {dt:.17e} {c01r:.17e} {mtype:.1f}\n")


def _run(exe, text):
    r = subprocess.run([str(exe)], input=text, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"}, timeout=30)
    assert r.returncode == 0, f"{exe} failed: {r.stderr}"
    return np.array([float(x) for x in r.stdout.split()])


@pytest.fixture(scope="module")
def exes():
    """Build the package standalone, and the same routine from the repo."""
    tmp = Path(tempfile.mkdtemp())

    pkg = tmp / "pkg"
    subprocess.run([sys.executable, str(_MAKE), "-o", str(pkg)], check=True,
                   capture_output=True)
    # Compiled with ONLY what the package ships. If the extraction had pulled
    # in a dependency, this link fails.
    pkg_exe = tmp / "pkg_wrap"
    subprocess.run([_FC, "-ffixed-line-length-132",
                    str(pkg / "wrapper_driver.f"),
                    str(pkg / "biofilm_material_v01.f"),
                    str(pkg / "biofilm_stress_core.f"),
                    "-o", str(pkg_exe)], check=True, cwd=tmp)

    # The in-repo route, which does need the bridge module and the C shim.
    o = {n: tmp / f"{n}.o" for n in ("hook", "core", "shim")}
    subprocess.run([_FC, "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(_AU / "coupling" / "usermat_py_hook.f"),
                    "-o", str(o["hook"])], check=True, cwd=tmp)
    subprocess.run([_FC, "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_AU / "usermat_biofilm.f"), "-o", str(o["core"])],
                   check=True, cwd=tmp)
    subprocess.run([_CC, "-c", "-fPIC",
                    str(_AU / "coupling" / "biofilm_py_eval.c"),
                    "-o", str(o["shim"])], check=True)
    repo_exe = tmp / "repo_wrap"
    subprocess.run([_FC, "-ffixed-line-length-132", "-I", str(tmp),
                    str(_AU / "crosscheck" / "wrapper_driver.f"),
                    str(_AU / "biofilm_material_v01.f"),
                    str(o["core"]), str(o["hook"]), str(o["shim"]),
                    "-o", str(repo_exe)], check=True)
    return pkg_exe, repo_exe, pkg


def test_package_builds_with_nothing_but_its_own_files(exes):
    """The fixture proves it by linking; this states it as a test so the
    intent survives a refactor of the fixture."""
    pkg_exe, _, pkg = exes
    assert pkg_exe.exists()
    shipped = {p.name for p in pkg.iterdir() if p.suffix == ".f"}
    assert shipped == {"biofilm_material_v01.f", "biofilm_stress_core.f",
                       "wrapper_driver.f"}, shipped


@pytest.mark.parametrize("case", CASES)
def test_package_answers_identically_to_the_repo(exes, case):
    pkg_exe, repo_exe, _ = exes
    biofilm, growth, eta, dt, c01r, mtype = case
    text = _stdin(F_TEST, I3, 1000.0, 1.0, 0.30, 0.30,
                  biofilm, growth, eta, dt, c01r, mtype)
    np.testing.assert_allclose(_run(pkg_exe, text), _run(repo_exe, text),
                               rtol=0.0, atol=0.0)


def test_package_honours_the_cutback_contract(exes):
    pkg_exe, repo_exe, _ = exes
    text = _stdin(np.diag([1e-3, 1e-3, 1e-3]), I3, 1000.0, 1.0, 0.30, 0.30,
                  1.0, 2.0, 5.0, 1.0e-4, 0.15, 1.0)
    out = _run(pkg_exe, text)
    assert int(out[15]) == 1
    np.testing.assert_allclose(out, _run(repo_exe, text), rtol=0.0, atol=0.0)


def test_the_extracted_core_is_verbatim(exes):
    """Extracted, not rewritten -- so every line must appear in the source."""
    _, _, pkg = exes
    src = (_AU / "usermat_biofilm.f").read_text()
    body = pkg.joinpath("biofilm_stress_core.f").read_text()
    body = body[body.index("      subroutine BIOFILM_STRESS_CORE"):]
    assert body.rstrip() in src, "the extracted core is not verbatim"

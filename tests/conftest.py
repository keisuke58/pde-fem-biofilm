"""Shared fixtures of the test suite.

The tests of the user element (test_uel_*.py) run in one process and share two things that outlive a test module:

- the environment variable BIOFILM_ECO_CASE, the file the Fortran point model reads. Other tests set it for their
  own case and leave it set (ecology_native_ctypes.init, used by test_numerics_nov4.py and the coupled Python code,
  points it at a two-species case), and most element tests only set it when it is empty;
- the library p3/uel/libunified_uel.so, whose point model is set once per process: from the file at the first call
  (eco_native_init), or by a test through driver.setup (test_uel_kinematics.py sets theta = 0), and whose switches
  (start rule of a new Gauss point, implicit gate, kinematics of the C entry points) keep their last value.
  test_uel_init_rule.py loads a second copy of it and leaves the Python modules bound to that copy in sys.modules.

So the result of an element test depended on which tests ran before it. The fixture below gives every test_uel_*
module the same start: the element's own point-model file (make_column_inp.eco_file) in BIOFILM_ECO_CASE, the
library built if it is missing or older than its sources, the file read again and the switches back at their
defaults, in every loaded copy. The variable is put back afterwards.
"""
import ctypes as ct
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
UEL = ROOT / "p3" / "uel"
ECO_SRC = ROOT / "ansys_usermat" / "coupling" / "ecology_native.f"


def _build_uel_library():
    """the build of the `lib` fixtures of test_uel_*.py, done here first so the library is never rebuilt after it
    has been loaded"""
    so = UEL / "libunified_uel.so"
    srcs = (UEL / "unified_uel.f90", ECO_SRC)
    if not so.exists() or max(x.stat().st_mtime for x in srcs) > so.stat().st_mtime:
        subprocess.run(["gfortran", "-O2", "-fPIC", "-c", str(ECO_SRC), "-o", "ecology_native.o"], cwd=UEL, check=True)
        subprocess.run(["gfortran", "-O2", "-fPIC", "-ffree-line-length-none", "-c", "unified_uel.f90", "-o", "unified_uel.o"],
                       cwd=UEL, check=True)
        subprocess.run(["gfortran", "-shared", "-o", "libunified_uel.so", "unified_uel.o", "ecology_native.o"], cwd=UEL, check=True)
    return so


def _reset_uel_library(lib):
    """the library as a new process has it: point model from BIOFILM_ECO_CASE, switches at their defaults"""
    # the module procedure that reads the file (gfortran's name for biofilm_py_bridge::eco_native_init); the file
    # of the element tests sets everything a test can have changed (theta, c*, alpha*, eta, phi_init, newton)
    lib.__biofilm_py_bridge_MOD_eco_native_init()
    lib.uel_set_init_rule(ct.c_int(0), ct.c_double(0.0))        # default start of a new Gauss point, registry cleared
    lib.uel_set_gate_law(ct.c_double(0.0), ct.c_int(0))         # no implicit gate, gate = 1
    lib.uel_set_finite(ct.c_int(1), ct.c_int(0))                # finite strain, no mean dilatation
    lib.uel_set_current(ct.c_int(0))                            # nutrient in the reference configuration


@pytest.fixture(scope="module", autouse=True)
def uel_isolated(request, tmp_path_factory):
    name = request.module.__name__.rpartition(".")[2]           # tests is a package: tests.test_uel_...
    if not name.startswith("test_uel_") or shutil.which("gfortran") is None:
        yield
        return
    before = os.environ.get("BIOFILM_ECO_CASE")
    so = _build_uel_library()
    if str(UEL) not in sys.path:
        sys.path.insert(0, str(UEL))
    from make_column_inp import eco_file        # imports driver.py, which loads the library: after the build
    case = tmp_path_factory.mktemp("uel_eco") / "uel_eco.txt"
    eco_file(case)
    os.environ["BIOFILM_ECO_CASE"] = str(case)
    # the library of p3/uel, and the copy test_uel_init_rule.py builds in a scratch directory: its fixture leaves
    # its own driver and emulate_abaqus (bound to that copy) in sys.modules for the modules that run after it
    libs = {str(so): ct.CDLL(str(so))}
    for mod in ("driver", "emulate_abaqus"):
        lib = getattr(sys.modules.get(mod), "lib", None)
        if lib is not None:
            libs.setdefault(lib._name, lib)
    for lib in libs.values():
        _reset_uel_library(lib)
    yield
    if before is None:
        os.environ.pop("BIOFILM_ECO_CASE", None)
    else:
        os.environ["BIOFILM_ECO_CASE"] = before

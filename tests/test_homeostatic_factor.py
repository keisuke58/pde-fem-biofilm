"""BIOFILM_HOMEOSTATIC_FACTOR (growth_from_phi.f): the three forms of the
compression factor on the growth increment, called from a small driver."""
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

_GFP = Path(__file__).resolve().parents[1] / "ansys_usermat" / "growth_from_phi.f"
_FC = shutil.which("gfortran")
pytestmark = pytest.mark.skipif(_FC is None, reason="gfortran unavailable")

DRIVER = """\
      program drive
      implicit none
      double precision p, pref, phi, f, g
      integer iform
      read(*,*) p, pref, iform, phi, f
      call BIOFILM_HOMEOSTATIC_FACTOR(p, pref, iform, phi, f, g)
      print '(ES24.16E3)', g
      end
"""


@pytest.fixture(scope="module")
def exe():
    tmp = Path(tempfile.mkdtemp())
    (tmp / "drive.f").write_text(DRIVER)
    out = tmp / "drive"
    subprocess.run([_FC, "-ffixed-form", str(tmp / "drive.f"), str(_GFP),
                    "-o", str(out)], check=True)
    return out


def g(exe, p, pref, iform, phi=0.5, f=1e-3):
    r = subprocess.run([str(exe)], input=f"{p} {pref} {iform} {phi} {f}\n",
                       capture_output=True, text=True, check=True)
    return float(r.stdout)


def test_off_is_one(exe):
    for form in (0, 1, 2):
        assert g(exe, 5.0, 0.0, form) == 1.0


def test_linear_forms(exe):
    assert abs(g(exe, 0.5, 2.0, 0) - 0.75) < 1e-15
    assert g(exe, 3.0, 2.0, 0) == 0.0
    assert g(exe, -1.0, 2.0, 0) == 1.0
    ph = 2.0 * (0.25 + 1e-3)
    assert abs(g(exe, 0.1, 2.0, 1) - (1.0 - 0.1 / ph)) < 1e-15


def test_hyperbolic_form(exe):
    assert abs(g(exe, 2.4e-10, 2.4e-10, 2) - 0.5) < 1e-15
    assert abs(g(exe, 2.4e-9, 2.4e-10, 2) - 1.0 / 11.0) < 1e-15
    assert g(exe, -1.0, 2.4e-10, 2) == 1.0          # tension: unchanged

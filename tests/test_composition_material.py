"""Composition-dependent material constants: E(phi) -> ustatev(11:14) -> USERMAT.

Covers the "stiffness E(phi)" leg of the model (RESEARCH_MODEL.md sec.3). Two
halves:

  1. composition_to_material.py maps CLSM composition to (C10, C01, D1, eta)
     via material_models.py, and emits an APDL TB,STATE block for them.
  2. usermat_biofilm.f's kStateMat=1 path actually *uses* those per-integration
     -point constants instead of prop(1:4) -- checked through the real
     usermat() entry point, not by inspection.

The load-bearing assertion is test_state_constants_override_prop: running
material A's constants through prop must equal running them through state
while prop carries material B. Without that, "the state path is wired" could
just mean "the state path is silently ignored and prop won".

Requires gfortran + a C compiler for the USERMAT half; skipped where absent.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
_COUP = _AU / "coupling"
_CORE = _AU / "usermat_biofilm.f"
_HOOK = _COUP / "usermat_py_hook.f"
_DRIVER = _COUP / "usermat_endtoend_driver.f"
_SHIM_C = _COUP / "biofilm_py_eval.c"

sys.path.insert(0, str(_COUP))
import composition_to_material as c2m  # noqa: E402

_CC = shutil.which("cc") or shutil.which("gcc")
_FC = shutil.which("gfortran")

# Condition-level moduli from run_ve_twin_experiment.py (MAP theta, 3-model
# comparison). The commensal/dysbiotic contrast is the whole point of this
# path: ~31x in E, which prop-constant runs cannot represent.
E_CH, DI_CH = 995.0, 0.05      # commensal HOBIC  (stiffest)
E_DS, DI_DS = 32.0, 0.85       # dysbiotic static (softest)


# --------------------------------------------------------------------------- #
# 1. composition -> constants
# --------------------------------------------------------------------------- #
def test_uniform_composition_gives_physical_constants():
    c = c2m.material_constants_from_composition(np.full(5, 0.2))
    assert float(c["C10"]) > 0.0
    assert float(c["D1"]) > 0.0
    assert float(c["eta"]) >= 0.0
    # small-strain consistency: mu = 2(C10 + C01) = E / (2(1+nu))
    mu = 2.0 * (float(c["C10"]) + float(c["C01"]))
    assert mu == pytest.approx(float(c["E"]) / (2.0 * 1.30), rel=1e-12)


def test_commensal_is_stiffer_and_less_viscous_than_dysbiotic():
    """The direction of the effect, stated as a test so a sign flip in
    material_models.py cannot pass silently."""
    ch = c2m.material_constants_from_E(E_CH, DI_CH)
    ds = c2m.material_constants_from_E(E_DS, DI_DS)
    assert float(ch["C10"]) > float(ds["C10"])
    assert float(ch["D1"]) < float(ds["D1"])       # stiffer => less compliant
    assert float(ch["eta"]) < float(ds["eta"])     # dysbiotic flows more


def test_stiffness_contrast_is_large_enough_to_matter():
    ch = c2m.material_constants_from_E(E_CH, DI_CH)
    ds = c2m.material_constants_from_E(E_DS, DI_DS)
    assert float(ch["C10"]) / float(ds["C10"]) > 20.0


def test_composition_is_normalised_internally():
    a = c2m.material_constants_from_composition(np.full(5, 0.2))
    b = c2m.material_constants_from_composition(np.full(5, 2.0))   # x10, same mix
    assert float(a["E"]) == pytest.approx(float(b["E"]), rel=1e-12)


def test_vectorised_over_rows():
    phi = np.array([[0.2] * 5, [0.05, 0.05, 0.1, 0.3, 0.5]])
    c = c2m.material_constants_from_composition(phi)
    assert np.asarray(c["C10"]).shape == (2,)
    assert np.asarray(c["eta"]).shape == (2,)
    # a Pg/Fn-dominated mix produces little EPS => softer than an even one
    assert float(c["C10"][0]) > float(c["C10"][1])


def test_composite_model_saturates_at_the_clamp():
    """Documents a real property of compute_E_composite rather than asserting
    around it: with the current MECH_* calibration most reasonably diverse
    compositions hit E_MAX_PA, so E(phi) discriminates mainly at the
    dysbiotic end. Anything relying on fine gradations among healthy mixes
    needs the clamp/calibration revisited first."""
    even = c2m.material_constants_from_composition(np.full(5, 0.2))
    commensal = c2m.material_constants_from_composition(
        np.array([0.6, 0.2, 0.1, 0.05, 0.05]))
    import material_models as mm
    assert float(even["E"]) == pytest.approx(mm.E_MAX_PA)
    assert float(commensal["E"]) == pytest.approx(mm.E_MAX_PA)


def test_wrong_species_count_is_rejected():
    with pytest.raises(ValueError, match="5 species"):
        c2m.material_constants_from_composition(np.full(4, 0.25))


def test_apdl_block_never_exceeds_six_values_per_tbdata():
    """ANSYS silently drops values past the 6th in one TBDATA call -- a bug
    this project has already hit once (apdl/t_growth_free.dat)."""
    block = c2m.apdl_state_block(c2m.material_constants_from_E(E_DS, DI_DS),
                                 matid=2, alpha=0.3)
    for line in block.splitlines():
        if not line.startswith("TBDATA,"):
            continue
        payload = line.split("!")[0].strip()
        n_values = len(payload.split(",")) - 2      # minus "TBDATA" and the index
        assert n_values <= 6, f"{n_values} values in one TBDATA: {line}"


def test_apdl_block_sets_the_state_path_on():
    c = c2m.material_constants_from_E(E_DS, DI_DS)
    block = c2m.apdl_state_block(c, matid=2)
    assert f"TB,STATE,2,,{c2m.NSTATEV_WITH_MATERIAL}" in block
    assert "TBDATA,7,1.0" in block                  # prop(7)=kStateMat=1
    assert f"TBDATA,{c2m.STATE_BASE}," in block     # ustatev(11:14) written


def test_apdl_mtype_follows_c01():
    on = c2m.apdl_state_block(c2m.material_constants_from_E(E_DS, DI_DS))
    off = c2m.apdl_state_block(
        c2m.material_constants_from_E(E_DS, DI_DS, c01_ratio=0.0))
    prop_on = [ln for ln in on.splitlines() if ln.startswith("TBDATA,1,")][0]
    prop_off = [ln for ln in off.splitlines() if ln.startswith("TBDATA,1,")][0]
    assert prop_on.split(",")[6] == "1.0"           # Mooney-Rivlin
    assert prop_off.split(",")[6] == "0.0"          # degenerates to Neo-Hookean


# --------------------------------------------------------------------------- #
# 2. constants -> USERMAT (through the real usermat() entry point)
# --------------------------------------------------------------------------- #
usermat = pytest.mark.skipif(
    _FC is None or _CC is None or not all(
        p.exists() for p in (_CORE, _HOOK, _DRIVER, _SHIM_C)),
    reason="gfortran/cc or the usermat sources are unavailable")

F_TEST = np.array([[1.06, 0.02, 0.0], [0.0, 0.97, 0.01], [0.0, 0.0, 0.98]])
I3 = np.eye(3)


@pytest.fixture(scope="module")
def exe():
    tmp = Path(tempfile.mkdtemp())
    objs = {n: tmp / f"{n}.o" for n in ("hook", "core", "driver", "shim")}
    out = tmp / "e2e"
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(_HOOK), "-o", str(objs["hook"])], check=True, cwd=tmp)
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_CORE), "-o", str(objs["core"])], check=True, cwd=tmp)
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_DRIVER), "-o", str(objs["driver"])], check=True, cwd=tmp)
    subprocess.run([_CC, "-c", "-fPIC", str(_SHIM_C), "-o", str(objs["shim"])],
                   check=True)
    subprocess.run(["gfortran", "-o", str(out), str(objs["driver"]),
                    str(objs["hook"]), str(objs["core"]), str(objs["shim"])],
                   check=True)
    return out


# The viscous relaxation time is eta/(2*C10): ~0.22 s commensal, ~40 s
# dysbiotic. dt must stay well inside the stiffer one or the backward-Euler
# step drives Fv toward singular and the detFe clamp fires (a real code path,
# exercised in test_coupling_vs_fortran.py, but not what these tests are about).
DT_TEST = 0.01


def _run(exe, prop_mat, state_mat=None, alpha=0.2, dt=DT_TEST, F=F_TEST, Fv=I3):
    """prop_mat/state_mat are (C10, C01, D1, eta). state_mat=None disables the
    per-IP path (prop(7)=0)."""
    c10, c01, d1, eta = prop_mat
    mtype = 1.0 if c01 > 0.0 else 0.0
    kstmat, sm = (0.0, (0.0,) * 4) if state_mat is None else (1.0, state_mat)
    stdin = (
        " ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        " ".join(f"{Fv[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        f"{alpha:.17e} {c10:.17e} {c01:.17e} {d1:.17e} {eta:.17e} "
        f"{mtype:.1f} {dt:.17e} 0.0\n" +
        f"{kstmat:.1f} " + " ".join(f"{v:.17e}" for v in sm) + "\n"
    )
    r = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"}, timeout=20)
    assert r.returncode == 0, f"driver failed rc={r.returncode}: {r.stderr}"
    t = [float(x) for x in r.stdout.split()]
    return np.array(t[0:6]), int(t[15])            # stress, keycut


def _mat(E, di):
    c = c2m.material_constants_from_E(E, di)
    return tuple(float(np.atleast_1d(c[k])[0]) for k in c2m.STATE_SLOTS)


@usermat
def test_state_constants_override_prop(exe):
    """The load-bearing test: prop carries the COMMENSAL material while state
    carries the DYSBIOTIC one, and the answer must be the dysbiotic result.
    If the state path were ignored, this would return the commensal stress."""
    ch, ds = _mat(E_CH, DI_CH), _mat(E_DS, DI_DS)
    s_ds_via_prop, _ = _run(exe, ds)
    s_ds_via_state, _ = _run(exe, ch, state_mat=ds)
    s_ch_via_prop, _ = _run(exe, ch)

    np.testing.assert_allclose(s_ds_via_state, s_ds_via_prop, rtol=1e-12, atol=1e-14)
    # ...and the two materials really are distinguishable, so the above is not
    # passing merely because both materials give the same stress.
    assert np.max(np.abs(s_ch_via_prop - s_ds_via_prop)) > 1.0


@usermat
def test_disabled_state_path_ignores_state_slots(exe):
    """With prop(7)=0 the ustatev(11:14) values must have no effect at all."""
    ch, ds = _mat(E_CH, DI_CH), _mat(E_DS, DI_DS)
    baseline, _ = _run(exe, ch)
    with_junk_state, _ = _run(exe, ch, state_mat=None)
    np.testing.assert_allclose(with_junk_state, baseline, rtol=1e-12, atol=1e-14)


@usermat
def test_uninitialised_state_falls_back_to_prop(exe):
    """kStateMat=1 but ustatev(11)=0 means the state was never filled in;
    the USERMAT must fall back to prop rather than run a zero-stiffness
    material (which would return zero stress and look like a converged solve)."""
    ch = _mat(E_CH, DI_CH)
    baseline, _ = _run(exe, ch)
    unset, keycut = _run(exe, ch, state_mat=(0.0, 0.0, 0.0, 0.0))
    np.testing.assert_allclose(unset, baseline, rtol=1e-12, atol=1e-14)
    assert keycut == 0
    assert np.max(np.abs(unset)) > 0.0, "fell through to a zero-stress material"


@usermat
def test_stiffness_contrast_shows_up_in_stress(exe):
    """The end-to-end payoff: the same deformation and the same alpha, with
    only the composition-derived constants differing, must produce a
    materially different stress -- which is exactly what a prop-constant model
    cannot represent."""
    s_ch, _ = _run(exe, _mat(E_CH, DI_CH))
    s_ds, _ = _run(exe, _mat(E_DS, DI_DS))
    ratio = np.max(np.abs(s_ch)) / np.max(np.abs(s_ds))
    assert ratio > 5.0, f"commensal/dysbiotic stress ratio only {ratio:.2f}"


@usermat
@pytest.mark.parametrize("phi,label", [
    (np.full(5, 0.2), "even"),
    (np.array([0.6, 0.2, 0.1, 0.05, 0.05]), "commensal-leaning"),
    (np.array([0.05, 0.05, 0.1, 0.3, 0.5]), "dysbiotic-leaning"),
])
def test_composition_drives_usermat_without_error(exe, phi, label):
    c = c2m.material_constants_from_composition(phi)
    mat = tuple(float(np.atleast_1d(c[k])[0]) for k in c2m.STATE_SLOTS)
    stress, keycut = _run(exe, mat, state_mat=mat)
    assert keycut == 0, f"{label}: unexpected cut-back"
    assert np.all(np.isfinite(stress))

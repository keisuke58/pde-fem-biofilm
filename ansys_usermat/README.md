# ANSYS USERMAT port — biofilm growth / viscoelastic law

**English** | [日本語](README.ja.md)

**✅ Verified on ANSYS MAPDL 2022 R2 (v222).** The `USERMAT` interface arguments (including `var0`, `var1..var8`, `tsstif`, `epsZZ`), the `keycut`/`cutFactor` adaptive-substepping controls, and the `dsdePl` material Jacobian were all validated in-solver via a 3-D solid-element benchmark (`SOLID185`, `NLGEOM,ON`, uniaxial tension) that runs and converges.

`usermat_biofilm.f` is an **ANSYS Mechanical APDL `USERMAT`** port of the verified
Abaqus UMAT (`umat_biofilm_visco.f` / `umat_biofilm_visco_phase2.f`). It lets
Felix/IKM's existing ANSYS FE model call the *same* growth/viscoelastic
constitutive law (`F = Fe·Fv·Fg`, `Fg=(1+α)I`) at each Gauss point in place of a
phenomenological material — the starting point for the proposed thesis.

The constitutive **algebra is a line-by-line mirror** of the verified Abaqus
core (Neo-Hookean deviator + `D1` pressure, backward-Euler viscous update,
F-perturbation consistent tangent). Only the **interface** is ANSYS-specific.

## Abaqus UMAT ↔ ANSYS USERMAT mapping

| Abaqus | ANSYS `usermat` | note |
|---|---|---|
| `DFGRD1` / `DFGRD0` | `defGrad` / `defGrad_t` | 3×3 deformation gradient |
| `STRESS(NTENS)` | `stress(ncomp)` | Cauchy stress |
| `DDSDDE` | `dsdePl(ncomp,ncomp)` | material Jacobian |
| `STATEV` | `ustatev(nStatev)` | state |
| `PROPS` | `prop(nProp)` | properties |
| `DTIME` | `dTime` | increment |
| `PNEWDT < 1` | `keycut = 1` (+ `cutFactor`) | cut-back request |
| `SSE` / `SPD` | `sedEl` / `sedPl` | energies |
| order `11,22,33,12,13,23` | order **`11,22,33,12,23,13`** | ⚠️ shear slots 5↔6 swapped |

The stress-component order difference is the classic porting trap; it is handled
here by the `VI/VJ` Voigt map (`data VI /1,2,3,1,2,1/`, `VJ /1,2,3,2,3,3/`).

## Properties & state

```
prop(1)=C10  prop(2)=C01  prop(3)=D1  prop(4)=eta  prop(5)=mtype  prop(6)=kUsePy
ustatev(1:9)=Fv (row-major 3×3)   ustatev(10)=alpha (growth driver)
```

- **Growth driver `alpha`** comes from the JAXFEM α-field mapped to each
  integration point (initialised via `TB,STATE` / a user field, or evolved).
- `kUsePy=1` selects the **Python material hook** (see below) instead of the
  inline Fortran law.

## Python material hook (per Gauss point)

The thesis' core deliverable — calling the paper's calibrated **Python** material
model at each Gauss point — is **wired and confirmed working end to end**
through the `PYTHON MATERIAL HOOK` in the source. Mechanism: an
`ISO_C_BINDING` / local-socket bridge that ships `(defGrad, Fv_old, alpha,
dTime, prop)` to Python and receives `(stress, Fv_new, dsdePl)` back into
`usermat()`'s own `stress`/`ustatev`/`dsdePl` outputs, with the
Abaqus↔ANSYS Voigt reindex (`MAP6`) applied on the way in. `kUsePy=1` runs
this path; if the Python server is unreachable or returns something invalid,
`usermat()` falls back to the verified inline core rather than failing the
solve. (Architecture: `ch5_flow/flow_python_material_hook`.)

The bridge lives in [`coupling/`](coupling/README.md): the Python side
(`material_server.py`, NumPy core + F-perturbation tangent + socket server),
the wire protocol, and the Fortran hook (`usermat_py_hook.f`). It is exercised
through the *real* `usermat()` entry point (not a bypass driver) by
`tests/test_usermat_kusepy_e2e.py`, which compiles `usermat_biofilm.f` +
`usermat_py_hook.f` + `biofilm_py_eval.c` into a standalone driver and checks
`kUsePy=1` against `kUsePy=0` across elastic/viscous/Mooney-Rivlin cases —
stress and the updated viscous state agree to numerical precision, and the
consistent tangent (`dsdePl`) agrees to floating-point-noise precision (both
sides use the same F-perturbation scheme, `PERT=1e-7`). That test also covers
the fallback path (server unreachable → `PYOK=.false.` → inline core, not a
crash).

This closes a real bug found while adding that dsdePl comparison: the wire
carries the Python side's 6×6 Jacobian as a row-major (NumPy/C-order) flatten,
but Fortran's `RESHAPE` fills column-major, so a plain `reshape(d36,[6,6])`
in `usermat_py_hook.f` silently returned the tangent's **transpose**. It was
invisible for near-symmetric elastic cases and only showed up (large,
sign-flipping discrepancies) once viscous/Mooney-Rivlin cases were checked —
fixed with an explicit `transpose(...)`.

## Build & use in ANSYS (outline)

```apdl
! compile & link with the ANSYS user-programmable-features toolchain
! (ANSUSERSHARED / usermat build), then in the model:
TB, USER, 1, 1, 6         ! 6 properties
TBDATA, 1, C10, C01, D1, eta, mtype, kUsePy
TB, STATE, 1, , 10        ! 10 state variables (Fv 1:9, alpha 10)
```

For local syntax checking without ANSYS:

```bash
gfortran -c -fsyntax-only -ffixed-line-length-132 usermat_biofilm.f
```

## Verification status

- ✅ **Compiles** clean with `gfortran` (`-fsyntax-only`, no warnings).
- ✅ **Bit-identical to the verified Abaqus UMAT** across 20 deformation states
  (`|Δσ|=|ΔFv|=|ΔJe|=0`) — see `crosscheck/` (compiles both real Fortran cores
  and compares). Also: isotropic growth patch; consistent tangent vs central
  difference **2.97e-8**; ANSYS shear ordering (`s12,s23,s13`) confirmed.
- ✅ **Runs and converges in ANSYS MAPDL 2022 R2 (v222).** The interface
  arguments (`var0`, `var1..var8`, `tsstif`, `epsZZ`), the `keycut`/`cutFactor`
  adaptive-substepping path, and the `dsdePl` material Jacobian were validated
  in-solver on a `SOLID185` uniaxial-tension benchmark with `NLGEOM,ON`.
  **The argument list is release-specific** — re-check it first when moving to
  another ANSYS version (`var0` sits between `coords` and `defGrad_t`;
  `var1..var8` trail `cutFactor`).
- The Abaqus core it mirrors is verified (tangent vs FD ~2.4e-8; patch tests
  13/13 in `phase2_patch_test.py`).

## Caveats / next steps

1. ~~Confirm the exact `usermat` argument list for the target ANSYS release.~~
   **Done for 2022 R2 (v222).** Redo this check for any other release.
2. ~~Confirm ANSYS expects the `dsdePl` convention produced by the
   F-perturbation.~~ **Confirmed by the converging `SOLID185` benchmark** — a
   wrong Jacobian convention shows up as poor or failed Newton convergence, and
   the run converges.
3. ~~Growth verification in-solver.~~ **Done 2026-08-19/20.** All four
   constrained-single-element cases (elastic/viscous × α=0.05/0.20, `F = I`,
   no FE solve needed to predict the answer) matched the closed form, plus a
   `KEYOPT(1,2)` ∈ {0,1,3} sweep (B-bar/enhanced/simplified enhanced strain)
   with no volumetric locking detected. Details and measured values:
   [`apdl/README.md`](apdl/README.md), [`apdl/RUNBOOK.md`](apdl/RUNBOOK.md).
   - Side finding: a preliminary cylinder-shell (two-layer, tooth-surface
     proxy) constrained-growth check converges at α=0.01 but not α=0.015, and
     even at the converged α the outer surface displacement is **not**
     uniform — a **two-lobe pattern**
     ([`assets/growth_cylinder_bulge.png`](../assets/growth_cylinder_bulge.png),
     `apdl/extract_cylinder_bulge.py`). Consistent with early buckling, not
     confirmed against a true eigenvalue analysis — flagged as interesting,
     not yet a settled conclusion.
4. ~~Wire the `PYTHON MATERIAL HOOK` (ISO_C_BINDING/socket).~~ **Done** — see
   above; verified end to end through the real `usermat()` entry point via
   `tests/test_usermat_kusepy_e2e.py`. The inline core stays as the
   verification reference and the automatic fallback. The wire protocol and
   per-call payload are diagrammed in
   [`ch5_flow/flow_python_material_hook`](../ch5_flow/README.md).
5. The Python-side material model is still the NumPy mirror of the verified
   Fortran law (`material_server.py`'s `stress_core`), not yet the calibrated
   JAX model (`JAXFEM/`/`material_models.py`) — swapping that in is a local
   change behind the same interface, not a wiring change.

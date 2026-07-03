# ANSYS USERMAT port — biofilm growth / viscoelastic law

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
model at each Gauss point — has an explicit extension point marked
`PYTHON MATERIAL HOOK` in the source. Intended mechanism: an
`ISO_C_BINDING` / local-socket bridge that ships `(defGrad, Fv_old, alpha,
dTime, prop)` to Python and receives `(stress, Fv_new, dsdePl)`. The inline
Fortran core is the reference/fallback used for verification. (Architecture:
`ch5_flow/flow_impl_architecture`.)

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
- ⚠️ **Not yet run inside ANSYS.** The interface (arg list, `keycut`/`cutFactor`,
  `dsdePl` convention) follows the standard `usermat` spec but should be checked
  against the target ANSYS version; the per-version argument list can vary
  slightly (e.g. `var1..var8`, `tsstif`, `epsZZ`).
- The Abaqus core it mirrors is verified (tangent vs FD ~2.4e-8; patch tests
  13/13 in `phase2_patch_test.py`).

## Caveats / next steps

1. Confirm the exact `usermat` argument list for the target ANSYS release.
2. `dsdePl` here is `∂σ/∂ε`-style from an F-perturbation; confirm ANSYS expects
   the Cauchy/Jaumann material Jacobian in the same convention (may need
   symmetrisation or the `NLGEOM` rotation `rotateM`).
3. Wire the `PYTHON MATERIAL HOOK` (ISO_C_BINDING/socket) — the actual thesis
   work; the inline core stays as the verification reference.

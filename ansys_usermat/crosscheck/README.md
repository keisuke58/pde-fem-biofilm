# Abaqus UMAT ↔ ANSYS USERMAT cross-check

Proves the ANSYS `USERMAT` constitutive core reproduces the **verified Abaqus
UMAT** law to **machine precision** — the confidence needed before handing the
port to Felix and for the thesis verification chapter.

It compiles the two **actual Fortran cores** (not re-implementations) and drives
them over a battery of deformation states, reconstructing the full 3×3 Cauchy
stress tensor from each (accounting for the different Voigt orderings) and
comparing stress, viscous state `Fv`, and `Je`.

| file | role |
|---|---|
| `xcheck_driver_abq.f` | reads `F, Fv, (α,C10,C01,D1,η,mtype,dt)` from stdin, calls `BIOFILM_STRESS_CORE` from `umat_biofilm_visco.f` (Abaqus order 11,22,33,12,13,23) |
| `xcheck_driver_ans.f` | same, calls the core from `usermat_biofilm.f` (ANSYS order 11,22,33,12,23,13) |
| `crosscheck.py` | builds both (gfortran), runs 20 cases, reconstructs tensors, asserts equal |

## Run

```bash
python ansys_usermat/crosscheck/crosscheck.py        # report
python -m pytest ansys_usermat/crosscheck/crosscheck.py
```

Requires `gfortran` + `numpy`. No Abaqus/ANSYS needed.

## Result

**20/20 cases bit-identical (`|Δσ|=|ΔFv|=|ΔJe|=0`)** — identity+growth, uniaxial,
shear, elastic (`η=0`), frozen (`η→∞`), large growth, Neo-Hookean, Mooney-Rivlin,
and 12 random finite-strain states.

> This harness caught a real bug: the initial ANSYS port applied the
> Mooney-Rivlin `C01` term only to the final Cauchy stress, not to the viscous
> flow driver (as the Abaqus reference does), so `Fv` evolved differently for
> `mtype=1, η>0`. Fixed — now bit-identical. (Since the top-level tangent is a
> finite-difference of this same core, core equivalence implies tangent
> equivalence.)

# What the viscous update actually is

Short version: the repository calls the `Fv` update **"backward Euler"** in
several places. The flow increment is evaluated at the **old** state, so the
scheme is explicit in the flow direction. The label matters because "backward
Euler" implies unconditional stability, and this update has a step limit —
which is exactly the constraint recorded in `biofilm_material_v01.f` and handed
to the partner group.

## The update

Both implementations of the growth law do, per Gauss point:

```
  Fe_trial = F·Fg⁻¹·Fv_n⁻¹          (trial elastic state, from Fv at t_n)
  τ        = τ_dev(Fe_trial)         (flow driver, evaluated at the OLD state)
  Fv_n+1   = ( I + Δt/(2·η·J_e)·τ )·Fv_n
  σ        = σ(F·Fg⁻¹·Fv_n+1⁻¹)      (stress re-evaluated at the NEW state)
```

The re-evaluation of the stress at `Fv_n+1` is presumably what earned it the
name: the *stress* is taken at the end of the step. But `Fv_n+1` itself is
obtained from a driver computed entirely at `t_n`, with no iteration and no
residual. That is a forward-Euler flow increment, not a backward-Euler solve.

A genuinely implicit update would require `τ` at `Fv_n+1`, i.e. solving

```
  Fv_n+1 = ( I + Δt/(2·η·J_e(Fv_n+1))·τ(Fv_n+1) )·Fv_n
```

which nothing here does.

## What follows from it

The step must resolve the relaxation time `τ_relax = η/(2·C10)`. Measured on
the ANSYS core at `η = 5`, `C10 ≈ 167` (`τ_relax ≈ 0.015 s`):

| Δt | Δt/τ_relax | σ₁₁ |
|---|---|---|
| 1e−4 | 0.007 | −513 |
| 5e−3 | 0.33 | −134 |
| 1e−2 | 0.67 | **+347** |
| 2e−2 | 1.34 | +3288 |

The stress crosses zero near `Δt/τ ≈ 0.5` and diverges past `Δt/τ ≈ 1` — the
signature of an explicit increment taken past the timescale it is resolving, and
not something an unconditionally stable scheme does. Pinned by
`tests/test_material_wrapper.py::test_the_viscous_step_must_resolve_the_relaxation_time`.

## Where the label appears

**Corrected in place** (the growth law on the thesis critical path):

- `umat_biofilm_visco.f` — the Abaqus UMAT
- `ansys_usermat/usermat_biofilm.f` — the ANSYS USERMAT, 0 ULP against it
- `ch5_flow/flow_python_material_hook.tex` — feeds a thesis figure

**Deliberately left alone:**

- `rigor_audit_growth_2026-06-26.md` — a dated audit record. Editing what a
  past audit said would be rewriting history; this note is the correction.
- Everything about the **φ/c diffusion solve** (`tooth_pde3d.py`,
  `fem_2d_extension.py`, `fem_report.tex`, `JAXFEM/`). Those genuinely are
  backward Euler — an implicit linear solve per step — and are unaffected.

**Unverified, and not on this thesis's critical path:**

- `umat_biofilm_visco_phase2.f` and `test_umat_2ch.py` implement a *different*
  two-channel Prony model. Its `test_large_dt_stability` asserts a coarse step
  (`Δt/τ₁ = 10`) stays finite and reaches the same long-time value as a fine
  one. That is consistent with the transient being wrong while the long-time
  limit is set by the equilibrium branch and therefore path-independent — but
  this has not been checked, so nothing here should be read as a claim about
  that model either way.

## For the write-up

Say what the scheme is rather than naming it: *a single-step update whose flow
increment is evaluated at the previous viscous state, with the stress
re-evaluated at the updated one*. That is accurate, and it makes the step
restriction follow rather than contradict.

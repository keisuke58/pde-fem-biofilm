# UserElement (USER300) — API notes and how it plugs into our USERMAT

Preparatory research for the open question in
[`THESIS_ASSIGNMENT.md`](../THESIS_ASSIGNMENT.md) §4.1: whether the spatial
biofilm-ecology field is (a) extra DOFs inside a `UserElem`, or (b) precomputed
externally and passed in as a state variable. Written before Felix's actual
implementation arrived, from the local ANSYS 2022 R2 install
(`C:\Program Files\ANSYS Inc\v222\ansys\customize\user\UserElem.F` — ships with
every install as a documented, working example; no external doc lookup needed).

## The element type

`USER300`, configured per-element via two APDL commands:

- `USRELEM` — declares shape: `NDIM`, `NNODES`, `NINTPNTS`, `NREAL`,
  `NSAVEVARS`, `KEYSHAPE`, and crucially **`KEYANSMAT`** (see below).
- `USRDOF` — declares which DOFs live on which nodes (not every node need
  carry every DOF — e.g. a mixed u-p element can put pressure only on corner
  nodes).

## The subroutine contract (`UserElem.F`)

One Fortran subroutine, called once per element per solve pass. Inputs include
nodal coordinates (`xRef`/`xCur`), current/incremental/iterative DOF values
(`TotValDofs`/`IncValDofs`/`ItrValDofs`), element key options, and a
persistent `saveVars` array (`nSaveVars` long — this is the UserElem's own
state storage, separate from and in addition to whatever a called material
routine keeps). It must return, depending on what `keyMtx(:)` requests: the
stiffness matrix `eStiff`, mass `eMass`, damping `eDamp`, internal/external
force vectors `fInt`/`fExt`, plus result quantities for `PRESOL`.

## The part that matters here: it already knows how to call our USERMAT

`KEYANSMAT` on `USRELEM` selects one of two modes inside the subroutine:

- `keyAnsMat = 1` — **"use standard ANSYS material"**: after computing the
  strain increment from the B-matrix (`CALL maxv(BMat, IncValDofs, IncStrain,
  ...)`), it calls `ElemGetMat(elId, matId, ..., IncStrain, defG0, defG, cMat,
  MatProp, Stress, Strain, ...)` — literally commented `USERMAT is called from
  here` (`UserElem.F:505`). `ElemGetMat` dispatches to whatever material is
  attached to `matId` via `TB` commands, `TB,STATE`/`TB,USER` included. **This
  is our existing, verified `usermat_biofilm.f`, unmodified.**
- `keyAnsMat = 0` — "make up your own material": the element inlines its own
  stress law instead of calling out.

So a `UserElem` that owns extra field DOFs (φ, α, ψ transported as nodal
unknowns — option (a) in the open question) does **not** require reimplementing
the constitutive law inside the element. The standard, sample-code-documented
pattern is:

1. `UserElem` computes the mechanical strain from the displacement DOFs as
   usual, and calls `ElemGetMat` → our `usermat_biofilm.f` for `Stress`/`cMat`
   exactly as today.
2. `UserElem` *separately* assembles the residual/stiffness contribution for
   the extra field DOFs (diffusion/reaction of the ecology state) — this part
   is new code, but it is additive to the element, not a rewrite of the
   material law.
3. `saveVars` on the `UserElem` side and the material's own state-variable
   array (`TB,STATE`, `nStatev`) are two independent stores; the element owns
   the field history, the material still owns `Fg`/`Fv`/`α` per Gauss point
   exactly as it does now.

## What this means for the open question

This weakens the "materially larger and more interesting piece of work"
framing in §4.1 somewhat: option (a) does not mean throwing away the
constitutive verification work (0 ULP vs Abaqus, `[HANDOFF.md](HANDOFF.md)`) —
that asset carries over unchanged as the thing `ElemGetMat` calls into. The
actual new work under (a) is the field-DOF residual/stiffness assembly (steady
diffusion-reaction on the element, roughly), which is a bounded, separately
testable piece — it can be prototyped and unit-tested against a manufactured
solution independent of the mechanical/material coupling.

**Not yet resolved / needs Oliver's version to confirm:** whether Felix's
actual `UserElem` uses this same `keyAnsMat=1` dispatch pattern, or does
something else (e.g. calls the material routine directly, bypassing
`ElemGetMat`, or stores the field unknowns differently). This note describes
the *documented sample* pattern, not Felix's code — diff against his source
once it arrives per §4.2.

## What the shipped example does NOT cover

Read the full 725 lines end to end (2026-08-19), not just the `ElemGetMat`
call. Correcting the framing above: this is **not** a template for the
diffusion-reaction / extra-field-DOF part of our problem. What's actually in
the file, section by section:

1. **Header/argument doc** (lines 1–211) — the full contract, as above.
2. **Setup** (~213–300) — declares work arrays, zeroes `eStiff`/`eMass`/
   `fExt`/`fInt` per what `keyMtx(:)` requests.
3. **Shape functions & B-matrix** (~300–500) — `ElemShpFn` gives isoparametric
   shape derivatives; these get assembled into a standard small-strain
   `BMat` (strain–displacement matrix) via the element Jacobian. This part is
   pure solid-mechanics kinematics — displacement DOFs only, `nUsrDof` here
   equals `nDim` per node, nothing else.
4. **Material call** (~500–533) — the `ElemGetMat`/`keyAnsMat` dispatch
   documented above.
5. **Stiffness/mass assembly** (~535–601) — `matba(BMat, cMat, eStiff, ...)`
   (the standard $\mathbf{B}^T\mathbf{C}\mathbf{B}$ integration), `ElemMass`
   for a consistent mass matrix. Nothing beyond textbook FE assembly.
6. **Result output** (~603–725) — extrapolate Gauss-point stress/strain to
   corner nodes for `PRESOL`/`OUTPR`; no physics here, just post-processing.

**Not present anywhere in this file:**
- Any non-mechanical field DOF (concentration/volume-fraction/whatever the
  ecology state would be) — `nUsrDof` in the whole example is sized purely
  off displacement components. `USRDOF` (the APDL command that would attach
  extra DOFs to nodes) is configured outside this file entirely and never
  referenced from it.
- Any diffusion or reaction term, any residual assembly beyond
  $\mathbf{f}_{int} = \int \mathbf{B}^T\boldsymbol{\sigma}\,dV$, and no
  boundary-condition handling — BCs (`D`, `F`, `SF` commands) are applied by
  ANSYS's standard DOF-constraint machinery *outside* `UserElem`; the element
  never sees a BC as such, only the DOF values ANSYS hands it each iteration
  (`TotValDofs`/`IncValDofs`).

**Conclusion for the open question (§4.1 of THESIS_ASSIGNMENT.md):** ANSYS's
own shipped example gives a solid, reusable template for "keep doing
mechanics as before, call out to the verified USERMAT for stress" — that part
is close to free. It gives **zero** template for the actual new piece under
option (a): a coupled residual that adds a diffusion-reaction equation for
the ecology field over the same element, with its own weak form, its own
contribution to `eStiff`/`fInt` (added to, not replacing, the mechanical
block), and its own entries in `saveVars`/`nSaveVars` for field-history state.
That has to be derived and coded from scratch (weak form → Galerkin residual
→ Newton tangent, the standard reaction-diffusion FE pattern used already in
`JAXFEM/`, just re-expressed in this Fortran contract) — it is **not** a
matter of finding the right ANSYS example to copy. Felix's version, once it
arrives, is the reference for how he actually did this same step; that is the
single most valuable thing to look at in his code, more than the material
dispatch (which we've already independently confirmed matches the documented
pattern).

## Source

`C:\Program Files\ANSYS Inc\v222\ansys\customize\user\UserElem.F` (ships with
every ANSYS 2022 R2 install; header comments + a complete worked example, 725
lines total, read in full). Read on IKMHIWI03, 2026-08-19.

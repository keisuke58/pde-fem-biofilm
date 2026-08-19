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

## Source

`C:\Program Files\ANSYS Inc\v222\ansys\customize\user\UserElem.F` (ships with
every ANSYS 2022 R2 install; header comments + a complete worked example are
in the file itself, ~700 lines). Read on IKMHIWI03, 2026-08-19.

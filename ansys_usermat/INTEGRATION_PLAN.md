# Integration plan — contributing the material law to Oliver's framework

Decided 2026-09-01. Background and evidence:
[`OLIVER_MODEL_NOTES.md`](OLIVER_MODEL_NOTES.md) ([日本語](OLIVER_MODEL_NOTES.ja.md)).

## The split

**Agreed with M. Soleimani:** their NEM part solves the biofilm field; **we
contribute the material law**. This is option (A) of the two laid out in the
notes, and it is settled — not an open question.

```
Oliver's framework                          this repo
─────────────────────────                   ──────────────────────────
USolBeg   parameters, NEM setup
usermat   per Gauss point ──── CALL ───▶    the growth + viscoelastic law
USSFin    field solve (PARDISO)             (Fg=(1+α)I, Mooney-Rivlin/D1,
                                             backward-Euler Fv)
```

**Platform: local port to ANSYS v222 on IKMHIWI03** — option (B). No cluster
account is being requested; steps in
[`apdl/V222_PORT_INSTRUCTIONS.md`](apdl/V222_PORT_INSTRUCTIONS.md).

## Why this combination is safer than it first looks

The obvious worry about "contribute into their framework, but build on a
different ANSYS release than they use" is that the deliverable ends up tied to
the wrong version. It does not, because of where the release-dependence
actually sits.

`AceGenNeoHookV04` — the call site our law would occupy — is **completely
ANSYS-independent**:

```fortran
      SUBROUTINE AceGenNeoHookV04(v, mDefGrad, vCauchy, mTangCC,
     &   sYoung, sYoungL, sNu, sNuL, sBiofilm, sAlpha, sElasticWork, sID)
      include 'sms.h'        ! AceGen's own runtime header, nothing else
```

No `usercm` common block, no `parevl`, no UPF calls: deformation gradient in,
Cauchy stress and tangent out. So:

- **the routine we deliver is release-independent.** Develop and verify it on
  v222 here, and it drops into their 2024 R2 build unchanged;
- **only the glue is release-specific** (`Usermat_*.F`, 41 vs 42 arguments) —
  and that is their file, not our deliverable.

## What this removes from the critical path

| dependency | status |
|---|---|
| Cluster account | **not needed** — building locally on v222 |
| AceGen Mathematica notebook | **probably not needed** — see below |
| Oliver's agreement on the interface | **the one real dependency**, and a small concrete ask |

On the notebook: their constitutive routines are AceGen-generated, so the
first assumption was that adding `Fg` meant editing a Mathematica notebook we
do not have. But `BIOFILM_STRESS_CORE` in [`usermat_biofilm.f`](usermat_biofilm.f)
is already hand-written and verified 0-ULP against the Abaqus UMAT
(`crosscheck/`). Wrapping it to match their calling convention needs no
AceGen at all — AceGen is their authoring tool, not a requirement of the
interface. Worth confirming with Oliver that a hand-written routine in the
pool is acceptable to them, which is in the draft to him.

## Work items

1. **Port the pool to v222 locally** so there is a test bed.
   `apdl/patch_usermat_to_v222.py` does the signature retarget; the rest is in
   `apdl/V222_PORT_INSTRUCTIONS.md`, including the pre-flight findings
   (`/fpp` on `userdata_*.f`, the integer-width question).
2. **Wrap `BIOFILM_STRESS_CORE` to their calling convention** — the same shape
   as `AceGenNeoHookV04`, plus arguments for the growth variable and the
   viscous state (`Vdp_Cv_n` already reserves a slot for the latter).
3. **Verify the wrapper** against the existing crosscheck battery, so the
   0-ULP Abaqus equivalence travels with it.
4. **Hand it over** and let them wire it at the `AceGenNeoHookV04` call site.

Steps 2 and 3 do not depend on step 1 succeeding — the routine is plain
Fortran and can be tested with the existing gfortran harness. Step 1 buys a
realistic test bed, not a prerequisite.

## Timeline note

Thesis submission is early December — about 13 weeks from this decision. The
plan is scoped so that the deliverable (steps 2–3) is under our own control
and testable here; only step 4 depends on someone else's schedule.

## Later, not now

Option (D) from the notes — borrowing the NEM idea to strengthen `JAXFEM/`'s
own derivative operator rather than using their code — is the intended
direction after returning to Keio. It is deliberately out of scope here: it
duplicates machinery that already exists on their side and would compete with
steps 1–4 for the same weeks.

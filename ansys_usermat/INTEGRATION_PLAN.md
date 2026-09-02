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

1. **Port the pool to v222 locally** so there is a test bed. **Done, 2026-09-02:**
   all 11 pool source files compile clean under v222 (see
   `apdl/V222_PORT_INSTRUCTIONS.md` §1.6) — `apdl/patch_usermat_to_v222.py`
   does the signature retarget; the rest is in
   `apdl/V222_PORT_INSTRUCTIONS.md`, including the pre-flight findings
   (`/fpp` on `userdata_*.f`, the integer-width question). Linking (a full
   custom `ANSYS.exe`) is not done yet — needs `ANSCUST.BAT`, which is
   interactive and needs a human at the console.

   **Also confirmed, same day: the actual deliverable compiles too, not just
   the stock pool.** `biofilm_material_v01.f` + `usermat_biofilm.f` (as
   `BIOFILM_STRESS_CORE`) + `usermat_py_hook.f` (the `biofilm_py_bridge`
   module it `use`s) all compile with zero errors under the identical v222
   flag set — the "release-independent, drops in unchanged" claim two
   paragraphs up is now backed by a real compile, not just an argument from
   the call signature. One real fixup needed along the way:
   `usermat_py_hook.f` is fixed-form Fortran with lines past column 72
   (the `bind(C, name=...)` declarations), which needs
   `/extend-source:132` — the same flag class as gfortran's
   `-ffixed-line-length-132` used elsewhere in this repo's syntax checks,
   just not needed for the other pool files, which stay under 72 columns.
2. **Wrap `BIOFILM_STRESS_CORE` to their calling convention** — the same shape
   as `AceGenNeoHookV04`, plus arguments for the growth variable and the
   viscous state (`Vdp_Cv_n` already reserves a slot for the latter).
   **Done:** [`biofilm_material_v01.f`](biofilm_material_v01.f),
   `BIOFILM_GROWTH_VISCO_V01`.
3. **Verify the wrapper** against the existing crosscheck battery, so the
   0-ULP Abaqus equivalence travels with it. **Done:**
   [`tests/test_material_wrapper.py`](../tests/test_material_wrapper.py), 13
   tests, driven through [`crosscheck/wrapper_driver.f`](crosscheck/wrapper_driver.f).
   The load-bearing one is `test_wrapper_is_only_an_adapter`: fed (E, ν), the
   wrapper's stress and viscous update must equal the core's fed the
   (C10, C01, D1) they map to, at `rtol=0, atol=0`. It does, so the core's
   0-ULP Abaqus equivalence applies to the wrapper unchanged.
4. **Hand it over** and let them wire it at the `AceGenNeoHookV04` call site.

### Two things step 3 pinned down, worth passing on with the routine

- **`sDt` must resolve the viscous relaxation time `η/(2·C10)`.** The core
  advances `Fv` with the flow increment evaluated at the old state, so the step
  is only accurate while `dt ≪ τ`. Measured at η=5, C10≈167 (τ≈0.015 s): σ₁₁
  goes from −513 Pa at dt=1e−4 to +347 Pa at dt=1e−2 — it crosses zero near
  `dt/τ ≈ 0.5` and diverges past `dt/τ ≈ 1`. This is a property of the existing
  verified core, not of the wrapper, and it is not something to fix here: the
  core is locked 0 ULP to the Abaqus UMAT, so changing its integrator is a
  separate decision carrying its own re-verification. It matters because in the
  handover case *their* framework picks the time step. `sEta = 0` selects the
  elastic path and removes the limit.
- **Cut-backs return defined outputs.** On `sKeyCut = 1` the stress and tangent
  are zeroed and `mFvN1` is restored to `mFvN` rather than left at the core's
  update off a collapsed configuration. Their framework passes work arrays that
  are not zeroed between calls, so "left untouched" would have meant returning
  uninitialised memory.

Note on terminology: the repo calls this update "backward Euler" throughout
(`umat_biofilm_visco.f`, `rigor_audit_growth_2026-06-26.md`), but the flow
increment is evaluated at `Fv_n`, which is what creates the step limit above.
The label was left alone rather than renamed across the six files that use it.

Steps 2 and 3 do not depend on step 1 succeeding — the routine is plain
Fortran and can be tested with the existing gfortran harness. Step 1 buys a
realistic test bed, not a prerequisite.

## Timeline note

**Submission is November 2026, defence December 2026** — 12 weeks from this
decision. (An earlier draft of this section said "early December submission";
that was wrong.) Full schedule: [`../ROADMAP_2026.md`](../ROADMAP_2026.md).

The plan is scoped so that the deliverable (steps 2–3) is under our own control
and testable here; only step 4 depends on someone else's schedule. **With the
November date, that makes step 1 more important than it looks above.** A
working local v222 build lets us produce results with the wrapper without step
4 landing at all, so it is the move that takes Oliver's calendar off the
critical path — see `ROADMAP_2026.md` §3.

## Later, not now

Option (D) from the notes — borrowing the NEM idea to strengthen `JAXFEM/`'s
own derivative operator rather than using their code — is the intended
direction after returning to Keio. It is deliberately out of scope here: it
duplicates machinery that already exists on their side and would compete with
steps 1–4 for the same weeks.

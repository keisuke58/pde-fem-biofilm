# n=2 -> n=3/4/5 growth-term trial (Oliver's ANSYS/NEM pipeline)

2026-09-07. Follow-up to `V222_PORT_INSTRUCTIONS.md` and the
`oliver_n2_vs_n5_species_gap` memory, which established that generalizing
Oliver's ANSYS pipeline from n=2 to n=5 species means rewriting `Ussfin`'s
PDE assembly, not just usermat glue. This is the first concrete, small,
incremental step in that direction: extend the explicit biofilm-growth
update (not the PARDISO-solved fields) to a 3rd species, mirroring
species 1/2's existing pattern exactly.

All edits are local to `F:\biofilm_upf_wired` (Oliver's source, never
committed to this repo). This file and `n3_growth_reference.py` are the
repo-side write-up/verification harness, matching the convention already
used for the ecology-bridge work (`ecology_4region_reference.py`, etc.).

> **⚠️ Correction, 2026-09-08 — the interaction term this whole file
> generalizes is not the published paper's equation.** Everything below
> treats `Interaction12/21` (and the `InteractionIJ` generalization) as "the
> paper's interaction scheme," per `OLIVER_MODEL_NOTES.md`'s 2026-09-01
> reading. That identification is wrong, found by reading the actual
> equations (Eq. 10, 16–18) of Klempt/Geisler/Soleimani/Junker,
> *A continuum multi-species biofilm model with a novel interaction scheme*
> (arXiv:2509.01274 / AAM 96, 164 (2026)) directly:
>
> - **The paper's scheme is additive/bilinear**: `Ia_i = Σ_j a_ij·φ̄_j =
>   (A·φ̄)_i` (φ̄=φψ) enters a dissipation-derived residual as
>   `-c*·ψ_i·Ia_i` (Eq. 16). `A` is stated symmetric, with diagonal
>   self-terms `a_ii` included, and the nutrient `c*` is a single scalar
>   entering linearly.
> - **Oliver's scheme is multiplicative**: `(MaxGrowth_I + Interaction_IJ·Bio_J)`
>   shifts the coefficient inside a per-nutrient Monod saturation curve, two
>   independent nutrient fields, no diagonal self-term, no symmetry
>   constraint (`Interaction12`/`21` are separate named constants).
>
> These are two different functional forms, not a notational difference —
> see `OLIVER_MODEL_NOTES.md`'s matching 2026-09-08 correction for the full
> comparison and for why this repo's own `hamilton_ode_jax.py`/
> `ecology_jax.py` (not this file's `InteractionIJ` work) is the one that
> actually matches the paper term-by-term.
>
> **Consequence for this file's own results, re-verified 2026-09-08**: the
> critical-coupling bisection below (0.0005 stable / 0.0007 diverges,
> Oliver's baseline `MaxGrowth=100/dt=0.1/Penalty1=5/HalfVelo=0.1`) was
> re-run on real ANSYS on 2026-09-08 and reproduced exactly (0 errors at
> 0.0005; `L-2 norm of the residual force overflowed` at 0.0007) — **it is a
> real, reproducible property**, but it characterizes Oliver's own
> non-paper-literal growth-law extension, not "the paper's critical coupling
> strength." Read every "the paper's interaction scheme" phrase below with
> that caveat; the bisection/numerics themselves are unaffected by the
> correction, only their scientific interpretation.

## Why this was tractable without touching PARDISO

Confirmed by reading `Ussfin_P21-V21_Conection_Test.F` in full: Bio1/Bio2
are **not** part of the PARDISO-solved system (only T/Nut1/Nut2 are). They
are updated by a plain explicit (forward-Euler-style) per-Gauss-point
formula, computed and written directly in the same `!$OMP DO` loop that
computes gradients/Laplacians for everything else. Adding a 3rd species
means duplicating that explicit block, not touching the sparse assembly at
all.

The `userdata` arena (`InitVals`/`GetVals`/`SetVals`) is also not a fixed-
size struct — `USolBeg` computes a running `ofs` offset counter and calls
`InitVals(ofs-1)` once at the end, so adding 4 new slots (`bio3_n`,
`bio3_n1`, `bioloc3_n`, `bioloc3_n1`) is just 4 more `ofs = ofs + nTot`
lines, no array-dimension change anywhere.

## What was added (all mirroring species 1's own structure exactly)

- `usercm.inc`: `sGdp_Beta3`, `sGdp_KLocal3`, `sGdp_MaxGrowth13/23`,
  `sGdp_HalfVelo13/23`, `sGdp_OriWeight13/23`, `sGdp_Bio3start`; offsets
  `ofs_dp_bio3_n/n1`, `ofs_dp_bioloc3_n/n1`.
- `USolBeg`: `parevl` reads for the 9 new APDL scalar parameters (all
  plain-text deck constants, e.g. `MAX_GROWTH13 = 100`); offset
  registration; IC seeding that **reuses species 1's own element region**
  (`Si_Bio1_ElemList`/`sGi_Bio1_ElemCnt`) rather than requiring a new
  `ELEM_LIST_BIO3` deck parameter — a deliberate trial-only shortcut, not
  a real per-species IC mechanism; default Bioloc3 = 1.0 fill, matching
  Bioloc1/2.
- `Ussfin`: full mirror of the Bio1/Bio2 explicit-update block — gradient,
  Laplacian, dot-products/norms, `OriBio3`, `GrowthBio3` (Monod kinetics,
  **no** `Interaction13/31`/`Interaction23/32` term — see below), the
  `Bio3_n` explicit update including `KLocal3*Bioloc3`, the `Bioloc3_n`
  local-retention update, all OMP `PRIVATE`/`SHARED` declarations, the
  n/n1 swap blocks, and the final `SetVals` write-back.

**Deliberately not done in this trial**: `Interaction13/31`/`23/32`
cross-species coupling (species 3 grows independently of species 1/2,
uncoupled) and a real per-species IC deck mechanism (reuses species 1's
element list instead of its own `ELEM_LIST_BIO3`). Both are natural next
increments, not attempted here.

## Bioloc3: tried, found necessary, hit a self-inflicted bug, fixed, now in

First draft of this trial dropped the `+ KLocal3*Bioloc3` term entirely
("just growth + diffusion + penalty, keep it minimal"). Running
`n3_growth_reference.py` (a pure-Python mirror of the exact Fortran
formula, no ANSYS/ifort/MKL dependency) against the n3trial deck's real
constants immediately showed `Bio3_n` going to -3.5 after a single
substep from an IC of 1.0. A second Python pass showed this negative dip
is actually a property of Oliver's own growth-formula/`DELTIM,0.1,0.1,0.1`
combination at these constants -- **species 3 computes the bit-identical
trajectory species 1 would compute under the same inputs** (same
`MaxGrowth`, `HalfVelo`, `OriWeight`, `KLocal`, `Beta`, deliberately
mirrored 1:1) -- so it is not a defect the n=3 extension introduces, and
`KLocal*Bioloc` alone doesn't fix it anyway (its magnitude,
`dt*KLocal3*Bioloc3 ~= 1e-3`, is three orders of magnitude below
`dt*OriBio3*GrowthBio3 ~= 4.5`). Wired it in anyway for structural parity
with species 1/2 (mirroring Bioloc1/2's allocate/GetVals/update/SetVals
chain in full).

First attempt at this triggered `Ussfin`'s existing NEM D-Matrix
neighbor-search computation to return `NaN` at `sID 57` -- reproducible
even on the completely unmodified `ds_oliver_wired_baseline.dat`,
deterministically. Root cause (confirmed by comparing `.obj` timestamps):
**not a real Oliver bug** -- after editing `usercm.inc` to add the 2 new
`ofs_dp_bioloc3_n/n1` offsets, only `USolBeg` and `Ussfin` were
recompiled; `Usermat_P21-V21_v222.obj`, `NEM_UserData_P21_V05.obj`, and
`userdata_P21-V21_Conection_Test.obj` were left stale from 4:45pm,
computed against the OLD (pre-Bioloc3) `/usercm/` layout, while the
freshly-rebuilt `USolBeg`/`Ussfin` used the NEW (shifted) layout --
exactly the same category of mistake as the earlier "forgot to recompile
USolBeg" SIG$SEGV a few steps before this. Recompiling all 5 files that
`#include "usercm.inc"` together (`USolBeg`, `Ussfin`, `Usermat_v222`,
`NEM_UserData`, `userdata_P21-V21_Conection_Test`) resolved it
completely -- baseline and n3trial (now WITH Bioloc3) both run 0 errors,
reproduced twice each.

**Lesson, worth repeating for n=4/n=5**: any edit to `usercm.inc` requires
recompiling ALL FIVE files that include it, together, in the same
`link_v222.ps1 -Sources` call -- never just the file(s) you think you
changed. `grep -l "usercm.inc" *.F *.f` in the WorkDir to get the full
list before any rebuild that touches the common block.

## Verification status

- **Structural/mechanical correctness**: verified via
  `n3_growth_reference.py` (species 3 reproduces species 1 exactly under
  identical inputs; species 3 does not appear in `GrowthBio1`/`GrowthBio2`
  at all, confirmed by inspection -- no unintended coupling).
- **Real ANSYS execution**: **verified**. `ds_oliver_wired_n3trial.dat`,
  now including the full Bioloc3 term, completes with 0 errors under
  `-smp -np 1` (see "DMP hang" below for why that flag is required),
  reproduced twice in a row for determinism. This is the actual
  deliverable: n=3 runs end-to-end in real ANSYS today, physically
  complete (matching species 1/2's own KLocal*Bioloc structure).

## Does this let us check correctness against the repo's own (n=5) ecology model?

Investigated 2026-09-07, before starting n=4/5: no -- not as a literal
numerical comparison. Oliver's `GrowthBio1/2` (and the mirrored
`GrowthBio3/4/5`) and this repo's own n=5 ecology model
(`ecology_jax.py` / `jax_hamilton_0d_5species_demo.py`) are structurally
different formulations, not two encodings of the same equation:

- Oliver's growth term is **substrate-explicit Monod kinetics**:
  `GrowthBio1 = sqrt(Lap^2) * [(MaxGrowth11 + Interaction12*Bio2)*Nut1/
  (HalfVelo11+Nut1) + ...]` -- two nutrient concentration fields drive
  growth, and the interaction term `Interaction12*Bio2` **multiplicatively
  shifts the max growth rate** applied to the Monod saturation curve.
- This repo's n=5 model has **no nutrient field at all** (CLSM-derived,
  closed volume-fraction simplex `phi_i`, `sum_i phi_i = 1`) and its
  interaction term `Ia = A @ (phi*psi)` enters the (implicitly,
  6-Newton-iteration-solved) residual **additively**:
  `-(c/Eta_i)*psi_i*Ia_i`, a bilinear coupling including a diagonal
  self-term `A_ii` (logistic-type self-growth) that has no Oliver-side
  analogue at all.

`OLIVER_MODEL_NOTES.md`'s "`Interaction12/21` are the off-diagonal of
matrix A" is a **structural analogy** (both are pairwise
growth-modulating terms), not an equation-level identity -- multiplicative
vs. additive, nutrient-driven vs. nutrient-free, explicit vs. implicit
time-stepping. Feeding identical inputs into both formulas would not
produce a meaningful pass/fail signal; there is no shared unit/variable
space to compare in. This is consistent with the "different scopes" table
already in `OLIVER_MODEL_NOTES.md` for the constitutive-law side of this
same integration.

**What actually is a meaningful correctness check** (and what n=3's
`n3_growth_reference.py` already did): confirming a new species block is
a faithful, non-interfering mirror of Oliver's *own* established n=2
pattern -- same formula shape, no accidental coupling into
`GrowthBio1/2`. That check applies unchanged to n=4 and n=5 below.

## Toolchain blocker found and resolved: DMP mode hangs on any rebuilt Ussfin

This is the first time in this repo's ANSYS-side work that
`Ussfin_P21-V21_Conection_Test.F` itself needed to be recompiled (every
prior session's wiring -- growth-law, ecology bridge -- touched only
`Usermat`, always reusing the untouched 2026-09-02 `Ussfin.obj`). Getting
it to compile at all, and then getting the resulting build to actually
run, surfaced two separate, real, pre-existing environment problems --
neither one caused by the n=3 biology code:

1. Raw `#INCLUDE 'mkl_spblas.fi'`/`'mkl_pardiso.fi'` doesn't parse under
   this fixed-form (`.F`) ifort compilation -- those two files use
   free-form-style trailing-`&` continuation. (`mkl_sparse_handle.fi` and
   `mkl_service.fi` have no such continuation and compile fine as-is.)
2. First fix attempt: `!DEC$ FREEFORM`/`NOFREEFORM` directives around just
   those two includes. This **compiles clean but produces a build that
   hangs at runtime**, reproducibly, on every deck tried -- including the
   completely unmodified `ds_oliver_wired_baseline.dat` with zero n=3
   content. The hang is silent (no error, no new output, zero further
   disk I/O) right at the boundary between `USolBeg` finishing and
   `Ussfin`'s first real substep (where PARDISO would first run for
   T/Nut1/Nut2). CPU keeps climbing (confirmed to 2+ hours on one run)
   with no progress.
3. Second fix attempt: the *correct* fix -- MKL ships proper free-form
   module sources (`mkl_spblas.f90`, `mkl_pardiso.f90`, declaring
   `MODULE MKL_SPBLAS`/`MODULE MKL_PARDISO`) alongside the raw `.fi`
   headers. Compiled those separately into `.mod`/`.obj` in the WorkDir
   and replaced the raw `#INCLUDE`s with `USE MKL_SPBLAS` / `USE
   MKL_PARDISO` (placed before `#include "impcom.inc"`, matching the
   USE-before-IMPLICIT ordering rule already established for the
   ecology-bridge's `biofilm_py_bridge` module). This compiles clean too
   -- **and still hangs, identically**, on the same unmodified baseline
   deck, clean single-job, no resource contention.

So neither MKL-header workaround was the compile-side cause -- both
compile clean. The real cause was found by testing serial (`-smp -np 1`,
no MPI) against the identical unmodified baseline deck: **it completed in
17 seconds, 0 errors**, immediately after a build that had just hung for
2+ hours under the default DMP (distributed memory parallel) launch. The
hang is specific to running a freshly-recompiled Ussfin.obj under DMP,
not to the MKL header fix technique, not ILP64/LP64 (that hypothesis was
wrong -- both the raw-`.fi`-with-`FREEFORM` build and the proper
`USE MKL_SPBLAS`/`USE MKL_PARDISO` module build hang identically under
DMP and both run fine under `-smp -np 1`), and not the n=3 code (the
unmodified baseline hangs under DMP too). Likely an MPI/OMP interaction
specific to this machine's build of the freshly-compiled object that the
untouched 2026-09-02 `Ussfin.obj` never triggered -- not investigated
further since the workaround (`-smp -np 1`) is sufficient and free.

**Resolution: always invoke this custom `ANSYS.exe` with `-smp -np 1`**
(not the plain `-custom .\ANSYS.exe` invocation `V222_PORT_INSTRUCTIONS.md`
documents elsewhere) whenever `Ussfin` has been recompiled. Confirmed
working, reproducibly, for both the plain baseline and the n=3 trial
deck.

The two isolation builds that established this (both hang under DMP,
both run clean under `-smp -np 1`, against the identical unmodified
`ds_oliver_wired_baseline.dat`):
- Pristine `Ussfin.F` + `FREEFORM`/`NOFREEFORM` directive fix for the MKL
  headers (no n=3 code at all).
- Pristine `Ussfin.F` + proper `USE MKL_SPBLAS`/`USE MKL_PARDISO` module
  fix (no n=3 code at all) -- this is the one kept in the final build,
  since it's the more correct fix even though both run equally well
  under `-smp -np 1`.

## n=4 / n=5 extension: done and verified in real ANSYS (2026-09-07)

n=3 is complete and running (growth + diffusion + penalty + KLocal*Bioloc,
uncoupled from species 1/2). Species 4 and 5 were then added the same
day, each a mechanical repeat of species 3's template with `3`->`4`/`5`
and fresh constant names -- no new design work, exactly as anticipated
below. **Both run cleanly in real ANSYS:**

- `ds_oliver_wired_n4trial.dat`: 0 errors, `-smp -np 1`, reproduced twice
  (`n4run`, `n4run2`).
- `ds_oliver_wired_n5trial.dat`: 0 errors, `-smp -np 1`, reproduced twice
  (`n5run`, `n5run2`).
- Regression checks after each new build: `ds_oliver_wired_n3trial.dat`
  still 0 errors after the n=4 build (`n3regress`); `n4trial` still
  0 errors after the n=5 build (`n4regress`) -- confirms each new
  species' additive-only edits didn't disturb the earlier ones.

All three (n=3/4/5) reuse species 1's element region for IC seeding and
stay uncoupled (no `InteractionN1/1N` terms) -- both "natural next
increments" flagged below were deliberately left for a future session,
consistent with keeping each step's edit surface small and verifiable.

The original per-species template this was built from, preserved as
reference for n=6+ (or for reproducing n=4/5 from scratch):

Concretely, per new species N:

1. **`usercm.inc`** (COMMON list + type-decl section, both places):
   `sGdp_BetaN`, `sGdp_KLocalN`, `sGdp_MaxGrowth1N/2N`,
   `sGdp_HalfVelo1N/2N`, `sGdp_OriWeight1N/2N`, `sGdp_BioNstart`;
   offsets `ofs_dp_bioN_n/n1`, `ofs_dp_biolocN_n/n1`.
2. **`USolBeg`**: 9 `parevl` reads (mirror species 3's block exactly,
   `MAX_GROWTH1N`/`MAX_GROWTH2N`/etc. as the APDL parameter names);
   2 offset-registration lines + 2 Bioloc offset lines; the "Fill local
   biofilmN = 1.0" SetVals block; the biofilmN IC-seed loop (decide then
   whether to keep reusing species 1's element region, or finally add a
   real per-species `ELEM_LIST_BIOn` deck mechanism -- flagged as a
   "natural next increment" back when n=3 was built, still true).
3. **`Ussfin`**: declare `vGdp_BioN_n/n1`, `vGdp_BiolocN_n/n1`; allocate;
   zero-init; GetVals; gradient+Laplacian block; dot-product/norm block;
   `OriBioN`; `GrowthBioN` (decide here whether to finally add the
   `InteractionN1/1N` etc. cross-species terms -- n=3 deliberately left
   these at zero); `PenForceN`; the explicit `BioN_n` update including
   `KLocalN*BiolocN`; the `BiolocN_n` update; OMP `PRIVATE`/`SHARED`
   additions (both parallel regions that touch species fields -- n=3
   needed 2 `PRIVATE` blocks and 4 `SHARED` line edits, not more); n/n1
   swap blocks (2 places); final `SetVals` write-back.
4. **Deck**: add the 9 new APDL scalar lines (mirror
   `ds_oliver_wired_n3trial.dat`'s inserted block) to a new
   `ds_oliver_wired_nNtrial.dat` copy of the baseline.
5. **Build**: `link_v222.ps1 -Sources` must list **all 5 files that
   `#include "usercm.inc"`** every time step 1 touches it (`USolBeg`,
   `Ussfin`, `Usermat_P21-V21_v222`, `NEM_UserData_P21_V05`,
   `userdata_P21-V21_Conection_Test`) -- this was the actual cause of
   both real bugs hit while building n=3 (a SIG$SEGV from missing
   `USolBeg`, then a NEM D-Matrix NaN from missing the other three).
6. **Run**: always `-smp -np 1`, never the plain DMP default, for any
   deck built against a freshly-recompiled `Ussfin` -- confirmed the
   default DMP launch hangs silently for hours on ANY deck (including
   the untouched baseline) once Ussfin has been rebuilt at all.
7. **Verify**: write a `nN_growth_reference.py` (copy `n3_growth_reference.py`,
   extend to N species) BEFORE running in ANSYS -- it caught a real
   physical-completeness bug (missing Bioloc term) for free, with no
   ANSYS/ifort/MKL dependency, well before any compile was needed.

Doing species 4 and 5 together (once resumed) is probably more efficient
than one at a time, since the mechanical edits are identical in shape --
but verify each species' Python reference independently before touching
Fortran, the way n=3's was.

## Automating steps 5/6 (build + run) after this session

Steps 5 and 6 above are exactly where both real bugs this session came
from (forgetting one of the 5 `usercm.inc`-dependent files; forgetting
`-smp -np 1` on a freshly-rebuilt `Ussfin`) -- entirely mechanical
mistakes, not modeling ones. `run_species_trial.ps1` (this directory)
wraps both into one command so a future species addition (n=6+) or a
re-run of n=3/4/5 cannot hit either mistake: it always recompiles the
fixed list of 5 files via `link_v222.ps1`, always launches with
`-smp -np 1`, clears stale `.lock` files first, and greps the output for
`NUMBER OF ERROR   MESSAGES ENCOUNTERED` to report PASS/FAIL directly
rather than requiring a manual log read. It does not automate steps 1-4
(the actual Fortran/deck edits for a new species) -- those still require
reading and editing Oliver's source by hand, the same as n=3/4/5 were.

## Full n=5 pairwise interaction terms: implemented, and a real instability found (2026-09-07)

Generalized Oliver's own `Interaction12`/`Interaction21` mechanism
(species 1/2 only) to all 20 directed pairs among species 1-5.
`InteractionIJ` appears in species I's growth equation, multiplying
`Bio_J` -- exactly mirroring the existing two-species pattern, just
summed over all `J != I` instead of a single partner:

```
InterSumI = sum_{J != I} InteractionIJ * Bio_J
GrowthBioI = sqrt(LapBioI**2) * (
    (MaxGrowth_nut1_I + InterSumI) * Nut1 / (HalfVelo_nut1_I + Nut1)
  + (MaxGrowth_nut2_I + InterSumI) * Nut2 / (HalfVelo_nut2_I + Nut2)
)
```

Added: 18 new `sGdp_InteractionIJ` scalars in `usercm.inc` (both the
COMMON list and the type-decl section), 18 new `parevl` reads in
`USolBeg`, `Sdp_InterSum1..5` (new PRIVATE scalars) plus the rewritten
`GrowthBio1..5` in `Ussfin`, and matching `PRIVATE`/`SHARED` OMP clause
additions. `interaction_terms_reference.py` (this directory) is the
Python reference, verified before any Fortran/ANSYS work: an all-zero
interaction matrix reproduces the pre-existing uncoupled formula exactly
(bit-for-bit, since `0 * Bio_J = 0`), and a modest test coupling
(`0.05` each) stays a small, bounded perturbation on `MaxGrowth=100` in
the 0D/well-mixed case the Python check represents.

**Real-ANSYS execution surfaced a genuine numerical instability, not a
code defect -- and, after narrowing it down, a real, calibratable
critical coupling threshold rather than "any coupling breaks it."**
`ds_oliver_wired_n5trial_inter0.dat` (all 20 interaction constants
explicitly `0`, exercising the exact same new code path as every other
case below) runs the complete 11 substeps cleanly -- 0 errors, force
convergence dropping to `1e-7`/`1e-8`, confirming the refactor is
behavior-preserving when interaction is off.

First pass made it look binary: `0.05`, `0.005`, even `-0.05`
(competitive sign), and even a single pair alone (`Interaction34`/
`Interaction43`, species 3/4 -- which don't touch the mechanical DOF
system at all, ruling out a stiffness-feedback explanation) all failed
at the load step's *final* substep with `The L-2 norm of the residual
force overflowed`. Bisecting the magnitude (`interaction_terms
0.05 -> 0.005 -> 0.001 -> 0.0007` fail, `0.0005 -> 0.0001 -> 0.00001 ->
0.000001` pass) found the real answer: **there is a genuine critical
coupling strength for this species-3<->4 pair, between `0.0005` (stable)
and `0.0007` (diverges), under Oliver's baseline `MaxGrowth=100`/
`dt=0.1`/`Penalty1=5`/`HalfVelo=0.1`.** This is a real nonlinear
stability boundary, not a memory/indexing bug -- a bug would not produce
a clean, monotonic, reproducible pass/fail split at a specific
magnitude.

This also explains why none of `dt` (10x smaller), `Penalty1` (10x
larger), or `MaxGrowth13/14` (10x smaller) individually rescued the
`0.05` test case (each tried once, at `0.05`, all still failed at
essentially the same ~10-12th substep regardless): `0.05` is roughly
70-100x above the ~`0.0006` threshold found for the baseline parameter
set, and each of those single 10x changes reshapes the threshold rather
than moving `0.05` below whatever the new one is. None of the three were
retested *combined*, or at a coupling magnitude actually below their own
(different) threshold -- that would be the natural next step if this is
picked back up.

A temporary debug write (`WRITE` gated by `ID.EQ.1`, added and reverted
in the same session -- not left in the tree) traced `GrowthBio3/4` for
the species-3/4-only `0.05` case (well above the threshold found later)
substep by substep: `0.159` (substep 1) -> `12024` (substep 2, a
~75,000x jump) -> `1.46e7` (substep 3) -> `-2.6e14` (substep 4) -> ... ->
`-Infinity` by substep 9 -- once past the threshold, the blow-up is fast
and total, not a gentle drift. The underlying reason a threshold this
low (~`0.0006`) exists at all traces back to the same thing
`n3_growth_reference.py`'s own trace flagged when n=3 was built:
**species 3's own *uncoupled* trajectory is already a knife-edge case**
("this dips negative and only slowly re-stabilizes via PenForce") --
i.e. a single species already sits close to the boundary between a
bounded oscillation and unbounded divergence under Oliver's own
`MaxGrowth=100`/`dt=0.1`/`Penalty1=5`, so it doesn't take much
additional push from a coupling partner to cross it -- just more than
this specific, now-measured amount.

**Conclusion: the interaction mechanism itself is implemented correctly
and is not the bug** -- `n5trial_inter0`'s clean 11-substep run through
the identical code path (just multiplying by literal `0.0`), and the
clean monotonic pass/fail split found by bisecting the magnitude, are
both the proof. **It is also not simply "unusable" -- coupling strengths
up to ~`0.0005` (a ~7% shift on `MaxGrowth=100` from a partner at
`Bio~1.0`, since `0.0005/HalfVelo=0.1` sets the scale) run cleanly, and
only strengths above roughly `0.0006` diverge.** Anyone using this
mechanism for a real (e.g. TMCMC-calibrated) multi-species study would
need to know where their intended coupling magnitudes sit relative to
this threshold *for their own parameter set* -- it will differ from
`~0.0006` for different `MaxGrowth`/`dt`/`Penalty1` choices, since single-
parameter 10x changes to each of those (tried once each, only at the
already-far-above-threshold `0.05`) didn't rescue that case, meaning they
reshape the threshold rather than trivially disable the effect. Mapping
that threshold's dependence on the other constants is the natural next
step if this is picked back up, not attempted this session.

**Update, 2026-09-08 (real ANSYS, last access day before ~1 month
cloud-only): this mapping done — see `threshold_param_sweep_README.md`
(same directory) for the full table and reasoning.** Headline: `MaxGrowth13/14` 10x
smaller raises the threshold by at least ~10x (0.0006 → ≥0.005, not
further bisected); `Penalty1` 10x larger barely moves it (same
`[0.0005, 0.0006)` bracket as baseline); `dt` 10x smaller makes it
*worse* — diverges even at 0.0005 where baseline is stable, most likely
because `TIME` was held fixed (so 10x smaller `dt` means 10x more
explicit growth-update substeps over the same elapsed duration, not the
same physical scenario resolved more finely) rather than a genuine
numerical-stiffness effect — flagged, not resolved, in that file. Also
narrowed the baseline threshold itself while at it: `0.0006` diverges
too, so baseline's true threshold sits in `(0.0005, 0.0006)`, tighter
than the `(0.0005, 0.0007)` bracket above.

**Current state left in `F:\biofilm_upf_wired`** (local only, see backup
note below): the full 20-constant interaction mechanism is wired into
usercm.inc/USolBeg/Ussfin and verified to preserve prior behavior with
all constants at their deck default of `0` (`n5trial_inter0`,
reproduced clean; `n3trial`/`n4trial` regression-checked clean against
the same rebuild). Test decks with nonzero coupling are kept alongside as documented
evidence, not deleted: the initial failing cases
(`n5trial_coupled(.dat/_small.dat)`, `n5trial_inter34only.dat`,
`n5trial_inter34neg.dat`, plus the `dt01`/`pen50`/`lowgrowth` single-
parameter retries, all still failing), and the bisection series that
found the actual threshold (`n5trial_inter_0p00001.dat` through
`n5trial_inter_0p001.dat`, passing below ~`0.0006` and failing above).

A full local backup of `F:\biofilm_upf_wired` (515 files, ~1.7 GB,
including the built `ANSYS.exe` and every deck used across the n=3/4/5
and interaction-term work) was made to `F:\biofilm_upf_wired_backup_20260907`
before this final round of edits, in case ANSYS access ends before
further work resumes -- not committed here (Oliver's source), but noted
so a future session on this machine knows where to look.

### Root cause narrowed further: this is not about species count at all

Natural follow-up question: did Oliver's own **unmodified n=2** deck
(`ds_oliver_wired_baseline.dat`, `Interaction12`/`Interaction21` set
nonzero, no n=3/4/5 code touched at all) also diverge? Tested directly
-- **it did not**: `ds_oliver_wired_baseline_inter12.dat`
(`INTERACTION12 = 0.05`, `INTERACTION21 = 0.05`, otherwise identical to
the delivered baseline) completes all 11 substeps with 0 errors.

The reason is visible directly in the deck, not a coincidence:
`MY_BIOSTART2 = 0.0` (species 2 starts completely unseeded, zero
everywhere) and `MAX_GROWTH12 = 0.1` (species 2's own growth rate is
1000x smaller than species 1's `MAX_GROWTH11 = 100`). **Species 2 in
Oliver's own delivered deck is a dormant placeholder for a genuinely
"mono-species" test, not a second active species** -- with `Bio2`
staying near zero throughout, `Interaction12 * Bio2` stays near zero
regardless of the constant's value, so turning `Interaction12/21` on
changes essentially nothing. Oliver's own delivery, as far as this deck
shows, has never actually exercised two *actively growing* species
coupled together.

Species 3/4/5, by contrast, were deliberately built with `MaxGrowth=100`
and `MY_BIOSTART=1.0` (mirroring species 1's real, active numbers, per
this doc's own "What was added" section above) specifically so the
n=3/4/5 extension would be a meaningful test of the growth mechanism --
not a dormant copy. Coupling species 3 and 4 was therefore the **first
time this deck's parameter scale has ever been asked to run two
genuinely active, growing species coupled together**, and that is
exactly where the instability appears (this baseline+`0.05` test itself
is also ~80x above the ~`0.0006` critical coupling strength found by the
bisection above -- consistent, not a separate phenomenon). **This
reframes the finding: it is not "n=5 destabilizes something n=2 handled
fine" -- it is "nobody, including Oliver, has verified this growth-rate
scale (`MaxGrowth=100`, `dt=0.1`, `Penalty1=5`) under real two-active-
species coupling before, and its actual safe coupling range turns out
to be much smaller than a value like `0.05` would suggest."** Worth
raising with Oliver directly, independent of anything n=3/4/5-specific.

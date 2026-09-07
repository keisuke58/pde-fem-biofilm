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

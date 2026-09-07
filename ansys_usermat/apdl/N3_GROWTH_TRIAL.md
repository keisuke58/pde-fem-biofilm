# n=2 -> n=3 growth-term trial (Oliver's ANSYS/NEM pipeline)

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
  `ofs_dp_bio3_n/n1`.
- `USolBeg`: `parevl` reads for the 9 new APDL scalar parameters (all
  plain-text deck constants, e.g. `MAX_GROWTH13 = 100`); offset
  registration; IC seeding that **reuses species 1's own element region**
  (`Si_Bio1_ElemList`/`sGi_Bio1_ElemCnt`) rather than requiring a new
  `ELEM_LIST_BIO3` deck parameter — a deliberate trial-only shortcut, not
  a real per-species IC mechanism.
- `Ussfin`: full mirror of the Bio1/Bio2 explicit-update block — gradient,
  Laplacian, dot-products/norms, `OriBio3`, `GrowthBio3` (Monod kinetics,
  **no** `Interaction13/31`/`Interaction23/32` term — see below), the
  `Bio3_n` explicit update (growth + diffusion + penalty, **no**
  `KLocal3*Bioloc3` term — see below), all OMP `PRIVATE`/`SHARED`
  declarations, the n/n1 swap blocks, and the final `SetVals` write-back.

**Deliberately not done in this trial**: `Interaction13/31`/`23/32`
cross-species coupling (species 3 grows independently of species 1/2,
uncoupled) and a real per-species IC deck mechanism (reuses species 1's
element list instead of its own `ELEM_LIST_BIO3`). Both are natural next
increments, not attempted here.

## Bioloc3: tried, found necessary for physical correctness, then reverted (real ANSYS bug found)

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
`dt*OriBio3*GrowthBio3 ~= 4.5`).

Wiring Bioloc3 in fully anyway (mirroring Bioloc1/2's allocate/GetVals/
update/SetVals chain, for structural parity with species 1/2) surfaced a
**real, separate ANSYS runtime bug**: it made `Ussfin`'s existing NEM
D-Matrix neighbor-search computation return `NaN` at `sID 57` --
reproducible even by rerunning the completely unmodified
`ds_oliver_wired_baseline.dat` against the Bioloc3-including build,
deterministically across repeated runs. This has nothing to do with
Bio3/Bioloc3's own math (NEM's D-Matrix is a geometry/neighbor-weighting
computation, computed once at initialization, that never reads Bio3 or
Bioloc3 state) -- most likely a latent, pre-existing arena-size or memory-
layout sensitivity in Oliver's own code that the extra 2*nTot-sized
Bioloc3 allocation happened to expose. **Reverted** back to the
growth-only version (Bioloc3 fully removed again) once this was
identified, restoring a clean, verified-good state -- fixing this latent
bug is out of scope for today's trial and doesn't need an ANSYS license,
so it's a good candidate for later.

## Verification status

- **Structural/mechanical correctness**: verified via
  `n3_growth_reference.py` (species 3 reproduces species 1 exactly under
  identical inputs; species 3 does not appear in `GrowthBio1`/`GrowthBio2`
  at all, confirmed by inspection -- no unintended coupling).
- **Real ANSYS execution**: **verified**. `ds_oliver_wired_n3trial.dat`
  (growth-only version, no Bioloc3) completes with 0 errors under
  `-smp -np 1` (see "DMP hang" below for why that flag is required),
  reproduced twice in a row for determinism. This is the actual
  deliverable: n=3 runs end-to-end in real ANSYS today.

## Toolchain blocker found and resolved: DMP mode hangs on any rebuilt Ussfin

This is the first time in this repo's ANSYS-side work that
`Ussfin_P21-V21_Conection_Test.F` itself needed to be recompiled (every
prior session's wiring -- growth-law, ecology bridge -- touched only
`Usermat`, always reusing the untouched 2026-09-02 `Ussfin.obj`). Getting
it to compile at all, and then getting the resulting build to actually
run, surfaced two separate, real, pre-existing environment problems --
neither one caused by the n=3 biology code:

This is the first time in this repo's ANSYS-side work that
`Ussfin_P21-V21_Conection_Test.F` itself needed to be recompiled (every
prior session's wiring -- growth-law, ecology bridge -- touched only
`Usermat`, always reusing the untouched 2026-09-02 `Ussfin.obj`). Getting
it to compile at all surfaced a real, pre-existing environment problem:

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

## Recommended next steps (not yet started)

1. Fix the Bioloc3 / NEM D-Matrix NaN bug found above -- license-free,
   can happen anytime. Needed before species 3 can be physically
   complete (matching species 1/2's own KLocal*Bioloc structure).
2. If further n-generalization is wanted: add
   `Interaction13/31`/`23/32`, then a real per-species IC mechanism, then
   repeat this same incremental pattern for species 4 and 5 -- explicitly
   deferred by the user for now ("n=3はできるようにしてほしい 4,5はまだ",
   2026-09-07).

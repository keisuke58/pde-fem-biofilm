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
  `Bio3_n` explicit update (now including `KLocal3*Bioloc3`, see below),
  the `Bioloc3_n` local-retention update, all OMP `PRIVATE`/`SHARED`
  declarations, the n/n1 swap blocks, and the final `SetVals` write-back.

**Deliberately not done in this trial**: `Interaction13/31`/`23/32`
cross-species coupling (species 3 grows independently of species 1/2,
uncoupled) and a real per-species IC deck mechanism (reuses species 1's
element list instead of its own `ELEM_LIST_BIO3`). Both are natural next
increments, not attempted here.

## Bioloc3 was not optional (caught by the Python reference check)

First draft of this trial dropped the `+ KLocal3*Bioloc3` term entirely
("just growth + diffusion + penalty, keep it minimal"). Running
`n3_growth_reference.py` (a pure-Python mirror of the exact Fortran
formula, no ANSYS/ifort/MKL dependency) against the n3trial deck's real
constants immediately showed `Bio3_n` going to -3.5 after a single
substep from an IC of 1.0 — i.e. the simplification wasn't a safe
minimal cut, `KLocal*Bioloc` is a load-bearing counterweight in Oliver's
own formula, not a separable extra. Fixed by wiring Bioloc3 in fully
(mirroring Bioloc1/Bioloc2's own allocate/GetVals/update/SetVals chain).

**Important nuance, found on the second Python pass**: adding
`KLocal3*Bioloc3` back does not, by itself, stop the negative dip — its
magnitude (`dt*KLocal3*Bioloc3 ~= 1e-3`) is three orders of magnitude
below `dt*OriBio3*GrowthBio3 ~= 4.5` at this deck's `DELTIM,0.1,0.1,0.1`.
The critical check that resolves whether this is a real bug: the
sanity-check block in `n3_growth_reference.py` confirms **species 3
computes the bit-identical trajectory species 1 would compute under the
same inputs** (same `MaxGrowth`, `HalfVelo`, `OriWeight`, `KLocal`, `Beta`
constants — deliberately mirrored 1:1 for this comparison). So the
negative dip is a property of Oliver's own growth-formula/`DELTIM`
combination at these constants, reproduced faithfully by species 3, not a
defect the n=3 extension introduced. Whether Oliver's real species-1
field also transiently goes negative at seed-region boundaries in
practice was not separately verified here (out of scope for this trial);
if it matters, the fix would be the same one already flagged for the
ecology bridge -- a substep size matched to the growth ODE's own
stiffness, not `DELTIM,0.1`.

## Verification status

- **Structural/mechanical correctness**: verified via
  `n3_growth_reference.py` (species 3 reproduces species 1 exactly under
  identical inputs; species 3 does not appear in `GrowthBio1`/`GrowthBio2`
  at all, confirmed by inspection -- no unintended coupling).
- **Real ANSYS execution**: **not verified**. Blocked by a toolchain issue
  unrelated to the n=3 code itself -- see below.

## Toolchain blocker found: Ussfin cannot currently be recompiled cleanly here

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

So neither MKL-header workaround is the actual cause. The common
denominator across every reproduction is simply "Ussfin was recompiled at
all in this environment" -- something about how this machine's
ifort 2025.3 + this MKL 2026.1 payload + this link_v222.ps1 flag set
produces a `Ussfin.obj` that doesn't behave like the working
2026-09-02 one, even for code paths (PARDISO for T/Nut1/Nut2) the n=3
trial never touches.

**Leading unconfirmed hypothesis**: an MKL integer-size (LP64 vs ILP64)
ABI mismatch between the PARDISO interface as compiled here and
`mkl_rt.lib` as linked (`.mod` files exist under both
`mkl\intel64\lp64\` and `mkl\intel64\ilp64\` in the payload -- not yet
tried explicitly selecting one). Not chased further today: this is a
toolchain problem, independent of the n=3 biology/code, and doesn't need
an ANSYS license to investigate -- unlike everything else in this
session, it can be worked on during the upcoming ANSYS-inaccessible
period.

**Two clean isolation tests support "any Ussfin recompile" as the
trigger**, both run resource-uncontended (only one ANSYS job at a time
after killing prior stuck jobs):
- Pristine `Ussfin.F` + `FREEFORM`/`NOFREEFORM` directive only (no n=3
  code at all) -> hangs identically.
- Pristine `Ussfin.F` + proper `USE MKL_SPBLAS`/`USE MKL_PARDISO` module
  fix (no n=3 code at all) -> hangs identically.

Both used the exact same unmodified `ds_oliver_wired_baseline.dat`.

## Recommended next steps (not yet started)

1. Resolve the Ussfin-recompile hang (ILP64/LP64 hypothesis first) --
   license-free work, can happen anytime.
2. Once Ussfin can be rebuilt and actually run again, rerun
   `ds_oliver_wired_n3trial.dat` (already created, includes the 9 new
   deck parameters) and confirm 0 errors + the species-3 field actually
   changes over the run.
3. If further n-generalization is wanted: add
   `Interaction13/31`/`23/32`, then a real per-species IC mechanism, then
   repeat this same incremental pattern for species 4 and 5.

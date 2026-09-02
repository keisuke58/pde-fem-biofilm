# Progress update — 2 September 2026

First of the fortnightly updates set out in
[`ROADMAP_2026.md`](../ROADMAP_2026.md) §4. Kept to one page by design; the
links go to the evidence rather than reproducing it.

---

Dear Prof. Soleimani,

A short update on the ANSYS side, and three things I would like to settle
early rather than late.

## Done

**The material law is packaged for Oliver's framework and verified.**
`BIOFILM_GROWTH_VISCO_V01` takes the same argument shape as the
`AceGenNeoHookV04` call site in their `Usermat_P21-V21_*.F`, and adds growth
kinematics and the viscous state over their purely elastic routine. It has no
ANSYS includes, no common block and no UPF calls, so — like theirs — it is
independent of the ANSYS release it is linked into. That matters because it
means I can develop on our v222 here and it still drops into their 2024 R2
build unchanged.

It is deliberately a thin adapter around the core that is already verified 0
ULP against the Abaqus UMAT, rather than a reimplementation, and a test holds
that at zero tolerance. So the verification travels with the routine instead of
having to be redone.

**A self-contained handover package for Oliver is ready** — two Fortran files,
no dependencies, generated from the sources under test so it cannot drift from
them.

**The routine now runs inside ANSYS.** The custom v222 executable builds and
solves here, so the law has been exercised by a commercial solver and not only
by a compiler — which closes the gap a reviewer would reasonably ask about,
since everything before this was checked outside any solver. Four things
came out of that, each one a claim a reviewer would otherwise have to take on
trust:

- Free and constrained growth, solved in ANSYS, against the closed form — and
  the constrained answer is reproduced by the ideal term *plus* the spherical
  discrepancy below, predicted rather than observed.
- Newton convergence on a 12,240-element two-layer shell. An inexact tangent
  shows up as iteration count, not as error; every substep converged in one to
  four iterations and automatic stepping never cut back.
- The state slots read back out of the solve. Growth has to travel in a state
  slot here, which fails silently if the write is short; all 176,800 samples
  came back as written.
- The same shell solved with the delivered routine and with the verified core,
  agreeing to every digit reported.

The decks, the solver listings and the scripts that read them are in
[`ansys_usermat/apdl/`](../ansys_usermat/apdl/); each figure is regenerated
from the listing rather than transcribed.

## Found, and I would value your view

Cross-implementation agreement establishes that two ports compute the same
expression identically. It does not establish that the expression is right, and
in checking that distinction I found that it is not: the isochoric split in
both implementations applies `J^(-2/3)` to the subtracted trace but not to the
tensor beside it. The coded and correct forms agree exactly at `J = 1` and
nowhere else — and growth is precisely what moves `J_e` away from 1.

The effect is bounded, and I want to be clear about that rather than alarming:
**the discrepancy is purely spherical**, so it is a pressure error and the von
Mises stress we report is unaffected to machine precision. It reaches the
deviator only weakly, through the viscous flow driver (≤0.01% over 40 steps).

My proposal is **not to change it before submission** — it does not move the
reported results, and correcting it would invalidate the whole reference chain
— but to write it into the verification chapter as a characterised limitation.
I think it reads as a stronger verification story than an unqualified "0 ULP",
since it shows the checks were examined rather than trusted. I would welcome
your view, and in particular whether you know if Felix's published formulation
defines the split as implemented; if it does, changing it would be a
disagreement with the paper rather than a correction, which is a different
decision.

Details, with a script that reproduces every number from source:
[`DEVIATOR_SCALING_FINDING.md`](../DEVIATOR_SCALING_FINDING.md).

## Next

Chapter 5 is drafted as far as it can go without results, with the
verification above written into it. The remaining piece is the results
section, which needs runs on the real geometry rather than the verification
cases — that is what I am moving to now. Producing them here rather than
waiting on Oliver's side to wire the routine in is what the local build was
for, so that dependency is gone.

## Three things to settle

1. **The exact submission date.** I am working to late November, aiming at the
   week of the 23rd so the last week is buffer, with the defence in December.
   Could you confirm?

2. **Chapter 5 scope.** I propose Chapter 5 is the ANSYS material-law
   contribution described above, and that the JAXFEM material I had drafted is
   held back for the continuation at Keio. My reasoning is that the assignment
   asks for integration into an *existing, working* FEM framework — which is
   Oliver's, not one of ours — and that the ANSYS side is the part that is
   verified. I would rather have your agreement on this now than surprise you
   with it in November.

3. **The handover to Oliver — and specifically what Chapter 5 should claim
   about it.** The package is ready to send whenever you think best. Since the
   local build works I no longer need his side wired in to produce results, so
   my inclination is to send it with the finished chapter rather than now, and
   to have Chapter 5 claim the law is *packaged for* their framework, which is
   provable from the code alone. If you would rather it claimed the law
   *running inside* their framework, that does depend on his calendar and I
   should send it immediately — so it is really a decision about the claim,
   and I would follow your preference.

With best regards,
Keisuke Nishioka

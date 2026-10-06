# Felix Klempt's implementation against the partner's element: what differs (6 Oct 2026)

Felix Klempt sent his current USERMAT / USOLBEG / USSFIN pool (the version
that runs on the LUH Linux server, two species and two nutrients; one of each
set to zero is the 2024 paper's model). It is confidential (CLAUDE.md, "Felix
Klempt's code and dissertation chapter"): it is kept outside git on IKMHIWI03
(`F:\felix_private\`), and this note says in words what differs. No code or
derivation of his is reproduced here.

Compared: Felix's pool against the partner's pool as received on 1 Sep 2026
(`F:\biofilm_upf_oliver`, before any of this repository's fragments were
pasted in). File by file:

| file | result |
|---|---|
| AceGen elements (phase/viscous, stress, three neo-Hooke versions, elasto-air), NEM user data, MySubroutines, userdata, `sms.h` | byte-identical |
| common block include (`usercm*.inc`) | differs: two new penalty parameters (species 2, species sum) |
| USolBeg (parameter reading) | differs: reading slip fixed, the two new penalty parameters read |
| Ussfin (field update) | differs: front term, bounding of phi, time step, surface points |
| Usermat (material) | differs: one line, the species sum |

So the NEM operators, the AceGen stress and phase routines and the nutrient
solve are the same code in both. Everything that differs is in the biofilm
update and in the parameter reading.

## The six points asked

| # | point | same / different | how |
|---|---|---|---|
| 1 | Why the front term is inactive in the partner's version | **different: fixed in Felix's** | In the partner's USolBeg, `ORI_WEIGHT12` and `ORI_WEIGHT22` are both stored in the species-1 weight, so the last one read (`ORI_WEIGHT22 = 0`) switches species 1's front term off (`apdl/FRONT_TERM_FIX.md`). Felix's version stores each in its own variable. The way Ussfin uses the weights is the same in both. |
| 2 | Consumption and the nutrient solve | **same** | Both put consumption parameter times phi (zero order, g phi, no factor c) on the right-hand side of a quasi-static diffusion solve (no time-derivative term in the nutrient matrix), with the same treatment of Dirichlet and Neumann rows. This matches Felix's mail of 5 Oct ("I think g phi"). |
| 3 | Keeping phi in [0, 1] | **different** | Partner: one explicit penalty force per species, pulling phi back when it leaves [0, 1]; species 2's force uses species 1's penalty coefficient. Felix: that penalty force is still computed (with the same species-1 coefficient for species 2) but no longer used. In its place, a semi-implicit update with three parts: an upper bound on the sum of the species (new parameter `PENALTYSUM`, split between the species by their shares), front growth that fades out as the sum approaches 1, and a per-species lower bound with its own coefficient (`PENALTY1`, new `PENALTY2`). Where the summed front contribution is negative, Felix's version adds no front growth at all instead of removing phi. |
| 4 | Form of the front term | **different** | Direction (normalised gradient product, the same lower bound 1e-14 on the denominator, the same weights) and the Monod factor with the interaction terms are the same. The magnitude differs: partner \|Laplacian of phi\|, Felix \|gradient of phi\| as in Eq. 34. Together with the one-sided growth of point 3, Felix's version grows the colony towards the nutrient and does not erode it on the far side, which is what his mail describes. |
| 5 | beta and the time step | **different in the step** | The beta term and the local source term are the same and explicit in both. Ussfin has an inner time loop inside each ANSYS substep. Partner: the phi and local-growth updates inside that loop use the full ANSYS substep, so a substep that the inner loop splits into several steps advances phi several times too far. Felix: they use the inner step. The current week decks (`deltim 0.025`, `MY_USSFIN_DT_MIN = MY_USSFIN_DT_MAX = 0.1`) take one inner step per substep, so the partner's runs here are not affected. With Felix's version the same decks would advance phi by the inner step 0.1 per substep of 0.025, i.e. four times too fast: `MY_USSFIN_DT_MIN/MAX` must not exceed the substep there. |
| 6 | Other differences | see below | |

Other differences (6):
- **Surface points of the phase field** (the loop that sets material
  properties at the interface points): Felix's version holds phi and the
  local growth at their previous values there; in the partner's version the
  corresponding lines are commented out, so that loop leaves them unassigned.
- **Penalty guard value**: Felix sets it once before the parallel region
  (comment: it was previously set only inside the solid-phase branch). It is
  still declared private in the OpenMP region in both versions, so whether
  the copy inside the region is initialised before the normalisation uses it
  is worth checking in both.
- **OpenMP lists** extended for the new variables; an end-of-routine print.
- **Usermat**: the species sum used by the material adds species 1 twice in
  the partner's version; Felix's adds species 1 and 2. This repository's
  wired build (`F:\biofilm_upf_wired`) already has the corrected sum.
- The partner's version has nothing that Felix's lacks: all differences are
  additions or corrections on Felix's side.

## His Workbench project (6 Oct 2026, read on IKMHIWI03, kept outside git)

The shared folder also holds an ANSYS Workbench project (`F:\felix_private\share2`),
three static analyses of a 1 mm cube (14^3 and 20^3 SOLID185 elements, a
third tiny one that stopped with an error), /units,MPA, T = 1 with 100
fixed substeps. Compared in words with what is used here:
- **It is not one of the 2024 test cases.** Its name says what it does, density-based growth only:
  the orientation weights are not set, so the front term is not active,
  and the new penalty parameters (species 2, sum) are not set either. It
  therefore does not tell which consumption and front term the 2024
  figures were computed with.
- Consumption enters as in the code, zero order (g phi), with g / d about
  20 per mm^2 on the 1 mm cube, i.e. a strongly consuming setting, in the same
  direction as the 2024 "low" case rather than Table 2's 4.1 value.
- beta is the same 0.02 mm^2/T* as this work's conversion of Table 2 to
  the partner's cube (CLAUDE.md, decided 5 Oct), the Monod constant k = 1
  and the maximum growth rate 1 /T*, as in Table 2 after scaling.
- Young's modulus 10 MPa with nu = 0.4999 (void 1000 times softer):
  neither the paper's 10 Pa nor the 10 kPa he suggested by mail.
- k_alpha (local growth) is 0.1 /T*, a hundred times Table 2's 1e-3.
So the g question for the 2024 runs stays open; the AceGen file
(Dr. Soleimani) or Felix himself has to settle it.

## What this means here

- The wired build used by the ANSYS week chain still has the reading slip
  (front term off) and the |Laplacian| magnitude; neither matters for the
  week runs, which run without the front term. The sum fix in Usermat is
  already in.
- For the Klempt 2024 reproduction (`KLEMPT2024_REPRODUCTION.md`): Felix's
  implementation uses Eq. 34's |grad phi|, zero-order consumption, and a
  front term that only adds phi and fades out at full occupancy. The
  reproduction's "growth on every face" (w = 0.5) and first-order
  consumption are not in his code. The one-sided, saturating front term is
  the candidate to test next in the reproduction and in Abaqus, stated as
  "as in the first author's current implementation" (his current code,
  which may differ in detail from the 2024 runs).
- No build or run was made with Felix's files; compiling them on Windows is
  to be decided after this comparison.

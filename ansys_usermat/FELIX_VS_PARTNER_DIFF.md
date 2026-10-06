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
- **Growth is the average over the two species, in both versions** (checked
  by reading the code, 6 Oct): the stress routine receives the mean of the
  two species' local growth variables, each starting at 1. In a run with one
  species, the other stays at 1, so the growth seen by the stress is half of
  that species' own growth. This explains the factor 1/2 in the ratio 0.49
  of chapter 5, where the partner's own growth variable was compared with
  Eq. 36 for one species. It does not affect this work's ANSYS runs, which
  compute the growth from Eq. 36 at the Gauss point instead. To be confirmed
  with Felix before it is used in a document.

## Values in his Workbench project (6 Oct 2026, read on IKMHIWI03, kept outside git)

The shared folder also holds an ANSYS Workbench project, kept in
`F:\felix_private\share2`. Its three static analyses are:
- a 1 mm cube with 14^3 SOLID185 elements;
- the same cube with 20^3 elements;
- a small model (237 nodes), which stopped with an error.

The two cube analyses use the same parameters; only the mesh differs.

**It is not one of the 2024 test cases**, and nothing in it corresponds to
4.1, 4.2 "high" or 4.2 "low". The set-up, as the project name says, is
density-based growth only:
- the nutrient is held on one face of the cube;
- species 1 and species 2 each start in a half ball either side of the
  mid-plane;
- the orientation weights are not set, so the front term is not active.

So the project shows the parameters he uses now, not the ones behind the
2024 figures. The table converts each value to the same footing as Table
2 (the domain size L and T* = 1), because the units differ: µm in the paper
and in the Abaqus reproduction, mm in his project.

| item | his project (1 mm cube, mm, MPa) | Klempt 2024 Table 2 (20 µm cube) | Abaqus reproduction (`make_klempt_inp.py`, 20 µm) | agree? |
|---|---|---|---|---|
| 1 consumption form | g phi (zero order), as in his code | g phi (Eq. 35) | g phi; g phi c with `--first-order` | form yes |
| 1 consumption, g L^2/d (the only group the quasi-static Eq. 35 sees) | 10 x 1 / 0.5 = **20** | 1e8 x 400 / 1e10 = **4** (4.1, 4.2 "high"); **400** ("low", g = 1e10) | 4 / 400 | **no: 5 times Table 2's 4.1 value** |
| 2 beta / L^2 | 0.02 /T* | 2 / 400 = 0.005 /T* | 0.005 /T* | **no: 4 times** |
| 2 nutrient diffusion d | 0.5 mm^2/T*; the solve is quasi-static, so only g/d matters | 1e10 µm^2/T* | 1e10 (only g/d enters) | through g L^2/d |
| 3 k_alpha | 0.1 /T* | 1e-3 /T* | 1e-3 /T* | **no: 100 times** |
| 3 nutrient boundary value, initial values | c = 1 on the elements next to one face; phi = 1 in the two half balls; c starts at 1 | c = 1 at the corner (4.1) or on the bottom face (4.2) | as Table 2 | set-up differs (not a 2024 case) |
| 4 front term: r, k | r = 1 mm/T* (r/L = 1 /T*), k = 1, interaction 50 between the species; inactive (weights not set) | r = 100 µm/T* (r/L = 5 /T*), k = 1 | as Table 2 | k yes; r not used in the project |
| 5 E, nu, units | /units,MPA; E = 10, nu = 0.4999 set under the names `MY_YOUNG_BIOFILM`, `MY_NU_BIOFILM`, **which this code does not read** (it reads `YOUNG_BIO`, `POISSON_BIO`, not defined in the deck) | 10 Pa, nu = 0.49 | 0.01 MPa (10 kPa, his mail) in the runs of 6 Oct | **no**; the value the solver used is undetermined |
| 6 geometry, mesh | 1 mm cube, 14^3 or 20^3 elements (h = 0.071 / 0.05 mm) | 20 µm cube, 1 µm elements (20^3) | 20^3, 1 µm | 20^3 in one of his runs |
| 7 time | T = 1, 100 fixed substeps (dt = 0.01), automatic stepping off, large deformation on; the inner time step parameters are not set, so the inner step equals the substep | T* in [0, 1]; 1e3 substeps nominal | T = 1, dt = 1e-3 | length 1, no time scale factor |
| 8 phi penalty, normalisation floor | `PENALTY1` = 100; `PENALTY2`, `PENALTYSUM` not set (zero); the floor of the normalised direction is 1e-14 in the code | not given | penalty 100 (or 0); eps = 1e-8 in the direction | penalty yes; floor differs, without effect on these runs |

`TB,USER` carries no constant block (`TB,USER,1,1,1` with no `TBDATA`). All
material and model constants reach the code as APDL parameters, which
USolBeg reads by name, so there is no ordered list of constants to check.
Where the names differ from what the code reads, that is listed above:
E and nu. The deck also sets a third species and an antibiotic. The
current code does not read these: it has two species and two nutrients.
The deck therefore belongs to an earlier version of the code.

Answers to the three questions:
- **g**: the form is g phi, as in the code. In the only dimensionless group
  the quasi-static nutrient equation sees, g L^2/d is 20, against Table
  2's 4 for 4.1. In the 4.1 run with Table 2's g, the nutrient already ran
  out (`abaqus_composition/README.md`), and a five times larger g consumes
  more, so g phi with this g makes 4.1 worse, not better. A run of 4.1
  with g x 5 is queued after the current runs, to show this.
- **time**: length 1 with 100 substeps and no factor such as s = 5/3. The
  other rates differ from Table 2 (k_alpha x 100, beta x 4 relative to the
  domain), so the project as a whole does not use Table 2.
- **E**: the deck intends 10 MPa (E = 10 with /units,MPA), neither 10 Pa nor
  10 kPa. Because of the name mismatch, the solver may not have used it.

The g question for the 2024 runs stays open: the AceGen file (Dr.
Soleimani) or Felix himself has to settle it.

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

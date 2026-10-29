# Paper outline (draft, 21 Oct 2026): coupling a multi-species point model to a 3D growth field

Status: skeleton only. Numbers in [..] come from the "paper values only" runs (`references/PARAMETER_TABLES.md`, last
table; Python `results/2026-10-3d/README_coupled_living.md`, ANSYS `u3_` after the trip). Co-authored (the user, 21 Oct).

**27 Oct 2026 (the user): this paper and P3 of `KEIO_PLAN.ja.md` become one paper** (graduation March 2028, so no
split): the coupling, the structure of the coupled equations, the verification, the examples (cube, bars, maze, molar
crown), plus the homeostatic-pressure growth law and the stress on real tooth/implant geometry. Aim: accepted before
the Keio final presentation (December 2027).

## Working title
Coupling a published multi-species biofilm point model to a continuum growth model: consumption by the living
biomass, verification in two implementations, and what the composition does to growth and stress

## Authors (to agree with everyone before any writing is shared)
| name | role in this work | to ask |
|---|---|---|
| Keisuke Nishioka | coupling, both implementations, verification, examples, writing | first and corresponding author |
| Oliver [surname] (IKM) | the element (NEM, AceGen), the coupling scheme agreed on 20 Oct, sketch | co-author? |
| Felix Klempt (IKM) | the growth model (Klempt et al. 2024) and the point model (Klempt et al. 2026) | co-author? his unpublished two-species work stays out |
| Dr. Meisam Soleimani (IKM) | supervision, model (Klempt 2024 co-author) | co-author |
| Prof. Junker (IKM) | head of the institute, Klempt 2024/2026 co-author | ask Meisam whether to include |
| Assoc. Prof. Mayu Muramatsu (Keio) | Keio supervisor | co-author (the user, 21 Oct); consent asked on 2 Nov; matters for the Keio master (single-author rule of the report: a co-authored paper is cited, not used as the report basis) |

Proposed order (the user, 29 Oct; to agree with Meisam in Jan 2027): K. Nishioka, Oliver [surname], (F. Klempt), M. Muramatsu,
M. Soleimani, P. Junker. Felix only if he wants to be (his published models are cited either way); Prof. Junker as head of
the institute last, if Meisam agrees.

Open: the Keio 課題研究報告 says a report "based on a submitted or published single-author paper" must state it; a
co-authored paper is still fine to publish, but check with Assoc. Prof. Muramatsu how it counts for the Keio side.
Felix's confidential two-species theory and code are not used; only the published 2024 and 2026 models.

## Target journals (in order)
1. Biomechanics and Modeling in Mechanobiology (where Klempt et al. 2024 appeared)
2. Archive of Applied Mechanics (where Klempt et al. 2026 appeared)
3. PAMM (short version first, GAMM annual meeting 2027) as a quick, citable step

## Structure
1. **Introduction.** Continuum biofilm growth models (Klempt 2024; others: Wanner and Gujer 1986, Eberl et al. 2001,
   Feng et al. 2021, Soleimani et al. 2023) describe the amount of biomass; point models (Klempt 2026) describe
   composition, living fraction and antibiotic response but have no space. Gap: how the two talk to each other
   without changing either. Contribution in three points.
2. **Models.** 2.1 Growth field (Eqs. 34-36 of Klempt 2024, front term, misprint g phi -> g phi c, dc/dt kept).
   2.2 Point model (Klempt 2026, used as published, own constraint). 2.3 Coupling: c*, alpha* from the fields
   (linear in time over the inner steps), consumption by the living biomass g phi_3D (sum_i phibar_i) c, antibiotic
   as the second field; Eq. 36 unchanged. Algorithm (`thesis/chapters/ch4_algorithm_draft.tex`).
   2.4 Parameters: all from Klempt 2024 Table 2 (converted) and Fritsch case I / Klempt 2026 Table 2; assumptions
   (phi_min, start of new points, source locations) listed with sensitivities.
3. **Implementations and verification.** Python (finite volumes) and ANSYS (neighboured element method).
   Uniform-field replay against the point model alone (1e-15 two/four species; ANSYS 1.8e-12); hand-over of
   sum phibar (difference 0); convergence in mesh and time; Python against ANSYS [..]. Limitations found on the
   way: neighbour truncation of the element (asymmetry), inner-step use of phi (fixed by one inner step).
4. **Results.**
   4.1 Two species, separate and mixed (Oliver's examples) [..].
   4.2 The nutrient field decides where the antibiotic acts: threshold c* = alpha* b_i; alpha*_0 sweep [..].
   4.3 Start share decides the winner (two species), sweep of the share [..]; four species of case 3 forget their
       start (same composition from separate or mixed start: phi = 0.41/0.25/0.17/0.15 [16^3]).
   4.4 Where the regions meet: boundary set where the fronts meet and frozen (no transport of composition) [fig].
   4.5 What reaches the stress: with Eq. 36 unchanged the composition changes the stress by [..] % through the
       consumption and the front term; front term off [..].
   4.6 Inputs of Klempt 2026 case 3 (sine nutrient, antibiotic from the start), two and four species [..].
5. **Discussion.** What the coupling can and cannot show (frozen boundaries, no displacement; Eq. 36 not living;
   slow-diffusion sensitivity). Outlook: transport of species (species-resolved continuum), calibrated data
   (Heine et al.), growth from the living biomass (Wanner and Gujer), real geometry (Keio).
6. **Conclusions.**

## Figures (planned)
1. Schematic of the coupling (`fig_coupling_schematic_agreed`) and the algorithm.
2. Verification: replay errors, hand-over, mesh/time convergence, Python vs ANSYS.
3. Examples: mid-plane maps (biomass, sum phibar, share, c, a, stress) for the four examples.
4. Antibiotic threshold: sum phibar vs alpha*_0 and local c.
5. Start-share map (which species wins).
6. Contact zone (`fig_contact_zone`).
7. Case 3, two vs four species.

## Before writing
- Agree authorship and target with Meisam (and Oliver, Felix) after the thesis is submitted (14 Dec) or earlier.
- Finish the ANSYS `u3_` runs and the Python sweeps.
- Decide the antibiotic onset (Klempt 2026 case 4: one third) and phi_min.

## Timeline (the user, 21 Oct 2026: submit around June 2027, with margin)
| when | what |
|---|---|
| until 14 Dec 2026 | thesis submitted (the coupling results as in the thesis) |
| 16 Dec 2026 | colloquium |
| Jan-Feb 2027 | agree authors and journal with Meisam, Oliver, Felix |
| 12 Mar 2027 | Keio 課題研究報告 (thesis content, separate from the paper) |
| Mar-May 2027 | write; ANSYS u3_ runs and missing Python runs |
| **June 2027** | **submit** (BMMB first choice) |
| until Dec 2027 | review and revision; aim: accepted before the Keio final presentation (Dec 2027) |
Option: a short PAMM paper from the GAMM annual meeting in spring 2027 as an earlier citable step. Do not let the
submission slip past June: a major revision can take more than six months.

## Scope and authorship notes (the user, 21 Oct 2026)
- Scope (corrected, the user): include a non-cubic geometry, as Klempt et al. 2024 did with the maze and the grate
  (their Tables 6-8): a maze/grate example in Python first, and the tooth/implant geometry of the Keio continuation
  (Abaqus, Keio server) if it is ready in time.
- **Assoc. Prof. Muramatsu is to be a co-author (the user, 21 Oct 2026).** Ask her consent at the meeting of
  2 Nov 2026 and tell Meisam when authorship is agreed (Jan-Feb 2027); all authors must agree.

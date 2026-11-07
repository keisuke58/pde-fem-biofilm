# USNCCM19 abstract (draft, 4 Nov 2026; numbers updated 6 Nov)

Conference: 19th U.S. National Congress on Computational Mechanics, Pasadena, 25-29 July 2027. Abstract due 15 Jan 2027.
Check the word limit and the minisymposium (biomechanics / multiphysics / growth) on the conference site before submitting;
this draft is about 250 words. Authors as agreed with Meisam (OUTLINE.md: proposed six); to be settled on 1 Dec.
Presenter: K. Nishioka.

## Title
Coupling a multi-species point model to a continuum growth model of oral biofilms: what the composition does to growth and stress

## Abstract
Continuum growth models of biofilms describe where biomass grows and which stresses it carries, but treat the biofilm
as one material. Point models describe how the species composition and the living fraction evolve, including the
response to an antibiotic, but have no space. This work couples the two. A published multi-species point model, derived
like the growth model from the extended Hamilton principle, runs unchanged at every Gauss point of a three-dimensional
phase-field growth model with a front term. It receives the nutrient and an antibiotic from diffusion-reaction fields,
interpolated in time, and only its living biomass consumes them. The growth enters the mechanics through a
multiplicative split with an isotropic growth tensor.

The coupled model is solved in a finite-volume code and in a user element in ANSYS and verified with uniform fields
(agreement with the point model alone to round-off), manufactured solutions and grid series. The front term is
advection dominated, so the stress requires a second-order upwind treatment; with a compressive limiter the stress of
a mixed start converges with a grid convergence index of about 6 % instead of about 40 % (grids up to 256^3).

In two-species examples an antibiotic lowers the living volume fraction of the seed by 84 to 96 %, while growth and
stress change by a few per cent: the growth counts dead biomass as well, and the stress follows the
stiffness-weighted gradient of the growth (correlation 0.94 over 51 runs). The composition reaches the mechanics
only through the nutrient and the front speed. Where the biomass starts and which species it is change the converged seed stress by
5 to 13 %. The results show which paths a composition-dependent mechanics of biofilms would need.

## Notes
- "5 to 13 %": extrapolated 64-192^3 (8 Nov): asymmetric -6 / -13 % (GCI 13 / 18 %), B as A +4.5 % (GCI 0.2 %). On one grid the gaps looked larger (-31/-40 % on 64^3).
- "about 6 %": superbee, mixed with antibiotic, 128/192/256^3: extrapolated 3.95e-4 Pa, GCI 6.3 % (64/128/192^3: 4.6 %); minmod GCI 41 %.
- No Felix material; Abaqus not named (third code) to keep the abstract short.

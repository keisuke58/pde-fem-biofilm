# Section 5 draft: discussion (cloud, 25 Oct 2026)

Draft text for Section 5 of `OUTLINE.md`. Numbers from `sec4_results_draft.md` and `KEIO_PLAN.ja.md` §10 (16^3
option runs, not in the thesis). Wording neutral; author list open.

## 5.1 What the coupling shows

- The composition reaches the growth field only through the consumption g phi_3D (sum_i phibar_i) c. When the
  antibiotic kills the biomass, less nutrient is consumed, the nutrient rises in the seed (0.54 to 0.91 in the
  separate example) and the term -r c/(k+c) n_c . grad phi_3D moves the edge faster. The effect on growth is small
  (alpha - 1 changes by about 1 %), because the growth law dalpha/dt = k_alpha phi_3D counts all biomass, living or
  dead.
- The local nutrient decides where the antibiotic acts: the point model compares c* with alpha* b_i, so the same dose
  removes a species in one part of the biofilm and not in another. A model without space cannot show this.
- The stress depends on how sharp the edge of the biofilm is, and that edge is set by the advection term of Eq. (1)
  (cell Peclet number about 31 on 16^3). The coupling therefore needs a second-order treatment of that term; a
  constant added diffusion gives a quarter of the converged seed stress.

## 5.2 Limits

1. **Front speed.** The speed r c/(k+c) depends on the nutrient only, so dead biofilm keeps spreading after the
   antibiotic. Tying the speed to the living biomass is a change of the model. In a test on 16^3 (front speed
   multiplied by the local living fraction) the seed stress rose by 51 % in the separate example and about fourfold in
   the mixed one, so this choice matters for the stress.
2. **Growth of dead biomass.** Eq. 36 lets dead biomass grow at the same rate as living biomass. Growth from the
   living biomass (as in Wanner and Gujer 1986) would make the antibiotic act on the stress; it is left as outlook.
3. **No transport of composition.** The point model has no space; the composition moves only through the rule for
   newly reached points (they take the living composition of their neighbours). In the maze example this rule decides
   whether the biofilm recovers after the antibiotic (living fraction 0.35 against 0.04). A species-resolved continuum
   would replace the rule.
4. **Mechanics one way.** The term -phi mu (I : C_e - 3) of Eq. 30 of Klempt et al. (2024) is neglected, as in that
   paper. With the stiffness used here (E = 10 Pa) it acts only for coefficients far beyond the paper values (an effect
   on growth from about 1e5 on, 16^3 test).
5. **Stiffness.** Living and dead biomass, and both species, have the same stiffness. A species-dependent stiffness
   changes the seed stress by 9 to 89 % in a 16^3 test, depending on which species is stiffer.
6. **Data.** The examples use synthetic data with known parameters. Calibrated data (for example the five-species
   data of Heine et al.) and a measured stress are needed before the stresses can be validated, not only verified.

## 5.3 Outlook

- Species-resolved transport, so that the composition moves with the biomass.
- Front speed and growth tied to the living biomass, with the consequences for stress quantified above.
- Calibrated multi-species data in the element.
- A realistic tooth or implant geometry.

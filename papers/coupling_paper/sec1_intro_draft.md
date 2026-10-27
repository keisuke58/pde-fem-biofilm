# Section 1 draft: introduction (cloud, 25 Oct 2026)

Draft text for Section 1 of `OUTLINE.md`. Wording neutral ("this work"); author list open. Check each statement on the
other papers against the PDFs in `references/` before it goes into the manuscript (marked [check]).

## 1.1 Context

Oral biofilms on teeth and implants are communities of several bacterial species. Their composition changes over
days: early colonisers prepare the surface for later ones, and pathogenic species such as P. gingivalis appear only
when the community allows it. Treatments act on the composition (antibiotics, antiseptics), while the mechanical load
on the surface comes from the amount and the growth of the biomass. A model that predicts stresses for a given
treatment therefore has to describe both: where the biofilm grows, and what it is made of.

## 1.2 Two kinds of model

- **Continuum growth models** describe the amount of biomass in space. Wanner and Gujer (1986) introduced the
  one-dimensional multispecies biofilm; Eberl et al. (2001) a density-dependent diffusion model of the spreading
  [check]. Feng et al. (2021) model a two-species oral biofilm in which the spreading follows the summed growth of the
  species, and Soleimani et al. (2023) the co-aggregation of two species in an ANSYS user element. Klempt et al.
  (2024) derive the growth of a biofilm from the extended Hamilton principle, with a phase field phi_3D, a nutrient c
  and the growth variable alpha of a multiplicative split F = F_e F_g, and compute the resulting stresses. Their model
  has one biomass and no composition.
- **Point models** describe the composition of a community at one point. Klempt et al. (2026) derive a
  multi-species model from the same variational principle: volume fractions phi_i and living fractions psi_i of N
  species, an interaction matrix A, and the response of each species to an antibiotic. Fritsch et al. (2025)
  calibrate such models to data with Bayesian inference [check scope]. These models have no space: the nutrient c*
  and the antibiotic alpha* are given inputs.

## 1.3 Gap

The point model needs local values of c* and alpha*, which a growth field can supply. The growth field needs to know
how much of the biomass is alive, because only living biomass consumes nutrient. Connecting the two at every
material point is the natural step, but it raises questions that neither paper answers: which variables are passed,
whether either model has to be changed, which new assumptions are needed (for example, the composition of biofilm
that a point reaches later), and whether the composition reaches the growth and the stress at all.

## 1.4 This work

This work couples the point model of Klempt et al. (2026) to the growth model of Klempt et al. (2024) without
changing either set of equations:

1. The point model runs at every integration point as published. It receives c* = c*_0 c and alpha* = alpha*_0 a from
   the nutrient and antibiotic fields, and returns the living volume fraction sum_i phibar_i, which sets the
   consumption g phi_3D (sum_i phibar_i) c in the field equations.
2. The coupling is implemented twice, in a finite-volume code in Python and in a neighboured-element user element in
   ANSYS, and verified: replay of the point model alone, agreement of the two codes, and grid convergence with the
   grid convergence index. The advection term -r c/(k+c) n_c . grad phi_3D of the growth equation turns out to set
   the accuracy of the stress.
3. Examples with two species and an antibiotic (synthetic data of Fritsch et al. 2025, case I, and the constants of
   Klempt et al. 2024 and 2026) show what the composition changes and what it does not: the antibiotic removes most
   of the living biomass, while growth and stress change little, because the growth law depends on phi_3D only.

Section 2 states the models and the coupling, Section 3 the implementations and their verification, Section 4 the
examples, Section 5 the limits of the coupling, and Section 6 concludes.

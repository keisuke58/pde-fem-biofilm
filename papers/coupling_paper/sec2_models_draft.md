# Section 2 draft: models and coupling (cloud, 25 Oct 2026)

Draft text for Section 2 of `OUTLINE.md`, from Chapter 4 of the thesis (`thesis/chapters/ch4_ansys.tex`, Sections
"The point model at every Gauss point" and "Two species with nutrient and antibiotic"). Wording is neutral ("this
work"); the author list is open. Notation as in Klempt et al. 2024 and 2026 (CLAUDE.md).

## 2.1 Growth field

The biofilm is described by the phase field phi_3D in [0, 1], the nutrient c and the growth variable alpha of
Klempt et al. (2024), with the multiplicative split F = F_e F_g, F_g = alpha I, alpha(0) = 1; the growth is alpha - 1.
After division by the viscosities the evolution equations read (Eqs. 34 to 36 of Klempt et al. 2024)

    dphi_3D/dt = beta lap phi_3D + k_alpha alpha - r c/(k + c) n_c . grad phi_3D,   n_c = grad c / |grad c|,   (1)
    dc/dt      = d lap c - g phi_3D (sum_i phibar_i) c,                                                      (2)
    dalpha/dt  = k_alpha phi_3D,                                                                                (3)

with zero flux on the boundary except where the nutrient is supplied. Two changes against the equations as printed
are used throughout, both confirmed by the authors: the consumption of Eq. 35 is g phi c (the printed g phi is a
misprint), and the front direction is regularised, grad c / (|grad c|^2 + eps^2)^(1/2), eps^2 = 1e-6 (eps^2 = 1e-8
changes the seed stress by 0.1 %). The mechanical term phi mu (I : C_e - 3) of their Eq. 30, which they neglect
because mu is small against the other parameters, is neglected here as well, so the mechanics is coupled one way.
The stiffness is E(phi_3D) = E (phi_3D^2 + f) with E = 10 Pa, nu = 0.49 and a void floor f = 1e-3.

Factor sum_i phibar_i in Eq. (2) is the coupling term of this work (Section 2.3); in Klempt et al. (2024) the whole
biofilm consumes.

## 2.2 Point model

At a material point the multi-species model of Klempt et al. (2026) describes the volume fractions phi_i and the
living fractions psi_i of N species and the void phi_0, with sum_i phi_i + phi_0 = 1. It follows from the extended
Hamilton principle with the free energy Psi = -1/2 c* phibar^T A phibar, phibar_i = phi_i psi_i, a dissipation
with viscosities eta_i and the constraint enforced by the multiplier gamma; the antibiotic of concentration alpha*
lowers psi_i with sensitivities b_i. The model is used as published, with its own constraint and its implicit time
integration (steps of 1e-4 in its own time unit, Newton's method to tolerance). For the examples the synthetic
two-species data of Fritsch et al. (2025), case I, are used with their true parameters (b = (1, 2)).

## 2.3 Coupling

The two models are connected at every integration point without changing either:

1. **Fields to point model.** The point model sees c* = c*_0 c and alpha* = alpha*_0 a, with a the antibiotic,
   which is transported by a second field with the diffusivity and uptake of the nutrient. Both are interpolated
   linearly in time over the inner steps of the point model.
2. **Point model to fields.** Only the living biomass consumes: the uptake of nutrient and antibiotic is
   g phi_3D (sum_i phibar_i), where phi_3D sum_i phibar_i <= phi_3D is the living volume fraction in the field. A dead
   biofilm consumes nothing. The growth law (3) is unchanged: alpha measures the biomass produced, alive or dead.
3. **Time.** T* = 1 is the duration covered by the point model's data, so the point model advances s Delta T* per
   increment with s the length of that window in its own units (s = 0.1 for case I), not a fitted value.
4. **Start of the point model.** It starts at a point once phi_3D reaches phi_min = 0.001. Where a seed is placed it
   starts from the seed composition. A point that the spreading biofilm reaches later takes the living composition of
   its colonised neighbours (their amounts divided in the ratio of phibar_i) and starts alive (psi_i = 0.999,
   gamma = 0); the papers do not fix this rule (assumption, compared in Section 4.4).

The composition reaches the growth and the stress only through the nutrient field and the front term of Eq. (1):
the antibiotic lowers psi_i, less nutrient is consumed, and more reaches the front.

## 2.4 Parameters

All constants are those of the papers: Klempt et al. (2024), Table 2, converted to a 2 mm cube (beta = 0.02 mm^2/T*,
r = 10 mm/T*, k = 1, k_alpha = 1e-3 1/T*, d = 1e8 mm^2/T* so that the fields are quasi-steady within about 4e-8 T*,
g l^2/d = 1 with l = 1 mm), and the antibiotic of case 4 of Klempt et al. (2026) (alpha*_0 = 100 from T* = 1/3).
The conversion of beta is an assumption (the paper gives no length and time scale). Values not given by the papers:
phi_min (results change by less than 1 % below 0.001), the start rule of new points and eps.

## Open for the paper

- Algorithm box (thesis Algorithm 4.1).
- Whether to show the Keio options (mechanical term of Eq. 30, front speed tied to the living biomass, species
  stiffness; `KEIO_PLAN.ja.md` §10) as an outlook or leave them out.

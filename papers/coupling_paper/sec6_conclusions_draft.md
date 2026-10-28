# Section 6 draft: conclusions (cloud, 27 Oct 2026)

Draft text for Section 6 of `OUTLINE.md`. Numbers from `sec3_verification_draft.md`, `sec4_results_draft.md` and the
thesis Ch. 4 grid table (`tab:ch5-gci`). Wording neutral ("this work"); author list open.

---

This work coupled the multi-species point model of Klempt et al. (2026) to the continuum growth model of Klempt et al.
(2024) at every integration point, without changing the equations of either model. The point model receives the local
nutrient and antibiotic from the fields, c* = c*_0 c and alpha* = alpha*_0 a, and its living volume fraction sets the
uptake g phi_3D (sum_i phibar_i) c of both. The main findings are:

1. **Implementation and verification.** The coupling was implemented in a finite-volume code, in the neighboured-element
   user element of the original authors in ANSYS and in a Galerkin user element in Abaqus. The point model inside each
   code reproduces the point model alone to round-off in uniform fields, the finite-volume and Galerkin codes converge
   at second order for a manufactured solution, and the three codes agree within their grid-convergence bands. Three
   numerical details of the neighboured element had to be corrected for this: the values at the surface points, the
   conservation of the Laplacian and the choice of neighbours.
2. **The advection term sets the accuracy of the stress.** The term -r c/(k+c) n_c . grad phi_3D of the growth equation
   has a cell Peclet number of about 31 on a 16^3 grid. A constant added diffusion keeps it stable but gives a quarter of
   the converged seed stress; a second-order upwind treatment is needed. With it, the growth converges with a grid
   convergence index of 0.7 to 3.3 % in all examples, the seed stress with 0.7 to 8 % when the species start apart and
   26 to 29 % when they start mixed.
3. **What the composition changes.** An antibiotic removes most of the living biomass (from 0.96 to 0.15 of the seed
   when the species start apart, from 0.96 to 0.04 when they start mixed), while the growth changes by about 1 % and the converged seed stress by less than its numerical uncertainty (6.55
   against 6.70e-4 Pa extrapolated when the species start apart, GCI 8 %; on one grid the difference reaches 10 %
   because the two cases converge from opposite sides).
   The local nutrient decides where the antibiotic acts, because the point model weighs c* against alpha* b_i.
4. **What the coupling needs beyond the papers.** The composition of biofilm that reaches a new point is not fixed by
   either model. Taking the living composition of the neighbours or the start state gives the same cube results within
   a few per cent, but decides in a maze whether the biofilm recovers after an antibiotic (living fraction 0.35 against
   0.04).
5. **Geometry.** Rigid obstacles restrain the swelling and raise the seed stress about tenfold when the species start
   apart; a maze makes the biofilm follow its corridors and approach the nutrient source.

The composition reaches growth and stress only through the nutrient and the front speed, because the growth law of
Klempt et al. (2024) counts all biomass. The next steps are therefore models in which the front speed, the growth or the
stiffness depend on the living biomass or on the species, and growth laws with a mechanical feedback such as a
homeostatic pressure. Each of them can be added to the coupling as presented here, and the verification procedure of
this work, several codes, a manufactured solution and grid convergence indices, carries over.

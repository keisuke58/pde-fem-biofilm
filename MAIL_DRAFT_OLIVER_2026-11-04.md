# Draft (not sent): results of the agreed coupling, to Oliver (to send around 4 Nov 2026)

Replaces `archive/MAIL_DRAFT_NEIGHBOURS_2026-10-21.md` (numbers of that draft are superseded). Plan (the user, 21 Oct): a
report of the results, not a list of questions. Recheck every number against the final thesis Ch. 4 (4 Nov) and the
grid series at the small step (IKMHIWI03 day7-day9). Attach: the progress deck `slides_progress_1024.pdf` (or the
pptx) and the schematic. No internal run names in the mail. Oliver's address: to fill in (not in the repository).
Cc Dr. Soleimani? (the user decides). Felix only if the user decides so, and nothing of his code or theory.

---

Subject: Coupling as agreed on 20 October: results and three changes in the element

Hi Oliver,

here are the results of the coupling we agreed on 20 October. I implemented it in your element in ANSYS and, as an
independent check, in a finite-volume code in Python.

1. What is implemented
- The point model of Klempt et al. (2026) runs at every Gauss point as published, with its own constraint
  sum phi_i + phi_0 = 1.
- c* comes from the nutrient field and alpha* from the second nutrient field, which carries the antibiotic; both
  are interpolated linearly in time over the inner steps of the point model.
- Eq. 35 with consumption by the living biomass, g phi_3D (sum_i phi_i psi_i) c, and the same for the antibiotic.
  Eq. 34 with the front term, Eq. 36 unchanged. All constants from Klempt et al. (2024), Table 2, converted to the
  2 mm cube; antibiotic from T* = 1/3 with alpha*_0 = 100 (case 4 of the 2026 paper).
- A point that the biofilm reaches later takes the living composition of its neighbours and starts alive.

2. Checks
- Every traced step of the point model in ANSYS equals the point model run on its own to 1.8e-12.
- Python and the element agree on 8^3 within 1 to 5 % once the changes below are in.
- Growth converges on refined grids (GCI 1.5 to 3 %). The seed stress needs a second-order treatment of the front
  term and converges more slowly: about 6.7e-4 Pa in the separate example (GCI 8 %), 4.0e-4 Pa in the mixed one
  (GCI 26 %). Your element with upwind gradients extrapolates to 5.8e-4 Pa (GCI 16.5 %), inside that band.

3. The examples you suggested (seed means at T* = 1)

| example | living volume fraction | share of A | nutrient c |
|---|---|---|---|
| A in the seed, B at the edges | 0.965 | 0.978 | 0.536 |
| the same, antibiotic | 0.152 | 0.951 | 0.911 |
| A and B mixed, equal amounts | 0.961 | 0.859 | 0.608 |
| the same, antibiotic | 0.038 | 0.610 | 0.973 |

- The antibiotic removes B first; whether A survives depends on the local nutrient, since the point model compares
  c* with alpha* b_i.
- Growth hardly changes with the antibiotic (about 1 %), since it lowers the living fraction but not phi_3D; the
  seed stress changes by about 10 %, through the nutrient and the front.
- In the maze of the 2024 paper the biofilm recovers after the antibiotic when new points take the neighbours'
  living composition (living fraction 0.35 at the end) and does not when they start from the start state (0.04).
- The front speed of Eq. 34 does not depend on the living biomass, so dead biofilm keeps spreading after the
  antibiotic.

4. Three changes in the element (patches available)
- Surface points: phi there was an inverse-distance mean that also took inner points, so at a wall with phi = 1 it
  came out low and phi leaked through the wall. Taking the nearest Gauss point reduced the leak to a quarter.
- Laplacian: the rows of the NEM operator sum to zero, the columns do not, so it does not conserve phi. A
  conservative form removes the rest of the leak.
- Neighbours: the search keeps the first points in element order; on 16^3 and 24^3 the stencils became asymmetric
  inside the domain. Selecting by distance makes them symmetric to 1e-15.

5. The front term
With r from Table 2 the cell Peclet number is about 31 on 16^3, and with a central gradient phi oscillates and the
run stops after two increments. A constant added diffusion r h/4 is stable but gives a quarter of the converged seed
stress on 16^3. A Lax-Wendroff term with an extremum sensor matches Python on 16^3 but is unstable on 24^3; upwind
gradients are stable from 8^3 to 24^3. I would be glad to hear how you or Felix stabilised this term.

I attach the slides with the figures and the convergence plots. Shall we talk about it in a short meeting?

Best regards,
Keisuke

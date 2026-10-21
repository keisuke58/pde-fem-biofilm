# Draft (not sent): results of the agreed coupling, to Oliver

Plan (the user, 21 Oct 2026): send as a report of the results, about two weeks after the meeting (around 4 Nov),
not as a list of questions. Update every number from the latest runs before sending (A-28..A-31 on IKMHIWI03, the
Python queue in the cloud). Attach: the schematic (`coupling_schematic.py --agreed`), the table of the two examples,
one animation frame figure (`assets/frames_u1_16_mixeq_ab.pdf`). Cc Felix only for item 3 if the user decides so.

---

Subject: Coupling as agreed on 20 Oct: first results

Hi Oliver,

here are the first results of the coupling we agreed on 20 October. I have implemented it in your element in
ANSYS and, as a check, in a separate Python version.

1. What is implemented
- The point model of Klempt et al. (2026) runs at every Gauss point as published, with its own constraint
  sum phi_i + phi_0 = 1. Its state is carried from one increment to the next.
- c* comes from the nutrient field and alpha* from the second nutrient field, which now carries the antibiotic;
  both are interpolated linearly in time over the inner steps of the point model.
- Eq. 35 with consumption by the living biomass only, g phi_3D (sum_i phi_i psi_i) c, and the same for the
  antibiotic. Eq. 34 with the front term (a tenth of the converted Table 2 value of r), Eq. 36 unchanged.
- Synthetic data: case I of Fritsch et al. (2025), two species, b = (1, 2), alpha*_0 = 50.

2. Checks
- Every traced step of the point model in ANSYS equals the point model run on its own to 1.8e-12.
- The living fraction used in the field equations equals the one the point model returned in the same
  increment (difference 0), and the constraint holds to 2e-16.
- ANSYS and Python agree on 16^3 to within 0.02 in the living volume fraction and the nutrient of the seed;
  the seed stress is about 12 % higher in ANSYS. Halving the outer increment changes the stress by 0.3 %.

3. The two examples you suggested (16^3, T* = 1, seed)   [QS: quasi-static ANSYS u1_ numbers; replace with the
   time-derivative results (Python _transc1 / ANSYS u2_) before sending]

| example | living volume fraction | nutrient c | share of A | stress [mPa] |
|---|---|---|---|---|
| A in the seed, B in the corners | 0.962 | 0.448 | 0.98 | 1.16 |
| the same, antibiotic from T* = 0.5 | 0.810 | 0.528 | 0.98 | 1.21 |
| A and B mixed, equal amounts | 0.955 | 0.456 | 0.79 | 1.18 |
| the same, antibiotic from T* = 0.5 | 0.226 | 0.836 | 0.77 | 1.27 |

- The antibiotic kills species B wherever it arrives; whether A survives depends on the local nutrient
  (c* = 100 c against alpha* b_1 = 50), so the nutrient field decides where the antibiotic acts.
- When the biomass dies, the nutrient is no longer consumed and rises (0.46 to 0.84 in the mixed seed).
- The stress changes little (+8 %), since Eq. 36 is unchanged: the point model reaches the stress only through
  the consumption, the nutrient field and the front term. [add: front term off, A-30]
- The mixed example depends on the start share: with Fritsch's start (0.25, 0.30) species B prevails, with
  equal amounts A survives. [add: share sweep, A-31]

Where A and B meet (your open question): the regions grow until the fronts meet and the boundary then stays in
place (position on the diagonal 0.795 mm at T* = 1, 1.5 and 2; the share of A changes from above 0.95 to below
0.05 within one cell). The composition is not transported with the biomass in this coupling, so neither species
invades the other; points reached between the fronts start mixed. Displacement would need the species in the field
equations. The same with the antibiotic: B in the corners dies (living fraction 0.84 -> 0.37 from T* = 0.5 to 2),
but A does not move into that space (share of A there 0.02-0.03, boundary still at 0.796 mm). A sharp boundary
is what experiments show; that it cannot move is a limit of the coupling. [figure: fig_contact_zone]

4. Two things I noticed in the element (no change to the code for these runs)
- Neighbour selection: GetNhgb enlarges its radius until there are enough points and then takes exactly
  NEIGHBOR_CNT of them in element order, so the last distance shell is cut by the numbering. Results of
  mirror-symmetric set-ups are therefore slightly asymmetric (phi up to 0.0085 on 8^3 with diffusion only,
  0.9 % in x in my earlier 16^3 runs). A distance-sorted selection or more neighbours reduce it; the seed mean
  changes by 0.5 %. I keep NEIGHBOR_CNT = 30 and describe it as a limitation.
- phi and the inner USSFin step: phi is advanced with the outer increment at every inner step. With one inner
  step per increment (MY_USSFIN_DT = 0.1, as in all your example decks and my Chapter 4 runs) this is exact;
  with a smaller inner step phi runs too fast. I had used 0.0025 for the first coupled runs and have rerun all
  of them with one inner step.
Would a distance-sorted neighbour selection be in line with how the operator is meant to work?

Best regards,
Keisuke

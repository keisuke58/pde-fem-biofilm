# Section 4 draft: results (cloud, 25 Oct 2026)

Draft text for Section 4 of `OUTLINE.md`. Numbers: Python TVD, inherit rule, seed means at T* = 1, from
`ansys_usermat/apdl/results/2026-10-3d/*inhliv*.json` (checked by `scripts/audit_thesis_numbers.py`, item 11) and
the IKMHIWI03 reports A-71, A-73, A-82, A-84, A-87, A-89. Replace the 64^3 values by the grid series at the small step
(32^3/64^3/128^3, outer step 0.0025) when all examples are in. Wording neutral; author list open.

## 4.1 Set-up

A 2 mm cube with the seed of the IKM example model at its centre. Separate example: species A in the seed, species B
in four small blocks along the edges. Mixed example: both species in the seed, in equal amounts (phi_1 = phi_2 =
0.275) or in the proportions of Fritsch et al. (0.25, 0.30). Nutrient held on the face y = -1 mm, antibiotic on all
outer faces from T* = 1/3. Front term on.

## 4.2 Composition and living biomass

| example | sum phibar | phi_1/(phi_1+phi_2) | c | alpha - 1 | von Mises [Pa] |
|---|---|---|---|---|---|
| separate | 0.965 | 0.978 | 0.536 | 9.06e-4 | 6.56e-4 |
| separate, antibiotic | 0.152 | 0.951 | 0.911 | 8.95e-4 | 5.68e-4 |
| mixed, equal amounts | 0.961 | 0.859 | 0.608 | 8.33e-4 | 3.50e-4 |
| mixed, equal amounts, antibiotic | 0.038 | 0.610 | 0.973 | 8.30e-4 | 3.23e-4 |
| mixed (0.25, 0.30), antibiotic | 0.037 | 0.510 | 0.973 | 8.30e-4 | 3.22e-4 |

(64^3, outer step 0.025; growth and stress converged values in Section 3.)

- The antibiotic lowers the living volume fraction of the seed from 0.965 to 0.152 when the species start separately
  and from 0.961 to 0.038 when they start mixed. B is removed first: its nutrient term c* A_22/eta_2 is half of A's, while the antibiotic term per viscosity, alpha* b_i/eta_i, is the same for both.
- **The nutrient field decides where the antibiotic acts.** The point model compares c* = c*_0 c with alpha* b_i:
  where the nutrient is high, A survives the antibiotic; where it is low, A dies as well. When the biomass dies the
  nutrient is no longer consumed and rises in the seed (0.54 to 0.91 in the separate example).
- Composition and living fraction change by less than 3 % with grid (32^3 to 128^3) and outer step (0.02 to 0.0025).

## 4.3 Antibiotic over time

With the same total dose given constantly, in pulses or as a sine from T* = 1/3 (64^3, A-71), alpha - 1 and the seed
stress change by less than 1 % (alpha - 1 by 0.1 %, von Mises by 0.8 %), while the living fraction follows the dose:
in the pauses of the pulsed dose the separate seed recovers from 0.17 to 0.30 (end value 0.298 against 0.152 with the
constant dose); the sine gives a lower minimum (0.110) and a higher end value (0.186) than the
constant dose. How the dose is given acts on composition and living biomass, hardly on the mechanics.

## 4.4 Which composition does new biofilm get?

The rule for newly reached points (Section 2.3, item 4) is not fixed by the papers. With the start values of the
nearer seed instead of the neighbours' living composition, the seed of the cube examples changes by a few percent
(sum phibar 0.147 against 0.152, separate example with antibiotic). In the maze of Klempt et al. (2024), Fig. 11,
where the biofilm keeps spreading after the antibiotic, the rule decides the outcome: at T* = 60 the living volume
fraction is 0.35 with the neighbours' composition and 0.04 with the start state. [figure: assets/anim_maze_rules_end.png]

## 4.5 What reaches the stress

- alpha - 1 hardly depends on the composition (Eq. 36 unchanged; the antibiotic lowers psi_i, not phi_3D).
- The mixed seed carries about half the stress of the separate one: its edge is less sharp because the front carries
  intermediate values of phi_3D.
- Stress needs a second-order treatment of the front term (Section 3); the extrapolated seed stress is about
  6.7e-4 Pa (separate, GCI 8 %; 6.55e-4 Pa and 0.7 % without antibiotic) and 4.0e-4 Pa (mixed, GCI 26 %; 4.5e-4 Pa
  and 29 % without antibiotic) at the small outer step. With the Fritsch start fractions the mixed stress is not yet
  in the asymptotic range on 128^3; the growth converges in every example (GCI 0.7 to 3.3 %).
- Limitation: the front speed of Eq. 34 depends on the nutrient only, so dead biofilm keeps spreading after the
  antibiotic (two thirds of the cells reached after the antibiotic in the maze border only dead biofilm). Tying the
  front speed to the living biomass changes the stress strongly (16^3: +51 % separate, about 4 times mixed;
  `KEIO_PLAN.ja.md` §10); this changes the model and is left as outlook.

## 4.6 Geometry: rigid obstacles and the maze

- Four rigid bars after Fig. 9 of Klempt et al. (2024) (positions read off the sketch), no biofilm or nutrient
  inside, frictionless contact: the biofilm grows around the bars and closes behind them; the bars restrain the
  swelling and raise the seed stress about tenfold. [figure: assets/fig_obstacle_columns_3d.png]
- The maze of Fig. 11 (560 um, 2D, Table 6 constants, h = 4 um): the biofilm follows the corridors, also into dead
  ends, and comes within about 40 um of the nutrient source by T* = 60; with the antibiotic held at the entrance from
  T* = 20 it stays about 130 um away. [figure: papers/coupling_paper/fig_maze_k24_h4.pdf]

## Still to add

- Grid series at the small step: done (thesis Ch. 4); 192^3 for the mixed example with the Fritsch fractions.
- 4.3 of the outline (start share decides the winner) and 4.4 (where the regions meet): re-run with the inherit
  rule before using the earlier numbers.
- Klempt 2026 case 3 inputs (sine nutrient), two and four species (`results/2026-10-3d-case3/`).

## 4.7 A clinical molar crown

To show the coupling on a real geometry, the same model and constants (Section 4.1, all in mm) are run around a
clinical molar crown (Open-Full-Jaw dataset, Gholamalizadeh et al. 2022; the derived surface in the repository,
CC BY-NC-SA). The crown (6.7 x 6.8 x 5.5 mm) stands with its cervical base on the bottom face of a 10 mm box and is
rigid for the fields (no biofilm, no nutrient inside, zero flux through its surface) and for the stress (u = 0).
Nutrient and, from T* = 1/3, antibiotic are held at 1 on the five faces around the crown (the saliva side); the
bottom face is zero flux (the root side). The biofilm starts as a cervical ring 0.25 mm thick and 1 mm high, species
A on the half x < 0 and B on x > 0. These choices are assumptions of this work. Grid 80^3 (h = 0.125 mm; 40^3 for
comparison), outer step 0.025, film = cells with phi_3D > 0.5, stress on the field grid with the film surface free
(elements with phi_3D <= 0.05 dropped), 17 min per run.

With the paper's front speed the biofilm spreads along the bottom face and outwards and reaches the side faces by
T* = 1 (every open cell is reached; the film holds 298 mm^3). The composition is carried with it: the A and B halves
stay separate (film mean share 0.50; 0.956 in the A ring, 0.051 in the B ring). The antibiotic enters with the
nutrient (a = c in the film, 0.73 on average, 0.24 at the ring, lowest value 0.14) and acts first on B: the living
fraction of the B ring falls from 0.71 at the onset to 0.13 at T* = 1, that of the A ring from 0.81 to 0.57; the
film mean is 0.37 (0.95 without antibiotic, where the ring reaches c = 0.13 and keeps 0.86 and 0.82 alive). Growth
and stress hardly depend on the antibiotic (film mean alpha - 1 4.3e-4 with and without; von Mises 3.0 and 3.2 mPa,
95th percentile 9.7 and 10.5 mPa, largest values in the corners between the film, the crown and the bottom face).
Between 40^3 and 80^3 the film means change by 3 % (living fraction 0.360 to 0.370, share 0.513 to 0.501, c 0.736
to 0.730, alpha - 1 4.1e-4 to 4.3e-4) and the film stress by 7 % (2.79 to 2.99 mPa; 2.80 mPa on the 80^3 field with
0.25 mm elements, so the difference is mostly the stress mesh). The example demonstrates the coupling on a real
surface; the stress of the cervical region with a growth law is not treated here.

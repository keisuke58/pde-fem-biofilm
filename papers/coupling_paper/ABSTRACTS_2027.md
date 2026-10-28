# Conference abstracts 2027 (drafts, 25 Oct 2026)

Not to be submitted before Dr. Soleimani and Oliver agree (order: permission from Meisam, then tell Assoc. Prof.
Muramatsu; `OUTLINE.md` for the author list). Numbers from `sec3_verification_draft.md` and the progress deck of
October 2026; recheck them against the final thesis (4 Nov) before sending. Deadlines are expected in January 2027
(USNCCM: mid-January in earlier editions; ESB: late January); confirm on the official sites.

Authors (to agree): K. Nishioka, Oliver [surname], M. Soleimani, (P. Junker), (M. Muramatsu).

---

## A. USNCCM19 (computational mechanics, about 300 words)

**Coupling a multi-species point model of oral biofilms to a continuum growth model: implementation in a
neighboured element method and verification**

Continuum models of biofilm growth describe how much biomass forms and where, but not which species it consists
of or how much of it is alive. Point models describe composition, living fractions and the response to an
antibiotic, but have no space. This work couples the multi-species point model of Klempt et al. (2026), used as
published, to the continuum growth model of Klempt et al. (2024), in which a phase field φ, a nutrient c and a
growth variable α with F_g = α I drive growth-induced stresses.

The point model runs at every integration point. It receives the local nutrient and antibiotic from the 3D fields,
interpolated linearly in time over its inner steps, and returns the living volume fraction Σᵢ φ̄ᵢ, which scales
the consumption of nutrient and antibiotic. The growth equations themselves are not changed. The coupling is
implemented in a user element of a neighboured element method in ANSYS and, independently, in a finite-volume code
in Python with up to 128³ cells.

Comparing the two implementations exposed four implementation issues in the element, among them a
non-conservative discrete Laplacian and the treatment of surface points; after the corrections both codes agree
within 1 to 5 %. The front term of the growth equation is advection-dominated, with a cell Péclet number of about
31 on a 16³ mesh. A first-order artificial diffusion keeps the solution stable but gives about a quarter of the
converged seed stress on that mesh; with second-order schemes (Lax–Wendroff with an extremum sensor, TVD) the seed
stress of the separate example converges to about 6 × 10⁻⁴ Pa (three-grid GCI below 10 %), while in the mixed
example it still changes between 64³ and 128³. The composition of the biofilm changes by less than 1 % between 16³ and 128³.

Examples with two species and an antibiotic show how the composition and the living fraction follow the nutrient
field, and examples with rigid obstacles and a maze geometry reproduce the growth patterns of the original model.

---

## B. ESB 2027 (biomechanics, about 300 words)

**Growth and stress of a multi-species biofilm: coupling a published multi-species point model to a continuum model with
nutrient and antibiotic fields**

Biofilms on dental implants consist of several species whose composition decides whether the biofilm is harmless
or pathogenic. Mechanical models of biofilm growth, such as the Hamilton-principle model of Klempt et al. (2024),
predict where biomass grows and the resulting stresses, but treat the biofilm as one material. The multi-species
point model of Klempt et al. (2026) describes composition, living fractions and antibiotic action, but not space.

This work connects the two without changing either model. At every integration point the point model receives the
local nutrient and antibiotic concentrations from the 3D fields and returns the living volume fraction, which
controls how much nutrient and antibiotic the biofilm consumes. Growth and stress follow the continuum model.
The coupling was implemented in two independent codes, a user element in ANSYS and a finite-volume code in Python,
which agree within a few percent after four implementation issues were found and corrected.

In a seed of two species with nutrient supplied from one side, the antibiotic, applied from the outer faces after
one third of the simulated time, lowers the living fraction of the seed from 0.96 to 0.15 when the species start
separately and from 0.96 to 0.04 when they start mixed. Where the nutrient is high, the antibiotic acts less,
since the point model compares the local nutrient with the antibiotic threshold. The growth variable and the seed
stress depend only weakly on the composition, because composition reaches the mechanics through the nutrient field
and the spreading of the biofilm. Stresses need second-order schemes for the front term and converge more slowly when the species start mixed; a
first-order scheme underestimates them by a factor of about four on a 16³ mesh.

Examples with rigid obstacles and a maze show growth around and through confined geometries, relevant for
biofilms in the gaps of implant–abutment connections.

---

Check before sending:
- the 0.96/0.15 and 0.96/0.04 values against the final Ch. 4 (A-85 base case, new points inherit the living
  composition: 0.965/0.152 and 0.961/0.038, unchanged at two digits);
- the examples use the synthetic case I of Fritsch et al. (2025), not the calibrated five-species model; say so in
  the talk;
- the implant–abutment sentence is a motivation, not a result.

Checked on IKMHIWI03 (25 Oct 2026) against the runs of 24-25 Oct:
- "three-grid GCI below 10 %" for the separate example: holds, but only just (9.8 % with the position rule, 8.4 % with
  the inherit rule; 128^3: 6.05e-4 and 5.90e-4 Pa). Safer: "converges to about 6 x 10^-4 Pa (three-grid GCI about
  10 %)". The ANSYS element with upwind gradients extrapolates to 5.8e-4 Pa (8^3/16^3/24^3, GCI 16.5 %), inside that
  band; this can be said for the element.
- "composition changes by less than 1 % between 16^3 and 128^3": true for the shares (0.949 to 0.951); the living
  fraction of the seed changes by about 10 % (0.133 to 0.147). Write "the shares of the species change by less than
  1 %", or add the living fraction separately.
- "Where the nutrient is high, the antibiotic acts less ..." (ESB): not checked by a run here; either show it (seed
  living fraction against the local c) or drop the sentence.
- New since the drafts, worth a sentence: the rule for newly reached points (they inherit the living composition of
  their neighbours) decides whether the biofilm recovers after the antibiotic (maze: 0.35 against 0.04), the same
  in Python and in a third, Galerkin implementation (Abaqus user element, SUPG); and the stress of the mixed example
  converges slowly (front with a nearly discontinuous edge).
- A-89 (25 Oct): the coupling step matters for the stress. At Delta T* = 0.025 the seed stress is about 6 % low
  (separate) and 18 % high (mixed) against Delta T* -> 0; the grid series are being redone at 0.0025. Recheck "about
  6 x 10^-4 Pa" and the mixed-example sentence ("still changes between 64^3 and 128^3") against them.
- The Abaqus element may appear in the paper and the abstracts (it is only kept out of the LUH thesis).

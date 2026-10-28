# Literature values for P3 (condition (f) of KEIO_PLAN.ja.md §1)

Written 17 Oct 2026 from a first web search. **Not yet checked against the original papers**: every value below
comes from a search summary or an abstract, and several of the sources are theses or secondary pages. Before a
value goes into a paper, read the original and note the page.

What P3 needs to place against the literature: (i) the thickness and growth of plaque (the model's seed 0.125 mm,
growth alpha - 1 up to 0.05-0.2); (ii) the stiffness E of a biofilm (the model works in vM/E and p/E; absolute
stresses scale with E); (iii) the homeostatic pressure p_h at which confined growth stalls, to set p_h/E.

| quantity | value found | source (as found) | status |
|---|---|---|---|
| In vivo plaque thickness on splints, 3 days | about 40-141 um (live smears), 65-129 um (dead smears); plaque on the splint 101-653 um | CLSM/SEM study, 5 participants, acrylic splints (search result: "Preliminary Study of In Vivo Formed Dental Plaque Using Confocal Microscopy and Scanning Electron Microscopy") | small study; read the paper |
| Enamel slabs in situ | 14.9 +/- 5.0 um after 1 day, 49.3 +/- 11.6 um after 7 days | in situ appliance study (PubMed 17446294 in the search) | read the paper |
| Early removable-disk study | wide variation between individuals; thickness set mainly by the number of microorganisms | PubMed 1063761 (1976) | old, qualitative |
| S. mutans biofilm elastic modulus | about 0.28 kPa (+/- 0.35, n = 30) in a group paper; 380 Pa in the thesis of the same group (conflict) | Southampton thesis and paper (eprints.soton.ac.uk 359755, 157185) | **values conflict; read both** |
| S. epidermidis biofilm | storage modulus 11 +/- 3 Pa, loss modulus 1.9 +/- 0.5 Pa at 1 Hz, yield stress about 20 Pa | rheometer study (Soft Matter 2013, c2sm27005f) | not a plaque organism |
| Plaque-simulating dextran gel | 17 kPa | same search | surrogate, not plaque |
| Confined E. coli (Chu et al. 2018, Nat. Commun. 9, 4087, doi 10.1038/s41467-018-06552-z; checked in the PMC full text, 17 Oct) | all tested strains generated pressures greater than about 1.5 psi (about 10 kPa) on the chamber roof (Fig. 1c); the stress-response reporter rpoH rose at about 0.75 psi (about 5 kPa) (Fig. 1d-f); turgor pressure of cells estimated at 30-100 kPa; **the paper gives no pressure at which growth stalls** (it says growth arrest happens at relatively low pressures at the boundary and attributes it to a possible adaptive response) | PMC6173693 | checked |
| Regrowth in an in vitro pocket model, 48 h | control 408 um and 290 um (two surfaces); treated 34 and 33 um | manufacturer-sponsored in vitro study | weak source; do not use |

## What the values say about the model

- **Thickness:** the seed layer of 0.125 mm lies inside the 14-650 um range found for plaque of one to several days.
  It is a reasonable thickness for a few days of undisturbed growth, not a measured one for the cervical margin.
- **Stiffness:** reported moduli of biofilm and plaque-like material are tens to a few hundred Pa. The Klempt 2024
  value E = 10 Pa lies in this range; the partner's example stiffness (E = 1000 MPa) does not and is only a bring-up
  value. Stresses in the model scale with E, so with E = 10-300 Pa the cervical band stresses are of order
  0.04 E = 0.4-12 Pa.
- **p_h:** with p_h/E = 0.05 and E = 10-300 Pa, p_h = 0.5-15 Pa. This is far below the ~10 kPa that confined
  colonies were found to exert (Chu et al. 2018: pressure generated, not a stall pressure), which concern single-species
  colonies in microfluidic chambers at much higher stiffness. If the model is
  read with E of the order of kPa (as the confined-colony work suggests for the matrix under pressure), p_h/E = 0.05
  gives p_h of the order of 100 Pa; to reach 10 kPa one needs p_h/E of about 1 at E = 10 kPa. The ratio p_h/E is
  therefore the quantity to vary (0.02-1), and P3 should say that its value cannot be set from the literature yet.

## Open items (to read before writing P3)

1. ~~Chu et al. 2018~~ checked on 17 Oct: it gives about 10 kPa generated and about 5 kPa for the stress response, **not a
   homeostatic (stall) pressure**. Earlier statements in this repository (keio_wp2/README.ja.md, the 2 Nov slides) that
   "growth stops at about 10 kPa (Chu et al.)" were corrected. A stall pressure for a biofilm still has to be found.
2. A review of dental plaque mechanics (modulus, yield stress, cohesion): not found in this search.
3. In vivo plaque growth after cleaning (toothbrushing-cessation studies) for the growth rate k_alpha and the time scale.
4. Whether any study gives the stress or deformation at the cervical margin (gingival margin) of plaque.

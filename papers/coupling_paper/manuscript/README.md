# Coupling paper: manuscript

First complete draft (27 Oct 2026) of the paper outlined in `../OUTLINE.md`, assembled from the section drafts
`../sec1_intro_draft.md` to `../sec6_conclusions_draft.md`, the figure plan `../FIGURES.md` and Chapter 4 of the thesis.
Target journal: Biomechanics and Modeling in Mechanobiology (Springer).

## Restructured 27 Oct 2026 (evening): weight on the modelling results

- Section 3 (implementations and verification) shortened to about one page: the three codes in one paragraph, the
  verification steps in one paragraph, the front term as a methods remark (cell Peclet number about 31, constant added
  diffusion gives a quarter of the converged seed stress, second-order treatment needed), the grid series with
  Figure 2 (`fig_grid_examples.pdf`). Moved to Appendix B: the manufactured-solution numbers, the four front-term
  treatments with the code comparison (`fig_convergence.pdf`), the time-step, regularisation and stability study, the
  GCI table. Moved to the supplementary list: the MMS error plot (`fig_mms.pdf`), the discretisation table of the three
  codes, the element's time-step series and stability analysis.
- New Section 2.4 "What the structure of the coupled equations implies", derived from the published equations only:
  the single path living fraction -> c -> front term -> phi_3D -> alpha -> stress (Eqs. 34-36), the bound
  alpha - 1 <= k_alpha T* (seed within 10 %) and the factor-two bound of the front speed; the quasi-steady nutrient
  as the linear problem lap c = (Lambda^2/l^2) phi_3D (sum phibar) c with Lambda^2 = g l^2/d = 1, depletion length
  l/(Lambda sqrt(kappa)) = 1 mm (living) and 2.6 mm (sum phibar = 0.15), slab estimate 0.41 / 0.82 at the centre;
  the antibiotic balance from Eq. 17 of Klempt et al. 2026, relative rate c* A_ii/eta_i - alpha* b_i/(eta_i phi_i^2)
  (decay rate grows as 1/phi_i^2; a balance with a graded outcome, no sharp threshold); the time link s equivalent to
  the viscosities eta_i/s. Sections 4.2, 4.4, 4.5, 5 and 6 read the examples against these statements (rates 2.7 and
  21 per T* in the separate and the mixed seed; the 1/phi_i^2 mechanism behind the new-point rule in the maze; the
  growth-from-living estimate 40 to 65 % from the seed time series). All numbers of the examples unchanged.
- Section 4.7 (molar crown) filled from `../sec4_results_draft.md` §4.7 with Figure 7; Figure 4 caption updated to
  the 2x2 layout; CC BY-NC-SA licence of the dataset in the data availability statement.
- 11 pt instead of 12 pt to keep 20 pages with the added content (the template is applied at submission anyway).

## Revised 27 Oct 2026 (night): review points

- Main message reframed as a result with a consequence (abstract, first contribution of the Introduction, lead of
  Section 2.4, new Discussion paragraph "What a composition-dependent mechanics needs", first item of the Conclusions):
  in the coupled equations the composition reaches growth and stress only through the nutrient and the front speed,
  so a model with composition-dependent mechanics needs growth or front speed tied to the living biomass or to the
  species, or a species-dependent stiffness; the estimates (growth from the living biomass 40 to 65 %, front speed
  51 % to fourfold, stiffness 9 to 89 %) are collected there and marked as estimates or 16^3 tests.
- Novelty sharpened in the Introduction: first coupling of a variational multi-species point model to a continuum
  growth model with stress, without changing either model; the derived structure; three independent implementations
  with verification; the start state of new points as an open modelling choice with its effect.
- New Discussion paragraph "Comparison with experimental scales" (plaque thickness 14.9 / 49.3 um after 1 / 7 days,
  Al-Ahmad et al. 2007; S. mutans modulus 0.28 kPa, Rmaile et al. 2013; growth arrest at about 10 kPa or below, Chu et
  al. 2018; 6 to 16 kPa, Helmlinger et al. 1997; all from `../../../p3/uel/LITERATURE_P3_checked.md`, DOIs checked
  against Crossref); no validation claimed, what a validation would need is stated. Four entries added to `refs.bib`.
- New Section 4.6 "A slower front" with Table 3: the four cube examples at r = 1 mm/T* (one tenth of the paper's
  speed; `coupled_living_consumption.py --front 0.1`, 64^3, outer step 0.025, JSONs
  `ansys_usermat/apdl/results/2026-10-3d/coupled_living_consumption_nf64_front0.1_phimin0.001_v2_fortran_ansysbc_inhliv*_tvd_wgmirror_fast_*.json`,
  colonised fraction from the `figdata/*.npz` fields, phi_3D >= 0.5). One sentence in Section 4.8 (the film around the
  crown is a consequence of the paper's front speed).
- To stay at 20 pages: the antibiotic-schedule figure (`fig_ab_schedule.pdf`) moved to the supplementary list (the
  numbers stay in Section 4.3), Appendix B.3 shortened, several figures narrower, float parameters relaxed, `\bibsep` 2 pt.

## Build

    pdflatex coupling_paper && bibtex coupling_paper && pdflatex coupling_paper && pdflatex coupling_paper

Offline, pdflatex + bibtex, no latexmk needed. The Springer Nature class `sn-jnl.cls` is not installed here, so the
draft uses `article` (11 pt since 27 Oct evening, A4, single column, `natbib` author-year, `lineno`); the template is applied at
submission. Figures are referenced in place (`../fig_*.pdf`, `../../../assets/*`), not copied. A figure that does not
exist yet is replaced by a box with its file name (`\figorbox`).

## Still open (visible as orange TODO boxes in the PDF)

- Author list, affiliations, contributions, funding, acknowledgements; corresponding author.
- Section 4.7 (clinical molar crown): run, numbers and the figure `../fig_real_tooth_coupled.pdf` are being produced
  separately; the dataset licence (CC BY-NC-SA) goes into the data availability statement.
- Fourth grid (192^3, computed on another machine) for the mixed example with the start fractions of Fritsch et al. (Table 3, Figure 4).
- Data availability: which parts of the element code may be published; release (DOI) of the repository.
- Springer Nature template at submission.
- Supplementary material (time-step series of the elements, stability analysis, sensitivities, maze animation) is
  announced at the end but not yet assembled.
- Figure 2 (verification across schemes) uses the runs that start newly reached points from the seed composition;
  redraw if that series is repeated with the inherit rule.
- All numbers come from the drafts, thesis Ch. 4 and the result JSONs; recheck against the final thesis before
  submission.

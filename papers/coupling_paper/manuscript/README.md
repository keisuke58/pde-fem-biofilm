# Coupling paper: manuscript

First complete draft (27 Oct 2026) of the paper outlined in `../OUTLINE.md`, assembled from the section drafts
`../sec1_intro_draft.md` to `../sec6_conclusions_draft.md`, the figure plan `../FIGURES.md` and Chapter 4 of the thesis.
Target journal: Biomechanics and Modeling in Mechanobiology (Springer).

## Build

    pdflatex main && bibtex main && pdflatex main && pdflatex main

Offline, pdflatex + bibtex, no latexmk needed. The Springer Nature class `sn-jnl.cls` is not installed here, so the
draft uses `article` (12 pt, A4, single column, `natbib` author-year, `lineno`); the template is applied at
submission. Figures are referenced in place (`../fig_*.pdf`, `../../../assets/*`), not copied. A figure that does not
exist yet is replaced by a box with its file name (`\figorbox`).

## Still open (visible as orange TODO boxes in the PDF)

- Author list, affiliations, contributions, funding, acknowledgements; corresponding author.
- Section 4.7 (clinical molar crown): run, numbers and the figure `../fig_real_tooth_coupled.pdf` are being produced
  separately; the dataset licence (CC BY-NC-SA) goes into the data availability statement.
- Fourth grid (192^3) for the mixed example with the start fractions of Fritsch et al. (Table 3, Figure 4).
- Data availability: which parts of the element code may be published; release (DOI) of the repository.
- Springer Nature template at submission.
- Supplementary material (time-step series of the elements, stability analysis, sensitivities, maze animation) is
  announced at the end but not yet assembled.
- Figure 2 (verification across schemes) uses the runs that start newly reached points from the seed composition;
  redraw if that series is repeated with the inherit rule.
- All numbers come from the drafts, thesis Ch. 4 and the result JSONs; recheck against the final thesis before
  submission.

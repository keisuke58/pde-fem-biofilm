# Third-party content bundled in this repository

This repository is public, so anything redistributed here needs its own licence
recorded. Everything else in the tree is covered by [`LICENSE`](LICENSE) (MIT).

---

## `references/Klempt2024_Hamilton_biofilm_growth_BMMB.pdf` — Klempt et al. (2024), the growth model this work builds on

> **A Hamilton principle-based model for diffusion-driven biofilm growth**
> Felix Klempt, Meisam Soleimani, Peter Wriggers, Philipp Junker
> Institute of Continuum Mechanics, Leibniz University Hannover
> *Biomechanics and Modeling in Mechanobiology* (2024) **23**:2091–2113
> doi:[10.1007/s10237-024-01883-x](https://doi.org/10.1007/s10237-024-01883-x)
> © The Author(s) 2024

**Licence: [Creative Commons Attribution 4.0 International (CC BY 4.0)](http://creativecommons.org/licenses/by/4.0/)**
— Open Access funding enabled and organized by Projekt DEAL. The licence permits
redistribution in any medium or format provided the authors and source are
credited, the licence is linked, and any changes are indicated. **The PDF is
bundled verbatim; no changes have been made to it.** Next to it,
`references/Klempt2024_Hamilton_biofilm_growth_BMMB.txt` is a plain-text
extraction (`pdftotext -layout`) of the same article, added so the text can be
searched; it is a change of format only, under the same licence and credit.

The file was renamed on 2026-10-02 from the publisher's DOI suffix
(`felix_s10237-024-01883-x.pdf`) to a searchable name. Where the code refers to "Klempt 2024", "Eq. 34–36" or "Table 2",
this is the document meant — see `JAXFEM/klempt2024_quantitative.py`,
`JAXFEM/klempt2024_sensitivity.py` and `klempt2024_gap_analysis.md`.

Note that this is the **single-species, interface-growth** PDE. The multi-species
model with the interaction scheme that this repository's own ecology follows is
the later paper (Klempt, Geisler, Soleimani, Junker, *Arch. Appl. Mech.* **96**,
164, 2026; arXiv:2509.01274), which is **not** bundled here — its fourteen
numerical examples are reproduced in `JAXFEM/klempt2026_reproduction.py`.

---

## `data/heine_species_distribution_biofilm.xlsx`

Unpublished experimental data provided by the experimental side (Heine) for the
companion Nishioka–Heine paper, not third-party licensed content. Provenance and
sheet layout are recorded in [`data/README.md`](data/README.md).

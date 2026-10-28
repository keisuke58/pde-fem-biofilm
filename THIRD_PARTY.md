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


## `references/Klempt2026_multispecies_interaction_AAM.pdf` — Klempt et al. (2026), the point model

> **A continuum multi-species bacterial growth model with a novel interaction scheme**
> Felix Klempt, Hendrik Geisler, Meisam Soleimani, Philipp Junker
> *Archive of Applied Mechanics* (2026) **96**:164
> doi:[10.1007/s00419-026-03160-y](https://doi.org/10.1007/s00419-026-03160-y)
> © The Author(s) 2026

**Licence: [Creative Commons Attribution 4.0 International (CC BY 4.0)](http://creativecommons.org/licenses/by/4.0/)**
— Open Access funding enabled and organized by Projekt DEAL. The PDF is the published version, bundled verbatim,
renamed on 2026-10-20 from the publisher's file name (`s00419-026-03160-y.pdf`). The `.txt` next to it is a plain-text
extraction (`pdftotext -layout`), a change of format only. The arXiv preprint (2509.01274) is not bundled: its arXiv
licence does not permit redistribution. The cases of this paper are transcribed in `JAXFEM/klempt2026_cases.py`.
## `references/Soleimani2023_coaggregation_SciRep.pdf` — Soleimani et al. (2023), two-species co-aggregation model in an ANSYS user element

> **Numerical and experimental investigation of multi-species bacterial co-aggregation**
> Meisam Soleimani, Szymon P. Szafranski, Taoran Qu, Rumjhum Mukherjee, Meike Stiesch, Peter Wriggers, Philipp Junker
> *Scientific Reports* (2023) **13**:11839
> doi:[10.1038/s41598-023-38806-2](https://doi.org/10.1038/s41598-023-38806-2)
> © The Author(s) 2023

Licensed under the Creative Commons Attribution 4.0 International License
(<http://creativecommons.org/licenses/by/4.0/>), which permits use and
redistribution in any medium or format provided the authors and source are
credited, the licence is linked, and any changes are indicated. **The PDF is
bundled verbatim; no changes have been made to it.** Next to it,
`references/Soleimani2023_coaggregation_SciRep.txt` is a plain-text extraction
(`pdftotext -layout`) of the same article, added so the text can be searched;
it is a change of format only, under the same licence and credit. Added
2026-10-04. Reading notes: `SOLEIMANI2023_NOTES.md`.

## `references/Feng2021_symbiotic_biofilm_BMB.pdf` — Feng et al. (2021), two-species oral biofilm, continuum advection-reaction model

> **Modeling of Symbiotic Bacterial Biofilm Growth with an Example of the Streptococcus–Veillonella sp. System**
> Dianlei Feng, Insa Neuweiler, Regina Nogueira, Udo Nackenhorst
> *Bulletin of Mathematical Biology* (2021) **83**:48
> doi:[10.1007/s11538-021-00888-2](https://doi.org/10.1007/s11538-021-00888-2)
> © The Author(s) 2021

Licensed under the Creative Commons Attribution 4.0 International License
(<http://creativecommons.org/licenses/by/4.0/>), which permits use and
redistribution in any medium or format provided the authors and source are
credited, the licence is linked, and any changes are indicated. **The PDF is
bundled verbatim; no changes have been made to it.**
`references/Feng2021_symbiotic_biofilm_BMB.txt` is a plain-text extraction
(`pdftotext -layout`), a change of format only, under the same licence and
credit. Added 2026-10-04. Reading notes: `ROADMAP_TWO_WAY.md`, literature section.

## `references/Rudolf2025_NEM_efficiency_FEAD.pdf` — Rudolf et al. (2025), the Neighbored Element Method used in the partner element

> **Computational efficiency and accuracy of the Neighbored Element Method**
> Tobias Rudolf, Felix Klempt, Hüray Ilayda Kök, Meisam Soleimani, Dustin Roman Jantos, Philipp Junker
> *Finite Elements in Analysis and Design* (2025) **249**:104353
> doi:[10.1016/j.finel.2025.104353](https://doi.org/10.1016/j.finel.2025.104353)
> © 2025 The Authors

Open access under the Creative Commons Attribution 4.0 International License
(<http://creativecommons.org/licenses/by/4.0/>), which permits use and
redistribution in any medium or format provided the authors and source are
credited, the licence is linked, and any changes are indicated. **The PDF is
bundled verbatim; no changes have been made to it.**
`references/Rudolf2025_NEM_efficiency_FEAD.txt` is a plain-text extraction
(`pdftotext -layout`), a change of format only, under the same licence and
credit. Added 2026-10-06. Reading notes: `KEIO_PLAN.ja.md` §8.

## `references/Vogel2020_adaptive_gradient_damage_IJNME.pdf` — Vogel and Junker (2020), adaptive reference for the gradient-enhanced damage model

> **Adaptive and highly accurate numerical treatment for a gradient-enhanced brittle damage model**
> Andreas Vogel, Philipp Junker
> *International Journal for Numerical Methods in Engineering* (2020) **121**:3108–3131
> doi:[10.1002/nme.6349](https://doi.org/10.1002/nme.6349)
> © 2020 The Authors

Open access under the Creative Commons Attribution License
(<http://creativecommons.org/licenses/by/4.0/>), which permits use,
distribution and reproduction in any medium provided the original work is
properly cited. **The PDF is bundled verbatim; no changes have been made to
it.** The `.txt` next to it is a plain-text extraction (`pdftotext -layout`),
a change of format only, under the same licence and credit. Added 2026-10-06.
Reading notes: `references/NEM_PAPERS.md`.

## `references/Junker2022_gradient_damage_large_deformation_IJNME.pdf` — Junker, Riesselmann and Balzani (2022), the Neighbored Element Method at large deformations

> **Efficient and robust numerical treatment of a gradient-enhanced damage model at large deformations**
> Philipp Junker, Johannes Riesselmann, Daniel Balzani
> *International Journal for Numerical Methods in Engineering* (2022) **123**:774–793
> doi:[10.1002/nme.6876](https://doi.org/10.1002/nme.6876)
> © 2021 The Authors

Open access under the Creative Commons Attribution-NonCommercial License
(<http://creativecommons.org/licenses/by-nc/4.0/>), which permits use,
distribution and reproduction in any medium provided the original work is
properly cited and is not used for commercial purposes. This repository is a
non-commercial academic project. **The PDF is bundled verbatim; no changes
have been made to it.** The `.txt` next to it is a plain-text extraction
(`pdftotext -layout`), a change of format only, under the same licence and
credit. Added 2026-10-06. Reading notes: `references/NEM_PAPERS.md`.

## `references/vonZabiensky2026_third_medium_contact_NEM_FEAD.pdf` — von Zabiensky, Jantos and Junker (2026), third-medium contact with the Neighbored Element Method

> **A fast and robust third medium contact approach using the neighbored element method**
> Max von Zabiensky, Dustin R. Jantos, Philipp Junker
> *Finite Elements in Analysis and Design* (2026) **255**:104489
> doi:[10.1016/j.finel.2025.104489](https://doi.org/10.1016/j.finel.2025.104489)
> © 2026 The Authors

Open access under the Creative Commons Attribution-NonCommercial-NoDerivatives
4.0 International License (<http://creativecommons.org/licenses/by-nc-nd/4.0/>),
which permits copying and redistribution in any medium for non-commercial
purposes with credit, without modification. This repository is a
non-commercial academic project. **The PDF is bundled verbatim; no changes
have been made to it.** The `.txt` next to it is a verbatim plain-text
extraction (`pdftotext -layout`), a change of format only, not an adaptation,
under the same licence and credit. Added 2026-10-06. Reading notes:
`references/NEM_PAPERS.md`.

Three further NEM papers that Felix Klempt sent (Jantos et al. 2019, Junker et
al. 2019, Blaszczyk et al. 2022) are not open access and are **not** bundled;
`references/NEM_PAPERS.md` gives their citations and a summary in my own words.

---
## `data/heine_species_distribution_biofilm.xlsx`

Unpublished experimental data provided by the experimental side (Heine) for the
companion Nishioka–Heine paper, not third-party licensed content. Provenance and
sheet layout are recorded in [`data/README.md`](data/README.md).

---
## `tier2b_real/_crown_real_watertight.stl` — clinical crown derived from the Open-Full-Jaw dataset

> **Open-Full-Jaw: An open-access dataset and pipeline for finite element models of human jaw**
> Torkan Gholamalizadeh et al.
> *Computer Methods and Programs in Biomedicine* (2022) **224**:107009
> Dataset and pipeline: <https://github.com/diku-dk/Open-Full-Jaw> (University of Copenhagen)

The file is derived from `Patient_1/Teeth/P1_Tooth_30.stl` of the OpenJaw dataset by
`tier2b_real/mesh_crown_real.py`: the clinical crown is cut out, uprighted, scaled and its cervical loop is
closed. It is a modified version. The Open-Full-Jaw repository is licensed under the Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International License
(<https://creativecommons.org/licenses/by-nc-sa/4.0/>); the derived file is shared under the same licence,
with the credit above, for non-commercial research use. The dataset archive asks for research use and
citation of the work above. Only this one derived crown is in the repository; the original OpenJaw STL files
and the meshes built from them are **not** redistributed (they stay on the Keio server). Provenance checked
on 16 Oct 2026 (`KEIO_CHAIN_LOG.md`).

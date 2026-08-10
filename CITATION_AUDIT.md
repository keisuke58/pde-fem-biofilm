# Citation audit (2026-08)

A check of every reference in the repository — the `.bib` files and the inline
citations in the docs — before anything is quoted in the thesis or shown to a
co-author. **Two findings need action; the primary reference is one of them.**

Scope: `biofilm_3tooth_refs.bib` (17), `jaw_biofilm_related.bib` (9),
`refs_openjaw.bib` (1), plus inline `Klempt …` citations across `*.md` / `*.tex`.

---

## 🔴 F1. The primary reference is cited under three different years — all wrong

The Klempt paper is the constitutive basis of the whole thesis. In the repo it
appears as:

| Form | Where |
|---|---|
| "Klempt 2024" | **20 files** (README, VERIFICATION…, FEM_README, BENCHMARK_PLAN, ch5_flow/*.tex, JAXFEM/algo_flow.tex, …) |
| "Klempt 2025" | **4 files** |
| "Klempt 2024/2025" | 2 occurrences |
| `@article{Klempt2025ContinuumBiofilm}` — **arXiv preprint 2509.01274 (2025)** | `biofilm_3tooth_refs.bib` (the only `.bib` a `.tex` actually loads) |

The paper is **published**, and the canonical record is neither 2024 nor 2025:

> Klempt, F., Geisler, H., Soleimani, M. et al. *A continuum multi-species
> bacterial growth model with a novel interaction scheme.*
> **Archive of Applied Mechanics 96, 164 (2026).**
> doi:[10.1007/s00419-026-03160-y](https://doi.org/10.1007/s00419-026-03160-y)

**Why this matters here:** Dr. M. Soleimani is a **co-author**. Citing his own
paper under the wrong year — three different wrong years — is the kind of thing
a co-author notices immediately.

The published entry has been added to `biofilm_3tooth_refs.bib` as
`Klempt2026ContinuumBacterialGrowth` (the preprint entry is retained and
cross-noted, since older text refers to the preprint's numbering).

### ⚠️ Two things to verify at the publisher page before submitting

Both are **unverifiable from this environment** — `doi.org`, `arxiv.org` and
`link.springer.com` are all blocked by the network egress proxy, and web-search
summaries for this paper proved unreliable (they returned unrelated hits). So:

1. **The title changed between preprint and publication.** The arXiv preprint
   (2509.01274, Sep 2025, confirmed via search) is *"A continuum multi-species
   **biofilm** model…"*; the published citation above reads *"A continuum
   multi-species **bacterial growth** model…"*. One search hit still showed
   "biofilm model" for the journal version. **Confirm the published title
   verbatim** and fix the `.bib` if needed.
2. **Equation and section numbers.** Text and figures cite specifics —
   "Klempt Eq. 34–36", "Klempt 2024, Sec. 2.1". Those numbers came from the
   version that was current when the text was written. **Re-check them against
   the published paper**; a right year pointing at wrong equation numbers is
   worse than the current state.

**Deliberately not done:** a blanket find-and-replace of the ~77 inline
"Klempt 2024" mentions. The year is only half the citation — the section and
equation numbers must be checked against the published PDF at the same time, and
that requires the paper itself. Fix them together, in one pass, with the PDF open.

---

## 🟡 F2. Five duplicate BibTeX keys with conflicting content

The same key is defined in two files with **different fields** — if both are ever
passed to BibTeX, one silently wins and the citation may lose its DOI:

| Key | Files | Conflict |
|---|---|---|
| `gholamalizadeh2022open` | `biofilm_3tooth_refs`, `refs_openjaw` | DOI present vs missing |
| `ODonnell2014InVitroPeriodontalBiofilm` | `biofilm_3tooth_refs`, `jaw_biofilm_related` | DOI present vs missing |
| `McCormack2017PDLFibres` | `biofilm_3tooth_refs`, `jaw_biofilm_related` | DOI present vs missing |
| `Asadzadeh2020PDLJawFEM` | `biofilm_3tooth_refs`, `jaw_biofilm_related` | fields differ (no DOI either side) |
| `Editorial2025DentalDigitalTwin` | `biofilm_3tooth_refs`, `jaw_biofilm_related` | fields differ (no DOI either side) |

Currently **latent, not live**: only `biofilm_3tooth_refs.bib` is loaded by a
`.tex` (`biofilm_3tooth_report.tex`); the other two are standalone reading lists.
It becomes a real bug the moment a document loads two of them.

**Recommended fix:** keep the DOI-bearing version as canonical (that is
`biofilm_3tooth_refs.bib` in 3 of the 5 cases), and have the other files either
drop the duplicate or use a distinct key.

---

## 🟢 F3. Hygiene items (no action needed before submission)

- **Missing DOIs** — 15 of 27 entries have no `doi` field; all of
  `jaw_biofilm_related.bib` and `refs_openjaw.bib`. Worth filling for the
  published ones.
- **Non-final references** that should be upgraded if they have since appeared:
  - `Fritsch2025BayesianMicrofilms` — `@unpublished`, "Manuscript" (the basis of
    the TMCMC pipeline). Co-authors include Klempt, Soleimani, Junker, Beer.
  - `Heine2025PeriImplant` — the companion experimental paper.
  - `Klempt2025ContinuumBiofilm` — superseded by the 2026 journal version (F1).
- `Abaqus2024` is a `@manual` — correct as-is; check the version string matches
  the release actually used.

---

## What was checked, and how

| Check | Method | Result |
|---|---|---|
| Every `.bib` entry parsed; year / type / journal / DOI | script over all `@…{}` blocks | 27 entries, table above |
| Duplicate keys across files, with content comparison | key→file map, field diff | 5 conflicts (F2) |
| Inline `Klempt` year usage | grep over `*.md` / `*.tex` | 2024 ×20 files, 2025 ×4 (F1) |
| Which `.bib` is actually loaded | `\bibliography{}` in `*.tex` | only `biofilm_3tooth_refs` |
| Published record of the primary reference | **not verifiable here** — publisher domains egress-blocked | taken from the author-supplied citation; see F1 caveats |

*Re-run the structural checks any time; they are pure text analysis over the
`.bib` files and docs.*

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

### ✅ Resolved 2026-08-20 — the user supplied the published PDF directly

The two items below are no longer blocked on network access; the full text of
`10.1007/s00419-026-03160-y` was read end-to-end. **This surfaced a bigger
problem than either item anticipated — see the new §F1b below.**

1. **Title, confirmed.** The published title is verbatim *"A continuum
   multi-species bacterial growth model with a novel interaction scheme"* —
   matches the `.bib` entry exactly. No fix needed; the "biofilm model" vs
   "bacterial growth model" wording concern was unfounded (that was the
   *preprint's* title only).
2. **Equation/section numbers — worse than "stale", actually wrong-paper.**
   See §F1b: the growth-kinematics equations these numbers were meant to
   point at (`Fg=(1+α)I`, `F=Fe·Fg`) **do not exist in this paper at all**.
   Any inline citation of the form "Klempt Eq. N" for growth kinematics needs
   to be redirected to a *different* Klempt paper (below), not merely
   renumbered against this one.

## 🔴 F1b. The growth-kinematics citation has been pointing at the wrong paper

The published AAM 2026 paper (`Klempt2026ContinuumBacterialGrowth` in the
`.bib`) is a **0-D material-point model for species population dynamics** —
state variables `φ` (volume fraction) and `ψ` (living fraction), derived via
Hamilton's principle, strong-form ODEs at Eqs. 16–18. It contains **no
deformation gradient, no `Fg`, no growth kinematics of any kind.**

This repo's UMAT constitutive basis — `F = Fe·Fg` — is **not in this paper**.
It is this paper's own reference **[11]**, now confirmed 2026-08-20 directly
against the published PDF (not just secondhand via the AAM paper's
bibliography):

> Klempt, F., Soleimani, M., Wriggers, P., Junker, P. *A Hamilton
> principle-based model for diffusion-driven biofilm growth.* Biomechanics
> and Modeling in Mechanobiology 23, 2091–2113 (2024).
> doi:[10.1007/s10237-024-01883-x](https://doi.org/10.1007/s10237-024-01883-x)

— a different journal, a different (though overlapping) author list
(Wriggers instead of Geisler), genuinely published in 2024 (not a
preprint-to-2026 story like the other one). **This paper was not in the
`.bib` files at all before today** — added as `Klempt2024DiffusionDrivenGrowth`
in `biofilm_3tooth_refs.bib`, with the `Klempt2026ContinuumBacterialGrowth`
entry's note corrected to say what it actually contains.

**Why this matters:** every doc/figure/UMAT comment in the repo that cites
"Klempt 2024" for growth kinematics (`RESEARCH_MODEL.md`, `HANDOFF.md`,
`README.md`, `VERIFICATION_SENSITIVITY_LIMITATIONS.md`, ~19 files touch this
— see the search below) was *coincidentally* using the right year, but if any
of them were ever repointed at the AAM 2026 DOI (the natural thing to do once
"the paper got published"), that would make it definitively wrong, not just
year-stale. Most docs already hedge with "see CITATION_AUDIT.md before
quoting equation numbers" (good practice already in place), so the exposure
is contained — but the fix is: **cite `Klempt2024DiffusionDrivenGrowth` for
growth kinematics, `Klempt2026ContinuumBacterialGrowth` only for the φ/ψ
species-interaction model.**

**Not yet done:** the ~77 inline "Klempt 2024" mentions across the repo
haven't been individually re-pointed at the correct one of the two papers —
that still needs a pass with both PDFs open (both now confirmed available).

## 🔴 F1c. New discrepancy found while confirming F1b: `Fg = αI`, not `Fg = (1+α)I`

Reading the confirmed 2024 paper directly (Sec. 2.1, kinematics of growth)
turned up a second, independent problem: **the paper's own definition is
`Fg = α·I`**, isotropic growth scaled directly by the local expansion
parameter — quote: *"The growth part is constructed with the local expansion
parameter α in the form of Fg = αI, where I is the identity tensor."*

This repo's UMAT and every doc (`RESEARCH_MODEL.md`, `HANDOFF.md`, README,
etc.) consistently uses **`Fg = (1+α)I`** instead. These are not the same
convention — under the paper's own `Fg = αI`, `α = 0` gives a singular,
non-invertible growth tensor (`Fg = 0`), which cannot be what's intended at
zero growth. Two explanations are possible: (a) this repo's `α` is defined
as "the paper's α minus 1" (a shifted variable, undocumented as such), or
(b) `Fg = (1+α)I` is a deliberate, sensible fix for well-posedness at zero
growth that was never written down as a departure from the cited source.

**Not fixed — flagged for a decision.** This is a modeling-convention
question, not a typo, and needs whoever owns the UMAT derivation to confirm
which explanation is correct (and if (a), that the shift is stated
explicitly wherever `Fg=(1+α)I` appears; if (b), that the departure from the
cited paper is documented as intentional). Surfaced to the user 2026-08-20;
unresolved.

---

## 🟢 F2. Five duplicate BibTeX keys with conflicting content — fixed 2026-08-20

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

**Done:** all 5 duplicates in `jaw_biofilm_related.bib` / `refs_openjaw.bib`
renamed with a `_jbr` / `_openjaw` suffix and a comment pointing at the
canonical entry (no content deleted — both files are informal reading lists,
not currently `\bibliography{}`-loaded by anything, so nothing was live-broken
either way). `gholamalizadeh2022open`'s `publisher` field, present only in the
`refs_openjaw.bib` copy, was folded into the canonical entry in
`biofilm_3tooth_refs.bib`.

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

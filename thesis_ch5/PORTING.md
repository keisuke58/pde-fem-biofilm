# Porting this into the thesis repository

A chapter 5 already exists there, so this is **material to merge into it**, not
a file to drop in whole. This session cannot see that repository (GitHub access
here is scoped to `pde-fem-biofilm`), so the skeleton is written to survive any
structure rather than to match a known one.

## What to copy

| From here | To |
|---|---|
| `thesis_ch5/ch5_ansys_contribution.tex` | the chapter file, merged section by section |
| `assets/flow_oliver_solution_loop.png` | wherever that tree keeps figures |
| the `Klempt2024DiffusionDrivenGrowth` and `Klempt2026ContinuumBacterialGrowth` entries from `biofilm_3tooth_refs.bib` | its `.bib`, if not already there |

## The one path to set

```latex
\providecommand{\ansysfigdir}{../assets/}   % <- point at the figure directory
```

Nothing else in the file refers to this repository.

## Why it should merge without a fight

- **No package it does not already have.** Only `amsmath`, `amssymb`, `bm`,
  `xcolor`, `graphicx`. Numbers are written out instead of `siunitx`, and the
  figure enters as a PNG, so **TikZ is not required** in the host tree. (If it
  does load TikZ with `shapes.geometric, arrows.meta, positioning, fit,
  backgrounds, calc`, there is a commented vector route in the figure
  environment — sharper in print.)
- **No label collisions.** Every label is `ch5-*`.
- **No macro collisions.** `\TierA` / `\TierB` / `\ansysfigdir` are
  `\providecommand`, so an existing definition wins rather than erroring.
  (Note for anyone editing: LaTeX macro names cannot contain digits — that is
  why it is `\ansysfigdir` and not `\ch5figdir`.)

## Merging against the chapter that is already there

Since the existing chapter 5 was written around JAXFEM and that material is
being held back for Keio ([`ROADMAP_2026.md`](../ROADMAP_2026.md) §1), expect
to **replace rather than interleave** in most sections. Two places need a real
decision rather than a paste:

1. **The model section (5.2).** If the existing chapter already states the
   growth kinematics, keep its statement and delete the one here — but check
   the citation: the kinematics are `Klempt2024DiffusionDrivenGrowth`, not the
   2026 paper. Citing the newer one here is the natural mistake, and
   `biofilm_3tooth_refs.bib` says so in the entry's own note.
2. **Anything the existing chapter claims about JAXFEM.** It is out of scope
   for this submission by decision, not by accident, so it should come out
   rather than sit alongside. If a sentence is needed to explain the absence,
   the honest one is that the field solution is contributed by the partner
   framework — which is §5.5 here.

## After merging

```
pdflatex <thesis> && bibtex <thesis> && pdflatex <thesis> && pdflatex <thesis>
```

and check that no `ch5-*` reference is undefined. `_build_check.tex` here does
the same thing in isolation, which is the faster place to catch a syntax error
before it reaches the thesis tree.

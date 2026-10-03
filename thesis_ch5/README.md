# Chapter 5 — skeleton and evidence map

Moved 3 Oct 2026: the chapter is now `../thesis/chapters/ch4_ansys.tex`, chapter 4
of the thesis (`../thesis/`). `_build_check.tex` builds it on its own.

`ch5_ansys_contribution.tex` (now `ch4_ansys.tex`) is a frame for the chapter with the facts already
pinned in place. It is not draft prose; the point is that writing becomes
filling in sentences rather than hunting for numbers and file paths.

Scope is fixed by [`ROADMAP_2026.md`](../ROADMAP_2026.md) §1: **Ch5 is the
ANSYS material-law contribution**, and the JAXFEM material is held back for
Keio.

**Porting:** the thesis lives in a different repository and already has a
chapter 5, so this is material to merge into it. See
[`PORTING.md`](PORTING.md) — the short version is that it needs no package the
host tree does not already have (no siunitx, no TikZ), all labels are
namespaced `ch5-*`, and there is exactly one path to set.

## Tier A / Tier B

Sections are marked in the margin. **A** = the claim is already true and
verifiable today, so it can be written now. **B** = waits on runs.

Everything except §5.6 (Results) is Tier A. That is the whole point of the
[Tier A/B split](../ROADMAP_2026.md#2-what-ch5-can-already-claim--today): the
methods-and-verification chapter is complete without a single new run, and
nothing in it waits on the partner framework's schedule.

## Evidence map (rewritten 2026-10-02)

Scope: one and two species, coupled into the partner element following Klempt
et al. 2024; composition from the Klempt et al. 2026 point model. Abaqus work
and more species are the Keio continuation and are not results here
(CLAUDE.md, document rules). **B** blocks hold numbers to re-run before
submission.

| Section | Claim | Where it is established |
|---|---|---|
| 5.1 | Scope, three contributions | scope decision 2026-10-01; `THESIS_ASSIGNMENT.md` §1 |
| 5.2 | Kinematics, `α = α_K − 1`; Klempt stiffness `E(φ) = (φ²+f)E` | `Klempt2024DiffusionDrivenGrowth` Eq. 20, Table 2; `biofilm_material_v01.f`, `tests/test_material_wrapper.py::test_klempt_*` |
| 5.3 | Adapter = core; tangent vs AD (648 cases) | `tests/test_material_wrapper.py`, `tests/test_tangent_quality.py` |
| 5.4 | Closed forms in ANSYS; curved shell; pressure term | `apdl/` decks, `closed_form_reference.py`, `DEVIATOR_SCALING_FINDING.md` |
| 5.5 | Eq. 36 at the Gauss points, once per increment; three framework properties | `apdl/ONE_SPECIES_COUPLING.md`, `apdl/callsite/`, `one_species_reference.py` |
| 5.6 | Point model, amount/composition split, s, φ_cap, φ_min | `JAXFEM/klempt2026_reproduction.py`, `JAXFEM/fritsch2025_cases.py`, `apdl/SPLIT_COUPLING.md`, `composition_reference.py`, `composition_s_sweep.py` |
| 5.7 | Exact solution, whole model, seeded-element stress (Klempt 2024 stiffness), two species; the earlier four-condition study is not in the thesis (decided 2026-10-02) | `apdl/figs_1005.py` and `assets/fig1005_*.png`; `apdl/SPLIT_COUPLING.md` |
| 5.8 | Limitations | this chapter |

The V&V figure is `ch5_flow/flow_vv_thesis.tex` (ANSYS evidence only), not
`flow_vv_hierarchy.tex`, which still lists the Abaqus comparison.

## Two things to get right

**Cite the right Klempt paper.** The growth kinematics
(`F = Fe·Fg`, `Fg = (1+α)I`) are `Klempt2024DiffusionDrivenGrowth`. The
5-species φ/ψ interaction dynamics are `Klempt2026ContinuumBacterialGrowth`.
Different papers; `biofilm_3tooth_refs.bib` says so in the entry's note.

**§5.4.3 is the section that earns the chapter its credibility.** It says
plainly that agreement between two implementations of the same expression is
not correctness, and then demonstrates it on this code. That is a stronger
verification story than an unqualified "0 ULP", not a weaker one.

## Checking it still builds

`_build_check.tex` is a throwaway wrapper — not part of the thesis — that
compiles the skeleton with the figure and bibliography so syntax errors and
broken references surface here rather than in the thesis tree:

```
cd thesis_ch5
pdflatex _build_check && bibtex _build_check && pdflatex _build_check && pdflatex _build_check
```

Currently: builds clean, zero undefined citations or references.

The skeleton needs only amsmath, amssymb, bm, xcolor and graphicx. Numbers are
written out rather than via siunitx, and the figure enters as a PNG rather than
TikZ source, so it adds no dependency to the tree it is merged into.

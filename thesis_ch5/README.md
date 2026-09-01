# Chapter 5 — skeleton and evidence map

`ch5_ansys_contribution.tex` is a frame for the chapter with the facts already
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

## Evidence map

| Section | Claim | Where it is established |
|---|---|---|
| 5.1 | The task and the division of labour | [`THESIS_ASSIGNMENT.md`](../THESIS_ASSIGNMENT.md) §1, [`INTEGRATION_PLAN.md`](../ansys_usermat/INTEGRATION_PLAN.md) |
| 5.2 | Growth kinematics, elastic and viscous branches | `Klempt2024DiffusionDrivenGrowth`; [`RESEARCH_MODEL.md`](../RESEARCH_MODEL.md) |
| 5.3 | USERMAT implementation, Voigt order, tangent | [`usermat_biofilm.f`](../ansys_usermat/usermat_biofilm.f), [`ansys_usermat/README.md`](../ansys_usermat/README.md) |
| 5.4.1 | 0 ULP over 8017 states | [`crosscheck/README.md`](../ansys_usermat/crosscheck/README.md), `crosscheck.py`, `adversarial.py` |
| 5.4.2 | Closed-form growth | [`apdl/`](../ansys_usermat/apdl/) decks + `RUNBOOK.md` |
| 5.4.3 | **What the checks cannot establish** | [`DEVIATOR_SCALING_FINDING.md`](../DEVIATOR_SCALING_FINDING.md) + `check_deviator_scaling.py` |
| 5.5.1 | The partner framework's loop | [`OLIVER_MODEL_NOTES.md`](../ansys_usermat/OLIVER_MODEL_NOTES.md); figure `ch5_flow/flow_oliver_solution_loop.tex` |
| 5.5.2 | The delivered routine is only an adapter | [`biofilm_material_v01.f`](../ansys_usermat/biofilm_material_v01.f), `tests/test_material_wrapper.py::test_wrapper_is_only_an_adapter` |
| 5.5.3 | `Fv` needs 9 slots; `dt` must resolve `η/(2·C10)` | same file's header notes 2–3; `test_the_viscous_step_must_resolve_the_relaxation_time` |
| 5.6 | Results | **Tier B** — awaiting runs |
| 5.7 | Limitations | [`VERIFICATION_SENSITIVITY_LIMITATIONS.md`](../VERIFICATION_SENSITIVITY_LIMITATIONS.md), `DEVIATOR_SCALING_FINDING.md` §7 |

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

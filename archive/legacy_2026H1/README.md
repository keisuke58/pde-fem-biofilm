# archive/legacy_2026H1 — root files of January to June 2026, moved here on 19 Oct 2026

Root-level scripts, reports and result files from the first half of 2026 (Abaqus tooth/jaw models, OpenJaw
assembly, phase-field and viscoelastic studies, early posterior pipelines, PINN and reaction-diffusion demos, the
June report). None of them is imported by the current code, run by the tests, or named in CLAUDE.md, README.md,
REPO_MAP.md or the thesis (checked when they were moved). They are kept for provenance; `git log --follow` shows
their history.

Scripts here import some modules that stayed in the root (e.g. `material_models`, `ecology_constants`); run them
from the repository root with `PYTHONPATH=.` if one is ever needed again.

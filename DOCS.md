# Documentation index

Index of the prose documents by topic. The guide to the whole tree, with the
current work first, is [REPO_MAP.md](REPO_MAP.md); this file lists the
documents behind it, most of them from before October 2026. Checked on
15 November 2026: every path below exists.

## Entry and reference

| Doc | What |
|---|---|
| [README.md](README.md) | Front page: what is here, the model, the machines |
| [REPO_MAP.md](REPO_MAP.md) | Guide to the tree: thesis, coupled model, papers, meetings, machines, Keio continuation |
| [RESEARCH_MODEL.md](RESEARCH_MODEL.md) | Model chain of the earlier tooth pipeline: equations, and what is measured, calibrated or assumed |
| [THESIS_ASSIGNMENT.md](THESIS_ASSIGNMENT.md) | Scope as set by the supervisor (August 2026): roles, deliverables, open questions |
| [THESIS_PLAYBOOK.md](THESIS_PLAYBOOK.md) | Process view as of August 2026 (inventory, write-up order, checks). It describes the earlier pipeline and its notation; the present state of the thesis is in [thesis/README.md](thesis/README.md) and [thesis/REVIEW_2026-11-13.md](thesis/REVIEW_2026-11-13.md) |
| [FEM_README.md](FEM_README.md) | Reference of the line that bridges a dysbiosis index to the stiffness |
| [methods_supplement_fem.md](methods_supplement_fem.md) | Methods supplement of that line (February 2026) |
| [CITATION_AUDIT.md](CITATION_AUDIT.md) | Audit of the references and citations (August 2026) |
| [THIRD_PARTY.md](THIRD_PARTY.md) | Third-party content in the repository and its licences |

## Verification of the earlier pipeline (June 2026)

| Doc | What |
|---|---|
| [VERIFICATION_SENSITIVITY_LIMITATIONS.md](VERIFICATION_SENSITIVITY_LIMITATIONS.md) | Consolidated verification, sensitivity and limitations; read this one first |
| [archive/rigor_audit_growth_2026-06-26.md](archive/rigor_audit_growth_2026-06-26.md) | Audit of the growth kinematics and numerics |
| [DS_composition_fix.md](DS_composition_fix.md) | The composition copy error of the dysbiotic static condition and its fix |
| [espec_sensitivity_findings.md](espec_sensitivity_findings.md) | Sensitivity of the stress ratio to the per-species stiffness |
| [archive/legacy_2026H1/clsm_composition_findings.md](archive/legacy_2026H1/clsm_composition_findings.md) | Where the compositions came from; recomputation from the raw CLSM data |
| [archive/legacy_2026H1/pde3d_rigor_findings.md](archive/legacy_2026H1/pde3d_rigor_findings.md) | Robustness of the 3D growth field against depth and nutrient |
| [archive/legacy_2026H1/eigenstrain_theory_roadmap.md](archive/legacy_2026H1/eigenstrain_theory_roadmap.md) | Consistency of the eigenstrain formulation |

## Plans and overviews

| Doc | What |
|---|---|
| [KEIO_PLAN.ja.md](KEIO_PLAN.ja.md) | The current plan (Keio continuation) |
| [TODO_2026-11.ja.md](TODO_2026-11.ja.md) | What is open in November 2026 |
| [PIPELINE.md](PIPELINE.md) | The config-driven pipeline entry point (`pipeline.py`) and its risk metric |
| [research_goals_1_2.md](research_goals_1_2.md) | Research goals in four levels, as of the first half of 2026 |
| [archive/PLAN_NEXT.md](archive/PLAN_NEXT.md) | Roadmap of July 2026 (superseded) |
| [archive/legacy_2026H1/klempt2024_gap_analysis.md](archive/legacy_2026H1/klempt2024_gap_analysis.md) | Gap analysis against Klempt et al. 2024 (superseded by [KLEMPT2024_REPRODUCTION.md](KLEMPT2024_REPRODUCTION.md)) |
| [archive/legacy_2026H1/BENCHMARK_PLAN.md](archive/legacy_2026H1/BENCHMARK_PLAN.md), [archive/legacy_2026H1/FEM_PLAN.md](archive/legacy_2026H1/FEM_PLAN.md) | Benchmark and FEM plans of the first half of 2026 |
| [archive/legacy_2026H1/related_work_jaw_biofilm.md](archive/legacy_2026H1/related_work_jaw_biofilm.md) | Literature on FEM of jaw and biofilm |
| [archive/legacy_2026H1/overview2602_en.md](archive/legacy_2026H1/overview2602_en.md), [archive/legacy_2026H1/overview2602.md](archive/legacy_2026H1/overview2602.md), [archive/legacy_2026H1/overview_tmcmc_fem_en.md](archive/legacy_2026H1/overview_tmcmc_fem_en.md) | Method overviews of February 2026 (English and Japanese) |

## Subsystem documents

| Doc | What |
|---|---|
| [thesis/README.md](thesis/README.md) | The thesis: structure, rules, build |
| [papers/coupling_paper/manuscript/README.md](papers/coupling_paper/manuscript/README.md) | The coupling paper manuscript |
| [JAXFEM/README.md](JAXFEM/README.md), [JAXFEM/RESULTS.md](JAXFEM/RESULTS.md) | JAX solvers for Klempt 2024 Eqs. 34 to 36 and their results |
| [ansys_usermat/README.md](ansys_usermat/README.md) | ANSYS USERMAT of the growth law |
| [ansys_usermat/apdl/RUNBOOK.md](ansys_usermat/apdl/RUNBOOK.md) | Building the custom ANSYS executable on Windows and running the closed-form checks |
| [ansys_usermat/apdl/README.md](ansys_usermat/apdl/README.md) | The verification decks: status, results table, how to read a failure |
| [ansys_usermat/USERELEM_NOTES.md](ansys_usermat/USERELEM_NOTES.md) | User element against user material in ANSYS |
| [ansys_usermat/crosscheck/README.md](ansys_usermat/crosscheck/README.md) | Comparison of the Abaqus and ANSYS material routines |
| [ansys_usermat/coupling/README.md](ansys_usermat/coupling/README.md) | The bridge between the Fortran material routine and Python |
| [abaqus_composition/README.md](abaqus_composition/README.md), [abaqus_living/README.md](abaqus_living/README.md) | The couplings in Abaqus (Keio continuation) |
| [keio_wp2/README.ja.md](keio_wp2/README.ja.md), [keio_p1/README.ja.md](keio_p1/README.ja.md) | Keio work packages: stress acting back on growth, NEM against Galerkin |

## Archived

Superseded notes, kept as records. Indexes: [archive/README.md](archive/README.md),
[archive/legacy_2026H1/README.md](archive/legacy_2026H1/README.md) and
[ansys_usermat/apdl/archive/README.md](ansys_usermat/apdl/archive/README.md).

| Doc | Superseded by |
|---|---|
| [archive/growth_analysis_doc.md](archive/growth_analysis_doc.md) | VERIFICATION_SENSITIVITY_LIMITATIONS.md |
| [archive/FEM_FILE_MANAGEMENT.md](archive/FEM_FILE_MANAGEMENT.md) | README.md and REPO_MAP.md |
| [archive/eeff_sensitivity_results.md](archive/eeff_sensitivity_results.md) | espec_sensitivity_findings.md |
| [archive/DI_Eeff_sensitivity.md](archive/DI_Eeff_sensitivity.md) | FEM_README.md |

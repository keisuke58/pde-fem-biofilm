# Documentation index

Complete catalogue of the prose documentation in this repository. Start from
[README.md](README.md); this file is the exhaustive map. Subsystem READMEs
([JAXFEM/README.md](JAXFEM/README.md)) document their own folders.

---

## Entry & reference

| Doc | What |
|---|---|
| [README.md](README.md) | Repository entry point: pipeline, verified-vs-assumed, two lineages, key scripts |
| [FEM_README.md](FEM_README.md) | DI-bridge FEM lineage reference (DI → E(DI), transverse isotropy; nominal GPa scale) |
| [methods_supplement_fem.md](methods_supplement_fem.md) | Methods supplement: DI timepoint, E(φ) physics, TMCMC → Monod coupling |

## Rigor audit (2026-06-26)

The authoritative verification of the Klempt growth-stress pipeline. Read the
consolidated document first; the others are the underlying investigations.

| Doc | What |
|---|---|
| [VERIFICATION_SENSITIVITY_LIMITATIONS.md](VERIFICATION_SENSITIVITY_LIMITATIONS.md) | **Consolidated** Verification / Sensitivity / Limitations — read first |
| [rigor_audit_growth_2026-06-26.md](rigor_audit_growth_2026-06-26.md) | Growth kinematics & numerics audit |
| [DS_composition_fix.md](DS_composition_fix.md) | Dysbiotic-static composition copy-bug (σ_DS ~8× too high) — fix status |
| [clsm_composition_findings.md](clsm_composition_findings.md) | Composition provenance; raw-CLSM recomputation of σ_CH/σ_DH |
| [espec_sensitivity_findings.md](espec_sensitivity_findings.md) | `E_SPEC` sensitivity of σ_CH/σ_DH (3.7–12×) |
| [pde3d_rigor_findings.md](pde3d_rigor_findings.md) | 3D growth-PDE depth/nutrient robustness (5.3–6.6×) |
| [eigenstrain_theory_roadmap.md](eigenstrain_theory_roadmap.md) | α-eigenstrain theory consistency |

## Theory & planning

| Doc | What |
|---|---|
| [PLAN_NEXT.md](PLAN_NEXT.md) | **Consolidated next-steps roadmap (2026-07-02)** — T1 thesis / T2 first paper / T3 continuation, prioritized |
| [PIPELINE.md](PIPELINE.md) | Config-driven pipeline entry point (T2.1) + `P[σ>threshold]` risk metric (T2.3); isolated-checkout audit notes |
| [klempt2024_gap_analysis.md](klempt2024_gap_analysis.md) | Klempt et al. (2024) gap analysis & response |
| [BENCHMARK_PLAN.md](BENCHMARK_PLAN.md) | Klempt 2024 → 5-species benchmark roadmap |
| [FEM_PLAN.md](FEM_PLAN.md) | FEM pipeline plan — 3-tooth biofilm assembly |
| [research_goals_1_2.md](research_goals_1_2.md) | Research goals (levels 1–4) |
| [related_work_jaw_biofilm.md](related_work_jaw_biofilm.md) | Jaw-biofilm FEM / digital-twin literature |

## Overviews

| Doc | What |
|---|---|
| [overview2602_en.md](overview2602_en.md) | Method overview — big picture (EN) |
| [overview2602.md](overview2602.md) | Method overview (JA; rendered by `visualize_overview.py`) |
| [overview_tmcmc_fem_en.md](overview_tmcmc_fem_en.md) | TMCMC + FEM integrated overview |

## Subsystem docs

| Doc | What |
|---|---|
| [JAXFEM/README.md](JAXFEM/README.md) | JAX PDE reproduction suite (Klempt Eq. 34–36) |
| [JAXFEM/RESULTS.md](JAXFEM/RESULTS.md) | JAXFEM result notes |

## Archived

Historical working notes, kept for provenance and superseded by the docs above.
See [archive/README.md](archive/README.md).

| Doc | Superseded by |
|---|---|
| [archive/growth_analysis_doc.md](archive/growth_analysis_doc.md) | rigor audit + VERIFICATION doc |
| [archive/FEM_FILE_MANAGEMENT.md](archive/FEM_FILE_MANAGEMENT.md) | README + this index |
| [archive/eeff_sensitivity_results.md](archive/eeff_sensitivity_results.md) | espec_sensitivity_findings.md |
| [archive/DI_Eeff_sensitivity.md](archive/DI_Eeff_sensitivity.md) | FEM_README.md (DI lineage) |

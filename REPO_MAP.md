# Repository map

A navigation guide to the code — the top level is large, so start here.
For the documentation index see [`DOCS.md`](DOCS.md); for the rendered project
overview see the [project site](https://keisuke58.github.io/pde-fem-biofilm/).

## Entry points

| File | What |
|---|---|
| [`pipeline.py`](pipeline.py) | Config-driven pipeline entry point (see [`PIPELINE.md`](PIPELINE.md)) |
| [`JAXFEM/audit_all.py`](JAXFEM/audit_all.py) | All-in-one thesis-quality audit (`--quick` / `--strict` / `--strict-env`) |
| [`validate_composition.py`](validate_composition.py) | Model ↔ Heine experiment composition validation (figure + metrics) |

## Source papers, read and reconciled (2026-09-30)

Notes taken from the papers this work builds on, each recording what it settles
here rather than summarising the paper. [`CLAIMS_AND_EVIDENCE.md`](CLAIMS_AND_EVIDENCE.md)
is the consolidation; these are the detail behind it.

| File | What it settles |
|---|---|
| [`KLEMPT2024_REPRODUCTION.md`](KLEMPT2024_REPRODUCTION.md) | Why the 2024 PDE is **not** reproduced: the variants that fit are not the paper's, and its own combination is the worst fit |
| [`SOLEIMANI2023_NOTES.md`](SOLEIMANI2023_NOTES.md) | The nearest precedent, previously uncited. Our condition degeneracy belongs to the **simplex**, not to multi-species modelling; the UserElement question; the group's validation bar |
| [`FEI2020_NOTES.md`](FEI2020_NOTES.md) | Fei et al. 2020 (PNAS): growth-induced stress and wrinkling of V. cholerae biofilms, F = Fe·Fg with Monod growth, compared quantitatively with experiment; a validation target for the mechanics. PDF not bundled (PNAS licence) |
| [`SOLEIMANI2019_NOTES.md`](SOLEIMANI2019_NOTES.md) | The unconditionally stable viscous integrator (Eq. 32) our explicit `Fv` update needs; and `E = 10 Pa`, which qualifies the stress comparison |
| [`SOLEIMANI2021_NOTES.md`](SOLEIMANI2021_NOTES.md) | The Heaviside cap on `α` (Eq. 17), now implemented; resolves `CITATION_AUDIT.md` F1c; the group's own warning about the advection term |
| [`CHU2018_NOTES.md`](CHU2018_NOTES.md) | What the computed stress is *for*: the ~5 kPa threshold at which the bacterial stress response turns on |
| [`READING_GAPS.md`](READING_GAPS.md) | References Klempt 2024 rests on that this repository does not cite |

## Analysis lineages

- **Klempt growth-stress pipeline** (thesis headline) — `gen_tooth_klempt_umat_inp.py`,
  `umat_klempt_alpha.f`, the `JAXFEM/` PDE α-field, `run_tooth_klempt*.sh`.
- **DI-bridge FEM** (second lineage) — `material_models.py` (`compute_E_di*`),
  `fem_*_extension.py`, see [`FEM_README.md`](FEM_README.md).

## Constitutive model (UMAT / USERMAT)

| File | What |
|---|---|
| `umat_biofilm_visco.f`, `umat_biofilm_visco_2ch.f`, `umat_biofilm_visco_phase2.f` | Verified Abaqus viscoelastic UMATs |
| `usdfld_biofilm.f` | USDFLD growth-driver field routine |
| [`ansys_usermat/`](ansys_usermat/) | ANSYS USERMAT port + `crosscheck/` (dual-solver equivalence, 0 ULP) |
| [`abaqus_composition/`](abaqus_composition/) | The composition coupling of the ANSYS runs (point model through the Python bridge) in an Abaqus UMAT built from the same fragments; one-element checks in Abaqus 2024 (5 Oct 2026). Keio continuation, not thesis work |
| [`keio_report/`](keio_report/) | Keio 課題研究報告 (Japanese, 6 pages, JSCES two-column style): the thesis content (TMCMC, ANSYS growth law and composition coupling, mesh convergence) for the Keio side. Build: `platex` x3 + `pbibtex` + `dvipdfmx`. Master copy also in the private luh_summer_2026 repo |
| [`keio_wp2/`](keio_wp2/) | Keio WP2, stress into growth: the mechanical term of Klempt 2024 Eq. 30 kept, 2D plane-strain prototype in Python; the dimensionless number that decides whether it matters, and time integration (Japanese notes). Keio continuation, not thesis work |
| [`keio_p1/`](keio_p1/) | Keio P1 (numerics paper): the growth field solved with NEM (explicit, weighted-least-squares Laplacian at the integration points) and with Galerkin finite elements (implicit), mesh and time-step study in Python. Keio continuation, not thesis work |
| `material_models.py` | Python material model (E(φ), E(DI), viscoelastic) |

## PDE / ecology model (JAX)

`jax_hamilton_*_5species_demo.py`, `jax_*_reaction_diffusion_*.py`,
`multiscale_coupling_*.py`, [`JAXFEM/`](JAXFEM/) (Klempt Eq. 34–36 testbed) —
require `jax[cpu]` (not pinned in `requirements.txt`).

## Plotting & figures

- `plot_*.py`, `generate_*figure*.py` — matplotlib result/analysis figures → `assets/`.
- TikZ figure libraries: `umat_flow/`, `ch5_flow/`, `JAXFEM/algo_flow*.tex` → `assets/`.

## Data, tests, CI

| Path | What |
|---|---|
| [`data/`](data/) | Experimental data (Heine species-distribution workbook) |
| [`configs/`](configs/) | Pipeline configs |
| [`tests/`](tests/) | Unit suite (`pytest tests/`) |
| `pytest.ini`, `requirements.txt`, `.github/workflows/ci.yml` | Test scoping, pinned deps, CI |
| `.claude/hooks/session-start.sh` | Restores what a cloud container loses on restart: git identity, the two anti-attribution hooks, pytest deps, gfortran |

## Thesis, plan and handover (added 2026-09-01)

| Path | What it is |
|---|---|
| `ROADMAP_2026.md` ([日本語](ROADMAP_2026.ja.md)) | Submission Nov 2026, defence Dec. The Tier A/B split, the cadence with the supervisors, week by week |
| `ROADMAP_TWO_WAY.md` | From the thesis's one-way coupling (field → point model) to a two-way one: four steps ordered by literature support, the first two proposed for Keio |
| `QA_1005.md` | Likely questions and short answers (English, with Japanese notes) for the 5 Oct 2026 meeting with Dr. Soleimani and Assoc. Prof. Muramatsu |
| `MEETING_2026-10-05.md` | Summary of the 5 Oct 2026 meeting with Dr. Soleimani and Assoc. Prof. Muramatsu: journal suggestions, the two-way nutrient exchange the supervisor expects, submission format, contacts after the return to Japan, to-do list |
| `KEIO_PLAN.ja.md` | Plan for the Keio continuation (Japanese): what Felix Klempt's and Meisam Soleimani's mails of 5-6 Oct 2026 settled, the tools in hand (ANSYS, Abaqus, Fortran point model), the Klempt 2024 reproduction in Abaqus, and the work packages in order |
| `KEIO_SERVER_HANDOFF.ja.md` | Hand-off to the Keio Linux server (Japanese, 7 Oct 2026): how to run the Keio work packages with Abaqus there, and which steps were only checked on IKMHIWI03 |
| `KEIO_SERVER_ENV.ja.md` | The fifa environment as measured (Japanese, 8 Oct 2026): hardware, disk and scratch location, parallel mode, package versions against `requirements.txt`, and `scripts/run_chain_keio.py`'s guards for the shared server |
| `CLOUD_TO_FIFA.ja.md` | Instructions from the cloud session to the Keio server (fifa), newest first (Japanese): where to report, following master, how the next runs are queued |
| `KEIO_CHAIN_LOG.md` | Record of what the cloud session decided fifa should run next, and why, one entry per decision, newest at the bottom |
| `scripts/setup_keio_server.sh` | One-time setup of a clone on the Keio Linux server (fifa): git identity and hooks, checks for Python, numpy/scipy and the Abaqus command |
| `references/` | Klempt et al. 2024 (BMMB), Soleimani et al. 2023 (Sci. Rep., co-aggregation) and Feng et al. 2021 (Bull. Math. Biol., Streptococcus–Veillonella), all CC BY 4.0: PDFs and searchable text extractions; licence in `THIRD_PARTY.md`. NEM papers (Rudolf 2025, Vogel 2020, Junker 2022, von Zabiensky 2026; open access) with the reading guide `references/NEM_PAPERS.md`, which also summarises the three that are not open access |
| `PAPER_CHECK_KLEMPT2024.md` | Chapter 5 and the decks checked against Klempt et al. 2024: notation now as in the papers, contradictions found and fixed, what is still open (neighbours in tension?) and the Monday IKMHIWI03 commands |
| `COUPLING_STATUS.md` | One page: what is done and what is not, for one and for two species (ANSYS model, checks, open items, more species) |
| `thesis/` | **The thesis itself** (from 3 Oct 2026; copied from LUH_summer_2026). `main.tex`, `chapters/`, structure in `thesis/README.md`; the ANSYS chapter is `chapters/ch4_ansys.tex` |
| `thesis_ch5/` | Evidence map and a stand-alone build check (`_build_check.tex`) of the ANSYS chapter, which now lives in `thesis/chapters/ch4_ansys.tex` |
| `handover/` | The self-contained package for the partner group — generated from the sources under test by `make_handover.py`, so it cannot drift |
| `reports/` | Written progress updates to the supervisors, kept next to the work they describe |
| `DEVIATOR_SCALING_FINDING.md` | A mis-scaled isochoric split in the verified core: a pure pressure error, von Mises unaffected. Documented, not fixed — with the reasoning |
| `VISCOUS_UPDATE_SCHEME.md` | What the `Fv` update actually is, why "backward Euler" is the wrong name for it, and the step limit that follows |
| `ansys_usermat/FELIX_VS_PARTNER_DIFF.md` | Felix Klempt's current implementation against the partner's element, in words only (his code stays outside git): reading slip fixed, \|grad phi\| front term, one-sided saturating growth, inner time step |
| `ansys_usermat/apdl/ONE_SPECIES_COUPLING.md` | The one-species coupling into the partner's element (decided 2026-10-01): Klempt 2024 Eq. 36 at the Gauss point, the call-site recipe, the three questions that must be answered first, and how to verify the first run |
| `PDE_VERIFICATION_FINDINGS.md` | Code verification of the nutrient/species solvers against exact solutions (Thiele cosh, Klempt Eq. 35's parabola, a Neumann eigenfunction). The zero-flux wall cost one order in 1D and did not converge at all in 2D; fixed in all three places, and the effect on the reported 2D condition spread measured (~1e−6, so that result stands) |
| `E_SATURATION_FINDING.md` | The production φ→E bridge clips to [10, 1000] Pa, and at the current calibration that bound is active over much of the healthy composition space — so distinct conditions can report identical stiffness |
| `PINN_DESIGN.md` | A physics-informed surrogate, written up as a Keio design rather than started |
| `ansys_usermat/biofilm_material_v01.f` | `BIOFILM_GROWTH_VISCO_V01` — the routine handed over, an adapter around the verified core |
| `ansys_usermat/growth_law_verification.ipynb` | Executable walkthrough of what the verifications establish and what they do not |
| `ansys_usermat/apdl/closed_form_reference.py` | Closed form for the two growth cases, derived independently of the implementation |
| `ansys_usermat/apdl/check_deck.py` | Pre-flight for any deck whose stress will be reported. Catches only failures that are **silent in ANSYS**: a step too coarse for the viscous relaxation time, too few `TB,STATE` slots (which leaves α unread, so the solve runs purely elastic), α declared but never written, and an over-long `TBDATA` whose tail APDL drops. Handles Abaqus `.inp` too, where growth arrives as the temperature field instead |
| `ansys_usermat/apdl/make_layered_material.py` | Turns a depth-resolved α(x) field into ANSYS layered materials, and reports what the binning costs against the field — the ANSYS route takes α per material, so a spatial field has to be discretised |
| `RESEARCH_IDEAS.ja.md` | Research ideas after the thesis (10 Oct 2026), mapped to the Keio papers: a symmetric neighbour choice for the NEM stencil, the short-run stability artefact, implicit NEM against Galerkin, the calibrated five-species point model in the element, posterior propagation to stresses, measured stiffness and geometry |
| `ansys_usermat/apdl/nem_stability.py`, `nem_finite_cube.py`, `wp3_reference.py`, `wp3_stress_reference.py`, `wp3_stress_gauss.py`, `wp3_alpha_profile.py`, `wp3_figs.py` | WP3 seed-growth numerics of the partner element: the Bloch stability limit of the explicit NEM update (Appendix D of the thesis), a finite-volume reference solution of the same problem with the NEM error per mesh against it, the seed stress for the reference field with element and integration-point values, the mesh-independent measures of α − 1 and the two figures of Chapter 4 |
| `ansys_usermat/apdl/MESH_STUDY.md`, `make_mesh_levels.py` | The light mesh study: show the **ratio** Ch5 reports is mesh-stable, rather than converging absolute stress. The helper halves `ESIZE` and changes nothing else |
| `ansys_usermat/apdl/plot_cylinder_3d.py` | Draws the v222 curved-shell run from the solver's own listing — geometry, a section showing the two layers and the interface, and the SEQV distribution. Validated against the listing's own min/max/mean |
| [`ansys_usermat/coupling/`](ansys_usermat/coupling/) | Gauss-point bridge to Python: `usermat()` calls out over a socket per Gauss point instead of using the inline Fortran core, verified end-to-end (gfortran driver + real ANSYS) for both the material law (`kUsePy`) and the 0D Hamilton ecology ODE that drives α (`kUseEcology`). See `coupling/README.md` for the interface contract, verification status, and open next steps (TMCMC theta / CLSM-seeded state are still placeholders) |

## Pre-existing directories the tour had not listed

| Path | What it is |
|---|---|
| `tier2b_real/` | Abaqus coupon/implant job generation and real-geometry meshing (`implant_coupon.py`, `mesh_bone_region.py`, STL/mesh assets) |
| `umat_tangent_test/` | Single-element UMAT tangent and eigenstrain cross-check harnesses, including the dual-UMAT growth comparison in `xcheck_eigenstrain/XCHECK_RESULTS.md` |
| `ANSYS_ENVIRONMENT.md` | Full hardware, licence and product inventory for the IKMHIWI03 ANSYS machine |
| [`runs/`](runs/) | Per-run validation logs / env configs (provenance) |

## Documentation

`README.md`, [`DOCS.md`](DOCS.md) (full index),
[`VERIFICATION_SENSITIVITY_LIMITATIONS.md`](VERIFICATION_SENSITIVITY_LIMITATIONS.md)
(read first for what is verified vs assumed), `PLAN_NEXT.md`, `methods_supplement_fem.md`.
Historical notes live under [`archive/`](archive/).

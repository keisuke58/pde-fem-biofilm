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
| `thesis_ch5/` | Chapter 5 skeleton with an evidence map, plus `PORTING.md` for merging it into the thesis repository |
| `handover/` | The self-contained package for the partner group — generated from the sources under test by `make_handover.py`, so it cannot drift |
| `reports/` | Written progress updates to the supervisors, kept next to the work they describe |
| `DEVIATOR_SCALING_FINDING.md` | A mis-scaled isochoric split in the verified core: a pure pressure error, von Mises unaffected. Documented, not fixed — with the reasoning |
| `VISCOUS_UPDATE_SCHEME.md` | What the `Fv` update actually is, why "backward Euler" is the wrong name for it, and the step limit that follows |
| `E_SATURATION_FINDING.md` | The production φ→E bridge clips to [10, 1000] Pa, and at the current calibration that bound is active over much of the healthy composition space — so distinct conditions can report identical stiffness |
| `PINN_DESIGN.md` | A physics-informed surrogate, written up as a Keio design rather than started |
| `ansys_usermat/biofilm_material_v01.f` | `BIOFILM_GROWTH_VISCO_V01` — the routine handed over, an adapter around the verified core |
| `ansys_usermat/growth_law_verification.ipynb` | Executable walkthrough of what the verifications establish and what they do not |
| `ansys_usermat/apdl/closed_form_reference.py` | Closed form for the two growth cases, derived independently of the implementation |
| `ansys_usermat/apdl/check_deck.py` | Pre-flight for any deck whose stress will be reported. Catches only failures that are **silent in ANSYS**: a step too coarse for the viscous relaxation time, too few `TB,STATE` slots (which leaves α unread, so the solve runs purely elastic), α declared but never written, and an over-long `TBDATA` whose tail APDL drops. Handles Abaqus `.inp` too, where growth arrives as the temperature field instead |
| `ansys_usermat/apdl/make_layered_material.py` | Turns a depth-resolved α(x) field into ANSYS layered materials, and reports what the binning costs against the field — the ANSYS route takes α per material, so a spatial field has to be discretised |
| `ansys_usermat/apdl/MESH_STUDY.md`, `make_mesh_levels.py` | The light mesh study: show the **ratio** Ch5 reports is mesh-stable, rather than converging absolute stress. The helper halves `ESIZE` and changes nothing else |
| `ansys_usermat/apdl/plot_cylinder_3d.py` | Draws the v222 curved-shell run from the solver's own listing — geometry, a section showing the two layers and the interface, and the SEQV distribution. Validated against the listing's own min/max/mean |

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

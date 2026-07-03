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
| [`runs/`](runs/) | Per-run validation logs / env configs (provenance) |

## Documentation

`README.md`, [`DOCS.md`](DOCS.md) (full index),
[`VERIFICATION_SENSITIVITY_LIMITATIONS.md`](VERIFICATION_SENSITIVITY_LIMITATIONS.md)
(read first for what is verified vs assumed), `PLAN_NEXT.md`, `methods_supplement_fem.md`.
Historical notes live under [`archive/`](archive/).

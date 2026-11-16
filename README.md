# pde-fem-biofilm

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Growth and stress of oral biofilms, built on the continuum growth model of
Klempt et al. (2024), coupled at every point to the multi-species point model
of Klempt et al. (2026).**
口腔バイオフィルムの成長と応力：Klempt et al. (2024) の連続体成長モデルを、各点で Klempt et al. (2026) の多菌種の点モデルと連成する。

Companion code of the LUH / IKM master's thesis of Keisuke Nishioka (PDF to
the supervisor in December 2026, colloquium 16 December 2026), of the coupling
paper that follows from it, and of the continuation at Keio University
(Muramatsu lab, 2027). State as of **15 November 2026**.

---

## What is here

| Part | Where | State |
|---|---|---|
| **The thesis** (LaTeX, 89 pages) | [`thesis/`](thesis/README.md), check of 13 Nov in [`thesis/REVIEW_2026-11-13.md`](thesis/REVIEW_2026-11-13.md) | written; read-through and PDF to the supervisor in December 2026 |
| Species interactions from in-vitro data: TMCMC on the five-species data of Heine et al. (2025) | `thesis/chapters/ch3_heine.tex` (the runs are in a separate repository) | thesis chapter 3, final since 19 Oct |
| Klempt 2024 growth law in the partner's ANSYS user element (NEM), one and two species | [`ansys_usermat/`](ansys_usermat/), [`ansys_usermat/apdl/CONFIRMED_2026-10-22.ja.md`](ansys_usermat/apdl/CONFIRMED_2026-10-22.ja.md) | thesis chapter 4, verified on ANSYS 2022 R2 |
| The coupling agreed on 20 Oct 2026: growth field and point model exchange nutrient, antibiotic and living volume fraction | [`MEETING_2026-10-20.md`](MEETING_2026-10-20.md), Python code [`ansys_usermat/apdl/coupled_living_consumption.py`](ansys_usermat/apdl/coupled_living_consumption.py), results guide [`README_coupled_living.md`](ansys_usermat/apdl/results/2026-10-3d/README_coupled_living.md) | thesis chapter 4; grid series up to 256³ |
| The coupling paper | [`papers/coupling_paper/`](papers/coupling_paper/OUTLINE.md) (manuscript in `manuscript/`) | draft |
| The same coupling as an Abaqus user element, and the earlier composition coupling as an Abaqus UMAT | [`abaqus_living/`](abaqus_living/README.md), [`abaqus_composition/`](abaqus_composition/README.md) | Keio continuation |
| Reproduction of the Klempt 2024 test cases | [`KLEMPT2024_REPRODUCTION.md`](KLEMPT2024_REPRODUCTION.md) | case 4.1 reproduced with consumption g φ c |
| Stress acting back on growth: homeostatic pressure, and the mechanical term of Klempt 2024 Eq. 30 | [`keio_wp2/`](keio_wp2/README.ja.md), [`MECH_K_SWEEP.md`](ansys_usermat/apdl/results/2026-10-3d/MECH_K_SWEEP.md) | Keio continuation |
| NEM against Galerkin | [`keio_p1/`](keio_p1/README.ja.md), `ansys_usermat/apdl/RUN_WP3_IKMHIWI03.md` | Keio continuation |

The field model is that of Klempt et al. (2024), Eqs. 34 to 36, with the
consumption by the living biomass agreed on 20 Oct 2026:

```
φ̇_3D = β ∇²φ_3D + k_α α − r c/(k+c) n_c·∇φ_3D,    n_c = ∇c/|∇c|
ċ    = d ∇²c − g φ_3D (Σ_i φ̄_i) c
α̇    = k_α φ_3D,     F = F_e F_g,  F_g = α I,  α(0) = 1     (growth written α − 1)
```

The point model of Klempt et al. (2026) runs at every point, used as
published. It receives the nutrient and the antibiotic from the 3D fields, and
its living volume fraction Σ_i φ̄_i (φ̄_i = φ_i ψ_i) sets the consumption in
the nutrient equation. The printed equations are the reference; if this
summary and [`MEETING_2026-10-20.md`](MEETING_2026-10-20.md) differ, the
meeting note with its term-by-term check is right. Notation follows the papers.
Parameter values and their conversion to the 2 mm cube are in
[`references/PARAMETER_TABLES.md`](references/PARAMETER_TABLES.md); converted
values are assumptions and are reported with their sensitivity.

---

## Quickstart

```bash
pip install -r requirements.txt     # versions pinned to those on the Keio server (see the file)
pytest tests/                       # unit tests; run after code changes (CI runs by hand only)
cd thesis && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

The coupled model in Python (needs `jax[cpu]` and `scipy`; the main run of
the separate example with antibiotic on 64³ cells, as in the docstring of the
script, here with the point model in JAX; `--fortran` uses the compiled one):

```bash
python ansys_usermat/apdl/coupled_living_consumption.py --nf 64 --fast --front 1 --gc 1 --gab 1 \
    --alpha0 100 --tab 0.333333 --phi-min 0.001 --tvd --wall-grad mirror --v2 --ansys-bc \
    --inherit-living --transient-c 1e8 --only separate_ab
```

Without options the script uses older bring-up values, so start from this
command. The docstring explains every option; which result files are the main
runs is in
[`README_coupled_living.md`](ansys_usermat/apdl/results/2026-10-3d/README_coupled_living.md).

Abaqus on Linux (the Keio server):

```bash
bash scripts/setup_keio_server.sh   # once per clone: git identity, hooks, checks
abaqus_composition/run_comp.sh abaqus_composition/<job>.inp 2sp_case6 4
```

ANSYS (only on the LUH machine IKMHIWI03, until the end of December 2026): see
[`ansys_usermat/apdl/RUNBOOK.md`](ansys_usermat/apdl/RUNBOOK.md).

---

## Machines and where work is pushed

| Machine | What runs there | Read first |
|---|---|---|
| **IKMHIWI03** (LUH, Windows) | ANSYS with the partner element, Abaqus | `CLAUDE.md`, [`CONFIRMED_2026-10-22.ja.md`](ansys_usermat/apdl/CONFIRMED_2026-10-22.ja.md), [`CLOUD_TO_IKMHIWI03.ja.md`](ansys_usermat/apdl/CLOUD_TO_IKMHIWI03.ja.md) |
| **fifa** (Keio, Linux, shared) | Abaqus 2024, Python; no ANSYS | [`KEIO_SERVER_HANDOFF.ja.md`](KEIO_SERVER_HANDOFF.ja.md), [`KEIO_SERVER_ENV.ja.md`](KEIO_SERVER_ENV.ja.md), [`CLOUD_TO_FIFA.ja.md`](CLOUD_TO_FIFA.ja.md) |
| Cloud session | reading, writing, Python, merges into `master` | `CLAUDE.md` |

`master` here is the reference. Results from the machines arrive on their own
branches and are merged by the cloud session. The Muramatsu lab keeps a
private mirror (`mmc-research-group/nishioka-biofilm-fem`, branch `main`).

---

## Start here

| File | What |
|---|---|
| [`REPO_MAP.md`](REPO_MAP.md) | Guide to the whole tree; its first section says where to look for what |
| [`TODO_2026-11.ja.md`](TODO_2026-11.ja.md) | What is open this month (Japanese) |
| [`thesis/README.md`](thesis/README.md), [`thesis/REVIEW_2026-11-13.md`](thesis/REVIEW_2026-11-13.md) | Structure and rules of the thesis; the check of 13 Nov with the open points |
| [`MEETING_2026-10-20.md`](MEETING_2026-10-20.md) | The agreed coupling and the check of its equations |
| [`references/PARAMETER_TABLES.md`](references/PARAMETER_TABLES.md) | Parameter values of the papers and their conversion |
| [`papers/coupling_paper/OUTLINE.md`](papers/coupling_paper/OUTLINE.md) | The coupling paper |
| [`KEIO_PLAN.ja.md`](KEIO_PLAN.ja.md) | Keio continuation: goals and deadlines (§0), paper targets (§1), work packages (§3) |
| [`KLEMPT2024_REPRODUCTION.md`](KLEMPT2024_REPRODUCTION.md) | What reproduces from Klempt 2024 and with which departures |
| [`ansys_usermat/apdl/FRONT_TERM_FIX.md`](ansys_usermat/apdl/FRONT_TERM_FIX.md) | Findings in the partner element (surface points, the term −r c/(k+c) n_c·∇φ of Eq. 34) |
| [`DOCS.md`](DOCS.md) | Index of the older prose documents; superseded notes in [`archive/`](archive/) |

---

## Earlier work in this repository (not in the thesis)

Before October 2026 the repository carried a tooth pipeline in Abaqus that
compared von Mises stress across four culture conditions (commensal / dysbiotic,
static / HOBIC), and a second line that bridged a dysbiosis index to the
stiffness. Both are kept for reference; their numbers are not thesis results.
Their state and limitations are in
[`VERIFICATION_SENSITIVITY_LIMITATIONS.md`](VERIFICATION_SENSITIVITY_LIMITATIONS.md)
and [`FEM_README.md`](FEM_README.md); the method figures are in
[`umat_flow/`](umat_flow/README.md) and [`ch5_flow/`](ch5_flow/README.md).
`reproduce.sh` regenerates that line's figures and metrics.

---

## References

- F. Klempt, M. Soleimani, P. Wriggers, P. Junker: A Hamilton principle-based model for diffusion-driven biofilm growth. *Biomech. Model. Mechanobiol.* 23 (2024) 2091–2113. In [`references/`](references/) (CC BY 4.0).
- F. Klempt, H. Geisler, M. Soleimani, P. Junker: A continuum multi-species bacterial growth model with a novel interaction scheme. *Arch. Appl. Mech.* 96 (2026) 164. In [`references/`](references/) (CC BY 4.0).
- T. Rudolf, F. Klempt, H. I. Kök, M. Soleimani, D. R. Jantos, P. Junker: Computational efficiency and accuracy of the Neighbored Element Method. *Finite Elem. Anal. Des.* 249 (2025) 104353. In [`references/`](references/).
- N. Heine et al.: Influence of species composition and cultivation condition on peri-implant biofilm dysbiosis in vitro. *Front. Oral Health* (2025). The CLSM and viability data.
- M. Soleimani et al. (2023, *Sci. Rep.*) and D. Feng et al. (2021, *Bull. Math. Biol.*): two-species biofilm models, in [`references/`](references/).

Full bibliography: `thesis/references.bib`. Cite this repository via
[`CITATION.cff`](CITATION.cff); third-party content and its licences are listed
in [`THIRD_PARTY.md`](THIRD_PARTY.md).

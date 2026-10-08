# pde-fem-biofilm

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Growth and stress of oral biofilms, built on the continuum growth model of
Klempt et al. (2024), with the species composition from a point model at the
Gauss points.**
口腔バイオフィルムの成長と応力：Klempt et al. (2024) の連続体成長モデルに、Gauss 点ごとの点モデルで菌種の組成を与える。

Companion code of the LUH / IKM master's thesis of Keisuke Nishioka (submission
October 2026) and of its continuation at Keio University (Muramatsu lab,
2027). State as of **8 October 2026**.

---

## What is here

| Part | Where | State |
|---|---|---|
| **The thesis** (LaTeX, 75 pages) | [`thesis/`](thesis/README.md) | content frozen 12 Oct 2026, submission 20 Oct 2026 |
| Species interactions from in-vitro data: TMCMC on the five-species data of Heine et al. (2025) | `thesis/chapters/ch3_heine.tex`, [`JAXFEM/`](JAXFEM/README.md) | thesis chapter 3 |
| Klempt 2024 growth law in the partner's ANSYS user element (NEM), one and two species, composition from the Klempt et al. (2026) point model | [`ansys_usermat/`](ansys_usermat/), [`COUPLING_STATUS.md`](COUPLING_STATUS.md) | thesis chapter 4, verified on ANSYS 2022 R2 |
| The same model in Abaqus (UMAT, UMATHT, UEL, point model in Fortran) | [`abaqus_composition/`](abaqus_composition/README.md) | Keio continuation; identical on Windows and Linux |
| Reproduction of the Klempt 2024 test cases | [`KLEMPT2024_REPRODUCTION.md`](KLEMPT2024_REPRODUCTION.md) | case 4.1 reproduced with consumption g φ c |
| Stress-dependent growth (homeostatic pressure), implant collar and tooth crown | [`keio_wp2/`](keio_wp2/README.ja.md) | Keio WP2: Python prototype and 12 Abaqus runs agree |
| NEM against Galerkin | [`keio_p1/`](keio_p1/README.ja.md), `ansys_usermat/apdl/RUN_WP3_IKMHIWI03.md` | Keio WP3: runs in progress |

The model, as in Klempt et al. (2024):

```
φ̇ = β ∇²φ + k_α α (+ front term),     c: quasi-static, consumption g φ c,
α̇ = k_α φ,     F = F_e F_g,  F_g = α I,  α(0) = 1     (growth written α − 1)
```

The species composition φ₁/(φ₁+φ₂) at each Gauss point comes from the point
model of Klempt et al. (2026); the field decides the amount, the point model
the shares. Notation follows the papers (`CLAUDE.md`, "Notation as in the
papers"). β = 0.02 mm²/T* is Table 2 of Klempt 2024 converted to the 2 mm
cube; the conversion is an assumption and is reported with its sensitivity.

---

## Quickstart

```bash
pip install -r requirements.txt     # versions pinned to those on the Keio server (see the file)
pytest tests/                       # unit tests; run before every push (CI runs by hand only)
cd thesis && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

Abaqus on Linux (the Keio server):

```bash
bash scripts/setup_keio_server.sh   # once per clone: git identity, hooks, checks
abaqus_composition/run_comp.sh abaqus_composition/<job>.inp 2sp_case6 4
```

ANSYS (only on the LUH machine IKMHIWI03, until mid-December 2026): see
[`ansys_usermat/apdl/RUNBOOK.md`](ansys_usermat/apdl/RUNBOOK.md).

---

## Machines and where work is pushed

| Machine | What runs there | Read first |
|---|---|---|
| **IKMHIWI03** (LUH, Windows) | ANSYS with the partner element, Abaqus | `CLAUDE.md`, [`ansys_usermat/apdl/IKMHIWI03_TO_SUBMISSION.ja.md`](ansys_usermat/apdl/IKMHIWI03_TO_SUBMISSION.ja.md) |
| **fifa** (Keio, Linux, shared) | Abaqus 2024, Python; no ANSYS | [`KEIO_SERVER_HANDOFF.ja.md`](KEIO_SERVER_HANDOFF.ja.md), [`KEIO_SERVER_ENV.ja.md`](KEIO_SERVER_ENV.ja.md), [`CLOUD_TO_FIFA.ja.md`](CLOUD_TO_FIFA.ja.md) |
| Cloud session | reading, writing, Python, merges into `master` | `CLAUDE.md` |

`master` here is the reference. Results from the machines arrive on their own
branches and are merged by the cloud session. The Muramatsu lab keeps a
private mirror (`mmc-research-group/nishioka-biofilm-fem`, branch `main`).

---

## Start here

| File | What |
|---|---|
| [`thesis/README.md`](thesis/README.md) | Structure of the thesis and the submission timeline |
| [`REPO_MAP.md`](REPO_MAP.md) | Categorised guide to the whole tree |
| [`KEIO_PLAN.ja.md`](KEIO_PLAN.ja.md) | Keio continuation: paper targets (§0) and work packages |
| [`COUPLING_STATUS.md`](COUPLING_STATUS.md) | One page: what is done in the coupling and what is not |
| [`KLEMPT2024_REPRODUCTION.md`](KLEMPT2024_REPRODUCTION.md) | What reproduces from Klempt 2024 and with which departures |
| [`ansys_usermat/apdl/FRONT_TERM_FIX.md`](ansys_usermat/apdl/FRONT_TERM_FIX.md) | Findings in the partner element (surface points, front term) |
| [`CLAIMS_AND_EVIDENCE.md`](CLAIMS_AND_EVIDENCE.md) | Which numbers can be quoted, and what changed under them |
| [`RESEARCH_MODEL.md`](RESEARCH_MODEL.md) | Model chain, equations, and what is measured, calibrated or assumed |
| [`DOCS.md`](DOCS.md) | Complete documentation index; historical notes in [`archive/`](archive/) |

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
- F. Klempt, H. Geisler, M. Soleimani, P. Junker: A continuum multi-species bacterial growth model with a novel interaction scheme. *Arch. Appl. Mech.* 96 (2026) 164.
- T. Rudolf, F. Klempt, H. I. Kök, M. Soleimani, D. R. Jantos, P. Junker: Computational efficiency and accuracy of the Neighbored Element Method. *Finite Elem. Anal. Des.* 249 (2025) 104353. In [`references/`](references/).
- N. Heine et al.: Influence of species composition and cultivation condition on peri-implant biofilm dysbiosis in vitro. *Front. Oral Health* (2025). The CLSM and viability data.
- M. Soleimani et al. (2023, *Sci. Rep.*) and Z. Feng et al. (2021, *Bull. Math. Biol.*): two-species biofilm models, in [`references/`](references/).

Full bibliography: `thesis/references.bib`. Cite this repository via
[`CITATION.cff`](CITATION.cff); third-party content and its licences are listed
in [`THIRD_PARTY.md`](THIRD_PARTY.md).

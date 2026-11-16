# Repository map

A guide to the tree as of 15 November 2026. The repository holds three things
that overlap: the LUH master's thesis (written, colloquium on 16 Dec 2026), the
coupling paper that grew out of its Chapter 4, and the continuation at Keio.
Start with section 1; the later sections are for looking things up.
[`DOCS.md`](DOCS.md) lists the older prose documents by topic, and
[`archive/`](archive/) holds what is superseded.

`CLAUDE.md` is the working rulebook of the project (machines, decisions with
their dates, people, document rules). Where this map and `CLAUDE.md` differ on
a decision, `CLAUDE.md` is the one that is kept current.

## 1. Start here

| I want to | Read |
|---|---|
| see what is open this month | [`TODO_2026-11.ja.md`](TODO_2026-11.ja.md) (Japanese; dated items, replies awaited, work in progress) |
| read the thesis or check its state | [`thesis/main.pdf`](thesis/main.pdf) (89 pages), [`thesis/REVIEW_2026-11-13.md`](thesis/REVIEW_2026-11-13.md) (read-through of all chapters with the open points) |
| understand the coupled model | [`MEETING_2026-10-20.md`](MEETING_2026-10-20.md) (the coupling agreed with Oliver, with the equation check), [`references/PARAMETER_TABLES.md`](references/PARAMETER_TABLES.md) (every parameter value and its conversion) |
| run the coupled model in Python | [`ansys_usermat/apdl/coupled_living_consumption.py`](ansys_usermat/apdl/coupled_living_consumption.py), results guide [`ansys_usermat/apdl/results/2026-10-3d/README_coupled_living.md`](ansys_usermat/apdl/results/2026-10-3d/README_coupled_living.md) |
| know which ANSYS results and builds to use | [`ansys_usermat/apdl/CONFIRMED_2026-10-22.ja.md`](ansys_usermat/apdl/CONFIRMED_2026-10-22.ja.md), then [`ansys_usermat/apdl/CLOUD_TO_IKMHIWI03.ja.md`](ansys_usermat/apdl/CLOUD_TO_IKMHIWI03.ja.md) |
| read the paper draft | [`papers/coupling_paper/manuscript/`](papers/coupling_paper/manuscript/) and [`papers/coupling_paper/OUTLINE.md`](papers/coupling_paper/OUTLINE.md) |
| prepare the colloquium | [`slides_defense.tex`](slides_defense.tex), [`speaker_notes_defense.tex`](speaker_notes_defense.tex), [`QUIZ_COLLOQUIUM.ja.md`](QUIZ_COLLOQUIUM.ja.md) |
| see the plan for Keio | [`KEIO_PLAN.ja.md`](KEIO_PLAN.ja.md), instructions to the Keio server in [`CLOUD_TO_FIFA.ja.md`](CLOUD_TO_FIFA.ja.md) |

## 2. Thesis and colloquium

| Path | What it is |
|---|---|
| [`thesis/`](thesis/) | The thesis: `thesis/main.tex`, `thesis/chapters/`, `thesis/figures/`, `thesis/references.bib`. Chapters in use: `ch1_introduction`, `ch2_theory`, `ch3_heine` (calibration), `ch4_ansys` (coupling), `ch6_conclusion`, three appendices. Build: `pdflatex`, `bibtex`, `pdflatex` twice. [`thesis/README.md`](thesis/README.md) has the structure and the rules; its status table and timeline are from October |
| [`thesis/REVIEW_2026-11-13.md`](thesis/REVIEW_2026-11-13.md) | Check of 13 Nov 2026: every chapter read against the papers and the result files, findings graded A (must fix) to C (minor) |
| [`thesis/RIGOUR_CHECK_2026-10-27.md`](thesis/RIGOUR_CHECK_2026-10-27.md) | Earlier check of 27 Oct: the printed numbers of Chapter 4 against the result files (`scripts/audit_thesis_numbers.py`), grid convergence recomputed |
| [`thesis/DISCUSSION_COLLOQUIUM.md`](thesis/DISCUSSION_COLLOQUIUM.md) | The reasoning behind the main numbers, with a spoken explanation in English and Japanese for each |
| [`thesis/CODE_TOUR_COLLOQUIUM.ja.md`](thesis/CODE_TOUR_COLLOQUIUM.ja.md) | Which slide corresponds to which place in the code (Japanese) |
| [`thesis/NUMERICS_OVERVIEW.ja.md`](thesis/NUMERICS_OVERVIEW.ja.md) | How each equation of the coupled model is advanced in time and solved (Japanese) |
| [`thesis/HANDOFF_EML_2026-10-19.md`](thesis/HANDOFF_EML_2026-10-19.md), [`thesis/EML_numbers_2026-10-19.md`](thesis/EML_numbers_2026-10-19.md), [`thesis/EML_thesis_figs_2026-10-19.md`](thesis/EML_thesis_figs_2026-10-19.md) | Source of the Chapter 3 numbers and figures: the final calibration runs of the Tmcmc202601 repository (19 Oct) |
| [`thesis/CH4_REWRITE_PLAN.md`](thesis/CH4_REWRITE_PLAN.md), `thesis/chapters/ch4_coupling_draft.tex` | Plan and draft of the Chapter 4 rewrite for the agreed coupling (20 Oct); the text is now in `thesis/chapters/ch4_ansys.tex` |
| [`thesis/MAIL_DRAFT_MEISAM_SUBMISSION.md`](thesis/MAIL_DRAFT_MEISAM_SUBMISSION.md) | Draft of the mail that sends the thesis PDF (planned for December) |
| [`thesis/admin/`](thesis/admin/) | LUH examination rules (MPO 2024, grading); index in `thesis/admin/README.md` |
| [`thesis_ch5/`](thesis_ch5/) | Evidence map and a stand-alone build check (`thesis_ch5/_build_check.tex`) of the coupling chapter from before it moved into `thesis/` |
| [`slides_defense.tex`](slides_defense.tex), [`speaker_notes_defense.tex`](speaker_notes_defense.tex) | Colloquium deck and speaker notes (PDFs next to them) |
| [`QUIZ_COLLOQUIUM.ja.md`](QUIZ_COLLOQUIUM.ja.md) | Question bank for colloquium practice with short answers, and a log of what was asked |
| [`figdata/`](figdata/) | Cached data (`.npz`) of the thesis figures that need a simulation, so they can be redrawn without a rerun |
| [`THESIS_ASSIGNMENT.md`](THESIS_ASSIGNMENT.md) | The scope as set by the supervisor (August 2026) |
| [`data/HEINE_DATA_NOT_IN_GIT.md`](data/HEINE_DATA_NOT_IN_GIT.md) | Why the measured data of Heine et al. are not in git and where scripts expect a local copy |

Other decks at the top level are records of earlier meetings:
`slides_1001_oliver.tex`, `slides_1005.tex` (the template for the background
appendix), `slides_1020.tex`, `slides_coupling_1020.tex`,
`slides_progress_1024.tex`, `slides_1102_nem_ja.tex`, `slides_3d_demo.tex`
with [`speaker_notes_3d_demo.md`](speaker_notes_3d_demo.md), and the Japanese
teaching decks `slides_kiso_ja.tex`, `slides_qa_ja.tex`,
`slides_usermat_ja.tex`, `slides_volume_ja.tex`, `slides_muramatsu_ja.tex`.

## 3. The coupled model (agreed 20 Oct 2026)

The growth field of Klempt et al. (2024, Eqs. 34 to 36) with the point model
of Klempt et al. (2026) at every point: the field passes the nutrient and the
antibiotic to the point model, and the living volume fraction of the point
model sets the consumption in the field.

| Path | What it is |
|---|---|
| [`MEETING_2026-10-20.md`](MEETING_2026-10-20.md) | The meeting, the decisions after it, and the table that checks each equation term by term |
| [`ansys_usermat/apdl/coupled_living_consumption.py`](ansys_usermat/apdl/coupled_living_consumption.py) | The Python implementation (finite volumes on the 2 mm cube): examples with separate and mixed starts, rigid bars, maze, asymmetric starts, the cervical window on a molar crown. Main code of Chapter 4 (Part 3) and of the paper |
| [`ansys_usermat/apdl/coupled_amr.py`](ansys_usermat/apdl/coupled_amr.py) | Prototype of a two-level adaptive grid for the same solver; results in `ansys_usermat/apdl/results/2026-11-amr/README.ja.md` |
| `ansys_usermat/apdl/maze_coupled_2d.py`, `ansys_usermat/apdl/real_tooth_coupled.py` | The maze example in 2D and the coupled run on the real molar surface |
| `ansys_usermat/apdl/fig_*.py` | Figure scripts, one per figure, each with a docstring that names its runs. Set-up and 3D views: `fig_setup_separate_ab.py`, `fig_coupled_3d_panels.py`, `fig_uncut_views.py`, `fig_symmetry_sections.py`. Time histories: `fig_coupled_time_history.py`, `fig_stress_living_over_time.py`, `fig_four_species_shares.py`, `fig_four_species_time_history.py`. Examples: `fig_asymmetric_starts.py`, `fig_obstacle_columns.py`, `fig_maze.py`. Numerics: `fig_limiters_general.py`, `fig_edge_vs_stress.py`. Mechanical term: `fig_mech_K_sweep.py`, `fig_mech_feedback_section.py`. Molar: `fig_tooth_window_setup.py`, `fig_tooth_window_section.py`, `fig_tooth_window_thickness.py` |
| `ansys_usermat/apdl/anim_*.py` | Videos; `anim_coupled_tvd_3d.py` makes the five-panel video of a coupled run |
| [`ansys_usermat/figstyle.py`](ansys_usermat/figstyle.py) | The figure style (Times) that every figure script applies |
| [`ansys_usermat/apdl/results/2026-10-3d/README_coupled_living.md`](ansys_usermat/apdl/results/2026-10-3d/README_coupled_living.md) | Python results: which files are the main runs and which are superseded |
| [`ansys_usermat/apdl/results/2026-10-3d/MECH_K_SWEEP.md`](ansys_usermat/apdl/results/2026-10-3d/MECH_K_SWEEP.md) | The mechanical term of Klempt 2024 Eq. 30 added to Eq. 34: sweep of K and the checks (Keio continuation, not in the thesis) |
| [`ansys_usermat/apdl/results/2026-10-3d/GROWTH_LIVING.md`](ansys_usermat/apdl/results/2026-10-3d/GROWTH_LIVING.md) | Growth from the living biomass in place of Eq. 36 (Keio continuation): the two forms, seed means on three grids, what the antibiotic then does to growth and stress. Deck: `meetings/keio_living_growth/` |
| [`ansys_usermat/apdl/results/2026-10-3d/TOOTH_WINDOW.md`](ansys_usermat/apdl/results/2026-10-3d/TOOTH_WINDOW.md) | The cervical window on the molar crown at physical size: runs, biofilm thickness against measured plaque, stress rescaled to measured stiffness |
| `ansys_usermat/apdl/results/2026-10-3d/` | Further notes on the real-molar runs (`REAL_TOOTH*.md`), the wall treatment (`WALL_CONVERGENCE.md`) and the nutrient demonstration (`NUTRIENT_3D.md`) |
| [`ansys_usermat/apdl/results/2026-10-living/README.md`](ansys_usermat/apdl/results/2026-10-living/README.md) | ANSYS results of the coupling: which run generation is the main one |
| `ansys_usermat/apdl/make_living_deck.py`, `ansys_usermat/apdl/run_living.ps1`, `ansys_usermat/apdl/check_living_trace.py` | ANSYS side: deck generator, run script, and the replay check of an ANSYS trace against the Python point model |
| [`ansys_usermat/apdl/oliver_patches/`](ansys_usermat/apdl/oliver_patches/) | The switches added to Oliver's element routines, as patches |
| [`ansys_usermat/apdl/NOTE_OLIVER_2026-10-24.md`](ansys_usermat/apdl/NOTE_OLIVER_2026-10-24.md) | What I found in the element while checking it against the two other implementations (surface points, Laplacian, neighbours); draft note |
| [`abaqus_living/`](abaqus_living/) | The same coupling as an Abaqus user element (`abaqus_living/living_uel.f90`), deck generators for the cube and the maze, manufactured solution, run summaries. Keio continuation |
| `ansys_usermat/coupling_schematic.py` | The schematic of the coupling (`--agreed`) |
| [`references/PARAMETER_TABLES.md`](references/PARAMETER_TABLES.md) | Parameter tables of all the papers and the conversion to the 2 mm cube |
| [`ecology_constants.py`](ecology_constants.py) | The one place where c* and the Hill gate of the calibrated point model are set |

## 4. Papers

| Path | What it is |
|---|---|
| [`papers/coupling_paper/manuscript/`](papers/coupling_paper/manuscript/) | The coupling paper: `coupling_paper.tex`, `refs.bib`, and a README with the history of the restructuring |
| [`papers/coupling_paper/OUTLINE.md`](papers/coupling_paper/OUTLINE.md), [`papers/coupling_paper/FIGURES.md`](papers/coupling_paper/FIGURES.md) | Outline with the proposed authors, and the figure plan |
| `papers/coupling_paper/sec1_intro_draft.md` to `sec6_conclusions_draft.md` | Section drafts the manuscript was assembled from |
| [`papers/coupling_paper/REAL_TOOTH_WINDOW.ja.md`](papers/coupling_paper/REAL_TOOTH_WINDOW.ja.md) | Set-up of the real-geometry example: the cervical window (Japanese) |
| [`papers/coupling_paper/ASYMMETRIC.ja.md`](papers/coupling_paper/ASYMMETRIC.ja.md), [`papers/coupling_paper/EDGE_STRESS.ja.md`](papers/coupling_paper/EDGE_STRESS.ja.md), [`papers/coupling_paper/SINE_SCHEDULES.ja.md`](papers/coupling_paper/SINE_SCHEDULES.ja.md) | Notes on single examples: asymmetric starts, edge sharpness against stress, time-varying supply |
| [`papers/coupling_paper/VALIDATION_PLAN.ja.md`](papers/coupling_paper/VALIDATION_PLAN.ja.md), [`papers/coupling_paper/VALIDATION_RESULTS.ja.md`](papers/coupling_paper/VALIDATION_RESULTS.ja.md), [`papers/coupling_paper/EXPERIMENT_COMPARISON.ja.md`](papers/coupling_paper/EXPERIMENT_COMPARISON.ja.md), `papers/coupling_paper/validation_joint/` | Comparison with measurements: plan, results, the joint fit |
| [`papers/coupling_paper/USNCCM19_ABSTRACT.md`](papers/coupling_paper/USNCCM19_ABSTRACT.md), [`papers/coupling_paper/ABSTRACTS_2027.md`](papers/coupling_paper/ABSTRACTS_2027.md) | Conference abstract drafts for 2027 |
| [`papers/coupling_paper/RELEASE_PLAN.ja.md`](papers/coupling_paper/RELEASE_PLAN.ja.md), [`papers/coupling_paper/REVIEW_2026-10-29.md`](papers/coupling_paper/REVIEW_2026-10-29.md) | In which order the results are shown to the co-authors, and a review of the draft |
| [`papers/species_continuum/IDEA.ja.md`](papers/species_continuum/IDEA.ja.md) | Idea for a later paper: one volume-fraction field per species (not in the thesis) |
| [`papers/surrogate_pilot/PLAN.ja.md`](papers/surrogate_pilot/PLAN.ja.md) | Pilot of a surrogate model that carries parameter uncertainty to the stress (not in the thesis) |
| [`references/`](references/) | The papers this work follows, each with a searchable text extraction: Klempt et al. 2024 and 2026, Soleimani et al. 2023, Feng et al. 2021, the NEM papers (reading guide `references/NEM_PAPERS.md`), and two experimental papers. Licences in [`THIRD_PARTY.md`](THIRD_PARTY.md) |
| [`references/EXPERIMENTAL_DATA.md`](references/EXPERIMENTAL_DATA.md) | Measurements collected for comparison (stiffness, thickness, antimicrobial tests) |
| [`references/CHLORHEXIDINE_VS_POINT_MODEL.md`](references/CHLORHEXIDINE_VS_POINT_MODEL.md) | Gränicher et al. 2021 (chlorhexidine on a six-species oral biofilm) against the antibiotic term of the point model |

## 5. Meetings and decks

| Path | What it is |
|---|---|
| [`meetings/2026-11-02_muramatsu/`](meetings/2026-11-02_muramatsu/) | Meeting with Assoc. Prof. Muramatsu on 2 Nov: notes (`NOTES.ja.md`) and deck |
| [`meetings/2026-11_coupled_examples/`](meetings/2026-11_coupled_examples/) | Decks on the coupled examples (English and Japanese), built with `build.sh`; what was sent to whom is in its `README.ja.md` |
| [`meetings/2026-11-30_muramatsu/`](meetings/2026-11-30_muramatsu/) | Next meeting (30 Nov): draft abstract of the Keio report |
| [`meetings/2026-12-01_junker/`](meetings/2026-12-01_junker/) | Agenda for the meeting with Prof. Junker and Dr. Soleimani on 1 Dec |
| [`meetings/keio_deformation/`](meetings/keio_deformation/) | Japanese deck on the two-way coupling with the deformation (the mechanical term of Klempt 2024 Eq. 30); extended as results arrive |
| [`meetings/2027-03_kadaikenkyu/`](meetings/2027-03_kadaikenkyu/) | Slides for the five-minute recorded talk of the Keio report (due 12 Mar 2027) |
| [`keio_report/`](keio_report/) | The Keio report (課題研究報告, Japanese, two-column): `keio_report/kadaikenkyu_nishioka.tex`; the announcement with the rules is `keio_report/KLMS_ANNOUNCEMENT.ja.md` |
| [`MEETING_2026-10-20.md`](MEETING_2026-10-20.md), `archive/MEETING_2026-10-05.md` | Records of the meetings of 20 Oct (Oliver) and 5 Oct (Dr. Soleimani, Assoc. Prof. Muramatsu) |
| [`reports/`](reports/) | Written progress updates to the supervisors |
| [`MAIL_DRAFT_OLIVER_2026-11-04.md`](MAIL_DRAFT_OLIVER_2026-11-04.md) | Draft mail to Oliver with the results of the agreed coupling (not sent; the decks of 6 Nov and the meeting of 10 Nov took its place) |

## 6. Machines: instructions and working records

Three places do the work: IKMHIWI03 (LUH, Windows, ANSYS), fifa (Keio, Linux,
Abaqus and Python) and a cloud session that writes and merges. Each pair has
an instruction file and a report file.

| Path | What it is |
|---|---|
| [`CLOUD_TO_FIFA.ja.md`](CLOUD_TO_FIFA.ja.md) | Instructions from the cloud session to fifa, newest first |
| [`KEIO_CHAIN_LOG.md`](KEIO_CHAIN_LOG.md) | What fifa was asked to run next and why, one entry per decision, newest at the bottom |
| [`KEIO_SERVER_HANDOFF.ja.md`](KEIO_SERVER_HANDOFF.ja.md), [`KEIO_SERVER_ENV.ja.md`](KEIO_SERVER_ENV.ja.md) | Setting up a clone on fifa and running Abaqus there (7 Oct), and the environment as measured (8 Oct) |
| `scripts/setup_keio_server.sh`, `scripts/run_chain_keio.py`, `scripts/claude_remote.sh` | One-time setup on fifa, the run chain with its guards for the shared server, and the remote session |
| [`ansys_usermat/apdl/CONFIRMED_2026-10-22.ja.md`](ansys_usermat/apdl/CONFIRMED_2026-10-22.ja.md) | ANSYS: confirmed facts, the builds, and what not to use. Read first |
| [`ansys_usermat/apdl/CLOUD_TO_IKMHIWI03.ja.md`](ansys_usermat/apdl/CLOUD_TO_IKMHIWI03.ja.md) | Current instructions to IKMHIWI03 (older ones in `ansys_usermat/apdl/archive/`) |
| [`ansys_usermat/apdl/IKMHIWI03_TO_CLOUD.ja.md`](ansys_usermat/apdl/IKMHIWI03_TO_CLOUD.ja.md) | Reports from IKMHIWI03, numbered A-1 onwards |
| [`ANSYS_ENVIRONMENT.md`](ANSYS_ENVIRONMENT.md) | Hardware, licence and product inventory of IKMHIWI03 |
| [`ansys_usermat/apdl/RUNBOOK.md`](ansys_usermat/apdl/RUNBOOK.md), `ansys_usermat/apdl/link_v222.ps1`, `ansys_usermat/apdl/run_apdl.ps1`, `ansys_usermat/apdl/run_chain.ps1` | Building the custom ANSYS executable and running decks, single or as a detached chain |
| `dev-env.ps1`, `run_abaqus.ps1`, `run_tests.ps1`, `run_notebooks.ps1`, `build_slides.ps1`, `commit.ps1`, `push.ps1` | Helper scripts for IKMHIWI03 (described in `CLAUDE.md`) |
| [`BACKUP_PLAN_2026-12.ja.md`](BACKUP_PLAN_2026-12.ja.md), `scripts/backup_untracked_ikmhiwi03.ps1` | What to back up from IKMHIWI03 before the end of December |
| [`DAILY_REPORT_TRIP.ja.md`](DAILY_REPORT_TRIP.ja.md) | Daily reports of the cloud session during the trip of 21 to 29 Oct |
| `.claude/hooks/session-start.sh`, `scripts/pre-commit-no-ai-identity.sh`, `scripts/commit-msg-no-ai-trailer.sh` | Restores git identity, hooks and test dependencies in a fresh cloud container; the two git hooks |

## 7. Keio continuation

Not part of the LUH thesis.

| Path | What it is |
|---|---|
| [`KEIO_PLAN.ja.md`](KEIO_PLAN.ja.md) | The plan: goals and deadlines, paper targets, work packages, then a dated record |
| [`abaqus_composition/`](abaqus_composition/) | The composition coupling of the ANSYS runs in an Abaqus UMAT built from the same fragments, with the point model in Fortran |
| [`abaqus_living/`](abaqus_living/) | The agreed coupling as an Abaqus user element (section 3) |
| [`keio_wp2/`](keio_wp2/) | Stress acting back on growth (homeostatic pressure): Python prototype and Abaqus runs on the implant collar and the tooth crown |
| [`keio_p1/`](keio_p1/) | NEM against Galerkin finite elements for the growth field, mesh and time-step study in Python |
| [`p3/`](p3/), [`UNIFIED_MODEL_CONCEPT.ja.md`](UNIFIED_MODEL_CONCEPT.ja.md), [`P3_OUTLINE.ja.md`](P3_OUTLINE.ja.md), [`P3_STEP3_FEM_PLAN.ja.md`](P3_STEP3_FEM_PLAN.ja.md), [`P3_LITERATURE.ja.md`](P3_LITERATURE.ja.md) | The two-level model from one variational principle: working note (`p3/unified_model.tex`), concept, paper outline, finite-element plan, literature |
| [`uq/`](uq/) | Uncertainty of growth and stress over the final posterior of thesis Chapter 3 |
| [`ROADMAP_TWO_WAY.md`](ROADMAP_TWO_WAY.md) | From a one-way to a two-way coupling, in four steps (2 Oct, with an update note of 20 Oct) |
| [`RESEARCH_IDEAS.ja.md`](RESEARCH_IDEAS.ja.md) | Research ideas after the thesis, mapped to the planned papers |
| [`PINN_DESIGN.md`](PINN_DESIGN.md) | Design of a physics-informed surrogate, not started |
| [`EML_REQUEST_VOLUME.ja.md`](EML_REQUEST_VOLUME.ja.md) | Draft request to add the measured biofilm volume to the calibration likelihood (postponed on 15 Oct) |

## 8. The ANSYS element work (thesis Chapter 4, Parts 1 and 2)

| Path | What it is |
|---|---|
| [`ansys_usermat/`](ansys_usermat/) | ANSYS USERMAT of the Klempt growth law (`ansys_usermat/usermat_biofilm.f`), with `ansys_usermat/crosscheck/` (comparison with the Abaqus UMAT) and [`ansys_usermat/coupling/`](ansys_usermat/coupling/) (bridge to Python and the point model in Fortran, `ecology_native.f`) |
| [`ansys_usermat/apdl/README.md`](ansys_usermat/apdl/README.md) | The verification decks of the growth law, results table, how to read a failure |
| [`ansys_usermat/apdl/V222_PORT_INSTRUCTIONS.md`](ansys_usermat/apdl/V222_PORT_INSTRUCTIONS.md), [`ansys_usermat/OLIVER_MODEL_NOTES.md`](ansys_usermat/OLIVER_MODEL_NOTES.md), [`ansys_usermat/INTEGRATION_PLAN.md`](ansys_usermat/INTEGRATION_PLAN.md) | Running Oliver's routines on ANSYS 2022 R2, notes on his model, and how the material law was contributed to it (September) |
| [`ansys_usermat/apdl/ONE_SPECIES_COUPLING.md`](ansys_usermat/apdl/ONE_SPECIES_COUPLING.md), [`ansys_usermat/apdl/SPLIT_COUPLING.md`](ansys_usermat/apdl/SPLIT_COUPLING.md), [`ansys_usermat/apdl/NSUB_WIRING.md`](ansys_usermat/apdl/NSUB_WIRING.md), [`ansys_usermat/apdl/N3_GROWTH_TRIAL.md`](ansys_usermat/apdl/N3_GROWTH_TRIAL.md) | Design notes of the one-way coupling (October): Eq. 36 at the Gauss point, the point model for the composition, sub-stepping, more than two species. The agreed coupling of section 3 replaced the coupling they describe |
| [`ansys_usermat/apdl/FRONT_TERM_FIX.md`](ansys_usermat/apdl/FRONT_TERM_FIX.md) | Findings in the element: surface points and the term $-r\,\frac{c}{k+c}\,\mathbf n_c\cdot\nabla\phi$ of Eq. 34 |
| [`ansys_usermat/apdl/MESH_STUDY.md`](ansys_usermat/apdl/MESH_STUDY.md), [`ansys_usermat/apdl/BETA002_THESIS_MAP.md`](ansys_usermat/apdl/BETA002_THESIS_MAP.md) | The light mesh study, and the table that mapped the β = 0.02 mm²/T* runs to thesis tables and figures (6 Oct) |
| [`ansys_usermat/apdl/RUN_WP3_IKMHIWI03.md`](ansys_usermat/apdl/RUN_WP3_IKMHIWI03.md), [`ansys_usermat/apdl/RUN_WP2_IKMHIWI03.md`](ansys_usermat/apdl/RUN_WP2_IKMHIWI03.md) | Run sheets with results: seed growth with NEM (WP3) and the homeostatic growth law (WP2), 7 Oct |
| `ansys_usermat/apdl/nem_stability.py`, `ansys_usermat/apdl/nem_finite_cube.py`, `ansys_usermat/apdl/wp3_reference.py`, `ansys_usermat/apdl/wp3_stress_reference.py`, `ansys_usermat/apdl/wp3_stress_gauss.py`, `ansys_usermat/apdl/wp3_alpha_profile.py`, `ansys_usermat/apdl/wp3_figs.py` | Seed-growth numerics of the element: stability limit of the explicit NEM update (thesis appendix), finite-volume reference, seed stress, the figures of Chapter 4 |
| `ansys_usermat/apdl/closed_form_reference.py`, `ansys_usermat/apdl/check_deck.py`, `ansys_usermat/apdl/make_layered_material.py`, `ansys_usermat/apdl/make_mesh_levels.py`, `ansys_usermat/apdl/plot_cylinder_3d.py` | Closed form of the two growth cases, the pre-flight check of a deck, layered materials from a growth field, mesh levels, the curved-shell plot |
| [`ansys_usermat/apdl/threshold_param_sweep_README.md`](ansys_usermat/apdl/threshold_param_sweep_README.md) | The convergence threshold of the curved-shell run against its parameters |
| [`ansys_usermat/HANDOFF.md`](ansys_usermat/HANDOFF.md), [`handover/`](handover/) | The hand-over of the material law to the partner group (August and September), generated by `handover/make_handover.py` |
| `ansys_usermat/biofilm_material_v01.f`, `ansys_usermat/growth_law_verification.ipynb` | The routine handed over and the walkthrough of what its verification establishes |
| [`KLEMPT2024_REPRODUCTION.md`](KLEMPT2024_REPRODUCTION.md), [`PAPER_CHECK_KLEMPT2024.md`](PAPER_CHECK_KLEMPT2024.md) | What reproduces from Klempt 2024 and with which departures; the chapter and decks checked against the paper (2 Oct) |

## 9. Notes on the source papers

[`CLAIMS_AND_EVIDENCE.md`](CLAIMS_AND_EVIDENCE.md) is the consolidation (which
numbers can be quoted); the others are the detail behind it.

| File | What it settles |
|---|---|
| [`SOLEIMANI2023_NOTES.md`](SOLEIMANI2023_NOTES.md) | The nearest precedent (two species in an ANSYS user element) and the group's validation bar |
| [`SOLEIMANI2021_NOTES.md`](SOLEIMANI2021_NOTES.md) | The cap on α, and the group's warning about the advection term |
| [`SOLEIMANI2019_NOTES.md`](SOLEIMANI2019_NOTES.md) | The unconditionally stable viscous integrator, and E = 10 Pa |
| [`FEI2020_NOTES.md`](FEI2020_NOTES.md) | Fei et al. 2020: growth-induced stress and wrinkling compared with experiment, a validation target |
| [`CHU2018_NOTES.md`](CHU2018_NOTES.md) | The stress level at which the bacterial stress response turns on |
| [`READING_GAPS.md`](READING_GAPS.md), [`CITATION_AUDIT.md`](CITATION_AUDIT.md) | References Klempt 2024 rests on that were not yet cited, and the audit of the citations (August) |

## 10. Findings in the earlier code

| File | What it is |
|---|---|
| [`VERIFICATION_SENSITIVITY_LIMITATIONS.md`](VERIFICATION_SENSITIVITY_LIMITATIONS.md) | What is verified and what is assumed in the tooth pipeline of the first half of 2026 |
| [`DEVIATOR_SCALING_FINDING.md`](DEVIATOR_SCALING_FINDING.md) | A mis-scaled isochoric split in the material core: a pressure error, von Mises unaffected |
| [`VISCOUS_UPDATE_SCHEME.md`](VISCOUS_UPDATE_SCHEME.md) | What the viscous update is and the step limit that follows |
| [`PDE_VERIFICATION_FINDINGS.md`](PDE_VERIFICATION_FINDINGS.md) | The nutrient and species solvers against exact solutions; the zero-flux wall fixed |
| [`E_SATURATION_FINDING.md`](E_SATURATION_FINDING.md) | The stiffness bridge clips to [10, 1000] Pa, so distinct conditions can report the same stiffness |
| [`DS_composition_fix.md`](DS_composition_fix.md), [`espec_sensitivity_findings.md`](espec_sensitivity_findings.md) | The composition copy error of one condition and the sensitivity to the per-species stiffness (June) |
| [`RESEARCH_MODEL.md`](RESEARCH_MODEL.md) | The model chain of that pipeline: equations, and what is measured, calibrated or assumed |

## 11. Earlier lineages (before October 2026, not thesis results)

| Path | What it is |
|---|---|
| [`pipeline.py`](pipeline.py), [`PIPELINE.md`](PIPELINE.md), [`configs/`](configs/) | Config-driven pipeline of the tooth analysis and its risk metric |
| `gen_tooth_klempt_umat_inp.py`, `umat_biofilm_visco.f`, `umat_biofilm_visco_2ch.f`, `umat_biofilm_visco_phase2.f`, `usdfld_biofilm.f` | Abaqus tooth job generation and the viscoelastic UMATs |
| `material_models.py`, [`FEM_README.md`](FEM_README.md), [`methods_supplement_fem.md`](methods_supplement_fem.md) | The line that bridges a dysbiosis index to the stiffness |
| [`JAXFEM/`](JAXFEM/) | JAX solvers: the Klempt 2024 equations in 1D and 2D, reproductions of Klempt 2024 and 2026 (`JAXFEM/klempt2026_reproduction.py`), the synthetic cases of Fritsch et al. 2025 (`JAXFEM/fritsch2025_cases.py`), audits (`JAXFEM/audit_all.py`) |
| `jax_hamilton_0d_5species_demo.py`, `multiscale_coupling_1d.py`, `multiscale_coupling_2d.py`, [`validate_composition.py`](validate_composition.py) | Point-model and multiscale demos, and the comparison of the model composition with the Heine experiment |
| [`tier2b_real/`](tier2b_real/) | Abaqus coupon and implant jobs and the real-geometry meshes (the molar crown STL used in section 3 is here) |
| [`umat_tangent_test/`](umat_tangent_test/) | Single-element tangent and eigenstrain cross-checks of the UMATs |
| [`umat_flow/`](umat_flow/), [`ch5_flow/`](ch5_flow/) | TikZ method figures of the earlier chapter |
| [`klempt2024_results/`](klempt2024_results/), [`runs/`](runs/) | Outputs and run logs of that work |
| [`THESIS_PLAYBOOK.md`](THESIS_PLAYBOOK.md), [`research_goals_1_2.md`](research_goals_1_2.md) | Process view and research goals as of August 2026. Both describe the earlier pipeline; for the present state use section 1 |
| [`archive/legacy_2026H1/`](archive/legacy_2026H1/) | Root scripts, reports and results of the first half of 2026 that nothing uses any more |

## 12. Data, tests, site

| Path | What it is |
|---|---|
| [`data/`](data/) | Experimental data that may be shared, and the final posterior of Chapter 3 |
| [`tests/`](tests/) | Unit tests (`pytest tests/`); `tests/test_docs_index.py` fails when a top-level document or directory is missing from this map |
| `pytest.ini`, `requirements.txt`, `.github/workflows/ci.yml` | Test scoping, pinned versions, CI (run by hand) |
| [`assets/`](assets/) | Figures and videos written by the figure scripts |
| [`docs/`](docs/) | The project site (GitHub Pages) |
| [`travel/`](travel/) | Personal travel notes, not part of the research |

## 13. Archive

[`archive/`](archive/) holds superseded plans, status notes and hand-offs
(index in `archive/README.md`), and `ansys_usermat/apdl/archive/` the
instruction files to IKMHIWI03 of 13 to 22 Oct and older run sheets. They are
kept as records; their settings are mostly outdated.

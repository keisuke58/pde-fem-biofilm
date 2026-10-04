# Thesis (working copy)

The thesis is now written here. Copied on 3 Oct 2026 from
`keisuke58/LUH_summer_2026`, branch `claude/claude-md-docs-dc6ldk`, commit
`c3b8d67`, folder `30_Masterarbeit/` (`main.tex`, `references.bib`,
`chapters/`, and the figures that chapters 3, 4 and the appendix use). The
first commit of this folder is that copy unchanged, so every later change can
be read as a diff against it. That repository is no longer edited for the
thesis.

## Structure (agreed 3 Oct, to be confirmed with the supervisor on 5 Oct)

| ch. | content | from |
|---|---|---|
| 1 | introduction: clinical background, the question (point model into an FEM element), three contributions | old ch1 background, new objectives |
| 2 | background: Hamilton principle, Klempt 2024 growth field, Klempt 2026 point model, F = F_e F_g, TMCMC | old ch2 + `../thesis_ch5` + deck appendices |
| 3 | calibration of the point model: Heine five-species TMCMC; last section (3-5 pages) the metabolic sign prior (internship work) | old ch3 + old ch4, condensed |
| 4 | coupling into the ANSYS element: one and two species, verification, results, limits | `../thesis_ch5/ch5_ansys_contribution.tex` |
| 5 | conclusion and outlook (Keio: more species, two-way coupling) | old ch6, rewritten |

The old ch5 (spatial PDE, FISH, Abaqus-era stress results) is out of scope;
at most one paragraph moves to ch3 or the outlook. Abaqus results are not
thesis results (CLAUDE.md).

## Rules

The document rules of `../CLAUDE.md` apply (first person singular, notation
as in the papers, no AI-sounding wording, units, figures in Times New Roman).

## Build

```bash
cd thesis && pdflatex main && bibtex main && pdflatex main && pdflatex main
```
Needs `algorithm.sty` (Debian/Ubuntu: `texlive-science`).

## Status (3 Oct 2026)

| part | state |
|---|---|
| abstract | rewritten for the new structure (German draft kept, disabled) |
| ch1 introduction | rewritten |
| ch2 background | rewritten; checked against Klempt 2024 and the Klempt 2026 reproduction |
| ch3 calibration | old text kept; in-vivo sign-prior section added at the end (condensed old ch4); style pass still to do (first person, dashes) |
| ch4 ANSYS coupling | moved in from `../thesis_ch5`; drafting markers (A/B in the margin) still in |
| ch5 conclusion | rewritten (file `chapters/ch6_conclusion.tex`) |
| appendix | AGORA2 strains, TMCMC details, Klempt 2024 reproduction (new) |

Title page (3 Oct 2026) set exactly as registered on the form of 18.08.2026:
title "Hamilton-Principle Models of Oral Biofilm Dysbiosis: GPU-Accelerated
Inference, Metabolic Priors, and Spatial Stratification"; Erstpruefer Prof.
Junker (IKM), Zweitpruefer Dr.-Ing. M. Wangenheim (IDS), Betreuer Dr.-Ing. M.
Soleimani and Dr.-Ing. H. Geisler. The form itself (personal data) is not in
the repository.

Open:
- the registered subtitle names "Spatial Stratification" (the FISH work, now out
  of scope); the spatial part of the thesis is now the FEM coupling. Ask whether
  the subtitle can be changed, or describe the spatial composition of chapter 4
  in those terms;
- the old ch4 and ch5 files (`chapters/ch4_dieckow.tex`, `chapters/ch5_integration.tex`, `chapters/appendix_unused.tex`) are no longer input and can be deleted once nothing more is taken from them.
- ch3 numbers: the TMCMC results may still change (4 Oct 2026). Once they are
  final, check every number in the ch3 text against its table or figure (on
  4 Oct one pair of gLV entries had mixed CS and CH values, and the posterior
  was mislabelled NUTS; both fixed).

## Timeline to submission and colloquium (decided 4 Oct)

| when | what |
|---|---|
| 5 Oct | meeting 9:00; IKMHIWI03 run sheet `ansys_usermat/apdl/RUN_1005_IKMHIWI03.md` (last ANSYS day before 12 Oct) |
| 6-11 Oct | away; cloud: ch4 updated with the 5 Oct results, read-through |
| 12-16 Oct | TMCMC final, ch3 numbers; **full draft to Prof. Soleimani and Dr. Geisler about 16 Oct** |
| 19-28 Oct | supervisor comments |
| 29-30 Oct | proofreading, printing, binding |
| **2 Nov** | **submission** (date to confirm) |
| 3-10 Nov, 20-28 Nov | away: read the submitted thesis chapter by chapter, Q&A practice |
| 30 Nov-2 Dec | colloquium, before the December trip (3-11 Dec) |

Reading plan while away in November (about one hour a day):
1. ch1 and ch2 (model and notation) with the variables and equations sheet;
2. ch4 verification and coupling, redoing the closed-form stress by hand;
3. ch4 results, assumptions and sensitivity (s, phi_cap, d);
4. ch3 calibration and ch6;
5. `QA_1005.md` and new questions, answered aloud.

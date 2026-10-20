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

- **Figures (19 Oct 2026): all in Times**, drawn at the printed width (16 cm) with text of 7.5 pt or more, vector PDF
  (3D renderings PNG at 400 dpi), no title inside the figure, panel letters as in the caption. Body text Latin Modern.
  Figure scripts call `figstyle.apply()`; data of figures that need a simulation are cached in `figdata/`.

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

## Timeline to submission and colloquium (decided 4 Oct, revised 6 Oct and 19 Oct)

| when | what |
|---|---|
| 5 Oct | meeting 9:00; IKMHIWI03 run sheet `ansys_usermat/apdl/RUN_1005_IKMHIWI03.md` (last ANSYS day before 12 Oct) |
| 6-11 Oct | away; cloud: ch4 updated with the 5 Oct results, read-through |
| **12 Oct** | **content frozen** (decided 6 Oct): ch4 with the week runs that have arrived by then, ch3 TMCMC as it is; afterwards only format and typos. **Nearly complete draft to Dr. Soleimani** (format check), PDF renamed `Nishioka_MasterThesis_draft_2026-10-12.pdf` |
| 16 Oct | his comments on the format |
| 19 Oct | ch3 replaced by the final EML results; figures redrawn at text width (ch1, ch4); ch3 to Assoc. Prof. Muramatsu (Slack, scheduled) |
| 21-30 Oct | away (leaves Hannover 21 Oct about 15:00) |
| 31 Oct - 2 Nov | in Hannover: **test print** of the whole thesis, ask the copy shop for the binding time; meeting with Assoc. Prof. Muramatsu on 2 Nov (`slides_1102_nem_ja.tex`) |
| 3-10 Nov | away |
| **11-12 Nov** | **final print and binding** (two copies; calendar: City Copy Shop on 12 Nov, unless Dr. Soleimani says the IKM printer is fine), PDF for the electronic version (the printing date of 19 Oct belonged to the old 20 Oct submission) |
| **13 Nov** | **submission** (moved from 20 Oct on 13 Oct): two printed copies to Dr. Soleimani, plus the electronic version (PO § 15 (4)) |
| 20-28 Nov | away: read the submitted thesis chapter by chapter, Q&A practice |
| 16 Dec, 13:00 | colloquium, IKM room 320, in person; talk 20 min, questions 10-20 min (confirmed 13 Oct). Dr. Wangenheim (second examiner) wrote on 19 Oct that his attendance is not required |

Dr. Soleimani's answer (20 Oct): printing, ask Timo (Timo Rabing, IKM, rabing@ikm.uni-hannover.de); students usually print
at the university, which is cheap. Schedule: the PDF may be submitted up to 24 hours before the colloquium (so by 15 Dec,
13:00), the printed copies afterwards. Mail to Timo sent 20 Oct (where, how, binding, lead time; two copies, ~90 pages, colour,
printing in December). Plan: PDF ready around 13 Nov anyway; print in Hannover on 14-15 Dec and hand the copies over at the
colloquium (fallback 17-18 Dec, before the Christmas closures). Reply sent 20 Oct: asked whether the PDF counts as the
formal submission or the examination office / Prof. Junker must be told, and again for comments on the draft.

Open with Dr. Soleimani (19 Oct):
- his comments on the 9 Oct draft (asked for by 15 Oct) have not arrived;
- mail of 19 Oct asks how to print (copy shop or the IKM printer; sides, colour, binding), the form of the electronic
  version, and whether handing in on Fri 13 Nov is fine;
- follow-up: when he answers the printing mail, ask in the same reply whether he had time to look at the draft
  ("Any comments on the format before the final print would be very welcome"); if nothing has come by about 26 Oct,
  send one reminder that covers both. Do not send him a newer draft unless he asks (decided 9 Oct).

Colloquium practice during the trips (decided 19 Oct):
- Adriatic, 21-27 Oct: a little quiz between buses (before submission, so light);
- Portugal, 3-10 Nov: quiz (`QUIZ_COLLOQUIUM.ja.md`) and the reading plan below;
- Greece, 21-28 Nov: full run-through of the talk, Q&A in English;
- Egypt, 3-9 Dec: only the questions that went wrong;
- Hannover, 10-15 Dec: two or three full rehearsals aloud.

Reading plan while away in November (about one hour a day):
1. ch1 and ch2 (model and notation) with the variables and equations sheet;
2. ch4 verification and coupling, redoing the closed-form stress by hand;
3. ch4 results, assumptions and sensitivity (s, phi_cap, d);
4. ch3 calibration and ch6;
5. `QA_1005.md` and new questions, answered aloud.

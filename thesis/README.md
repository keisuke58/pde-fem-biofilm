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

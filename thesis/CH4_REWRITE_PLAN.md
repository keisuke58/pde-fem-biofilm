# Chapter 4 rewrite: the coupling agreed with Oliver (plan, 20 Oct 2026)

Decision: `MEETING_2026-10-20.md` (into the thesis, replacing coupling (a)). Internal deadline 13 Nov.
Draft of the new method text: `thesis/chapters/ch4_coupling_draft.tex` (not yet included in the build).

## What changes, section by section (`thesis/chapters/ch4_ansys.tex`)

| section | now | new |
|---|---|---|
| 4.1 Scope | point model gives the composition, one-way | point model used as published (unchanged); nutrient and antibiotic fields feed it, its living biomass sum phibar_i feeds the consumption back; growth unchanged |
| 4.2 Model | Eq. 34 without front term, Eq. 35 g phi c | Eq. 34 with the front term (as printed), Eq. 35 with consumption g phi L c; one sentence that g phi in the paper is a misprint; the second nutrient field carries the antibiotic |
| 4.3 Implementation | bridge, server, native | unchanged in structure; the call now passes (c_n, c_{n+1}, a_n, a_{n+1}) and returns L |
| 4.4 Verification | closed forms | unchanged |
| 4.5 Coupling to the field | framework | add: the front term (ORI weight, its slip of 2 Oct), upwinding / overshoot check |
| 4.6 Two species | scaled scheme main, unscaled scheme second, phi_cap, s = 0.15 | **unscaled scheme (mode 8) is the scheme**; the scaled scheme and phi_cap move to a short remark or are dropped; s from the data window; inputs interpolated in time; consumption by the living share; start values per region and for new points; algorithm rewritten |
| 4.7 Results | two species, packed seed; five species calibrated | A/B separate and mixed, with and without antibiotic (synthetic, Fritsch Case I); ANSYS if A-19 finishes, otherwise Python with the replay check; the calibrated five-species run kept as the one-way case (no antibiotic data, b = 0) |
| 4.8 Limitations | | add: new points' start values, phi_min, one-step lag of L, no transport of species between Gauss points (A and B do not mix), front term overshoot |
| 4.9 Outlook | steps 1-4 in Python | step 1 (local nutrient) is now done; keep 2-4 shortened |

Also: abstract, Ch. 1 contributions, Ch. 6 conclusion (one sentence each); the Chapter 3 link (calibrated model has b = 0).

## Numbers that the text will need (to fill in when the runs are in)

- replay check: coupled point model against the point model alone (Python: running; ANSYS: A-19 check 1)
- seed and edge columns: living share, share of A, c, a at T* = 0.5, 1 for the four runs (front r and r/10)
- phi range with the front term
- sensitivity: phi_min 0.01 / 0.1, outer step 0.025 / 0.0125, mesh 32^3 / 48^3 (Python)

## Order of work

1. Draft of 4.6 and the model changes in 4.2 (now, no numbers needed): done as `ch4_coupling_draft.tex`.
2. Numbers from the Python runs (tonight).
3. ANSYS (A-19, IKMHIWI03 on 21 Oct): replay check, then the four runs.
4. Results section and figures; limitations; abstract and conclusion.
5. Build, check the rendered pages, page count.

Wording: do not call the point model a "black box" in the thesis (the user, 20 Oct); write that it is used as
published, unchanged, with given inputs.

## Also to change with Chapter 4 (found 20 Oct)

- Ch. 2 `ch2_theory.tex` l. 116, 322-332: "I combine the two models in the same one-way manner" -> the new coupling.
- Ch. 6 `ch6_conclusion.tex` l. 73 (table: "yes (one-way)") and l. 137 (limitation "One-way coupling").
- Defence deck: Part 3 slides, Backup 5 (sensitivity to s and phi_cap), the "Inactive front term" limitation
  (the front term is now on), "Assumptions" item; speaker notes for these.
- `ansys_usermat/COUPLING_BASICS.ja.md`, `REPO_MAP.md`, `abaqus_composition/README.md`, `thesis_ch5/README.md`:
  one line each pointing to the new coupling.
- Keio report (`keio_report/`) after the thesis is fixed.

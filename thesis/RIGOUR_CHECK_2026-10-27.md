# Rigour check of the thesis, 27 October 2026

Scope: `thesis/main.tex` with the chapters it inputs (Ch. 1, 2, 3, 4 = `ch4_ansys.tex`, 6, appendices).
No new results; errors in what is stated or computed were looked for and, where the source is unambiguous,
fixed. Build: pdflatex, bibtex, pdflatex, pdflatex in `thesis/`.

## 1. Numbers (`scripts/audit_thesis_numbers.py`)

The script now checks 202 printed values of Chapter 4, the abstract and the conclusion against the result
files (0 FAIL, 1 WARN). Added in this check (sections 12 to 18 of the script):

- Grid table `tab:ch5-gci`: all 30 grid values; the three series without antibiotic or with the Fritsch
  start fractions come from the cloud result files (`ansys_usermat/apdl/results/2026-10-3d/*front1_dtc0.0025*inhliv*`,
  outer step 0.0025, 32^3/64^3/128^3), the two series with antibiotic and equal start amounts from report
  A-89 of `ansys_usermat/apdl/IKMHIWI03_TO_CLOUD.ja.md` (parsed from its table). The observed order p, the
  extrapolated value and the GCI were recomputed with the three-grid procedure of Celik et al. 2008
  (r = 2, Fs = 1.25): all 27 printed values agree; the "--" row (mixed (0.25, 0.30) with antibiotic, von
  Mises) is confirmed non-monotone (differences 0.95 and 1.05, not asymptotic). For the two A-89 rows the
  recomputation from the rounded grid values differs from the printed p/GCI by up to 0.05 and 0.2 %, which
  lies within the rounding of the printed grid values (checked by interval arithmetic); the printed values
  were computed on IKMHIWI03 from the unrounded numbers.
- Rigid bars (`*obstcolumns*dtc0.0025*inhliv*`): 7.1e-3 / 7.4e-4 / 9.5e-4 / 2.1e-4 Pa, "about ten times"
  (9.6), 9 % and 21 % lower growth, 14 % above the extrapolated value: all agree.
- Maze (`results/2026-10-maze/maze_h4_tend60*_k24fig*.json`): Table 6 constants, 42 and 133 um, living
  fraction 0.92 / 0.12 at T* = 24 / 0.35 at T* = 60 / 0.04 with the start-state rule: all agree.
- Antibiotic schedules (`*nf64*inhliv*abpulse*`, `*absine*`): growth and stress within 1 % of the constant
  dose (0.1 % and 0.8 %); the recovery in the pauses (see fixes).
- Time-step paragraph: the first-order limit 6.07e-4 Pa, "6 % below", "15 % above", 2.5 % and 12 % without
  antibiotic, the composition and growth changes (see fixes).
- Two-species composition across abstract, Ch. 4 and Ch. 6 (0.96 / 0.15 / 0.04, growth within 2 %).
- Largest phi_3D in the seed of the (0.25, 0.30) series (0.60, 0.69, 0.73).
- The five-species run's growth equals that of the earlier one-way two-species run (6.752e-4).
- Constants and conversions: Peclet number 31, diffusion length 0.14 mm, relaxation time 4e-8 T*,
  g l^2/d = 1, beta/d/r converted from Table 2 of Klempt et al. 2024 to the 2 mm cube, maze r = 200 um/T*.
- The molar crown is not mentioned in the thesis (nothing to check).

Not checkable from the repository (no result file here): the closed-form and shell numbers of Section 4.3
(ANSYS listings on IKMHIWI03), the example-input constants of the partner's deck, the step-2 Python values
of Section 4.7 (0.704, 0.735, 0.648, +6/-13 %, +5/-8 %), the Thiele number Lambda = 4 of step 1, and the
observed order 1.9 of the curved shell. These are unchanged.

## 2. Equations and units

- Ch. 2, Eqs. 34 to 36 as printed (Klempt et al. 2024, `references/*.txt` lines 500 to 516), Eq. 20 (energy),
  Eqs. 21 to 24 (growth and consumption terms), Table 2 and Table 6: agree with the thesis. Eq. (4.x)
  `eq:ch5-field` is the same equation written with the front term moved to the right-hand side
  (||grad phi|| n_phi = grad phi); the `\fterm` macro matches.
- Point model, Eqs. 16 to 18 of Klempt et al. 2026: the thesis's `eq:gov_phi`, `eq:gov_psi` and
  `eq:ch5-point` agree, including the derivatives of the dissipation function (checked by hand) and the
  +gamma of the printed Eq. 17, which the thesis drops and says so. Free energy and dissipation as in the paper.
- Coupling equations against the table "Equation check against Oliver's sketch" in `MEETING_2026-10-20.md`
  (E1 to E12): the thesis has dc/dt and da/dt (E2, E5), f(c) = c (E3), linear interpolation (E4, E6), Eq. 34
  with the front term (E7), the point model with its own constraint (E8, E9), the two examples (E10, E11)
  and Eq. 36 unchanged (E12). No change needed; the table stays as it is.
- Dimensions: every term of `eq:ch5-field`, `eq:ch5-eq36`, `eq:ch5-consumption`, `eq:ch5-cstar`,
  `eq:ch5-stiffness`, `eq:ch5-constrained` and `eq:ch5-flow` is consistent with beta, d in mm^2/T*,
  k_alpha, g in 1/T* and r in mm/T*. The one-element check: mu = 2 C10 = 4e-4 MPa, K = 2/D1 = 4e-4 MPa,
  nu = 0.125, closed form plus the spherical term as printed.
- Units of r: Table 2 of the paper prints r in 1/T*; the front term needs a speed. `PARAMETER_TABLES.md`
  already records this ("used as um/T* in this work"); the thesis said nothing. Fixed (see below).
- Notation rules: no `\varphi`, no `\alpha_K`, no `(1+\alpha)`, phi_3D in the field equations, phi_i /
  psi_i / phibar_i in the point model, k_alpha, no "L" for the living fraction. One `L` (the Thiele number
  in the list of symbols) was the cube half edge and is now `\ell`, as in Chapter 4.

## 3. Claims and consistency

Fixed (source unambiguous):

- Abstract and Ch. 6: "living volume fraction ... from 0.97 to 0.15" -> 0.96 (the no-antibiotic value is
  0.9646 / 0.961; the start value is 0.55). Ch. 4 already said "near 0.96".
- Ch. 4, two-species composition: "change by less than 3 % between the grids and the time steps" -> 4 %
  (separate example with antibiotic: 0.1515 at 64^3, step 0.025, against 0.1572 at 128^3, step 0.0025,
  a change of 3.8 %; A-89: 0.1508 -> 0.1562 between the steps 0.02 and 0.0025, 3.6 %).
- Ch. 4, time-step paragraph: "the composition and alpha-1 change by less than 3 and 0.3 %" held only for
  the separate example; in the mixed example alpha-1 changes by 3 % (8.30 -> 8.04e-4). Now: composition
  < 4 %, alpha-1 < 0.3 % (separate) and 3 % (mixed).
- Ch. 4, pulsed antibiotic: "recovers from 0.17 to 0.30" -> recovers from 0.17 to 0.34 in the last pause,
  0.30 at T* = 1 after the next pulse (result file).
- Ch. 4, five-species section: "the setting of the two-species runs" now names the earlier one-way
  computations (Section 4.7), since the two-species section now holds the agreed coupling; phi_cap, which
  was not defined anywhere, is written out as "the point model scaled to the field amount phi_3D up to a
  cap of 0.9".
- Ch. 6, mesh limitation: "the error of the seed stress is estimated at about 2 % from three meshes" ->
  "about 2 % below a finer solution of the same mechanics on up to 48^3 elements" (Ch. 4 states that a
  three-mesh fit would be misleading and uses the 48^3 reference).
- Ch. 6: "The ANSYS element solves the same equations and agrees within these bands" -> growth within 4 %,
  extrapolated stress within the band of the Python schemes (what Ch. 4 shows).
- Ch. 4, Table caption `tab:ch5-two-species`: "growth and stress converged values in Table ..." -> "grid
  series, extrapolated values and GCI" (the table has GCIs of 26 to 29 % and one non-asymptotic row).
- Ch. 4: "about a quarter of the converged value" -> "of the extrapolated Python value".
- Ch. 2: a fourth point on the unit of r (printed 1/T*, read as a speed; this reading is what the
  conversion to r = 10 mm/T* uses).

Checked and left: "verified / not validated" wording (consistent in abstract, Ch. 4 and Ch. 6); the GCI
statements of the abstract (1 to 29 %), Ch. 4 (0.7 to 3.3 %, 8 %, 26 to 29 %) and Ch. 6 (below 4 %, 8 %,
0.7 %, 26 to 29 %) agree with the table; the mesh numbers of Ch. 6 (1 %, 18 %, 0.7 %, 1.5 to 1.75) agree with
Ch. 4; 8.9e-4, 3e-11, 150 times, 15 coefficients agree across abstract, Ch. 4 and Ch. 6; "up to 8 %" in the
outlook matches "+5 and -8 %" of Section 4.7.

## 4. References and cross-references

- 66 cited keys, all in `references.bib`; the key entries (Klempt 2024 BMMB, Klempt 2026 Arch. Appl.
  Mech., Klempt 2023 PAMM, Fritsch 2025, Celik 2008 J. Fluids Eng., Rudolf 2025 FEAD, Blaszczyk 2022 CMAME,
  Feng 2021 Bull. Math. Biol., Soleimani 2023 Sci. Rep., Thiele 1939, Wanner and Gujer 1986, Pattem 2018,
  Bode 2025, Geisler and Junker 2025, Junker, Bode and Hackl 2025) have the expected year and journal.
- Every `\ref`/`\eqref` resolves (0 undefined in the log); bibtex gives 0 warnings.
- Figure captions checked against the figure scripts: `fig_coupled_3d_panels` (panels (a) to (e), T* = 0.325
  and 1, 64^3, one octant removed), `fig_maze_k24` ((a) without, (b) with antibiotic, T* = 12, 24, 36, 60),
  `fig_obstacle_columns_3d` ((a) set-up, (b) to (d) T* = 0.05, 0.1, 0.2, phi_3D = 0.5 and 0.2),
  `fig_thesis_cal5` ((a) shares, (b) living volume fractions). All match.

## 5. Build

86 pages, 0 errors, 0 undefined references or citations, 0 overfull boxes above 5 pt, audit 0 FAIL
(1 WARN: the 8^3 seed mean stress -2.515e-4 is printed as -2.52e-4, a half-unit rounding).

## Left for the user to decide

- The two GCI rows taken from report A-89 (separate and mixed with antibiotic, equal start amounts) have no
  result file in the repository; their p, extrapolated value and GCI are reproduced only within the rounding
  of the printed grid values. Exporting those five JSONs from IKMHIWI03 would make the check exact.
- The Thiele number of step 1 (Section 4.7, "Lambda = 4, l = 1 mm") is consistent with the symbol list
  (Lambda^2 = g l^2/d = 16) but could not be traced to a result file; the 20 Oct instructions to IKMHIWI03
  used Lambda^2 = 4 for the ANSYS deck of that step. Worth a look if the step-1 Python script is at hand.
- Observed order "1.9" of the curved shell (Section 4.3): the two errors quoted (2.7 % and 0.58 % against
  the finest mesh, one halving apart) give 2.2 by themselves; 1.9 presumably comes from the three-mesh fit.
  Not changed.
- Ch. 6 now says the element's growth agrees "within 4 %" (as Ch. 4); if the user prefers the earlier,
  shorter sentence, the wording is in this file's section 3.

## Update 28 Oct (192^3)

The 192^3 runs of the three antibiotic examples (A-92, A-94) are in `tab:ch5-gci`: GCI from 64^3/128^3/192^3. The
numbers above (8.3 %, 26 to 29 %, 1 to 29 %, the "--" row) are superseded: growth 0.5 to 4.3 %, separate stress 0.7 and
1.4 %, mixed 19 to 41 %, abstract 1 to 41 %. The antibiotic lowers the separate extrapolated stress by 2 % (6.55 to
6.42e-4 Pa), about the two bands together. Audit: 0 FAIL.

# Figures of the coupling paper (27 Oct 2026)

Plan for the manuscript, replacing the list "Figures (planned)" of `OUTLINE.md`. All figures in Times via
`ansys_usermat/figstyle.py` (TrueType in the PDF), drawn at the printed width, panel letters (a), (b) outside the axes,
no titles inside. 3D renderings as PNG at 400 dpi or full-resolution images inside the PDF.

| # | content | file | script | state |
|---|---|---|---|---|
| 1 | Coupling schematic: fields, point model at every point, mechanics | `assets/fig_coupling_schematic_agreed.png` | `ansys_usermat/coupling_schematic.py --agreed` | ready |
| 2 | Verification across codes and schemes (separate example with antibiotic; constant added diffusion against second-order schemes, ANSYS upwind) | `papers/coupling_paper/fig_convergence.pdf` | `fig_convergence.py` | ready; uses the start-position rule of A-52..A-64, redraw with the inherit-rule data if IKMHIWI03 repeats the series |
| 3 | Grid convergence of all two-species examples at the small outer step, with extrapolation | `papers/coupling_paper/fig_grid_examples.pdf` | `fig_grid_examples.py` | ready (27 Oct); add 192^3 of the mixed example when it comes |
| 4 | Manufactured solution: second-order convergence of the finite-volume and Galerkin codes (L2 errors of phi and c against h, slope-2 line) | `papers/coupling_paper/fig_mms.pdf` | `fig_mms.py` (Python from the `_mms_` JSONs of results/2026-10-3d, Abaqus from the table of A-50) | ready (27 Oct); observed orders 8-16-32: Python phi 1.85, 1.90, c 2.03, 2.01; Abaqus phi 1.68, 1.81, c 2.00, 2.05 |
| 5 | Separate and mixed example with antibiotic in 3D (share, sum phibar, c, a, von Mises; two times) | `assets/fig_coupled_3d_panels_{separate_ab,mixed_ab}.pdf` | `ansys_usermat/apdl/fig_coupled_3d_panels.py` | ready; one figure with both rows for the paper |
| 6 | Antibiotic over time (constant, pulsed, sine, same dose): held value a_b, seed sum phibar over time, seed alpha - 1 and von Mises at T* = 1 | `papers/coupling_paper/fig_ab_schedule.pdf` | `fig_ab_schedule.py` (runs `*nf64*inhliv*_abpulse0.2222_*`, `*_absine0.2222_*` of 27 Oct with the inherit rule, 64^3, about 6 min each in the cloud; constant = the main run) | ready (27 Oct); the history holds no alpha - 1 over time, so (c) shows the end values |
| 7 | Rule for new points: maze at T* = 60 with the antibiotic, neighbours' state (a) against start state (b), living fraction of the colonised cells over time (c) | `papers/coupling_paper/fig_maze_rules.pdf` | `fig_maze_rules.py` (maze_coupled_2d.py reruns of 27 Oct, `figdata/maze_h4_tend60_ab_k24fig{,_fresh}.npz` and the JSONs in results/2026-10-maze) | ready (27 Oct); the video frame `assets/anim_maze_rules_end.png` stays for talks |
| 8 | Geometry: rigid bars (3D surfaces) and the maze (phi_3D at four times) | `assets/fig_obstacle_columns_3d.pdf`, `assets/fig_maze_k24_h4.pdf` | `fig_obstacle_columns_3d.py`, `fig_maze_k24.py` | ready |

Order in the paper: 1 (Section 2), 2-4 (Section 3), 5-6 (Section 4.1-4.3), 7 (Section 4.4), 8 (Section 4.6).
Still to decide: whether the ANSYS element appears in Figure 2 only or also in Figure 5 (it gives the same fields within
the bands); whether the Keio options (homeostatic growth, front speed tied to the living biomass) get a figure in the
outlook or stay out.

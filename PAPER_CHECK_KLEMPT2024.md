# Chapter 5 and the decks checked against Klempt et al. 2024 (2026-10-02)

Paper: `references/Klempt2024_Hamilton_biofilm_growth_BMMB.pdf` (text in the
`.txt` next to it). Fixed in commits `2a1ba20` (notation) and `33507dc`
(content).

## Notation, now as in the papers

| Quantity | Paper | Was here |
|---|---|---|
| volume fraction | `\phi` (Klempt 2024, 2026, PAMM 2023) | `\varphi` |
| growth variable | `\alpha`, `F_g = \alpha I`, `\alpha(0) = 1`; growth is `\alpha - 1` | `\alpha_K`, own `\alpha = \alpha_K - 1`, `F_g = (1+\alpha) I` |
| growth rate | `k_\alpha` (Eq. 34/36, Table 2) | `k_a` |
| Monod constant | `k` (Eq. 34) | `K` |
| share of species i | no symbol in Klempt 2026: write `\phi_1/(\phi_1+\phi_2)` | `\chi_i` |

The UMAT's internal variable is still `\alpha - 1`; say so where code values
are quoted. Figure labels follow the same rule (`figstyle.py` for the font).

## Contradictions found and fixed

1. The paper's material is incompressible (`det F_e = 1`, Lagrange
   multiplier, H1P0); `\nu = 0.49` only defines `\mu`. Here it is
   approximated by `\nu = 0.49` in ANSYS (stated; 27 % change 0.45 -> 0.49).
2. The stiffness floor `f = 10^-3` is mine; the paper's void has no
   stiffness (stated as not from the paper).
3. Stress pattern: the paper reports pressure inside the biofilm, a ring of
   **tension** at its edge, zero outside. The seed's compression agrees;
   whether the neighbours are in tension is **not yet checked** (only their
   von Mises stress is known). The "same stress ring" wording was removed.
4. `\beta, r, k, d, g` in the ANSYS runs are the example input's, on a 2 mm
   cube (paper: 20 um), not the paper's. With the front term off, `r, k` and
   the nutrient do not act.
5. Time: the paper uses `10^3` substeps over `T* in [0, 1]`; the runs use
   `\Delta t = 0.1` up to `T* = 1.1`. Exact for the gradient-free growth;
   for the composition where the amount changes in time it is a limit
   (`assets/fig_coupling_convergence.png`: up to 0.07 / 0.29 in the share of
   species 1 during the transient at `\Delta t = 0.1`, first order).

## Inside the paper (noted, nothing to fix here)

- Table 1's local residual for `\alpha` differs in form from Eq. 36
  (`KLEMPT2024_REPRODUCTION.md`); this work follows Eq. 36.
- `\mu` is "in MPa" in the text, Pa in Table 2; Table 3 plots `p` in MPa.
  This work uses Pa.

## Open, needs IKMHIWI03 (Monday 8:00, about one hour, no new ANSYS run needed)

```powershell
git pull
$R = "ansys_usermat\apdl\results\2026-10-01_paper_values"
python ansys_usermat\apdl\plot_3d.py      "$R\all_stress_ds_pv_eq36.csv" --grid 8 --size 2.0 --out assets\stress3d_eq36.png
python ansys_usermat\apdl\plot_section.py "$R\all_stress_ds_pv_eq36.csv" --grid 8 --size 2.0 --out assets\section_eq36.png
python ansys_usermat\apdl\figs_1005.py    # redraw in Times New Roman with the new labels
```

- `plot_3d` must print OK (assumed element numbering); otherwise do not use
  the pictures.
- Read the neighbours' sign from the mean-stress panel: tension there would
  match the paper's ring (item 3).
  Faster, and as a number: `python ansys_usermat\apdl\neighbour_sign.py --csv
  "$R\all_stress_ds_pv_eq36.csv" --grid 8` prints p per layer around the seed. Python
  predicts (8^3, 3 Oct): seed -3.8e-5 Pa, first layer +5.8e-6 Pa with 56 % of
  its elements in tension, second layer +1.6e-6 Pa (79 %). `--grid 8` uses the
  same assumed numbering as plot_3d, so trust it only if plot_3d printed OK.
- If time allows: one run with `--post all` writes element centroids, so the
  pictures need no assumed numbering.
- Still to do on that machine: the ANSYS-vs-reference scatter over all Gauss
  points (from the existing trace CSVs) and the "macro unchanged with the
  point model on/off" comparison.
- Added 3 Oct: in Oliver's material routine, check how the stiffness adds
  the two species fields. The sources read on 2 Oct had `bio1 + bio1` where
  `bio1 + bio2` is expected. Oliver confirmed it is a typo (3 Oct, via the
  user): fix it to `bio1 + bio2` in the working copy. Then, before trusting
  the 5 Oct stress numbers:
  1. find whether the coupled runs pass through that line (one-species runs
     have bio2 = 0, so `bio1 + bio1` = 2 phi there; if the stiffness uses
     phi^2 the seed is 4x too stiff);
  2. if they do, rebuild and re-run the one-species paper-value decks
     (Eq. 36 and partner variants) and redraw `figs_1005.py`; Checks 2-3 may
     change in magnitude, the composition results do not (they do not depend
     on the stiffness);
  3. note the fix on the "Four things found" slide.
- Added 4 Oct: local nutrient in the point model (`ROADMAP_TWO_WAY.md`,
  "Step 1 tried in Python"). First only look, no new run needed: in the
  partner's material routine, find whether the nutrient at the Gauss point
  (`Nut1`, solved in `USSFin`) can be read in the material call (`GetVals`
  pool or similar). If yes, the next step is a fourth input to the bridge and
  c* = c*_0 c in the material server; then case 6 with a consumption strong
  enough for c to drop across the seed (Python: Thiele number >= 3,
  `ansys_usermat/composition_local_nutrient.py`).
- Added 3 Oct: run `ETLIST` in the deck (or read `ds.dat`) and check that
  type 1 (SOLID185) has KEYOPT(2) = 0, i.e. B-bar. With full integration the
  mean stress at nu = 0.49 would be 2-4x too large (volumetric locking,
  `ansys_usermat/apdl/locking_check.py`, MESH_STUDY.md).

## Open questions on the paper (3 Oct; sent to Felix by email on 3 Oct, also for Prof. Soleimani on 5 Oct)

Details in `KLEMPT2024_REPRODUCTION.md` sec. 12-13.
- Which form of Eq. 34/35 produced Fig. 4 and 7: growth on both faces of the
  colony, consumption g phi c? (Both are needed to reach the curves.)
- t_ref of each simulation: the curves fit Table 2 with a time scale of 1.5
  (4.1), 10 (4.2 high) and 4-5 (4.2 low).
- Is mu = 3.3557 in Pa or MPa? Table 3's pressure fits MPa (ANSYS micro-MKS).
- The nutrient source of 4.1: the edge line (Table 3) or the 5 um strip of
  Fig. 2?
- Added 3 Oct, not in the email: alpha by Eq. 36 (alpha-dot = k_alpha phi) or
  by Table 1 step 3a (driven by phi-dot)? Under Table 1 a seed that is full
  from the start does not grow (`KLEMPT2024_REPRODUCTION.md` sec. 4).

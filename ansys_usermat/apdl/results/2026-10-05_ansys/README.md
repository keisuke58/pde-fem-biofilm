# ANSYS runs of 5 Oct 2026 (IKMHIWI03), as JSON

One file per run, written by `ansys_usermat/apdl/export_runs_json.py` from the
git-ignored CSV outputs in `F:\biofilm_upf_wired`. Each file holds the deck's
first line (how it was generated), the seed elements (`seed_BIOFILM1`) and:

- `all_stress`: every element at the last result set: `elem, seqv, sx, sy, sz,
  alpha, cx, cy, cz` (stresses in MPa, deck units; `alpha` is the material
  routine's internal variable, i.e. alpha − 1 in the thesis notation; centroid
  in mm);
- `elem_stress`: one element over time (the deck's post element): `lstep, sbstep,
  time, sx … sxz, seqv, alpha, seqv_nbr_max`;
- `nut_field` (point-model runs only): `elem, nut1, phi1, phi2, alpha, cx, cy, cz`.
  Points where the point model never ran keep the start state phi1 = phi2 = 0.5
  (sum 1); points that ran have phi1 + phi2 ≤ phi_cap.

Element numbering is not x-fastest; use the centroid columns, not the element
number, for positions. What each run is: `ansys_usermat/apdl/RUN_1005_IKMHIWI03.md`.

| prefix | mesh |
|---|---|
| `ds8_`, `ds_` | 8³ (512 elements, seed 32) |
| `ds16_` | 16³ (4096, seed 256) |
| `ds24_` | 24³ (13824, seed 864) |

`ds8_ref` is the 8³ paper-value run (beta = 1e−4) re-run with the 5 Oct
executable; it equals `ds_kl_eq36` of 2 Oct to the last digit.

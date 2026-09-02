# Before the next ANSYS session

Written 2026-09-02, ahead of the first real-geometry runs. Everything here was
established without ANSYS, so none of it costs time at the machine. Three of
the four items are things that would have failed *at* the machine.

## 1. The layer generator could not read either real input. Fixed.

`make_layered_material.py` looked for a depth column in `x_mm, x_norm, x,
depth` and an alpha column in `alpha_x, alpha`. Neither producer writes those:

| Producer | depth | alpha |
|---|---|---|
| `JAXFEM/extract_alpha_field_1d.py` | `x_mm`, `x_norm` | `alpha_x_monod` |
| `_multiscale_results/macro_eigenstrain_<condition>_hybrid.csv` | `depth_mm`, `depth_norm` | `alpha_monod` |

The second is the per-condition field the results section needs. Both spellings
are now accepted, and the chosen columns are printed on the first output line.
Where a file carries more than one alpha column the generator refuses rather
than picking by list order, because `alpha_monod` and `alpha_x` are different
fields, not aliases; pass `--alpha-col` to resolve it.

The tests read both schemas out of the producers' own source, so a rename on
either side fails in CI rather than at the machine. The previous test invented
the schema and therefore agreed with the bug.

A second defect on the same path: the per-condition files carry a `#` comment
block, and `genfromtxt(names=True)` takes the first line as the header *after*
stripping `#` — so a comment containing a comma became the column names and the
read failed with a column-count error that named nothing relevant. The block is
now counted and skipped explicitly.

## 2. How many layers — check this first, it takes two seconds

```
python ansys_usermat/apdl/make_layered_material.py <field>.csv --convergence
```

It decides the deck structure, so run it before building geometry. On the 1D
extractor, swept across a 40x range of nutrient consumption:

| `--g-eff` | alpha range | 1 layer | 16 layers |
|---|---|---|---|
| 5 | 2.866–2.888e-4 | 0.16 % | 0.04 % |
| 20 | 2.802–2.888e-4 | 0.65 % | 0.16 % |
| 50 (default) | 2.669–2.888e-4 | 1.8 % | 0.42 % |
| 200 | 1.956–2.888e-4 | 10.1 % | 2.2 % |

On that model **one uniform layer is within a couple of percent of the field**,
and layering buys little — which would make the N-swept-volume deck unnecessary
and let the existing two-material shell stand.

**Do not carry that conclusion over to the per-condition fields.** It was
measured on the 1D extractor at its defaults, and those files are not in this
checkout, so their variation is unknown here. The `--convergence` run on the
real file is what settles it. The earlier table in `make_layered_material.py`'s
docstring (1 layer → 0.75) came from a deliberately steep synthetic profile and
is not representative either.

## 3. alpha magnitude collides with a known mesh limit

The 1D extractor at its defaults gives alpha ≈ 2.9e-4. The curved-shell deck
uses 0.01, and that number is **not a modelled value** — the deck's own header
records an alpha sweep finding its convergence threshold between 0.01 (clean)
and 0.015 (six errors, the same corner-distortion pattern as the earlier
`element highly distorted` failures).

So the mesh has a ceiling near 0.01 and the model, at that horizon, sits about
35x below it. alpha here is `k_alpha * integral(phi_total dt)`, so its size is
set by the integration horizon rather than being a constant: roughly, the
ceiling corresponds to ~35x the default `T*` if phi stays near its current
value. Two consequences worth having in mind rather than discovering:

- there is real headroom, so a longer horizon is available;
- if a per-condition field comes in near or above 0.01, the binding constraint
  is the mesh, and the fix last time was a boundary condition, not refinement
  (`MESH_STUDY.md`).

## 4. What is genuinely blocked

The condition comparison needs an alpha per condition, and this checkout cannot
produce one:

- `_multiscale_results/macro_eigenstrain_<condition>_hybrid.csv` is absent.
- `extract_alpha_field_1d.py` has **no condition argument** — conditions enter
  only through a TMCMC MAP theta (`--theta-json`) or the nutrient parameters,
  and no per-condition theta is in this checkout either.

So a session here can produce *one* alpha, not four. If those files exist on
IKMHIWI03, the path is short: `--convergence` on each, then the generator, then
`check_deck.py`. If they do not, regenerating them is the first task and it is
not a small one.

## Order of work

1. Locate `_multiscale_results/macro_eigenstrain_*_hybrid.csv`. Everything else
   depends on whether they exist.
2. `--convergence` on each. Decide layer count from the number.
3. Generate the material blocks; paste into a deck on the
   `t_growth_cylinder_shell_wrapper.dat` pattern.
4. `check_deck.py` before running — it catches the silent ones (alpha never
   written, too few `TB,STATE` slots, over-long `TBDATA`, dt/tau too coarse).
5. Run. Confirm dt is inside the stable range of the chapter's §5.5.3 for the
   actual job, and record it with the results.

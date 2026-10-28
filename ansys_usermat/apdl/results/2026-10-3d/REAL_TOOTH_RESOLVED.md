# Real molar: resolved stress of the nutrient and gate runs (20 Oct 2026)

`real_tooth_physics_resolved.py` re-solves the final state of every run of `REAL_TOOTH_PHYSICS.md` with the film submodel of
`real_tooth_submodel.py` (rigid tooth, crop box, void phi <= 0.05 dropped, wall = the field's own voxel staircase) at
h = 0.083 mm, six elements across the 0.5 mm film. The runs themselves used the 32^3 full-box stress (h = 0.5 mm), which
`REAL_TOOTH_CONVERGENCE.md` shows is not converged; their growth and nutrient fields are converged to about 2 %. The gate
runs were driven by the 32^3 pressure during the run; only their end state is re-solved here.

The runs are on the 96^3 field, whose tooth wall is a 0.167 mm staircase. Resolved below that staircase, the base run's film
mean is 3.66 mPa, against 3.26 to 3.29 mPa on the 192^3 field with its finer wall (`REAL_TOOTH_SUBMODEL.md`): the wall
geometry, not the mesh, makes the 12 % difference. All runs share the same wall, so I compare them by the ratio to the base
run, and take 3.3 mPa as the level of the base run.

| run | film mean vM (mPa, 96^3 wall) | ratio to base | band/core | 0.1 mm at the tooth (mPa) | 95th pct (mPa) | mean p (mPa) |
|---|---|---|---|---|---|---|
| no nutrient, no gate | 3.66 | 1.00 | 1.21 | 7.78 | 11.86 | 2.84 |
| nutrient (g/d = 10 /mm^2) | 0.51 | 0.14 | 1.12 | 0.94 | 1.29 | 0.38 |
| gate (p_h = 1.85 mPa on the 32^3 run) | 2.60 | 0.71 | 1.26 | 5.36 | 8.89 | 1.80 |
| nutrient + gate | 0.47 | 0.13 | 1.13 | 0.85 | 1.15 | 0.34 |
| nutrient, Table 2 ratio (1 /mm^2) | 2.47 | 0.68 | 1.20 | 5.20 | 7.86 | 1.92 |

What this gives, with the base run at 3.3 mPa:

- Nutrient supplied by the saliva and used up in the film (g/d = 10 /mm^2, chosen): film mean stress 0.14 of the base run, about 0.5 mPa.
- Nutrient at the Table 2 ratio: 0.68 of the base run, about 2.2 mPa.
- Pressure gate alone: 0.71 of the base run, about 2.3 mPa. On the 32^3 mesh the gate looked stronger (0.44): the coarse mesh
  underestimated the stress and the gate threshold was set from that coarse pressure, so this ratio depends on how p_h is chosen.
- Nutrient and gate together: 0.13; the nutrient dominates.
- The stress is highest in the first 0.1 mm next to the tooth (about twice the film mean) and slightly higher in the cervical band than
  above and below it (band/core 1.1 to 1.3). The pointwise maximum has no converged value (a corner of the staircase wall).

Limitations: the gate threshold and the gate's effect on growth come from the 32^3 pressure; a gate run with the submodel
stress in every frame would take about 1 min per frame (66 frames). The film volume does not respond to the nutrient or the
gate (no front term).

# Mesh study — the light version

Purpose: show that **the ratio Ch5 reports is mesh-stable**. Not a convergence
study of absolute stress; that is a bigger job and a different question (§3).

## Why the ratio and not the stress

§5.6 compares conditions — σ_CH/σ_DH and the like. Both runs use the same mesh,
so discretisation error largely cancels in the ratio. That is a much weaker
requirement than converged absolute stress, and it is what the claim actually
rests on.

The same logic is why §5.4 needs no mesh study at all: wrapper against core is
a build-to-build comparison on one identical mesh, so the error cancels exactly.

## The run

Three mesh levels, `ESIZE` halved each time — the current deck is **2 elements
through each layer's thickness**, which is why level 0 is a floor, not a
baseline:

| level | substrate / growth `ESIZE` | elements through thickness |
|---|---|---|
| 0 (current) | 0.15 / 0.05 | 2 / 2 |
| 1 | 0.075 / 0.025 | 4 / 4 |
| 2 | 0.0375 / 0.0125 | 8 / 8 |

`python ansys_usermat/apdl/make_mesh_levels.py <deck>.dat` writes the variants.
Run each condition at each level. Extract max and mean SEQV per run.

## Acceptance

**The ratio between conditions changes by less than 5 % from level 1 to level 2.**

If it does, report the level-2 ratio and state that it moved less than 5 % under
a halved element size. If it does not, the ratio is not yet mesh-independent and
the honest options are to refine further or to report the comparison
qualitatively — not to pick the level that looks best.

Absolute stresses will still be moving at these levels. Say so rather than
implying otherwise.

## Known hazard — expect level 2 to be the hard one

This is not hypothetical. The deck's own header records `element highly
distorted` failures at **`ESIZE = 0.033` and `ESIZE = 0.05`**, with different
elements failing each time (1216/1461/2, then 6061/1298/6232/6182), and
independent of α magnitude (0.02 and 0.05 both failed) and of substep count
(up to `NSUBST` 200 with `AUTOTS,ON`). The levels above put the growth layer at
0.025 and 0.0125 — inside and below that range.

So plan for level 2 to fail rather than being surprised by it. If it does:

- a level that will not solve is **a result about the mesh**, to be reported,
  not dropped quietly;
- level 0→1 alone still says something — if the ratio barely moves across one
  halving, that is weaker evidence than two, and should be described as such
  rather than presented as convergence;
- the fix that helped last time was a boundary condition, not a finer mesh
  (a symmetry-style `UZ=0` at the axial ends took the error count from 12 to
  3). Distortion here has been a BC problem as much as a resolution one, which
  is worth knowing before spending a session on element sizes.

## Not in scope

A convergence study of absolute stress on the curved-shell deck. Its own header
says it is a smoke test and that such a study "should precede using this for
anything quantitative" — whether that deck becomes the quantitative vehicle at
all is undecided, and settling it is not this study's job.

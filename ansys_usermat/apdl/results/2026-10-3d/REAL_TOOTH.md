# Seed model on a real molar surface (3D demonstration, 20 Oct 2026)

Script `ansys_usermat/apdl/anim_real_tooth_3d.py`; video `assets/anim_real_tooth_3d.mp4`
(2400 x 912, 16 fps, 90 frames, 10 MB), last frame `assets/anim_real_tooth_3d.png`,
numbers `real_tooth.json` (this directory), cache `figdata/anim_real_tooth_3d.npz`
(27 MB, git-ignored, rebuilt by the script when missing). This is a demonstration
of the Chapter 4 seed model on a realistic tooth shape, not a thesis result.

## Geometry: source and licence

I looked for an openly licensed tooth surface first. The marketplaces
(CGTrader, Cubebrush, Remeshy, MakerWorld) have molars only under royalty-free,
non-commercial or unstated terms; the DTU/3Shape FDI 16 crowns are CC BY-NC-SA
and open meshes without a root; the NIH 3D entry 3DPX-003003 is CC0 but a whole
dental arch on a base. The one that fits is the Mendeley Data set of Diaz, Ferro
Sanchez and Campo (Universidad Autonoma de Occidente), *Data of synthetic 3D
models of the human jaw, including teeth, ligaments, and bone structures*,
V1, 12 Sep 2024, doi:10.17632/xjsx7nfhj8.1, **CC BY 4.0**. It holds the 14 lower
teeth (FDI 31 to 47, STL/IGES/SpaceClaim), their periodontal ligaments and the
mandible. I use `Teeth/STL/D46.stl`, the lower right first molar, bundled
unchanged as `assets/geometry/mendeley_xjsx7nfhj8_D46.stl` (573 kB, closed
manifold surface, 11456 triangles, mm) with the attribution in
`assets/geometry/README.md`. Entry for `THIRD_PARTY.md` (not added there yet):

> `assets/geometry/mendeley_xjsx7nfhj8_D46.stl`: Diaz C., Ferro Sanchez C. A.,
> Campo O. (2024), Data of synthetic 3D models of the human jaw, including
> teeth, ligaments, and bone structures, Mendeley Data V1,
> doi:10.17632/xjsx7nfhj8.1, file Teeth/STL/D46.stl, CC BY 4.0, unchanged.

One point to be aware of: despite the word "synthetic" in the title, the
authors built the models from a CBCT scan of one adult (roots) and an intraoral
scan (crowns), and the dataset page says nothing about consent or anonymisation.
The authors released the data under CC BY 4.0 themselves, and a single tooth
surface does not identify anyone, so I consider it usable; the parametric
(idealised) molar remains the fallback if this is not wanted.

## Scale and placement

The tooth is used at real size (crown about 10 x 12 mm, 19 mm with the roots).
It is placed as oriented in the jaw (occlusal plane horizontal, z up), centred on
the crown, in a cube of 16 mm with the crown top at z = +5.5 mm; the bottom face
at z = -8 mm cuts the two roots about 4 mm below the furcation and is clamped,
which holds the tooth. The cemento-enamel junction is near z = -1.6 mm. 855 mm^3
of the tooth lie inside the cube. The gap between the crown and the walls is
at least 1.7 mm, so the film does not reach the box wall during the run (the
weakness of the half-cylinder case).

Diffusion coefficient: the thesis value beta = 0.02 mm^2/T* belongs to the 2 mm
cube, where the 0.25 mm seed spreads by sqrt(beta T) = 0.15 mm in T = 1.1 T*,
60 % of its thickness. The film here is 0.5 mm thick, so I scale beta by
(0.5/0.25)^2 = 4: **beta = 0.08 mm^2/T*** (spreading 0.3 mm). I tried the two
obvious alternatives and rejected them. The thesis value unchanged spreads
0.15 mm, below one cell of any affordable grid. Scaling beta with the cube edge
(x 8^2 = 1.28 mm^2/T*, the same dimensionless problem as the 2 mm cube with the
tooth shrunk 1:8) spreads 1.2 mm and dilutes the 0.5 mm film to about a quarter
within the run (max alpha - 1 = 3.3e-4 at T* = 0.9 instead of 8.8e-4), so the
film turns into a cloud. k_alpha is unchanged, so alpha - 1 and the stresses are
comparable with the 2 mm demonstrations (the small-strain mechanics has no
length scale). beta is the one assumption of this demonstration beyond the
Chapter 4 model.

## Settings

| | |
|---|---|
| field | d phi/dt = beta lap phi + k_alpha alpha - penalty, d alpha/dt = k_alpha phi; beta = 0.08 mm^2/T*, k_alpha = 1e-3 /T*, penalty 5, T = 1.1 T* |
| field grid | 96^3 cells, h = 0.167 mm, RK4, dt = 0.0167 T* (66 steps, 1 per frame); zero flux at the walls and at the tooth (masked Laplacian) |
| tooth | cells whose centre lies inside the STL (signed distance from `compute_implicit_distance`), 184601 cells; no growth, E = 1e4 Pa |
| seed | phi = 1 within 0.5 mm (3 cells) of the tooth surface, z = -3 .. -1 mm (cervical third of the crown, around the gum line), within 60 degrees of the -x flank; 2954 cells = 13.7 mm^3 |
| mechanics | small-strain hex8 B-bar (`mesh_study_seed.element_matrices`), 32^3 elements (h = 0.5 mm), E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain (alpha - 1) I, bottom face clamped; tooth elements by majority of the 3^3 field cells |
| solver | pyamg smoothed aggregation (rigid-body near-nullspace) + CG, tol 1e-8, 52 to 65 iterations, warm start; stress solved every third frame (23 solves) and interpolated linearly in between (no feedback on the growth) |
| video | 66 frames + 24 at the end, 16 fps; panels phi, alpha - 1, von Mises in the film; film = envelope of phi > 0.05, Taubin-smoothed, coloured by the maximum of the field through the film (sampled 1.5 mm along the inward normal); tooth = the STL, drawn half a cell inside its voxel set |

## Results (T* = 1.1)

| quantity | value |
|---|---|
| seed: mean phi | 0.54 (the film spreads normal to the surface and along it) |
| seed: mean alpha - 1 | 7.6e-4; maximum in the domain 1.01e-3 |
| film volume (sum phi h^3) | 13.7 -> 17.2 mm^3 |
| cells with phi > 0.5 / > 0.05 | 1857 / 10111 (seed 2954) |
| max von Mises in the film (phi > 0.3 on the stress grid, 130 elements) | **6.0 mPa** at (x, y, z) = (-3.75, 4.75, -2.25) mm, the end of the band on the +y side |
| 95th percentile of von Mises | 5.2 mPa (colour scale of the video) |
| max pressure -tr(sigma)/3 in the film | 16.5 mPa |
| wall time | signed distance 10 s, field 16 s, stress 655 s (23 solves, 28 to 33 s each), render 60 s; 12.5 min in all on 4 cores |

For comparison, the half-cylinder case of `anim_scenarios_3d.py` (0.25 mm film
on a 0.7 mm cylinder in the 2 mm cube, beta = 0.02) has the same model and the
same k_alpha; its stresses are of the same order (mPa).

## What the video shows

The film sits as a band around the cervical third of the crown on one flank
(a plaque location near the gum line). Over the run it thickens and spreads a
little along the tooth (0.3 mm), the growth alpha - 1 builds up uniformly in the
seed to about 1e-3, and the von Mises stress in the film grows with it, with the
highest values at the ends of the band and where the surface curves most. The
tooth does not deform visibly (E ratio 1000, displacement not magnified).

## Limitations

- The stress mesh is coarse: with h = 0.5 mm the film is one element thick, so
  the von Mises pattern in the video shows the element grid (the chequer
  pattern). A finer stress mesh (48^3, 330k dofs) costs about 100 s per solve
  with pyamg; the stress is solved only every third frame for the same reason.
- beta is scaled (x 4) as explained above; the thesis value at real size would
  need a 0.05 mm grid to resolve its front.
- Constant growth rate, no pressure law, no composition: the model is the
  Chapter 4 seed model only. The tooth is rigid and not grown; nothing below the
  gum line (ligament, bone, gingiva) is modelled; the space outside the film has
  the small floor stiffness, not a fluid.
- The STL surface is voxelised on the field grid (cell centre inside/outside);
  the no-flux condition at the tooth is therefore a staircase at 0.17 mm. The
  rendering hides this (smoothed envelope, STL drawn half a cell inside), the
  computation does not.
- The film colouring in the video is the maximum of the field through the film
  (within 1.5 mm below the surface), chosen so that the thin diffuse front at the
  envelope does not hide the film; it is stated on the frame.
- The dataset's crown comes from an intraoral scan and the roots from CBCT; the
  mesh has a few flat facets at the occlusal fissures (visible in the render).

## How to rerun

```
python ansys_usermat/apdl/anim_real_tooth_3d.py            # full run, about 13 min
python ansys_usermat/apdl/anim_real_tooth_3d.py --png-only # re-render the last frame from the cache
```
Needs pyvista (off screen), imageio-ffmpeg, scipy and pyamg (`pip install pyamg`).

# The agreed coupling as an Abaqus user element (A-46, 22 Oct 2026)

The coupling agreed with Oliver (`MEETING_2026-10-20.md`; Klempt et al. 2024 Eqs. 34-36 with the front term, the point
model of Klempt et al. 2026 as published at the Gauss points, consumption by the living biomass, the antibiotic field)
in a standard Galerkin finite element, as a third implementation next to the ANSYS element (NEM, `F:\biofilm_upf_living8`)
and the Python finite-volume reference (`ansys_usermat/apdl/coupled_living_consumption.py`). The equations are the same;
only the discretisation differs. In a Galerkin element the phi Laplacian is symmetric and conserves phi and the walls
need no treatment (A-38..A-45: both had to be fixed in the other two codes).

Abaqus is not part of the LUH thesis (CLAUDE.md); this is the Keio continuation and a cross-check.

| file | what |
|---|---|
| `living_uel.f90` | the element (free form): 8-node hexahedron, dofs u (1-3), phi (11), c (12), a (13); equations in its header |
| `make_abaqus_source.py` | `living_uel_abaqus.for` = `ecology_native.f` + the element in fixed form (the p3 converter) |
| `test_tangent.f90` | closed-form tangent against central differences of the residual (22 Oct: 1e-10) |
| `make_living_inp.py` | deck: n^3 cube, `sep` or `mixeq`, increment, antibiotic start; the set-up of the Python reference `--v2 --ansys-bc --start-position --front 1 --artdiff` |
| `summary_living.py` | means per time from the Gauss-point file (seed: wp3 seed cells) |
| `eco_u3_a100*.txt` | the point-model files of the ANSYS runs (Fritsch et al. 2025 case I, alpha*_0 = 100) |

Run on IKMHIWI03 (Abaqus 2024; not from a shell with `dev-env.ps1` dot-sourced):

```
python abaqus_living/make_abaqus_source.py
python abaqus_living/make_living_inp.py 8 sep abaqus_living/lv8_sep_ab.inp --dt 0.005
$env:BIOFILM_ECO_CASE = "<repo>\abaqus_living\eco_u3_a100.txt"; $env:BIOFILM_UEL_GP = "F:\abaqus_work\lv8_sep_ab\living_gp.txt"
.\run_abaqus.ps1 -InpFile abaqus_living/lv8_sep_ab.inp -UserFile abaqus_living/living_uel_abaqus.for -LibsDll
python abaqus_living/summary_living.py F:\abaqus_work\lv8_sep_ab\living_gp.txt 8
```

Procedure as p3/uel (A-11): a 1e-8 step prescribes phi, c, a at every node, then the growth step; transient
`*COUPLED TEMPERATURE-DISPLACEMENT` with fixed increments, `UNSYMM`, dummy C3D8T elements with 1e-9 material.
All `*BOUNDARY` of a step must use the same OP (Abaqus refuses OP=NEW mixed with OP=MOD).

"""write_eco_cfg.py -- the case constants of the point model for the Fortran
version (ecology_native.f, eco_native_init), taken from material_server.set_case
so the native run uses exactly what the server would.

    python abaqus_composition/write_eco_cfg.py CASE OUT.txt

OUT: n_active / c* alpha* / eta(5) / theta(20).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "ansys_usermat" / "coupling"), str(ROOT / "ansys_usermat")]
import material_server as ms  # noqa: E402


def main(case, out):
    ms.set_case(case)
    c = ms.ECOLOGY_CASE
    hp = c["hp"]
    f = lambda v: " ".join(repr(float(x)) for x in v)  # noqa: E731
    Path(out).write_text(f"{ms.ECOLOGY_ACTIVE}\n{f([hp['c'], hp['alpha']])}\n{f(hp['eta'])}\n{f(c['theta'])}\n")
    print(f"wrote {out}: {case}, {ms.ECOLOGY_ACTIVE} species, c* {hp['c']}, alpha* {hp['alpha']}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

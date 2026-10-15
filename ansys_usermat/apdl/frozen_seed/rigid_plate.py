"""rigid_plate.py -- request A-18: may the tooth be replaced by u = 0 on the
surface the biofilm sits on? The A-16 problem (frozen_seed.dat, 8^3) as the
reference, and the same cube on a plate of thickness 0.25 (x, y as the cube,
two elements through the thickness, SOLID185, linear elastic, nu = 0.3, no
growth), u = 0 on the plate's underside, for three plate stiffnesses.

The plate and the cube are glued as volumes (VGLUE) before meshing, so they
share the nodes of z = 0; the run checks this by the node and element counts.

    python rigid_plate.py [WORKDIR] [N]      (default F:\\biofilm_upf_frozen, 8)

Needs the A-16 exe in WORKDIR (build_and_run.ps1). Writes
results/2026-10-uel/rigid_plate_<N>.json: per run the measures of A-18 and
their relative difference to the reference.
"""
import csv
import itertools
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "p3" / "uel"))
from make_seed_inp import seed_phi  # noqa: E402

PLATES = [None, 1e10, 1e4, 1e2]        # None: the A-16 reference (u = 0 at z = 0)
THICK, NU_PLATE = 0.25, 0.3


def deck(n, e_plate):
    s = (HERE / "frozen_seed.dat").read_text().replace("\r\n", "\n")
    s = s.replace("\nN = 8\n", f"\nN = {n}\n", 1)
    if e_plate is None:
        return s.replace("\nFINISH\n/POST1\n", "\n*GET,NNMAX,NODE,0,NUM,MAX\nFINISH\n/POST1\n", 1)
    rep = [
        ("TB,STATE,1,,12\n",
         f"TB,STATE,1,,12\nMP,EX,2,{e_plate!r}\nMP,PRXY,2,{NU_PLATE}\n"),
        ("BLOCK,0,1,0,1,0,1\nESIZE,1/N\nMSHKEY,1\nVMESH,ALL\n",
         f"BLOCK,0,1,0,1,0,1\nBLOCK,0,1,0,1,-{THICK},0\nVGLUE,ALL\nESIZE,1/N\nMSHKEY,1\n"
         # VSEL,LOC with one value needs the centroid exactly: give ranges
         "VSEL,S,LOC,Z,0,1\nVATT,1,,1\nVSEL,S,LOC,Z,-1,0\nVATT,2,,1\nVSEL,ALL\nVMESH,ALL\n"
         "ESEL,S,MAT,,2\nESEL,ALL\n"),
        ("NSEL,S,LOC,Z,0\n", f"NSEL,S,LOC,Z,-{THICK}\n"),
    ]
    for old, new in rep:
        if s.count(old) != 1:
            raise RuntimeError(f"frozen_seed.dat changed: {old!r} found {s.count(old)} times")
        s = s.replace(old, new)
    return s


def counts(wd, job):
    """node and element counts from the solver output (the mesh summary)"""
    out = (wd / f"{job}.out").read_text(errors="replace")
    nn = ne = n2 = None
    for line in out.splitlines():
        if n2 is None and "SELECTED BY  ESEL" in line:
            n2 = int(line.split()[0])
        if "MAXIMUM ELEMENT NUMBER" in line:
            ne = int(line.split("=")[-1])
        if "MAXIMUM NODE NUMBER" in line:
            nn = int(line.split("=")[-1])
    return nn, ne, n2


def measures(wd, n):
    h = 1.0 / n
    rows = [{k.strip(): float(v) for k, v in r.items()} for r in csv.DictReader(open(wd / "frozen_elems.csv"))]
    cen = np.array([[r["cx"], r["cy"], r["cz"]] for r in rows])
    vm = np.array([r["seqv"] for r in rows])
    p = np.array([(r["sx"] + r["sy"] + r["sz"]) / 3 for r in rows])
    xi = np.array([r["alpha"] for r in rows])
    bio = cen[:, 2] > 0
    corners = np.array(list(itertools.product((-0.5, 0.5), repeat=3))) * h
    pe = np.array([seed_phi(c + corners).mean() for c in cen])
    seed = bio & (pe > 0.5)
    core = bio & (np.sqrt((cen[:, 0] - .5) ** 2 + (cen[:, 1] - .5) ** 2 + cen[:, 2] ** 2) < 0.2)
    lay = bio & np.isclose(cen[:, 2], h / 2)
    m = dict(seed_elements=int(seed.sum()), seed_alpha_minus_1=float(xi[seed].mean()), seed_vM=float(vm[seed].mean()),
             core_elements=int(core.sum()), core_vM=float(vm[core].mean()),
             layer1_elements=int(lay.sum()), layer1_vM=float(vm[lay].mean()), layer1_p=float(p[lay].mean()),
             max_vM_biofilm=float(vm[bio].max()), max_vM_centroid=cen[bio][vm[bio].argmax()].round(6).tolist(),
             biofilm_elements=int(bio.sum()), plate_elements=int((~bio).sum()))
    if (~bio).any():
        m["max_vM_plate"] = float(vm[~bio].max())
    return m


def run(wd, n, e_plate):
    job = f"rigid_plate_{n}_" + ("ref" if e_plate is None else f"E{e_plate:.0e}".replace("+", ""))
    for f in wd.glob(f"{job}.*"):
        f.unlink()
    for f in ("frozen_elems.csv", "frozen_gp.csv"):
        (wd / f).unlink(missing_ok=True)
    (wd / f"{job}.inp").write_text(deck(n, e_plate).replace("\n", "\r\n"), encoding="ascii", newline="")
    exe = Path(os.environ["AWP_ROOT222"]) / "ANSYS" / "bin" / "winx64" / "ANSYS222.exe"
    r = subprocess.run([str(exe), "-b", "-smp", "-np", "1", "-custom", r".\ANSYS.exe", "-j", job,
                        "-i", f"{job}.inp", "-o", f"{job}.out"], cwd=wd, capture_output=True)
    out = (wd / f"{job}.out").read_text(errors="replace")
    nerr = out.count("*** ERROR ***")
    nn, ne, n2 = counts(wd, job)
    m = measures(wd, n)
    m.update(run=job, E_plate=e_plate, exit=r.returncode, errors=nerr, nodes=nn, elements=ne)
    if e_plate is not None:
        m["plate_elements_mat2"] = n2
    print(job, "exit", r.returncode, "errors", nerr, "nodes", nn, "elements", ne, "material 2", n2)
    return m


def main(wd=r"F:\biofilm_upf_frozen", n=8):
    wd, n = Path(wd), int(n)
    res = [run(wd, n, e) for e in PLATES]
    ref = res[0]
    keys = ["seed_alpha_minus_1", "seed_vM", "core_vM", "layer1_vM", "layer1_p", "max_vM_biofilm"]
    for m in res[1:]:
        m["rel_diff"] = {k: m[k] / ref[k] - 1 for k in keys}
        # shared nodes at z = 0: the plate adds (n+1)^2 * 2 nodes and n^2 * 2 elements
        m["shared_interface"] = (m["nodes"] == ref["nodes"] + 2 * (n + 1) ** 2
                                 and m["elements"] == ref["elements"] + 2 * n * n)
        if m["plate_elements_mat2"] != 2 * n * n:
            raise RuntimeError(f"{m['run']}: {m['plate_elements_mat2']} elements of material 2, expected {2 * n * n}")
    out = dict(request="A-18", n=n, plate_thickness=THICK, plate_nu=NU_PLATE, plate_layers=2, E_biofilm=10.0,
               reference="frozen_seed.dat (A-16), u = 0 at z = 0",
               measures="seed = biofilm elements with mean nodal phi > 0.5; core = centroid within 0.2 of (0.5, 0.5, 0); "
                        "layer1 = biofilm elements with centroid z = h/2; p = (SX + SY + SZ)/3 (tension > 0)",
               runs=res)
    dst = REPO / "ansys_usermat" / "apdl" / "results" / "2026-10-uel" / f"rigid_plate_{n}.json"
    dst.write_text(json.dumps(out, indent=1))
    for m in res:
        print(m["run"], {k: (f"{m[k]:.6e}" if isinstance(m[k], float) else m[k]) for k in keys + ["max_vM_plate"] if k in m},
              {k: f"{v:+.2e}" for k, v in m.get("rel_diff", {}).items()}, "shared", m.get("shared_interface"))
    print("wrote", dst)


if __name__ == "__main__":
    main(*sys.argv[1:3])

"""make_living_deck.py -- run decks for the coupling agreed with the partner on 20 Oct 2026 (A-19,
MEETING_2026-10-20.md, CLOUD_TO_IKMHIWI03_2026-10-20.ja.md): the point model as a black box at every Gauss point
(usermat prop(28) = 11), consumption by the living biomass g phi L c in the element (LIVING_CONS = 1), the
front term of Klempt 2024 Eq. 34 with |grad phi| (FRONT_GRAD = 1), the second nutrient field as the antibiotic
(source on the face y = +1 mm from T* = NUT2_TSTART on). Same set-up as the Python version
coupled_living_consumption.py.

    python ansys_usermat/apdl/make_living_deck.py BASE.dat OUT.dat --layout sep|mix [--ab] [--front 1|0.1|0]
        [--l1] [--cons G]

  --layout sep : species A (point-model species 1) in the Chapter 4 seed (component BIOFILM1 of the base deck),
                 species B in the four vertical edge columns |x|, |y| > 0.75 mm; one field phi (bio1) = 1 in both
  --layout mix : A and B together in the seed, the file's start (Fritsch Case I: 0.25, 0.30)
  --ab         : antibiotic a = 1 on the elements of the layer y = +1 mm from T* = 0.5 (Dirichlet, as the nutrient)
  --front F    : r = F x 10 mm/T* (Klempt 2024 Table 2, r = 100 um/T* on the 20 um cube, converted to the 2 mm
                 cube), k = 1; 0 switches the front term off (ORI_WEIGHT11 = 0)
  --l1         : L = 1 (prop(57) = 1): the consumption without the living share, for check 2
  --cons G     : g for both fields (default 4: Thiele number g L^2/d = 4 with L = 1 mm, d = 1 mm^2/T*)
  --v2         : example setting v2 (CLOUD_TO_IKMHIWI03_2026-10-20.ja.md, addendum 2; coupled_living_consumption.py
                 --v2): phi_min = 0.1, species B in four small corner blocks (|x|, |y| > 0.75 mm and |z| within the
                 seed's height), antibiotic uptake g_a = 0.25 (Lambda_a^2 = 0.25) and a = 1 on the elements of all
                 outer faces from T* = 0.5. alpha*_0 = 50 comes from the point-model file (eco_fritsch_2sp_a50.txt).
  --start-step S : a first load step of one increment of length S (e.g. 0.0025, one USSFin step; shorter steps are
                 rounded up to MY_USSFIN_DT_MIN inside the element) in which only the fields are solved (prop(60) = 1:
                 the point model takes its start state but is not advanced), then the run from S to 1 + S in the
                 base deck's increments, so the point model starts from the solved nutrient field (addendum 4).

The base deck is one of the partner's decks as derived for Chapter 4 (base_w<n>_c6_g1_s015.dat, never
committed). The nutrient keeps the base deck's component NUTRIENT1 (the element layer at y = -1 mm, c = 1).
Point model: BIOFILM_ECO_CASE = abaqus_composition/eco_fritsch_2sp.txt (theta in prop(8:27) from it).
Post: post_fields_living.mac (fields at T* = 0.25, 0.5, 0.75, 1) after the base deck's own post macros.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from deck_components import centroids, cmblock, read_components, read_mesh  # noqa: E402

CASE = ROOT / "abaqus_composition" / "eco_fritsch_2sp.txt"
R_TABLE2 = 100.0e-3 * (2.0 / 20.0e-3)   # 10 mm/T*


def case_theta():
    lines = [l.split() for l in CASE.read_text().splitlines() if l.strip()]
    th = [float(x) for x in lines[3]]
    assert len(th) == 20
    return th


def set_scalar(text, name, value):
    rx = re.compile(rf"^({re.escape(name)}\s*=\s*)([^\s!]+)", re.M)
    hits = rx.findall(text)
    if len(hits) != 1:
        sys.exit(f"STOP: {name}: {len(hits)} assignments, expected 1")
    return rx.sub(lambda m: m.group(1) + value, text)


def replace_cmblock(text, name, ids):
    rx = re.compile(rf"^CMBLOCK,\s*{name}\s*,ELEM,\s*(\d+)\s*\r?\n\(8i10\)\r?\n((?:[ \d]+\r?\n)+)", re.M | re.I)
    m = rx.search(text)
    if not m or len(rx.findall(text)) != 1:
        sys.exit(f"STOP: component {name} not found once")
    nl = "\r\n" if "\r\n" in m.group(0) else "\n"
    return text[:m.start()] + cmblock(name, ids).replace("\n", nl) + nl + text[m.end():]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("out")
    ap.add_argument("--layout", choices=("sep", "mix"), required=True)
    ap.add_argument("--ab", action="store_true")
    ap.add_argument("--front", type=float, default=0.1)
    ap.add_argument("--l1", action="store_true")
    ap.add_argument("--cons", type=float, default=4.0)
    # inner USSFin step. 0.1 (>= the outer increment) = one inner step, as in all Chapter 4 decks. The element advances
    # phi by the outer increment at every inner step (A-25), so a smaller value makes phi run (outer / inner) times
    # too fast; the earlier default 0.0025 did exactly that (living runs before u1_, 21 Oct 2026).
    ap.add_argument("--dt-ussfin", type=float, default=0.1)
    ap.add_argument("--old-cons", action="store_true",
                    help="LIVING_CONS = 0: the element's own consumption (g phi, zero order) for comparison")
    ap.add_argument("--v2", action="store_true")
    ap.add_argument("--start-step", type=float, default=0.0)
    # sensitivities (addendum 5): each overrides one value of the setting
    ap.add_argument("--phi-min", type=float, help="prop(29), default 0.1 with --v2, else 0.01")
    ap.add_argument("--cons-ab", type=float, help="CONSUMPTION21, default 0.25 with --v2, else --cons")
    ap.add_argument("--deltim", type=float, help="outer increment of the run (with --start-step), default the base's")
    ap.add_argument("--neigh", type=int, help="NEIGHBOR_CNT, the element's maximum number of neighbours (base: 30); "
                    "the search keeps the first ones found in element order (A-23, the symmetry check)")
    ap.add_argument("--artdiff", type=int, metavar="N",
                    help="stabilise the front term (addendum 27): MY_BETA1 = beta + r h / 4 with h = 2 mm / N (N elements "
                         "per edge), i.e. the numerical diffusion v h / 2 of first-order upwinding at the largest front "
                         "speed v = r c/(k + c) = r / 2; the equations are unchanged as h -> 0")
    ap.add_argument("--uniform", action="store_true",
                    help="phi = 1 in every element (BIOFILM1 = all), for the uniform-field check S3 (addendum 15); "
                         "use with --cons 0 --front 0 --l1 and no --ab, so c = 1 stays the solution")
    ap.add_argument("--tend", type=float,
                    help="end time T* of the run (with --start-step; default the base's, 1); fields are then written at "
                         "T*/4, T*/2, 3T*/4, T* (post_fields_living_t<T*>.mac, e.g. _t2 for 2). Deck values only "
                         "(addendum 21); the antibiotic still starts at T* = 0.5")
    ap.add_argument("--t-ab", type=float, default=0.5,
                    help="with --ab: T* at which the antibiotic source starts (NUT2_TSTART, default 0.5; addendum 22)")
    ap.add_argument("--transient-c", type=float, metavar="D",
                    help="dc/dt in both fields (addendum 28, build F:\\biofilm_upf_living5): NUT_TRANSIENT = 1, "
                         "MY_DIFF1 = MY_DIFF2 = D and both consumptions times D, so g/d stays as given by --cons/--cons-ab "
                         "(the quasi-static solution depends on g/d only); D = 1e8 is Klempt 2024 Table 2")
    ap.add_argument("--nut1-init-all", action="store_true",
                    help="NUT1_INIT_ALL = 1 (A-38, build F:\\biofilm_upf_living6): the nutrient starts at c = 1 in every "
                         "element; otherwise only on the supply elements and 0 elsewhere, so the first increments "
                         "see c = 0 (no front term, no consumption)")
    ap.add_argument("--front-dump", type=float, metavar="T",
                    help="FRONT_DUMP = T (A-39, build living6): USSFin writes the terms of the phi update of every point "
                         "(phi, grad phi, grad c, c, lap phi, Ori, Growth, penalty) to front_dump.txt while the run time "
                         "is below T; for short 8^3 runs only")
    ap.add_argument("--surf-nearest", action="store_true",
                    help="SURF_NEAREST = 1 (A-39, build living6): each surface point takes phi of its nearest "
                         "integration points (zero normal gradient) instead of the inverse-square mean of all its "
                         "integration-point neighbours, which let phi leak through the walls")
    ap.add_argument("--nem-dump", action="store_true",
                    help="NEM_DUMP = 1 (A-42, build F:\\biofilm_upf_living7): USSFin writes the NEM operator (neighbours, "
                         "dx, dy, dz, Lap of every point) and the point coordinates to nem_dump.txt at its first call")
    ap.add_argument("--lap-file", action="store_true",
                    help="LAP_FILE = 1 (A-42, build living7): the phi Laplacian from lap_cons.txt in the work directory "
                         "(checks_a38_a39/build_lap_cons_a42.py from nem_dump.txt of the same mesh): the NEM Laplacian "
                         "made conservative; use with --surf-nearest")
    a = ap.parse_args()

    nodes, elems = read_mesh(a.base)
    cen = centroids(nodes, elems)
    comps = read_components(a.base)
    seed = set(comps["BIOFILM1"])
    h = 2.0 / round(len(elems) ** (1 / 3))
    cols = {e for e, c in cen.items() if abs(c[0]) > 0.75 and abs(c[1]) > 0.75}
    if a.v2:
        zmax = max(abs(cen[e][2]) for e in seed) + 1e-9
        cols = {e for e in cols if abs(cen[e][2]) <= zmax}
        top = {e for e, c in cen.items() if max(abs(c[0]), abs(c[1]), abs(c[2])) > 1.0 - h}
    else:
        top = {e for e, c in cen.items() if c[1] > 1.0 - h}
    bio1 = seed | cols if a.layout == "sep" else seed
    if a.uniform:
        bio1 = set(elems)

    th = case_theta()
    props = {28: 11, 29: 0.1 if a.v2 else 0.01, 31: 0.1, 37: 2, 56: 1 if a.layout == "sep" else 2,
             57: 1 if a.l1 else 0, 58: 0.55, 59: 0.5625}
    if a.start_step > 0:
        props[60] = 1
    if a.phi_min is not None:
        props[29] = a.phi_min
    if (a.deltim is not None or a.tend is not None) and a.start_step <= 0:
        sys.exit("STOP: --deltim and --tend only with --start-step")
    props.update({7 + i: v for i, v in enumerate(th, start=1)})
    tmp = Path(a.out).with_suffix(".tmp.dat")
    cmd = [sys.executable, str(HERE / "make_wired_deck.py"), a.base, str(tmp),
           "--props", ",".join(f"{k}={v!r}" for k, v in sorted(props.items())), "--post", "all"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(r.stdout + r.stderr)
    text = tmp.read_text(encoding="latin-1")
    tmp.unlink()

    r_front = a.front * R_TABLE2
    scal = {"MY_BIOSTART1": "1.0", "MY_BIOSTART2": "0.0", "MY_NUTSTART1": "1.0",
            "MY_NUTSTART2": "1.0" if a.ab else "0.0", "MY_DIFF1": "1", "MY_DIFF2": "1",
            "CONSUMPTION11": repr(a.cons), "CONSUMPTION12": "0",
            "CONSUMPTION21": repr(a.cons_ab if a.cons_ab is not None else 0.25 if a.v2 else a.cons),
            "CONSUMPTION22": "0", "MAX_GROWTH11": repr(r_front), "MAX_GROWTH21": "0", "HALF_VELO11": "1",
            "ORI_WEIGHT11": "1" if a.front > 0 else "0", "ORI_WEIGHT21": "0",
            "MY_USSFIN_DT_MIN": repr(a.dt_ussfin), "MY_USSFIN_DT_MAX": repr(a.dt_ussfin)}
    if a.neigh is not None:
        scal["NEIGHBOR_CNT"] = str(a.neigh)
    if a.transient_c:
        d = a.transient_c
        scal["MY_DIFF1"] = scal["MY_DIFF2"] = repr(d)
        for k in ("CONSUMPTION11", "CONSUMPTION21"):
            scal[k] = repr(float(scal[k]) * d)
    if a.artdiff:
        mb = re.search(r"^MY_BETA1\s*=\s*([\d.eE+-]+)", text, re.M)
        beta0 = float(mb.group(1))
        scal["MY_BETA1"] = repr(beta0 + r_front * (2.0 / a.artdiff) / 4.0)
        print(f"artdiff: MY_BETA1 {beta0:g} -> {float(scal['MY_BETA1']):g} (h = {2.0 / a.artdiff:g} mm, r = {r_front:g})")
    for k, v in scal.items():
        text = set_scalar(text, k, v)
    nl = "\r\n" if "\r\n" in text else "\n"
    m = re.search(r"^PENALTY1\s*=.*$", text, re.M)
    new = (f"{nl}! pde-fem-biofilm 20 Oct 2026 (A-19): living share, front term, antibiotic{nl}"
           f"LIVING_CONS = {0 if a.old_cons else 1}{nl}FRONT_GRAD = 1{nl}NUT2_TSTART = {a.t_ab if a.ab else 0.0}{nl}")
    if a.transient_c:
        new += f"! addendum 28: dc/dt in both fields (build living5){nl}NUT_TRANSIENT = 1{nl}"
    if a.nut1_init_all:
        new += f"! A-38: the nutrient starts at c = 1 in every element (build living6){nl}NUT1_INIT_ALL = 1{nl}"
    if a.surf_nearest:
        new += f"! A-39: surface points take phi of their nearest integration points (build living6){nl}SURF_NEAREST = 1{nl}"
    if a.lap_file:
        new += f"! A-42: conservative phi Laplacian from lap_cons.txt (build living7){nl}LAP_FILE = 1{nl}"
    if a.nem_dump:
        new += f"! A-42: write the NEM operator once (build living7){nl}NEM_DUMP = 1{nl}"
    if a.front_dump:
        new +=f"! A-39: dump the terms of the phi update (build living6){nl}FRONT_DUMP = {a.front_dump!r}{nl}"
    text = text[:m.end()] + new + text[m.end():]
    text = replace_cmblock(text, "BIOFILM1", bio1)
    text = replace_cmblock(text, "NUTRIENT2", top)
    if a.start_step > 0:
        s = a.start_step
        rx = re.compile(r"^deltim,([^,\s!]+)[^\r\n]*\r?\ntime,([^\s!]+)[^\r\n]*$", re.M | re.I)
        m = rx.search(text)
        if not m or len(rx.findall(text)) != 1:
            sys.exit("STOP: deltim/time pair not found once")
        dt, t_end = float(m.group(1)), float(m.group(2))
        if a.deltim is not None:
            dt = a.deltim
        if a.tend is not None:
            t_end = a.tend
        text =text[:m.start()] + f"deltim,{s!r},{s!r},{s!r}{nl}time,{s!r}" + text[m.end():]
        rs = re.compile(r"^OUTRES,ALL,ALL\r?\nsolve[ \t]*$", re.M)
        m = rs.search(text)
        if not m or len(rs.findall(text)) != 1:
            sys.exit("STOP: final solve not found once")
        text = (text[:m.end()] + f"{nl}! start step done: the run itself{nl}deltim,{dt!r},{dt!r},{dt!r}{nl}"
                f"time,{t_end + s!r}{nl}solve" + text[m.end():])
    head = (f"! make_living_deck.py: {'v2, ' if a.v2 else ''}{'uniform, ' if a.uniform else ''}"
            f"layout {a.layout}, antibiotic {a.ab}{f' from T* {a.t_ab:g}' if a.ab else ''}, "
            f"front r = {r_front:g} mm/T*, L = 1 {a.l1}, old consumption {a.old_cons}, g = {a.cons:g}, "
            f"g_a {scal['CONSUMPTION21']}, phi_min {props[29]:g}, start step {a.start_step:g}, "
            f"outer dt {a.deltim if a.deltim is not None else 'base'}, USSFin dt {a.dt_ussfin:g}, "
            f"end T* {a.tend if a.tend is not None else 'base'}, "
            f"{f'dc/dt with d = {a.transient_c:g}' if a.transient_c else 'quasi-static c'}; BIOFILM1 {len(bio1)} elements, "
            f"NUTRIENT2 ({'all outer faces' if a.v2 else 'y = +1 layer'}) {len(top)} elements{nl}")
    post = "post_fields_living.mac" if a.tend is None or a.tend == 1 else f"post_fields_living_t{a.tend:g}.mac"
    text = head + text.rstrip() + f"{nl}FINISH{nl}*USE,{post}{nl}FINISH{nl}"
    Path(a.out).write_text(text, encoding="latin-1", newline="")
    print(head.strip())


if __name__ == "__main__":
    main()

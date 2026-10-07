"""make_week_decks.py -- the unattended run set of 5 Oct 2026 (one week of ANSYS
while nobody is at IKMHIWI03): writes the decks into F:\\biofilm_upf_wired with
make_wired_deck.py and the run list (deck:case:minutes per line) for run_chain.ps1.

    python ansys_usermat/apdl/make_week_decks.py [--workdir F:\\biofilm_upf_wired]
        [--list _week_1005_runs.txt] [--check]

Every run uses beta = 0.02 mm^2/T* unless named otherwise, and the two-species
runs use the Klempt 2024 stiffness (E = 10 Pa, nu = 0.49), so their stresses can
be put next to the one-species runs. The set covers the open choices, so the
results are there whichever way they are decided:
  - one species at beta = 0.02 for ch4: the partner's growth variable,
    nu = 0.499, the example stiffness (cloud review, CLOUD_REVIEW_1005_ANSWERS.md);
  - two species, 8^3: consumption x clock s, T* = 2, phi_cap, start value,
    start value of late points (prop(35)), beta;
  - two species, 16^3: the consumption series and the main sensitivities;
  - one species: beta 0.005 ... 0.1 on 8^3 / 16^3 / 24^3, T* = 2, dt on 16^3;
  - then the idle days: two species on 24^3, one species to T* = 5, more 16^3.
Consumption 1 (with MY_DIFF1 = 1) is Klempt 2024 Table 2 (d = 1e10 um^2/T*,
g = 1e8 /T*) scaled to the 2 mm cube like beta (lengths x 100 from the paper's
20 um cube): the quasi-static Eq. 35 depends on g/d only, penetration length
sqrt(d/g) = half the cube in both. Consumption 0/2/4/6 are example inputs.
The phi update is explicit, so dt beta / h^2 is kept below about 0.1
(3-D limit 1/6): beta = 0.1 on 16^3 and beta = 0.05 on 24^3 use dt 0.0125.
--check only prints the list (no deck is written).
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
E10 = ["YOUNG_BIO=1e-5", "YOUNG_VOID=-1e-3", "POISSON_BIO=0.49"]
CASE = {"c3": "2sp_case3", "c6": "2sp_case6"}
BASE2 = {(8, "c3"): "ds_c3_nut_g6_b002", (8, "c6"): "ds_c6_nut_g6_b002",
         (16, "c3"): "ds16_c3_nut_g6_b002", (16, "c6"): "ds16_c6_nut_g6_b002",
         (24, "c6"): "ds24_c6_nut_g6_b002"}         # refine_deck.py --n 24 of the 8^3 deck
BASE1 = {8: "ds8_beta002_dt4", 16: "ds16_beta002_dt4", 24: "ds24_beta002_dt4"}


def num(x):
    return f"{x:g}".replace(".", "").replace("-", "m")


def runs():
    """(name, base, props, sets, deltim, time, case, minutes) in run order."""
    out = []

    def two(n, c, cons=6, s=0.15, cap=None, chi0=None, late=None, cref=None, beta=None, T=None,
            gw=None):
        name = f"w{n}_{c}_g{cons}_s{num(s)}"
        props = {31: s}
        if cap is not None:
            props[32] = cap; name += f"_cap{num(cap)}"
        if chi0 is not None:
            props[30] = chi0; name += f"_chi{num(chi0)}"
        if late is not None:
            props[35] = late; name += f"_late{num(late)}"
        if cref is not None:
            props[33] = cref; name += f"_cref{num(cref)}"
        if gw is not None:
            props[36] = gw; name += f"_gw{num(round(gw, 3))}"
        sets = E10 + [f"CONSUMPTION11={cons}"]
        if beta is not None:
            sets.append(f"MY_BETA1={beta}"); name += f"_b{num(beta)}"
        if T is not None:
            name += f"_T{num(T)}"
        mins = {8: 60, 16: 180, 24: 900}[n] * (2 if T else 1) * (2 if s >= 0.5 else 1)
        out.append((name, BASE2[(n, c)], props, sets, None, T, CASE[c], mins))

    def one(n, beta, dt=None, T=None, tag="", props=None, sets=()):
        name = (f"w{n}_1sp{tag}_b{num(beta)}" + (f"_dt{num(dt)}" if dt else "")
                + (f"_T{num(T)}" if T else ""))
        mins = {8: 20, 16: 90, 24: 400}[n] * (2 if dt else 1) * (int(T + 0.999) if T else 1)
        out.append((name, BASE1[n], props or {}, [f"MY_BETA1={beta}", *sets], dt, T, "", mins))

    # 0. one species for ch4 at beta = 0.02 (cloud review, 5 Oct): A the partner's
    #    own growth variable (prop(28) = 3), B nu = 0.499, C the framework's
    #    example stiffness (1000 MPa, nu = 0.3, linear blend)
    for n in (8, 16):
        one(n, 0.02, tag="partner", props={28: 3})
        one(n, 0.02, tag="nu0499", sets=("POISSON_BIO=0.499",))
        one(n, 0.02, tag="Eex", sets=("YOUNG_BIO=1000", "YOUNG_VOID=1", "POISSON_BIO=0.3"))
    # 0b. species-weighted growth law, prop(36) = d (cloud, 5 Oct evening;
    #     needs the rebuild of tonight, two_way_growth_rate.py): case 6 with
    #     consumption 1 (Table 2) and 6, case 3 with 6; d = 1/3, and 0.5 once
    for n in (8, 16):
        two(n, "c6", cons=1, gw=1 / 3); two(n, "c6", cons=6, gw=1 / 3)
        two(n, "c3", cons=6, gw=1 / 3); two(n, "c6", cons=0, gw=1 / 3)
    two(8, "c6", cons=6, gw=0.5)
    # 1. two species, 8^3 (about 3-15 min each); s = 1.0 only with consumption 4/6
    for c in ("c6", "c3"):
        for cons in (1, 0, 2, 4, 6):
            for s in (0.05, 0.1, 0.15, 0.25, 0.5, 1.0):
                if s == 1.0 and cons < 4:
                    continue
                two(8, c, cons=cons, s=s)
        for s in (0.05, 0.15, 0.5):
            for cons in (1, 0, 6):
                two(8, c, cons=cons, s=s, T=2.0)
        for cap in (0.85, 0.95):
            two(8, c, cap=cap)
        for chi0 in (0.2, 0.8):
            two(8, c, chi0=chi0)
        for late in (0.2, 0.5, 0.8):
            two(8, c, late=late)
        for beta in (0.01, 0.05):
            two(8, c, beta=beta)
    # 2. one species, 8^3 and 16^3 (beta, T* = 2, dt)
    for beta in (0.005, 0.01, 0.05):
        one(8, beta); one(16, beta)
    one(8, 0.1); one(16, 0.1, dt=0.0125)
    one(8, 0.02, T=2.0); one(16, 0.02, T=2.0)
    one(16, 0.02, dt=0.0125)
    # 3. two species, 16^3 (about 30-60 min each)
    for c in ("c6", "c3"):
        for cons in (1, 6, 0, 4, 2):
            two(16, c, cons=cons)
    two(16, "c6", s=0.05); two(16, "c6", s=0.5)
    two(16, "c6", late=0.2); two(16, "c6", late=0.8)
    two(16, "c6", beta=0.01); two(16, "c6", beta=0.05)
    two(16, "c6", T=2.0)
    # 4. one species, 24^3 (about 2-4 h each)
    one(24, 0.01); one(24, 0.05, dt=0.0125)
    # 5. the idle days (cloud review): two species on 24^3 (third mesh for the
    #    composition; memory not checked, a FATAL stop is not retried), one
    #    species to T* = 5, then a denser 16^3 s x consumption grid
    two(24, "c6", cons=1); two(24, "c6", cons=6); two(24, "c6", cons=0)
    one(8, 0.02, T=5.0); one(16, 0.02, T=5.0)
    for s in (0.1, 0.25):
        two(16, "c6", s=s)
    for s in (0.05, 0.5):
        two(16, "c6", cons=0, s=s); two(16, "c3", s=s)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=r"F:\biofilm_upf_wired")
    ap.add_argument("--list", default="_week_1005_runs.txt")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    W = Path(a.workdir)
    rs = runs()
    names = [r[0] for r in rs]
    assert len(names) == len(set(names)), "duplicate run names"
    lines = []
    for name, base, props, sets, dt, T, case, mins in rs:
        lines.append(f"{name}:{case}:{mins}")
        if a.check:
            continue
        b = (W / f"{base}.dat").read_text(errors="replace")
        pe = re.search(r"^\*USE,post_elem_stress\.mac,(\d+)", b, re.M)
        cmd = [sys.executable, str(HERE / "make_wired_deck.py"), str(W / f"{base}.dat"), str(W / f"{name}.dat"),
               "--post", "both", "--post-elem", pe.group(1) if pe else "220"]
        if props:
            cmd += ["--props", ",".join(f"{k}={v:g}" for k, v in props.items())]
        for s in sets:
            cmd += ["--set", s]
        if dt:
            cmd += ["--deltim", f"{dt:g}"]
        if T:
            cmd += ["--time", f"{T:g}"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"{name}: make_wired_deck failed\n{r.stdout}{r.stderr}")
    print("\n".join(lines))
    hours = sum(int(l.rsplit(":", 1)[1]) for l in lines) / 60
    print(f"# {len(lines)} runs, timeouts add up to {hours:.0f} h", file=sys.stderr)
    if not a.check:
        (W / a.list).write_text("\n".join(lines) + "\n")
        print(f"# decks and {W / a.list} written", file=sys.stderr)


if __name__ == "__main__":
    main()

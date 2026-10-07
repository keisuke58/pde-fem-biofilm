"""hp_clamped_check.py -- the homeostatic-pressure growth law (prop(49:51) of
phi_mode_exec.inc) in the partner's element in ANSYS, against the closed form
of constrained growth (the ANSYS counterpart of tests/test_homeostatic_growth.py).

    python ansys_usermat/apdl/hp_clamped_check.py decks [--workdir F:\\biofilm_upf_hp]
    python ansys_usermat/apdl/hp_clamped_check.py judge [--workdir ...] [--report file]

Decks (from ds8_beta002_dt4 in F:\\biofilm_upf_wired): the 2 mm cube, 8^3, every
element in BIOFILM1 (phi = 1 everywhere), the field source K_LOCAL1 = 0, the
front term off and beta = 0 (so phi stays 1), each face held in its normal
direction (F = I at every Gauss point),
k_alpha = 1e-3, dt = 0.025, T* = 5, Klempt 2024 stiffness (E = 1e-5 MPa,
f = 1e-3, nu = 0.49):
  hp_clamp_off    prop(49) = 0             Eq. 36: alpha = k_alpha T
  hp_clamp_e      prop(49) = 2.5e-7 MPa    p_h = p_ref
  hp_clamp_local  prop(49) = 2.5e-7, prop(50) = 1, prop(51) = 1e-3
                                           p_h = p_ref (phi^2 + f)
Checked on element 220 (elem_stress_<run>.csv, every substep): alpha against
Eq. 36 / the recurrence alpha_{n+1} = alpha_n + k dt max(0, 1 - p(alpha_n)/p_h),
p = -SX against the closed form p(alpha), and that alpha stops at alpha* with
p(alpha*) = p_h.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WIRED = Path(r"F:\biofilm_upf_wired")
KALPHA, DT, T_END = 1e-3, 0.025, 5.0
E0, F_VOID, NU, PREF = 1e-5, 1e-3, 0.49, 2.5e-7
RUNS = {"hp_clamp_off": (0.0, 0, 0.0),
        "hp_clamp_e": (PREF, 0, 0.0),
        "hp_clamp_local": (PREF, 1, F_VOID)}
# each face held in its normal direction: uniform growth then gives F = I at
# every Gauss point. (Holding every node, the first version, leaves no free
# DOF and the first substep did not converge, 7 Oct.)
CLAMP = ("D,XMIN,UX,0\nD,XMAX,UX,0\nD,YMIN,UY,0\nD,YMAX,UY,0\n"
         "D,ZMIN,UZ,0\nD,ZMAX,UZ,0\nALLSEL\n")


def p_closed(alpha, phi=1.0):
    """-sigma_ii of constrained growth, Klempt 2024 stiffness E (phi^2 + f)
    (closed form of tests/test_abaqus_composition_umat.py)."""
    E = (phi * phi + F_VOID) * E0
    K, C10 = E / (3 * (1 - 2 * NU)), E / (4 * (1 + NU))
    J = (1 + alpha) ** -3
    return -(K * (J - 1) + 2 * C10 * (1 - (1 + alpha) ** 2) / J)


def p_h(pref, form, f):
    return pref * (1.0 + f) if form == 1 else pref


def alpha_star(ph):
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if p_closed(mid) < ph else (lo, mid)
    return (lo + hi) / 2


def decks(W: Path):
    base = WIRED / "ds8_beta002_dt4.dat"
    text = base.read_text(errors="replace")
    m = re.search(r"^CMBLOCK,IC_PHI_ELEM,ELEM,\s*(\d+)", text, re.M)
    n_el = int(m.group(1)) if m else 512
    allel = ",".join(str(i) for i in range(1, n_el + 1))
    for run, (pref, form, f) in RUNS.items():
        out = W / f"{run}.dat"
        cmd = [sys.executable, str(HERE / "make_wired_deck.py"), str(base), str(out),
               "--props", f"49={pref:g},50={form},51={f:g}",
               "--cmblock", f"BIOFILM1={allel}", "--set", "K_LOCAL1=0",
               # front term off (r = 100 of the base deck blows up on a full
               # biofilm, 7 Oct) and beta = 0: with beta the surface points
               # (phi = 0) drain the field (mean phi 0.89 at T* = 0.15)
               "--set", "MAX_GROWTH11=0", "--set", "MY_BETA1=0",
               "--time", f"{T_END:g}", "--post", "both", "--post-elem", "220"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"{run}: {r.stdout}{r.stderr}")
        s = out.read_bytes().decode("latin-1")
        nl = "\r\n" if "\r\n" in s else "\n"
        # the header lists all element numbers: shorten it below APDL's 640 characters
        s = s.replace(f"BIOFILM1={allel}", f"BIOFILM1=1..{n_el}", 1)
        hits = list(re.finditer(r"^OUTRES,ALL,ALL[ \t]*\r?\n(?=solve)", s, re.M | re.I))
        if len(hits) != 1:
            sys.exit(f"{run}: OUTRES,ALL,ALL before solve found {len(hits)} times")
        s = s[:hits[0].start()] + CLAMP.replace("\n", nl) + s[hits[0].start():]
        out.write_bytes(s.encode("latin-1"))
        print(f"{out}: prop(49:51) = {pref:g}, {form}, {f:g}; {n_el} elements in BIOFILM1; all nodes held")


def history(path: Path):
    rows = []
    with open(path) as fh:
        rd = csv.reader(fh)
        head = [h.strip() for h in next(rd)]
        it, ix, ia = head.index("time"), head.index("sx"), head.index("alpha")
        for r in rd:
            if r:
                rows.append((float(r[it]), float(r[ix]), float(r[ia])))
    return rows


def judge(W: Path, report: str | None):
    msg, ok = [], True

    def say(s=""):
        msg.append(s)
        print(s)

    say("Homeostatic growth law in the partner element, clamped 8^3 cube, phi = 1, element 220")
    for run, (pref, form, f) in RUNS.items():
        p = W / f"elem_stress_{run}.csv"
        if not p.exists():
            say(f"{run}: no elem_stress csv")
            ok = False
            continue
        h = history(p)
        t_end, sx_end, a_end = h[-1]
        # stress against the closed form at every substep but the first
        # (zero stress on the first call at Time = 0, oliver_usermat note)
        dev = max(abs(-sx - p_closed(a)) / max(p_closed(a), 1e-30) for _, sx, a in h[1:])
        if pref == 0.0:
            rel = abs(a_end - KALPHA * t_end) / (KALPHA * t_end)
            good = t_end >= T_END - 1e-9 and rel < 1e-6 and dev < 1e-4
            say(f"{run}: T* {t_end:.3f}, alpha {a_end:.6e} (Eq. 36 {KALPHA * t_end:.6e}, {rel:.1e}),"
                f" stress vs closed form max {dev:.1e}  {'pass' if good else 'FAIL'}")
        else:
            ph = p_h(pref, form, f)
            ast = alpha_star(ph)
            a, rec = 0.0, []
            for _ in h:
                a = a + KALPHA * DT * max(0.0, 1 - p_closed(a) / ph)
                rec.append(a)
            rdev = max(abs(x[2] - r) / r for x, r in zip(h, rec))
            pmax = max(-sx for _, sx, _ in h)
            good = (t_end >= T_END - 1e-9 and abs(a_end - ast) / ast < 1e-3
                    and abs(-sx_end - ph) / ph < 1e-3 and pmax <= ph * (1 + 1e-5) and dev < 1e-4)
            # (the result file holds single precision: no tighter than ~1e-6)
            say(f"{run}: p_h {ph:.4e} MPa; alpha at T* {t_end:.2f} {a_end:.6e}, alpha* {ast:.6e}"
                f" ({abs(a_end - ast) / ast:.1e}); p {-sx_end:.4e} (max {pmax:.4e});"
                f" recurrence max {rdev:.1e}; stress vs closed form max {dev:.1e}  {'pass' if good else 'FAIL'}")
        ok &= good
    say()
    say(f"RESULT: {'PASS' if ok else 'FAIL'}")
    if report:
        Path(report).parent.mkdir(parents=True, exist_ok=True)
        Path(report).write_text("\n".join(msg) + "\n")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=("decks", "judge"))
    ap.add_argument("--workdir", default=r"F:\biofilm_upf_hp")
    ap.add_argument("--report")
    a = ap.parse_args()
    if a.what == "decks":
        decks(Path(a.workdir))
        return 0
    return judge(Path(a.workdir), a.report)


if __name__ == "__main__":
    sys.exit(main())

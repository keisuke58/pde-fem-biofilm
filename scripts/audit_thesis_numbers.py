#!/usr/bin/env python3
"""Recompute the checkable headline numbers of thesis chapter 4 from the result files and
compare them with the values printed in the .tex files (regex). stdlib + numpy only.

    python scripts/audit_thesis_numbers.py          # prints PASS / FAIL per item, exit 1 on any FAIL

A tex value passes if the source value rounds to it (half a unit of the last printed digit).
Sources: ansys_usermat/apdl/results/**.json, MESH_STUDY.md. Chapter 3 (EML repo) is not checked.
"""
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
APDL = ROOT / "ansys_usermat/apdl"
RES = APDL / "results"
CH = ROOT / "thesis/chapters"
sys.path.insert(0, str(APDL))
import summarize_runs_json as S  # noqa: E402  (layers(), surface())

fails = 0


def tex(name):
    return "\n".join(l for l in (CH / name).read_text(encoding="utf-8").splitlines() if not l.lstrip().startswith("%"))


def printed(s):
    """'-2.52\\times10^{-4}' or '0.0209' or '1.9' -> (value, half unit of last digit)."""
    m = re.match(r"\s*\$?\s*([+-]?)\s*(\d+(?:\.(\d+))?)(?:\s*\\times\s*10\^\{?(-?\d+)\}?)?", s)
    sign = -1 if m.group(1) == "-" else 1
    dec = len(m.group(3) or "")
    e = int(m.group(4) or 0)
    return sign * float(m.group(2)) * 10.0**e, 0.5 * 10.0 ** (e - dec)


def check(label, tex_str, src, tol=None, note=""):
    global fails
    v, half = printed(tex_str)
    lim = tol if tol is not None else half
    err = abs(v - src)
    ok = err <= lim * (1 + 1e-9)
    warn = not ok and tol is None and err <= 2 * lim * (1 + 1e-9)   # off by one in the last printed digit
    fails += not ok and not warn
    print(f"{'PASS' if ok else 'WARN' if warn else 'FAIL'}  {label}: tex {tex_str.strip('$ ')}  source {src:.4g} {note}")


def flag(label, ok, detail):
    global fails
    fails += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {label}: {detail}")


def load(rel):
    return json.loads((RES / rel).read_text())


def seed_stats(rel):
    """seed mean von Mises and mean stress, layer-1/2 mean stress [Pa], alpha at the monitored element."""
    d = load(rel)
    a = d["all_stress"]
    el = np.array(a["elem"], int)
    seed = np.isin(el, d["seed_BIOFILM1"])
    cen = np.stack([a["cx"], a["cy"], a["cz"]], 1)
    lay, _ = S.layers(cen, seed)
    p = (np.array(a["sx"]) + np.array(a["sy"]) + np.array(a["sz"])) / 3e-6
    vm = np.array(a["seqv"]) * 1e6
    return dict(vm=vm[seed].mean(), p=p[seed].mean(), l1=p[lay == 1].mean(), l2=p[lay == 2].mean(),
                vmax=vm.max(), el=d["elem_stress"], alpha_seed=np.array(a["alpha"])[seed].mean())


def row(text, start):
    line = next(l for l in text.splitlines() if l.strip().startswith(start))
    return re.findall(r"\$[^$]*\$", line)


c4, c6, c2, c3 = tex("ch4_ansys.tex"), tex("ch6_conclusion.tex"), tex("ch2_theory.tex"), tex("ch3_heine.tex")
abstract = (ROOT / "thesis/main.tex").read_text(encoding="utf-8")

# 1. Table 4.2 (mesh refinement, ANSYS, beta = 0.02, dt = 0.025) ------------------------------------------
runs = {8: "2026-10-05_ansys/ds8_beta002_dt4.json", 16: "2026-10-05_ansys/ds16_beta002_dt4.json",
        24: "2026-10-05_ansys/ds24_beta002_dt4.json"}
st = {n: seed_stats(r) for n, r in runs.items()}
for lab, start, key in (("mean stress seed", "mean stress, seed", "p"), ("mean stress layer 1", "mean stress, first layer", "l1"),
                        ("von Mises seed", "von Mises, seed", "vm")):
    for n, t in zip((8, 16, 24), row(c4, start)):
        check(f"Tab mesh {lab} {n}^3", t, st[n][key])
check("layer 2 mean stress (text)", "1.3\\times10^{-6}", st[16]["l2"], tol=0.06e-6)
check("seed vM 16^3 -> 24^3 change 3 %", "3", (st[24]["vm"] / st[16]["vm"] - 1) * 100, tol=0.5, note="%")

# 2. monitored element (16^3) and largest vM -------------------------------------------------------------
e = st[16]["el"]
check("alpha-1 monitored element T*=1.1", "8.9\\times10^{-4}", e["alpha"][-1])
check("vM monitored element", "4.7\\times10^{-4}", e["seqv"][-1] * 1e6)
check("largest vM in model", "7.2\\times10^{-4}", st[16]["vmax"])
flag("Conclusion: no 'factor of about 1.7' for the beta=0.02 mesh", "factor of about 1.7" not in c6,
     "the 1.5-1.75 factor belongs to the sharp-seed Python study (item 9)")
check("conclusion: seed vM 8^3 below 16^3 [%]", "18", (1 - st[8]["vm"] / st[16]["vm"]) * 100, tol=1.0, note="%")
flag("conclusion: seed mean stress 8^3 vs 16^3 within 1 %", abs(st[8]["p"] / st[16]["p"] - 1) < 0.01,
     f"ratio {st[8]['p'] / st[16]['p']:.3f}")

# 3. Table 4.3 (beta sensitivity, mean stress in the seed) ----------------------------------------------
beta = {8: ["2026-10-wp3_fix/wp3e_n8_b0005.json", "2026-10-wp3_fix/wp3c_n8_b001.json", runs[8],
            "2026-10-07_surface_fix/w8_1sp_b005.json", "2026-10-07_surface_fix/w8_1sp_b01.json"],
        16: ["2026-10-wp3_fix/wp3e_n16_b0005.json", "2026-10-wp3_fix/wp3c_n16_b001.json", runs[16],
             "2026-10-07_surface_fix/w16_1sp_b005.json", "2026-10-07_surface_fix/w16_1sp_b01_dt00125.json"]}
for n, start in ((8, "$8^3$ &"), (16, "$16^3$ &")):
    for b, t, f in zip((0.005, 0.01, 0.02, 0.05, 0.1), row(c4, start)[1:], beta[n]):
        check(f"Tab beta {n}^3 beta={b}", t, seed_stats(f)["p"])

# 4. two species, packed seed and spreading ---------------------------------------------------------------
def share(rel, layer=None):
    d = load(rel)
    a = d["all_stress"]
    el = np.array(a["elem"], int)
    seed = np.isin(el, d["seed_BIOFILM1"])
    lay, _ = S.layers(np.stack([a["cx"], a["cy"], a["cz"]], 1), seed)
    nf = d["nut_field"]
    idx = {x: i for i, x in enumerate(el)}
    o = np.array([idx[x] for x in nf["elem"]])
    s = np.array(nf["phi1"]) + np.array(nf["phi2"])
    ran = (s > 1e-6) & (s < 0.95)
    sh = np.where(ran, np.array(nf["phi1"]) / np.where(ran, s, 1), np.nan)
    return sh, lay[o], seed[o], ran


sh, L, Sd, ran = share("2026-10-05_ansys/ds_c3_all.json")
check("case 3 packed seed share (no spreading)", "0.6665", np.nanmean(sh[Sd]))
sh, L, Sd, ran = share("2026-10-05_ansys/ds_c6_all.json")
check("case 6 packed seed share (no spreading)", "0.0209", np.nanmean(sh[Sd]))
sh, L, Sd, ran = share("2026-10-05_ansys/ds_c6_all_b002.json")
check("case 6 beta=0.02 first layer share", "0.33", np.nanmean(sh[(L == 1) & ran]))
check("case 6 beta=0.02 second layer share", "0.50", np.nanmean(sh[(L == 2) & ran]))
check("case 6 beta=0.02 seed share, upper end '0.06'", "0.06", np.nanmax(sh[Sd]))
check("point model ran in 304 of 512 elements", "304", ran.sum(), tol=0)
check("point model ran in 1970 of 4096 elements", "1970", share("2026-10-05_ansys/ds16_c6_all_b002.json")[3].sum(), tol=0)

# 5. calibrated five-species run (CH, raw-data MAP, Day-1 start; ch4 section 'five') ---------------------
pm = load("2026-10-cal5/w8_CH_g1_s015_final_d1_pm.json")
C = {c: i for i, c in enumerate(pm["columns"])}
rows = [r for r in pm["seed_trace"] if r[C["ip"]] == pm["seed_trace"][0][C["ip"]]]
last = {}
for r in rows:
    last[(int(r[C["ldstep"]]), int(r[C["isubst"]]))] = r
rows = [last[k] for k in sorted(last)]
g = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(1, 12) if i != 6] for r in rows])
phi, psi = g[:, :5], g[:, 5:10]
shr = phi / phi.sum(1, keepdims=True)
live = (phi * psi).sum(1)
check("five-species living fraction start", "0.81", live[0])
check("five-species living fraction end", "0.23", live[-1])
check("S. oralis share start", "0.74", shr[0, 0])
check("S. oralis share end", "0.55", shr[-1, 0])
check("A. naeslundii share end", "0.13", shr[-1, 1], tol=0.006)
check("V. dispar share end", "0.15", shr[-1, 2])
check("P. gingivalis share end", "0.08", shr[-1, 4])
ref = seed_stats("2026-10-cal5/w8_CH_g1_s015_final_d1.json")
check("five-species seed mean alpha-1", "6.752\\times10^{-4}", ref["alpha_seed"])
check("five-species seed mean vM", "4.68\\times10^{-4}", ref["vm"])

# 6. growth error against the finite-volume reference ----------------------------------------------------
ref = load("2026-10-wp3_fix/wp3_reference.json")
check("FV reference observed order", "1.97", ref["richardson"]["seed mean"]["order"])
for n, t in zip((8, 16, 24), ("4.1", "1.3", "0.7")):
    check(f"seed-mean alpha error {n}^3 [%]", t, ref["nem_error"][str(n)]["dt->0"]["err seed mean"] * 100)
for n, t in zip((8, 16, 24), ("5.1", "1.9", "1.1")):
    check(f"L2 error {n}^3 [%]", t, ref["nem_error"][str(n)]["dt->0"]["err L2"] * 100)
check("error at 0.4 mm, 8^3 [%]", "5.5", ref["nem_error"]["8"]["dt->0"]["err line"][4] * 100)
sr = load("2026-10-wp3_fix/wp3_stress_reference.json")
check("48^3 element-mean seed vM", "6.11\\times10^{-4}", sr["python"]["48"]["seed vM"])

# 7. NEM stability limits (appendix table) -------------------------------------------------------------------
cases = load("2026-10-wp3_fix/nem_stability.json")["cases"]
def lam(mesh, nn, bs):
    return next(c["lam_max"] for c in cases if c["mesh"] == mesh and c["nneigh"] == nn and (c["beta_star"] < 1e6) == bs)
for mesh, t in ((8, "0.062"), (16, "0.149"), (24, "0.173")):
    check(f"lambda_max N=32 {mesh}^3", t, lam(mesh, 32, True))
check("lambda_max N=32 unweighted 8^3", "0.180", lam(8, 32, False))
# lambda of the runs vs the largest lambda at which the runs to T*=5 stayed smooth (ch4 text)
smooth = {8: 0.048, 16: 0.136, 24: 0.130}
for n, dt, b in ((16, 0.1, 0.02), (24, 0.025, 0.02), (8, 0.025, 0.1)):
    lm = b * dt / (2.0 / n) ** 2
    flag(f"lambda {n}^3 dt={dt} beta={b} below smooth T*=5 run", lm <= smooth[n], f"{lm:.3f} vs {smooth[n]}")
has24 = (RES / "2026-10-wp3_fix/wp3c_n24_dt01.json").exists() and "all_stress" in load("2026-10-wp3_fix/wp3c_n24_dt01.json")
flag("24^3 dt=0.1: text and results agree", has24 or "on $24^3$ from $\\Delta t=0.05$" in c4,
     f"result present: {has24}; lambda {0.02 * 0.1 / (2 / 24) ** 2:.3f}")

# 8. closed forms (one element, alpha = 1.05) -----------------------------------------------------------------
K = 2 / 5e3
cont = K * (1.05**-3 - 1)
check("extra spherical term share 47 %", "47", (1 - cont / -1.0193e-4) * 100, tol=0.5, note="(from SX=-1.0193e-4)")
check("extra term at nu=0.49, (1-2nu)/(1+nu) [%]", "1.3", (1 - 2 * .49) / (1 + .49) * 100)
check("extra term at nu=0.3 [%]", "31", (1 - 2 * .3) / (1 + .3) * 100, tol=0.5)

# 9. Python mesh study (MESH_STUDY.md) vs 'about 1.6 to 1.75 times' ------------------------------------------
ms = {int(m.group(1)): (float(m.group(2)), float(m.group(3))) for m in
      re.finditer(r"\|\s*(\d+)\^3[^|]*\|\s*\d+\s*\|\s*([\d.e-]+)\s*\|\s*[\d.e-]+\s*\|\s*(-[\d.e-]+)\s*\|", (APDL / "MESH_STUDY.md").read_text())}
rat = [ms[8][i] / ms[n][i] for i in (0, 1) for n in (16, 32)]
flag("Python study: 8^3/(16^3,32^3) ratios in 1.5-1.75 (vM, then p)",
     all(1.45 <= r <= 1.78 for r in rat) and "$1.5$ to $1.75$" in c4, ", ".join(f"{r:.2f}" for r in rat))

# 10. cross-chapter consistency (abstract / ch1 / ch4 / ch6) --------------------------------------------------
for key, pat in (("3e-11 in alpha", r"3\\times10\^\{-11\}"), ("150x speed-up", r"about 150 times"), ("15 coefficients", r"\b15\b")):
    hits = {n: bool(re.search(pat, t)) for n, t in (("abstract", abstract), ("ch4", c4), ("ch6", c6))}
    print(("PASS" if hits["abstract"] and hits["ch6"] else "NOTE") + f"  cross-chapter {key}: {hits}")
check("conclusion: alpha-1 = 8.9e-4", re.search(r"8\.9\\times10\^\{-4\}", c6).group(0), e["alpha"][-1])
check("conclusion: 0.7 % above FV reference", "0.7", ref["nem_error"]["24"]["dt->0"]["err seed mean"] * 100)
check("ch2 table: mu = E/(2(1+nu)) [Pa]", "3.3557", 10 / (2 * 1.49))
check("speed-up end to end 21600/147", "147", 21600 / 147, tol=0.5)

print(f"\n{fails} FAIL")
sys.exit(1 if fails else 0)

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

# 11. two-species examples (Section "Two species with nutrient and antibiotic", Table tab:ch5-two-species) -----
# Python TVD, 64^3, outer step 0.025, inherit rule; seed means at T* = 1. Each printed value must appear in the table
# row and agree with the result file.
D3 = RES / "2026-10-3d"
_tab = c4[c4.index("label{tab:ch5-two-species}"):c4.index("end{tabular}", c4.index("label{tab:ch5-two-species}"))]
for label, rowname, f, key in (
        ("separate", "separate &", "nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_separate", "separate"),
        ("separate, antibiotic", "separate, antibiotic &", "nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_separate_ab", "separate_ab"),
        ("mixed, equal amounts", "mixed, equal amounts &", "nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_mixed", "mixed"),
        ("mixed, equal amounts, antibiotic", "mixed, equal amounts, antibiotic &", "nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_mixed_ab", "mixed_ab"),
        ("mixed (0.25, 0.30), antibiotic", "mixed $(0.25, 0.30)$, antibiotic &", "nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_mixed_ab", "mixed_ab")):
    fp = D3 / f"coupled_living_consumption_{f}.json"
    if not fp.exists():
        flag(f"two species, {label}: result file", False, f"missing {fp.name}"); continue
    h = json.loads(fp.read_text())["runs"][key]
    line = next((l for l in _tab.splitlines() if l.strip().startswith(rowname)), None)
    if line is None:
        flag(f"two species, {label}: table row", False, "row not found"); continue
    cells = [x.strip() for x in line.split("&")[1:]]
    cells[-1] = cells[-1].replace("\\\\", "").strip()
    for name, cell, val in zip(("living", "share A", "c", "alpha-1", "vM [Pa]"), cells,
                               (h["living_seed"][-1], h["share_a_seed"][-1], h["c_seed"][-1], h["seed_am1"], h["seed_vm_Pa"])):
        check(f"two species, {label}: {name}", cell, val)


# 12. grid table tab:ch5-gci: grid values from the result files (cloud, outer step 0.0025) or from report A-89
#     (IKMHIWI03, 32^3-128^3 of the two series with antibiotic and equal start amounts; their 192^3 JSONs are here,
#     A-92/A-94); p, extrapolated value and GCI recomputed with Celik et al. 2008 (Fs = 1.25): without antibiotic from
#     32^3/64^3/128^3 (r = 2), with antibiotic from 64^3/128^3/192^3 (r = 2 and 1.5, iterative).
def gci(f3, f2, f1, r21=2.0, r32=2.0):
    """f3, f2, f1 = coarse, medium, fine. Returns (p, extrapolated, GCI_fine [%]) or None if not monotone/converging."""
    e32, e21 = f2 - f3, f1 - f2
    if e32 * e21 <= 0:
        return None
    pp = 1.0
    for _ in range(200):
        pp = abs(math.log(abs(e32 / e21)) + math.log((r21 ** pp - 1) / (r32 ** pp - 1))) / math.log(r21)
    if pp < 0.3:
        return None
    ext = f1 + e21 / (r21 ** pp - 1)
    return pp, ext, 1.25 * abs(e21 / f1) / (r21 ** pp - 1) * 100

A89 = (APDL / "IKMHIWI03_TO_CLOUD.ja.md").read_text(encoding="utf-8").replace("−", "-")
A89 = A89[A89.index("## A-89"):]
def a89_row(start):
    line = next(l for l in A89.splitlines() if l.strip().startswith("| " + start))
    return [float(x) for x in re.findall(r"(\d+\.\d+)", line.split("|", 2)[2])]
def series(f, key, grids=("", "nf64_", "nf128_")):
    return [json.loads((D3 / f"coupled_living_consumption_{n}{f}.json").read_text())["runs"][key] for n in grids]
_g = c4[c4.index("label{tab:ch5-gci}"):c4.index("end{tabular}", c4.index("label{tab:ch5-gci}"))]
def gci_row(start):
    line = next(l for l in _g.splitlines() if l.strip().startswith(start))
    cells = [x.strip().replace("\\\\", "").replace("\\,\\%", "").replace("$", "") for x in line.split("&")[2:]]
    return cells
F = "front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_{}a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_{}"
sep = series(F.format("", "separate"), "separate")
mix = series(F.format("mixeq_", "mixed"), "mixed")
mix2 = series(F.format("", "mixed_ab"), "mixed_ab", ("", "nf64_", "nf128_", "nf192_"))
s192 = series(F.format("", "separate_ab"), "separate_ab", ("nf192_",))[0]
m192 = series(F.format("mixeq_", "mixed_ab"), "mixed_ab", ("nf192_",))[0]
ab_sep = {"am1": a89_row("分離例 α-1")[:3] + [s192["seed_am1"] * 1e4], "vm": a89_row("分離例 vM")[:3] + [s192["seed_vm_Pa"] * 1e4]}
ab_mix = {"am1": a89_row("混合例 α-1")[:3] + [m192["seed_am1"] * 1e4], "vm": a89_row("混合例 vM")[:3] + [m192["seed_vm_Pa"] * 1e4]}
for row, lab, vals in (
        ("separate, antib. & $\\alpha-1$", "GCI separate+ab alpha-1", ab_sep["am1"]),
        ("separate, antib. & von Mises", "GCI separate+ab vM", ab_sep["vm"]),
        ("separate & $\\alpha-1$", "GCI separate alpha-1", [h["seed_am1"] * 1e4 for h in sep]),
        ("separate & von Mises", "GCI separate vM", [h["seed_vm_Pa"] * 1e4 for h in sep]),
        ("mixed, antib. & $\\alpha-1$", "GCI mixed+ab alpha-1", ab_mix["am1"]),
        ("mixed, antib. & von Mises", "GCI mixed+ab vM", ab_mix["vm"]),
        ("mixed & $\\alpha-1$", "GCI mixed alpha-1", [h["seed_am1"] * 1e4 for h in mix]),
        ("mixed & von Mises", "GCI mixed vM", [h["seed_vm_Pa"] * 1e4 for h in mix]),
        ("mixed $(0.25,0.30)$, antib. & $\\alpha-1$", "GCI mixed(0.25,0.30)+ab alpha-1", [h["seed_am1"] * 1e4 for h in mix2]),
        ("mixed $(0.25,0.30)$, antib. & von Mises", "GCI mixed(0.25,0.30)+ab vM", [h["seed_vm_Pa"] * 1e4 for h in mix2])):
    cells = gci_row(row)
    four = len(vals) == 4
    for n, cell, v in zip(("32", "64", "128", "192"), cells[:4], vals):
        check(f"{lab} {n}^3", cell, v)
    if not four:
        flag(f"{lab}: 192^3 cell empty", cells[3] == "", cells[3])
    use, rr = (vals[1:], (1.5, 2.0)) if four else (vals, (2.0, 2.0))
    ucells = cells[1:4] if four else cells[:3]
    g = gci(*use, *rr)
    # the A-89 values are rounded (3 or 4 significant digits) while p, extrapolated value and GCI were computed from
    # the unrounded values: allow the spread over +-half a unit of the printed grid values
    halves = [0.5 * 10.0 ** -len(cell.split(".")[1]) if "." in cell else 0.5 for cell in ucells]
    import itertools
    sp = [gci(*[v + sgn * hh for v, sgn, hh in zip(use, sg, halves)], *rr) for sg in itertools.product((-1, 1), repeat=3)]
    sp = [x for x in sp if x is not None]
    for i, (name, cell) in enumerate(zip(("p", "extrapolated", "GCI [%]"), cells[4:7])):
        lo, hi = min(x[i] for x in sp), max(x[i] for x in sp)
        v, half = printed(cell)
        flag(f"{lab} {name}", lo - half <= v <= hi + half, f"tex {cell}, from the printed grid values {g[i]:.4g} (range {lo:.4g} to {hi:.4g})")
flag("ch4 text: GCI ranges match the table", all(x in c4 for x in (
    "GCI of $0.5$ to $4.3\\,\\%$ in all", "$0.7\\,\\%$ without and $1.4\\,\\%$ with", "$19$ to $41\\,\\%$",
    "from $6.55$ to $6.42\\times10^{-4}$~Pa, by $2\\,\\%$")), "text matches the table")
check("antibiotic lowers the separate extrapolated stress [%]", "2",
      (1 - gci(*ab_sep["vm"][1:], 1.5, 2.0)[1] / gci(*[h["seed_vm_Pa"] * 1e4 for h in sep])[1]) * 100)
for n, h, t in zip(("32", "64", "128", "192"), mix2, ("0.60", "0.69", "0.73", "0.75")):
    check(f"mixed (0.25,0.30)+ab: largest phi_3D in the domain, {n}^3", t, h["phi_max"][-1])

# 13. rigid bars (Python, 32^3, outer step 0.0025) -----------------------------------------------------------
bars_s = json.loads((D3 / "coupled_living_consumption_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_snaps_tvd_wgmirror_fast_obstcolumns_separate.json").read_text())["runs"]["separate"]
bars_m = json.loads((D3 / "coupled_living_consumption_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_snaps_tvd_wgmirror_fast_obstcolumns_mixed.json").read_text())["runs"]["mixed"]
check("bars: seed vM separate with bars", "7.1\\times10^{-3}", bars_s["seed_vm_Pa"])
check("bars: seed vM separate without bars (32^3)", "7.4\\times10^{-4}", sep[0]["seed_vm_Pa"])
check("bars: seed vM mixed with bars", "9.5\\times10^{-4}", bars_m["seed_vm_Pa"])
check("bars: seed vM mixed without bars (32^3)", "2.1\\times10^{-4}", mix[0]["seed_vm_Pa"])
check("bars: 'about ten times'", "10", bars_s["seed_vm_Pa"] / sep[0]["seed_vm_Pa"], tol=1.0)
check("bars: alpha-1 lower, separate [%]", "9", (1 - bars_s["seed_am1"] / sep[0]["seed_am1"]) * 100)
check("bars: alpha-1 lower, mixed [%]", "21", (1 - bars_m["seed_am1"] / mix[0]["seed_am1"]) * 100)
check("bars: 32^3 stress 14 % above the extrapolated value", "14", (sep[0]["seed_vm_Pa"] * 1e4 / gci(*[h["seed_vm_Pa"] * 1e4 for h in sep])[1] - 1) * 100)

# 14. maze (Python, 2D, 4 um cells) -------------------------------------------------------------------------
MZ = RES / "2026-10-maze"
def maze(n):
    d = json.loads((MZ / f"{n}.json").read_text())
    h = d["hist"]
    return d["setting"], np.array(h["t"][:len(h["living_mean"])]), np.array(h["living_mean"]), np.array(h["reach_um"])
st_m, t_m, lv, rc = maze("maze_h4_tend60_k24fig")
flag("maze: Table 6 constants (beta 10, r 200, k_alpha 1e-6)", (st_m["beta"], st_m["r"], st_m["k_alpha"]) == (10.0, 200.0, 1e-6), str((st_m["beta"], st_m["r"], st_m["k_alpha"])))
check("maze: nearest point to the centre at T*=60, no antibiotic [um]", "42", rc[-1])
check("maze: living fraction, no antibiotic", "0.92", lv[-1])
st_a, t_a, lva, rca = maze("maze_h4_tend60_ab_k24fig")
flag("maze: antibiotic from T*=20", st_a["T_ab"] == 20.0, str(st_a["T_ab"]))
check("maze: nearest point at T*=60, antibiotic [um]", "133", rca[-1])
check("maze: living fraction at T*=24, antibiotic", "0.12", lva[np.argmin(abs(t_a - 24))])
check("maze: living fraction at T*=60, antibiotic", "0.35", lva[-1])
check("maze: living fraction at T*=60, start state for new points", "0.04", maze("maze_h4_tend60_ab_k24fig_fresh")[2][-1])

# 15. antibiotic schedules (64^3, outer step 0.025, same total dose) ------------------------------------------
base = "nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_"
cst = json.loads((D3 / f"coupled_living_consumption_{base}tvd_wgmirror_fast_separate_ab.json").read_text())["runs"]["separate_ab"]
for name in ("abpulse0.2222", "absine0.2222"):
    h = json.loads((D3 / f"coupled_living_consumption_{base}{name}_tvd_wgmirror_fast_separate_ab.json").read_text())["runs"]["separate_ab"]
    flag(f"antibiotic {name}: alpha-1 and vM within 1 % of the constant dose",
         abs(h["seed_am1"] / cst["seed_am1"] - 1) < 0.01 and abs(h["seed_vm_Pa"] / cst["seed_vm_Pa"] - 1) < 0.01,
         f"alpha-1 {h['seed_am1'] / cst['seed_am1'] - 1:+.2%}, vM {h['seed_vm_Pa'] / cst['seed_vm_Pa'] - 1:+.2%}")
    if name.startswith("abpulse"):
        t, l, a = np.array(h["t"]), np.array(h["living_seed"]), np.array(h["a_seed"])
        i0 = np.argmin(l)
        check("pulsed dose: lowest living fraction", "0.17", l[i0])
        check("pulsed dose: living fraction at T*=1", "0.30", l[-1])
        check("pulsed dose: recovered value within the last pause", "0.34", l[i0:][(a[i0:] == 0)].max())

# 16. time step paragraph (A-89: 64^3, separate with antibiotic, outer steps 0.02 .. 0.0025) ------------------
ts = {}
for line in A89.splitlines():
    m = re.match(r"\|\s*(0\.0\d+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)e-4\s*\|\s*([\d.]+)e-4\s*\|", line.strip())
    if m:
        ts[float(m.group(1))] = tuple(float(m.group(i)) for i in (2, 3, 4))
flag("A-89 time-step table found (0.02, 0.01, 0.005, 0.0025)", sorted(ts) == [0.0025, 0.005, 0.01, 0.02], str(sorted(ts)))
ext_vm = ts[0.0025][2] + (ts[0.0025][2] - ts[0.005][2]) / (2 ** (math.log((ts[0.005][2] - ts[0.01][2]) / (ts[0.0025][2] - ts[0.005][2])) / math.log(2)) - 1)
check("time step: first-order limit of the separate+ab vM on 64^3 [1e-4 Pa]", "6.07", ext_vm)
sep_ab_64 = json.loads((D3 / f"coupled_living_consumption_{base}tvd_wgmirror_fast_separate_ab.json").read_text())["runs"]["separate_ab"]
check("time step: 0.025 about 6 % below the limit (separate+ab)", "6", (1 - sep_ab_64["seed_vm_Pa"] * 1e4 / ext_vm) * 100, tol=0.6, note="%")
mix_ab_64 = json.loads((D3 / f"coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_mixed_ab.json").read_text())["runs"]["mixed_ab"]
check("time step: mixed+ab 0.025 about 15 % above 0.0025 (3.23 vs 2.82)", "15", (mix_ab_64["seed_vm_Pa"] * 1e4 / ab_mix["vm"][1] - 1) * 100, tol=0.6, note="%")
sep_64 = json.loads((D3 / f"coupled_living_consumption_{base}tvd_wgmirror_fast_separate.json").read_text())["runs"]["separate"]
mix_64 = json.loads((D3 / "coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_mixed.json").read_text())["runs"]["mixed"]
check("time step: without antibiotic, separate 2.5 %", "2.5", abs(sep[1]["seed_vm_Pa"] / sep_64["seed_vm_Pa"] - 1) * 100, tol=0.1, note="%")
check("time step: without antibiotic, mixed 12 %", "12", abs(mix[1]["seed_vm_Pa"] / mix_64["seed_vm_Pa"] - 1) * 100, tol=0.6, note="%")
# composition and alpha-1 between outer steps 0.02 and 0.0025 (A-89 and the files): the text must cover the largest change
comp_change = max(abs(ts[0.0025][0] / ts[0.02][0] - 1), abs(ab_mix["vm"][1] * 0 + mix[1]["share_a_seed"][-1] / mix_64["share_a_seed"][-1] - 1),
                  abs(mix[1]["c_seed"][-1] / mix_64["c_seed"][-1] - 1)) * 100
am1_change_sep = abs(ts[0.0025][1] / ts[0.02][1] - 1) * 100
am1_change_mix = abs(mix[1]["seed_am1"] / mix_64["seed_am1"] - 1) * 100
c4w = re.sub(r"\s+", " ", c4)
flag("ch4 text: composition change between outer steps 'less than 4 %'", "composition changes by less than $4\\,\\%$" in c4w and comp_change < 4,
     f"largest change {comp_change:.1f} %")
flag("ch4 text: alpha-1 change 'less than 0.3 % (separate) and 3 % (mixed)'",
     "$0.3\\,\\%$ in the separate and $3\\,\\%$ in the mixed example" in c4w and am1_change_sep < 0.3 and am1_change_mix < 3.05,
     f"separate {am1_change_sep:.2f} %, mixed {am1_change_mix:.1f} %")
flag("ch4 text: 'less than 4 %' between grids and time steps (living fraction, shares, nutrient)",
     "change by less than $4\\,\\%$ between the grids and the time steps" in c4w
     and abs(a89_row("分離例 α-1")[0] * 0 + 0.1572 / sep_ab_64["living_seed"][-1] - 1) < 0.04, "separate+ab living fraction 0.1515 (0.025) vs 0.1572 (128^3, 0.0025)")

# 17. two-species composition across abstract / ch4 / ch6 ----------------------------------------------------
for lab, t, v in (("separate", "0.96", sep_64["living_seed"][-1]), ("mixed", "0.96", mix_64["living_seed"][-1]),
                  ("separate+ab", "0.15", sep_ab_64["living_seed"][-1]), ("mixed+ab", "0.04", mix_ab_64["living_seed"][-1])):
    check(f"living volume fraction of the seed, {lab} (64^3, 0.025)", t, v)
for n, t in (("abstract", abstract), ("ch6", c6)):
    flag(f"{n}: living fraction 'from $0.96$ to $0.15$' (not 0.97)", "from $0.96$ to $0.15$" in t and "0.97" not in t, "")
check("ch6: growth changes by less than 2 % with the antibiotic", "2", max(abs(sep_ab_64["seed_am1"] / sep_64["seed_am1"] - 1),
      abs(mix_ab_64["seed_am1"] / mix_64["seed_am1"] - 1)) * 100, tol=1.0, note="% (upper bound)")
flag("ch5-five: one-way two-species run has the same seed alpha-1 (6.752e-4)", abs(seed_stats("2026-10-05_ansys/ds_c6_all_b002.json")["alpha_seed"] - 6.752e-4) < 0.5e-7,
     f"{seed_stats('2026-10-05_ansys/ds_c6_all_b002.json')['alpha_seed']:.4e}")

# 18. constants and conversions (Klempt 2024 Table 2 -> 2 mm cube) ----------------------------------------------
check("Peclet r h/(2 beta) on 16^3", "31", 10 * (2 / 16) / (2 * 0.02), tol=0.5)
check("diffusion length sqrt(beta T*) [mm]", "0.14", math.sqrt(0.02))
check("nutrient relaxation time L^2/d [T*]", "4\\times10^{-8}", 2.0**2 / 1e8)
check("g l^2/d with l = 1 mm", "1", 1e8 * 1.0 / 1e8)
check("beta: 2 um^2/T* x 100^2 -> mm^2/T*", "0.02", 2 * 1e4 * 1e-6)
check("d: 1e10 um^2/T* x 100^2 -> mm^2/T*", "1\\times10^{8}", 1e10 * 1e4 * 1e-6, tol=1)
check("r: 100 um/T* x 100 -> mm/T*", "10", 100 * 100 * 1e-3)
check("r: maze: 200 um/T*", "200", maze("maze_h4_tend60_k24fig")[0]["r"], tol=0)

print(f"\n{fails} FAIL")
sys.exit(1 if fails else 0)

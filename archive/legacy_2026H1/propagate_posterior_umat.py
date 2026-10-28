#!/usr/bin/env python3
"""
propagate_posterior_umat.py  --  Rigorous direct posterior propagation.

For each posterior phi sample (or the 4 MAP for validation):
  1. Solve the Klempt 2D PDE (run_one_condition) -> alpha_max, phi_max field.
  2. Build the mode-A tooth INP (E_gated from phi, alpha/phi_local from PDE).
  3. Run Abaqus with umat_klempt_voigt.f.
  4. Extract max von Mises (= the sigma used for the headline).
No surrogate. Replaces posterior_klempt_stress_ci.py's power-law for the CI.

Usage:
  python3 propagate_posterior_umat.py validate          # 4 MAP, check vs MAP_SIGMA
  python3 propagate_posterior_umat.py full   [n_par]     # all posterior samples
"""
import sys, os, json, subprocess, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "JAXFEM"))

import klempt_pde_multispecies as pde
import gen_tooth_klempt_umat_inp as g
pde.N_STEPS = 2000  # headline alpha used N_STEPS=2000 (code default 400 is stale)

OUT = HERE / "_propagate_posterior"
OUT.mkdir(exist_ok=True)
EXTRACT = HERE / "extract_mises.py"
UMAT_A = str(g.UMAT["A"])

MAP_SIGMA = {  # kPa, headline reference (tooth, mode A)
    "commensal_hobic": 13.72, "dysbiotic_hobic": 2.13,
    "commensal_static": 8.77, "dysbiotic_static": 13.63}

# reuse the (heavy) template parse once
_BLOCKS = g.parse_template(g.TMPL)
_DEPTH = np.array([g.depth_norm(n) for n in range(g.N_NODES)])


def make_inp(cond, phi_vec, tag):
    """PDE -> fields -> mode-A INP. Returns (inp_path, alpha_max, phi_max, k_eff)."""
    res = pde.run_one_condition(cond, np.asarray(phi_vec), pde.NX, pde.NY, pde.N_STEPS)
    a2d, p2d = res["alpha_final"], res["phi_final"]
    alpha_nodes = np.array([g.alpha_at_depth(d, a2d) for d in _DEPTH])
    phi_nodes = np.array([g.phi_local_at_depth(d, p2d) for d in _DEPTH])
    inp = OUT / ("p_%s.inp" % tag)
    g.write_mode_A(inp, _BLOCKS, alpha_nodes, phi_nodes, np.asarray(phi_vec), cond)
    return inp, res["alpha_max"], res["phi_max"], res["k_alpha_eff"]


def launch(inp, tag):
    job = "p_" + tag
    return subprocess.Popen(
        ["abaqus", "job=" + job, "input=" + inp.name, "user=" + UMAT_A,
         "cpus=1", "interactive", "ask_delete=OFF"],
        cwd=str(OUT), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def extract(tag):
    odb = OUT / ("p_%s.odb" % tag)
    if not odb.exists():
        return None
    r = subprocess.run(["abaqus", "python", str(EXTRACT), str(odb)],
                       cwd=str(OUT), capture_output=True, text=True)
    for ln in r.stdout.splitlines():
        if ln.startswith("MISES_MAX"):
            return float(ln.split()[1])  # MPa
    return None


def run_batch(items, n_par):
    """items: list of (cond, phi_vec, tag). Returns {tag: sigma_MPa}."""
    # 1) build all INPs (PDE solves -- cheap, serial)
    meta = {}
    for cond, phi, tag in items:
        inp, amax, pmax, keff = make_inp(cond, phi, tag)
        meta[tag] = {"cond": cond, "alpha_max": amax, "k_eff": keff, "inp": inp}
    # 2) run Abaqus with bounded concurrency
    sigma = {}
    running = {}  # tag -> Popen
    queue = list(meta.items())
    done = 0
    total = len(queue)
    while queue or running:
        while queue and len(running) < n_par:
            tag, m = queue.pop(0)
            running[tag] = launch(m["inp"], tag)
        time.sleep(2)
        for tag in list(running):
            if running[tag].poll() is not None:
                running.pop(tag)
                sigma[tag] = extract(tag)
                done += 1
                s = sigma[tag]
                print("  [%d/%d] %-26s sigma=%s kPa  (alpha_max=%.3f k_eff=%.3f)"
                      % (done, total, tag, ("%.2f" % (s * 1e3)) if s else "FAIL",
                         meta[tag]["alpha_max"], meta[tag]["k_eff"]), flush=True)
    return sigma, meta


def validate():
    tm = g.load_tmcmc()
    items = [(c, tm[c], "MAP_" + c) for c in MAP_SIGMA]
    print("=== VALIDATION: 4 MAP phi (must reproduce MAP_SIGMA) ===")
    sigma, _ = run_batch(items, n_par=4)
    print("\n  cond                       pipeline   MAP_ref   err%")
    ok = True
    for c in MAP_SIGMA:
        s = sigma.get("MAP_" + c)
        if s is None:
            print("  %-26s   FAIL" % c); ok = False; continue
        err = (s * 1e3 / MAP_SIGMA[c] - 1) * 100
        flag = "" if abs(err) < 5 else "  <-- MISMATCH"
        if abs(err) >= 5:
            ok = False
        print("  %-26s  %7.2f   %7.2f  %+5.1f%%%s" % (c, s * 1e3, MAP_SIGMA[c], err, flag))
    print("\nVALIDATION", "PASS" if ok else "FAIL")
    return ok


def full(n_par):
    sys.path.insert(0, str(HERE / "JAXFEM"))
    from posterior_klempt_stress_ci import load_phi_samples, COND_TO_CI0D
    items = []
    for cond in MAP_SIGMA:
        phis = load_phi_samples(cond)
        if phis is None:
            continue
        for i, phi in enumerate(phis):
            items.append((cond, phi, "%s_%03d" % (cond, i)))
    print("=== FULL PROPAGATION: %d samples, %d parallel ===" % (len(items), n_par))
    sigma, meta = run_batch(items, n_par)
    # aggregate per condition
    out = {}
    for cond in MAP_SIGMA:
        vals = [sigma[t] for (c, _, t) in items if c == cond and sigma.get(t)]
        vals = np.array([v for v in vals if v])
        if len(vals) == 0:
            continue
        out[cond] = {"n": len(vals), "sigma_kPa": (vals * 1e3).tolist(),
                     "p05": float(np.percentile(vals, 5) * 1e3),
                     "p50": float(np.percentile(vals, 50) * 1e3),
                     "p95": float(np.percentile(vals, 95) * 1e3),
                     "map_ref_kPa": MAP_SIGMA[cond]}
        print("  [%s] N=%d  sigma 90%%CI=[%.2f, %.2f] kPa  median=%.2f  MAP=%.2f"
              % (cond, len(vals), out[cond]["p05"], out[cond]["p95"],
                 out[cond]["p50"], MAP_SIGMA[cond]))
    # ratio CH/DH
    ch = np.array(out["commensal_hobic"]["sigma_kPa"])
    dh = np.array(out["dysbiotic_hobic"]["sigma_kPa"])
    rng = np.random.default_rng(42)
    rb = ch[rng.integers(0, len(ch), 5000)] / dh[rng.integers(0, len(dh), 5000)]
    out["ratio_ch_dh"] = {"map": MAP_SIGMA["commensal_hobic"] / MAP_SIGMA["dysbiotic_hobic"],
                          "p05": float(np.percentile(rb, 5)), "p50": float(np.percentile(rb, 50)),
                          "p95": float(np.percentile(rb, 95))}
    print("\n  sigma_CH/sigma_DH MAP=%.2fx  direct 90%%CI=[%.2f, %.2f]x"
          % (out["ratio_ch_dh"]["map"], out["ratio_ch_dh"]["p05"], out["ratio_ch_dh"]["p95"]))
    json.dump(out, open(OUT / "direct_propagation_ci.json", "w"), indent=2)
    print("saved", OUT / "direct_propagation_ci.json")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    if mode == "validate":
        validate()
    elif mode == "full":
        full(int(sys.argv[2]) if len(sys.argv) > 2 else 6)

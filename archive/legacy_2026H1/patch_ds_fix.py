"""Fix the dysbiotic_static (DS) composition bug: ref_0d had a commensal-like
So-dominant [0.944,...] (copy error). Replace with the raw CLSM dysbiotic-static
composition (V.dispar-dominant, ~D10), recompute di_0d = 1 - H/ln(5), and update
the CI-script headline values (MAP_PHI[DS], MAP_SIGMA[DS] = corrected Abaqus
1.55 kPa)."""
import io, json
import numpy as np

# raw CLSM dysbiotic-static D10 (V.parvula-dominant), normalized
phi_new = np.array([0.036, 0.057, 0.568, 0.129, 0.209]); phi_new = phi_new / phi_new.sum()
H = -np.sum(phi_new * np.log(phi_new + 1e-12))
di_new = float(1.0 - H / np.log(5.0))
print("corrected DS phi_final =", [round(x, 4) for x in phi_new])
print("corrected DS di_0d     = %.4f (was 0.8454)" % di_new)

# 1) ref_0d_dysbiotic_static.json
p = "_multiscale_2d_results/ref_0d_dysbiotic_static.json"
d = json.load(open(p))
old = d["phi_final"]
d["phi_final"] = [float(x) for x in phi_new]
d["di_0d"] = di_new
d["_fix_note"] = "2026-06-26: phi_final was commensal-like So-dominant [0.944,...] (copy bug); replaced with raw CLSM dysbiotic-static D10 (V.dispar-dominant). di_0d recomputed."
json.dump(d, open(p, "w"), indent=1)
print("patched", p, " (old phi_final So=%.3f -> new So=%.3f)" % (old[0], phi_new[0]))

# 2) CI script MAP_PHI[DS] + MAP_SIGMA[DS]
c = "JAXFEM/posterior_klempt_stress_ci.py"
s = io.open(c, encoding="utf-8").read()
a1 = '    "dysbiotic_static": np.array([0.944, 0.011, 0.011, 0.011, 0.011]),'
b1 = '    "dysbiotic_static": np.array([%.4f, %.4f, %.4f, %.4f, %.4f]),  # FIX 2026-06-26: was So-dom copy bug; raw CLSM DS D10 (Vd-dom)' % tuple(phi_new)
assert s.count(a1) == 1, ("MAP_PHI DS", s.count(a1)); s = s.replace(a1, b1)
a2 = '    "dysbiotic_static": 13.63e-3,'
b2 = '    "dysbiotic_static": 1.55e-3,  # FIX 2026-06-26: corrected Abaqus run w/ raw CLSM DS comp (was 13.63e-3 from So-dom bug)'
assert s.count(a2) == 1, ("MAP_SIGMA DS", s.count(a2)); s = s.replace(a2, b2)
io.open(c, "w", encoding="utf-8").write(s)
import ast; ast.parse(s)
print("patched", c, "(MAP_PHI[DS], MAP_SIGMA[DS]=13.63->1.55 kPa); parse OK")

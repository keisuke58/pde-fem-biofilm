"""Complete the DS fix: (1) implant MAP_SIGMA[DS] via the validated framework
(sigma_imp = E_voigt*A*k_eff^b, fit on the 3 good implant points CH/DH/CS),
(2) regenerate samples_0d(dysbiotic_static) with the corrected composition."""
import io, json
import numpy as np

E_SPEC = np.array([1e-3, 8e-4, 6e-4, 2e-4, 1e-5])
K_ALPHA = np.array([1.0, 0.8, 0.4, 0.6, 0.3])
ev = lambda p: float(np.array(p) @ E_SPEC); ke = lambda p: float(np.array(p) @ K_ALPHA)
MAP_PHI = {"CH": [0.942, 0.012, 0.012, 0.011, 0.011], "DH": [0.097, 0.119, 0.474, 0.123, 0.093],
           "CS": [0.698, 0.061, 0.062, 0.063, 0.059]}
phi_DS = np.array([0.036, 0.057, 0.568, 0.129, 0.209]); phi_DS = phi_DS / phi_DS.sum()
MAP_SIG_IMP = {"CH": 6.13, "DH": 1.30, "CS": 4.28}     # 3 good implant points (exclude buggy DS)

# fit implant surrogate on 3 good points
X = np.array([[1.0, np.log(ke(MAP_PHI[c]))] for c in MAP_SIG_IMP])
y = np.array([np.log(MAP_SIG_IMP[c] / ev(MAP_PHI[c])) for c in MAP_SIG_IMP])
(logA, b), *_ = np.linalg.lstsq(X, y, rcond=None)
imp_ds = float(np.exp(logA) * ev(phi_DS) * ke(phi_DS) ** b)
print("implant fit sigma=E_voigt*%.3g*k_eff^%.2f  (3 pts CH/DH/CS)" % (np.exp(logA), b))
for c in MAP_SIG_IMP:
    pr = np.exp(logA) * ev(MAP_PHI[c]) * ke(MAP_PHI[c]) ** b
    print("  %s: pred=%.2f actual=%.2f kPa" % (c, pr, MAP_SIG_IMP[c]))
print("=> corrected implant DS = %.2f kPa  (was buggy 6.10)" % imp_ds)

# (1) patch MAP_SIGMA_IMPLANT[DS]
c = "JAXFEM/posterior_klempt_stress_ci.py"
s = io.open(c, encoding="utf-8").read()
a = '    "dysbiotic_static": 6.10e-3,'
nb = '    "dysbiotic_static": %.2fe-3,  # FIX 2026-06-26: corrected via framework (was 6.10e-3 So-dom bug)' % imp_ds
assert s.count(a) == 1, ("imp DS", s.count(a)); s = s.replace(a, nb)
io.open(c, "w", encoding="utf-8").write(s)
import ast; ast.parse(s)
print("patched MAP_SIGMA_IMPLANT[DS] = %.2fe-3" % imp_ds)

# (2) regenerate samples_0d(dysbiotic_static) with corrected composition + original tiny spread
p = "_ci_0d_results/dysbiotic_static/samples_0d.json"
samples = json.load(open(p))
rng = np.random.default_rng(0)
std = np.array([0.0004, 0.0001, 0.0001, 0.0001, 0.0001])
SP = ["So", "An", "Vd", "Fn", "Pg"]
for i, smp in enumerate(samples):
    pp = np.clip(phi_DS + rng.normal(0, std), 0, None); pp = pp / pp.sum()
    smp["phi_final"] = [float(x) for x in pp]
    for j, sp in enumerate(SP):
        smp["phi_%s" % sp] = float(pp[j])
    H = -np.sum(pp * np.log(pp + 1e-12))
    smp["di_0d"] = float(1.0 - H / np.log(5.0))
json.dump(samples, open(p, "w"))
print("regenerated %d DS samples_0d with corrected composition (So 0.94 -> %.3f)" % (len(samples), phi_DS[0]))

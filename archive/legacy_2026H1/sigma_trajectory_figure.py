"""sigma(t) trajectory for all 4 conditions (fig3 raw CLSM, correct naming):
prints the DS composition check (MAP_PHI vs raw) and saves a 2-panel thesis
figure (sigma curves + CH/DH, CS/DS ratio curves)."""
import csv, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E_SPEC = np.array([1e-3, 8e-4, 6e-4, 2e-4, 1e-5])
K_ALPHA = np.array([1.0, 0.8, 0.4, 0.6, 0.3])
order = ["So", "An", "Vd", "Fn", "Pg"]
SPM = {"S. oralis": "So", "A. naeslundii": "An", "V. dispar": "Vd", "V. parvula": "Vd",
       "V. dispar/parvula": "Vd", "F. nucleatum": "Fn",
       "P. gingivalis": "Pg", "P. gingivalis_20709": "Pg", "P. gingivalis_W83": "Pg"}
MAP_PHI = {"CH": [0.942, 0.012, 0.012, 0.011, 0.011], "DH": [0.097, 0.119, 0.474, 0.123, 0.093],
           "CS": [0.698, 0.061, 0.062, 0.063, 0.059], "DS": [0.944, 0.011, 0.011, 0.011, 0.011]}
MAP_SIG = {"CH": 13.72, "DH": 2.13, "CS": 8.77, "DS": 13.63}
ev = lambda p: float(np.array(p) @ E_SPEC); ke = lambda p: float(np.array(p) @ K_ALPHA)
X = np.array([[1.0, np.log(ke(MAP_PHI[c]))] for c in MAP_PHI])
yv = np.array([np.log(MAP_SIG[c] / ev(MAP_PHI[c])) for c in MAP_PHI])
(logA, b), *_ = np.linalg.lstsq(X, yv, rcond=None)
sigma = lambda p: float(np.exp(logA) * ev(p) * ke(p) ** b)

rows = list(csv.DictReader(open("../Tmcmc202601/data_5species/experiment_data/fig3_species_distribution_summary.csv")))
g = collections.defaultdict(dict)
for r in rows:
    sp = SPM.get(r["species"])
    if sp:
        k = (r["condition"], r["cultivation"], int(r["day"])); g[k][sp] = g[k].get(sp, 0.0) + float(r["mean"])
def comp(cond, cult, d):
    v = np.array([g[(cond, cult, d)].get(s, 0.0) for s in order]); s = v.sum(); return v / s if s else v
NAME = {("Commensal", "HOBIC"): "CH", ("Dysbiotic", "HOBIC"): "DH",
        ("Commensal", "Static"): "CS", ("Dysbiotic", "Static"): "DS"}
days = sorted({d for (c, cu, d) in g})

# DS composition check
print("DS (Dysbiotic Static) composition: headline MAP_PHI vs raw CLSM")
print("  MAP_PHI[DS] =", MAP_PHI["DS"], "(So-dominant!)")
for d in days:
    print("  raw D%-3d   = %s" % (d, np.round(comp("Dysbiotic", "Static", d), 3)))

sig = {nm: [sigma(comp(c, cu, d)) for d in days] for (c, cu), nm in NAME.items()}
fig, ax = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
col = {"CH": "#1f77b4", "DH": "#d62728", "CS": "#2ca02c", "DS": "#ff7f0e"}
for nm in ["CH", "DH", "CS", "DS"]:
    ax[0].plot(days, sig[nm], "o-", color=col[nm], label=nm)
ax[0].set_ylabel("max von Mises stress  [kPa]"); ax[0].legend(); ax[0].grid(alpha=.3)
ax[0].set_title("Growth-induced stress vs biofilm age (raw CLSM composition)")
rH = [sig["CH"][i] / sig["DH"][i] for i in range(len(days))]
rS = [sig["CS"][i] / sig["DS"][i] for i in range(len(days))]
ax[1].plot(days, rH, "s-", color="#1f77b4", label="CH/DH (HOBIC)")
ax[1].plot(days, rS, "^-", color="#2ca02c", label="CS/DS (Static)")
ax[1].axhline(6.44, ls="--", color="grey", label="headline 6.44x")
ax[1].set_xlabel("cultivation day"); ax[1].set_ylabel("commensal/dysbiotic stress ratio")
ax[1].legend(); ax[1].grid(alpha=.3)
fig.tight_layout()
fig.savefig("sigma_trajectory_4cond.pdf"); fig.savefig("sigma_trajectory_4cond.png", dpi=130)
print("\nsaved sigma_trajectory_4cond.pdf / .png")

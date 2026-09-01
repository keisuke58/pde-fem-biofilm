"""Reproduces DEVIATOR_SCALING_FINDING.md.

The isochoric split in usermat_biofilm.f (and identically in
umat_biofilm_visco.f) applies J^(-2/3) to the subtracted trace but not to BE
beside it, so the "deviator" of an isotropic elastic state is non-zero and
spherical. It is exact only at J = 1 -- which growth, by construction, leaves.

Run directly:  python ansys_usermat/crosscheck/check_deviator_scaling.py
Needs gfortran + cc. Reports numbers only; changes nothing.
"""
import sys, numpy as np, tempfile, subprocess
from pathlib import Path
R = Path("/home/user/pde-fem-biofilm")
tmp = Path(tempfile.mkdtemp())
FC="gfortran"
subprocess.run([FC,"-c","-ffixed-line-length-132","-J",str(tmp),
    str(R/"ansys_usermat/coupling/usermat_py_hook.f"),"-o",str(tmp/"h.o")],check=True,cwd=tmp)
subprocess.run([FC,"-c","-ffixed-line-length-132","-I",str(tmp),
    str(R/"ansys_usermat/usermat_biofilm.f"),"-o",str(tmp/"c.o")],check=True,cwd=tmp)
subprocess.run(["cc","-c","-fPIC",str(R/"ansys_usermat/coupling/biofilm_py_eval.c"),"-o",str(tmp/"s.o")],check=True)
exe=tmp/"x"
subprocess.run([FC,"-ffixed-line-length-132","-I",str(tmp),
    str(R/"ansys_usermat/crosscheck/xcheck_driver_ans.f"),str(tmp/"c.o"),str(tmp/"h.o"),str(tmp/"s.o"),"-o",str(exe)],check=True)

I3=np.eye(3)
def run(F,Fv,alpha,c10,c01,d1,eta,mtype,dt):
    s=(" ".join(f"{F[i,j]:.17e}" for i in range(3) for j in range(3))+"\n"+
       " ".join(f"{Fv[i,j]:.17e}" for i in range(3) for j in range(3))+"\n"+
       f"{alpha:.17e} {c10:.17e} {c01:.17e} {d1:.17e} {eta:.17e} {mtype:.1f} {dt:.17e}\n")
    r=subprocess.run([str(exe)],input=s,capture_output=True,text=True,env={"PATH":"/usr/bin:/bin"},timeout=30)
    assert r.returncode==0, r.stderr
    v=[float(x) for x in r.stdout.split()]
    return np.array(v[:6]), np.array(v[6:15]).reshape(3,3)

# reference_values.json settings
C10,C01,D1,MTYPE,DT = 2e-4, 0.0, 5000.0, 0.0, 5.0
print("F = I, Fv = I, isotropic growth -> Fe is ISOTROPIC, so the deviatoric")
print("flow driver must vanish and Fv must stay I for ANY eta.")
print("(t_growth_free.dat states exactly this expectation.)\n")
for alpha in (0.05, 0.20):
    se,_  = run(I3,I3,alpha,C10,C01,D1,0.0,   MTYPE,DT)
    sv,fv = run(I3,I3,alpha,C10,C01,D1,8e-3,  MTYPE,DT)
    print(f"alpha={alpha}")
    print(f"  eta=0     s11={se[0]: .6e}")
    print(f"  eta=8e-3  s11={sv[0]: .6e}   <- differs by {abs(sv[0]-se[0])/abs(se[0])*100:.1f}%")
    print(f"  Fv moved off I by {np.max(np.abs(fv-I3)):.3e}   (must be 0)")
    print(f"  Fv deviatoric part = {np.max(np.abs(fv-np.trace(fv)/3*I3)):.3e}  -> purely spherical\n")

# what the deviator SHOULD be, analytically, for this state
for alpha in (0.05,0.20):
    b = 1.0/(1.0+alpha)**2
    Be = b*I3
    detFe = 1.0/(1.0+alpha)**3
    t1 = detFe**(-2.0/3.0)
    I1B = t1*np.trace(Be)
    tau_code    = 2*C10*t1*(Be - (I1B/3.0)*I3)
    tau_correct = 2*C10*t1*(Be - (np.trace(Be)/3.0)*I3)
    print(f"alpha={alpha}: tau_code(1,1)={tau_code[0,0]: .6e}  tau_correct(1,1)={tau_correct[0,0]: .6e}")

print("\n" + "="*66)
print("Is the discrepancy specifically activated by GROWTH (Je != 1)?")
print("="*66)
for alpha in (0.05, 0.20):
    # traction-free growth: F = (1+a)I  ->  Fe = I exactly, Je = 1
    Ffree = (1.0+alpha)*I3
    sf,_ = run(Ffree,I3,alpha,C10,C01,D1,8e-3,MTYPE,DT)
    # fully constrained: F = I -> Je = 1/(1+a)^3 != 1
    sc,_ = run(I3,I3,alpha,C10,C01,D1,8e-3,MTYPE,DT)
    print(f"alpha={alpha}: free-growth (Je=1) peak|s|={np.max(np.abs(sf)):.3e}"
          f"   constrained (Je={1/(1+alpha)**3:.4f}) s11={sc[0]:.6e}")

print("\nWhat the CONSTRAINED elastic stress should be (eta=0), analytically:")
for alpha in (0.05,0.20):
    J = 1.0/(1.0+alpha)**3
    vol = (2.0/D1)*(J-1.0)              # correct: deviator vanishes, only volumetric
    se,_ = run(I3,I3,alpha,C10,C01,D1,0.0,MTYPE,DT)
    print(f"  alpha={alpha}: correct s11={vol: .6e}   code s11={se[0]: .6e}"
          f"   ratio={se[0]/vol:.3f}x")

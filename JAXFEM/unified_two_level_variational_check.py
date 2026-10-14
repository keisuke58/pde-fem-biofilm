"""Machine checks of the variational structure of the two-level model (unified_two_level_1d.py).

1. The residual used by the solver is the JAX gradient of the discrete functional
   E_h + dt D_h + constraint (functional()), so it is the discrete Euler-Lagrange equation by
   construction.  Here it is compared with the earlier hand-written residual (residual_hand):
   every composition row agrees to round-off.  The phase row differs on purpose: the hand
   version approximated the derivative of the composition-gradient energy (relative 9e-5), and
   it let the point model's barrier terms act on the phase, which the model no longer does
   (they are a regularised constraint, not stored energy; see functional()).
2. With one node, no growth, no gradients and lam = 1, the composition rows equal the point
   model as coded in this repository (core_hamilton_1d.residual) to round-off: the point model
   of Klempt et al. 2026 is itself the stationarity condition of the energy Psi_pm with the
   dissipation D_pm used here.

Run:  python JAXFEM/unified_two_level_variational_check.py
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / 'JAXFEM', ROOT / 'ansys_usermat/coupling'):
    sys.path.insert(0, str(p))
import numpy as np, jax, jax.numpy as jnp
jax.config.update("jax_enable_x64", True)
import unified_two_level_1d as t2, core_hamilton_1d as ch, ecology_jax as ej
from ecology_constants import C_STAR
rng=np.random.default_rng(1); N=7
th=rng.normal(size=20); A,b=ch.theta_to_matrices(jnp.asarray(th))
P={"N":N,"dt":1e-3,"dx":1/(N-1),"c":C_STAR,"A":A,"b":jnp.ones(5),"ast":0.3,"lam":0.7,"ka":1.3,"beta":0.02,"beta_s":0.05,"eta_phi":1.0,"eps":1e-6}
X=np.zeros((N,14));Xn=np.zeros((N,14))
for Y in (X,Xn):
  Y[:,0]=rng.uniform(0.2,0.9,N); s=rng.uniform(0.05,0.3,(N,5)); Y[:,1:6]=s; Y[:,6]=1-s.sum(1)
  Y[:,7:12]=rng.uniform(0.3,0.95,(N,5)); Y[:,12]=rng.normal(size=N); Y[:,13]=1+rng.uniform(0,0.1,N)
U,Un=jnp.asarray(X.ravel()),jnp.asarray(Xn.ravel())
Ra=np.asarray(t2.residual(U,Un,P)).reshape(N,14); Rh=np.asarray(t2.residual_hand(U,Un,P)).reshape(N,14)
print('autodiff vs hand, per column rel:',np.round(np.abs(Ra-Rh).max(0)/(np.abs(Rh).max(0)+1e-30),12))
# point model rows vs the repository's point model (core_hamilton_1d.residual), one node, no growth, no gradients
P2=dict(P,ka=0.0,beta_s=0.0,ast=0.0,lam=1.0,N=1)
U1=jnp.asarray(X[0]); U1n=jnp.asarray(Xn[0])
R1=np.asarray(t2.residual(U1,U1n,P2))
p=jax.tree_util.tree_map(jnp.asarray,ej.default_hparams(1e-3)); p['A']=A; p['b_diag']=b
g=jnp.concatenate([U1[1:12],U1[12:13]]); gn=jnp.concatenate([U1n[1:12],U1n[12:13]])
Rc=np.asarray(ch.residual(g,gn,p))
print('composition rows vs core point model:',np.round(np.abs(R1[1:13]-Rc)/(np.abs(Rc)+1e-30),12))

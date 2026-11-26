"""Produces Table 2: normalized error q^{3/2}|E(q,p)| in Theorem 3."""
import numpy as np, math
from family_prime import *
H=np.load("H_10000000.npy").astype(np.float64)
def J(nu,x):
    th=np.linspace(0,np.pi,2001); return np.trapezoid(np.cos(nu*th-x*np.sin(th)),th)/np.pi
P=primes_upto(6000)
for k in [2,3,4]:
  for q in [101,211,401,809]:
    D=family_dim(q,k); errs=[]
    for p in P:
        p=int(p)
        if p>=3*q: break
        if p==q: continue
        v=family_sum_conj(q,p,k,H).real/D
        errs.append(abs(v-12/(math.pi*math.sqrt(q))*(1+1/p)*J(k-1,4*math.pi*math.sqrt(p/q))))
    e=np.array(errs)
    print(f"k={k} q={q}  max err*q^1.5={e.max()*q**1.5:.2f}  mean err*q^1.5={e.mean()*q**1.5:.2f}   max err*q^0.5 (rel. to main)={e.max()*q**.5:.4f}")

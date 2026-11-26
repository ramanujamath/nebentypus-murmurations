"""Section 6.1: shows that the q^{-3/2} error in Theorem 3 has a smooth part (binned curves agree across q)."""
import numpy as np, math
from family_prime import *
H=np.load("H_10000000.npy").astype(np.float64)
def J(nu,x):
    th=np.linspace(0,np.pi,2001); return np.trapezoid(np.cos(nu*th-x*np.sin(th)),th)/np.pi
P=primes_upto(8000)
bins=np.linspace(0,3,16)
for k in [2,3]:
  rows=[]
  for q in [401,601,809]:
    D=family_dim(q,k); ys=[];es=[]
    for p in P:
        p=int(p)
        if p>=3*q: break
        if p==q: continue
        v=family_sum_conj(q,p,k,H).real/D; y=p/q
        main=12/(math.pi*math.sqrt(q))*(1+1/p)*J(k-1,4*math.pi*math.sqrt(y))
        ys.append(y); es.append((v-main)*q**1.5)
    ys=np.array(ys); es=np.array(es); idx=np.digitize(ys,bins)-1
    rows.append([es[idx==b].mean() for b in range(15)])
    # within-bin scatter vs between-bin signal
    within=np.sqrt(np.mean([es[idx==b].var() for b in range(15)]))
    print(f"k={k} q={q}: within-bin sd {within:.1f}")
  R=np.array(rows); print("binned means (rows q=401,601,809):"); print(np.round(R,1))
  print("corr(401,601)=%.2f corr(601,809)=%.2f"%(np.corrcoef(R[0],R[1])[0,1],np.corrcoef(R[1],R[2])[0,1]))

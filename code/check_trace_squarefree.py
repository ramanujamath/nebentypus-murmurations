"""Checks the squarefree-level trace formula: Hecke multiplicativity and integrality of dimensions."""
import numpy as np, math, cmath
from trace_prime import hurwitz_table
from trace_squarefree import *
H=hurwitz_table(3000); hw=hw_table(H,3000)
print("hw check", [round(hw[m],3) for m in [3,4,7,8,11,12,15,16,20,23,27,28]])
# (1) prime level q=5,k=5: a(19) now consistent with multiplicativity tests -> check eta-product style: Ramanujan bound & integrality
for (N,k) in [(5,5),(7,3)]:
    for tab,prim in chars_mod(N):
        if not prim or abs(tab[N-1]-(-1)**k)>1e-9: continue
        d=round(trace(1,N,tab,k,H,hw).real)
        if d!=1: continue
        a={n:trace(n,N,tab,k,H,hw) for n in range(1,41)}
        # check multiplicativity a(mn)=a(m)a(n) for coprime and Hecke at p^2
        errs=[abs(a[m*n]-a[m]*a[n]) for m in range(2,7) for n in range(2,7) if math.gcd(m,n)==1 and m*n<=40]
        p2=[abs(a[p*p]-(a[p]**2-tab[p%N]*p**(k-1))) for p in [2,3] if p*p<=40 and N%p]
        print(N,k,"mult err",max(errs),"p^2 err",max(p2) if p2 else None)
# (2) dims for N=15,21,33,35
for (N,k) in [(15,2),(15,3),(21,2),(21,3),(35,2)]:
    ds=[]
    for tab,prim in chars_mod(N):
        if not prim or abs(tab[N-1]-(-1)**k)>1e-9: continue
        ds.append(np.round(trace(1,N,tab,k,H,hw),4))
    print(N,k,"dims",ds)

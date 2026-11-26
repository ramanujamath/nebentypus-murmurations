"""Checks the factorization S_chi(p,q;qr) = chi(p rbar) tau(chibar) S(p qbar,1;r) behind the single-character formula (6)."""
import numpy as np, math, cmath
from trace_prime import char_table
def Sx(m,n,c,chi,q):
    s=0
    for d in range(c):
        if math.gcd(d,c)!=1: continue
        s+=chi[d%q]*cmath.exp(2j*math.pi*(m*pow(d,-1,c)+n*d)/c)
    return s
def S(m,n,c):
    return sum(cmath.exp(2j*math.pi*(m*pow(d,-1,c)+n*d)/c) for d in range(c) if math.gcd(d,c)==1) if c>1 else 1
q=11
for j in [1,2,5]:
    chi=char_table(q,j); tau_bar=sum(np.conj(chi[x])*cmath.exp(2j*math.pi*x/q) for x in range(q))
    for p in [7,13]:
        for r in [2,3,4,5,6,9]:
            lhs=Sx(p,q,q*r,chi,q)
            rhs=chi[(p*pow(r,-1,q))%q]*tau_bar*S(p*pow(q,-1,r)%r,1,r)
            print(j,p,r,np.round(lhs,6),np.round(rhs,6))

"""Compares the closed-form family sum of eps_f * conj(lambda_f(p)) (Proposition 7) with direct eigenform computations."""
from trace_prime import hurwitz_table, char_table, trace
from family_prime import *
import cmath, math
H=hurwitz_table(4000)
for (q,k) in [(11,3),(7,4),(13,2)]:
    tot={p:0 for p in [2,3,5,7,13,17,23] if p!=q}
    for j in range(1,q-1):
        chi=char_table(q,j)
        if abs(chi[q-1]-(-1)**k)>1e-9: continue
        d=round(trace(1,q,chi,k,H).real)
        if d==0: continue
        a=[0]+[trace(n,q,chi,k,H) for n in range(1,31)]
        f=lambda z,c: sum(c[n]*cmath.exp(2j*math.pi*n*z) for n in range(1,31))
        ac=[x.conjugate() for x in a]
        z=cmath.exp(1.3j)/math.sqrt(q); w=-1/(q*z)
        eps=(1j)**k*q**(k/2)*(q*z)**(-k)*f(w,a)/f(z,ac)
        for p in tot: tot[p]+=eps*a[p].conjugate()/p**((k-1)/2)
    for p in tot: print(q,k,"p",p,np.round(tot[p],6),np.round(family_sum_conj(q,p,k,H),6))

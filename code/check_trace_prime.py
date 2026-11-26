"""Checks the prime-level trace formula: Hurwitz values, Hecke multiplicativity of traces, and the eta-product
eta(z)^3 eta(7z)^3 in S_3(7, (./7))."""
from trace_prime import *
H=hurwitz_table(5000)
print({n:round(H[n],4) for n in [3,4,7,8,11,12,15,16,19,20,23,24,27,28]})
# eta(z)^3 eta(7z)^3 coefficients
M=60
c=[0]*M; c[1]=1  # q^1 * prod
import itertools
poly=[0]*M; poly[0]=1
def mul(p,e,step):
    # multiply by (1-q^{step*m})^e
    for m in range(1,M):
        s=step*m
        if s>=M: break
        for _ in range(e):
            for i in range(M-1,s-1,-1):
                p[i]-=p[i-s]
mul(poly,3,1); mul(poly,3,7)
a=[0]+poly[:M-1]
print("eta coeffs a(1..14):",a[1:15])
q,k=7,3
# quadratic char mod 7 = exponent j=(q-1)/2
chi=char_table(q,3)
print("chi values",[round(chi[x].real) for x in range(7)])
for n in [1,2,3,4,5,6,8,9,10,11,12,14,21]:
    print(n, round(trace(n,q,chi,k,H).real,6), round(trace(n,q,chi,k,H).imag,6), "eta:",a[n])
# dims for all odd chars mod 7 weight 3 and 5, mod 11, mod 13
for (q,k) in [(5,3),(7,3),(7,5),(11,2),(11,3),(13,2),(13,4)]:
    out=[]
    for j in range(1,q-1):
        chi=char_table(q,j)
        if abs(chi[q-1]-(-1)**k)>1e-9: continue
        out.append(round(trace(1,q,chi,k,H).real,6))
    print(q,k,"dims",out)

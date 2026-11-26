"""Checks the steps of the proof of Theorem 3: Kronecker-Hurwitz, smallness of Chebyshev modes, twisted sum."""
import numpy as np, math
from family_prime import *
H=np.load("H_10000000.npy").astype(np.float64)
def J(nu,x):
    th=np.linspace(0,np.pi,4001); return np.trapezoid(np.cos(nu*th-x*np.sin(th)),th)/np.pi
def U(m,u):
    u0,u1=np.ones_like(u),2*u
    if m==0: return u0
    for _ in range(m-1): u0,u1=u1,2*u*u1-u0
    return u1
# Step 1: Kronecker-Hurwitz  sum_t H(4n-t^2) = 2sigma(n) - 2 - 2min(p,q), n=pq
for (p,q) in [(1193,797),(53,457),(911,467)]:
    n=p*q; r=math.isqrt(4*n-1); t=np.arange(-r,r+1)
    print("KH",p,q, H[4*n-t*t].sum(), 2*(p+1)*(q+1)-2-2*min(p,q))
# Step 2: Chebyshev modes m>=2 are O(m tau(n) sqrt n)
p,q=1193,797; n=p*q; r=math.isqrt(4*n-1); t=np.arange(-r,r+1); u=t/(2*math.sqrt(n)); h=H[4*n-t*t]
for m in [0,2,4,6,10,20]:
    print("mode",m, round((U(m,u)*h).sum(),1), " sqrt(n)=",round(math.sqrt(n)))
# Step 3: twisted sum vs c0 * KH
for k in [2,3,4]:
    for p in [79,397,1193,1993]:
        n=p*q; r=math.isqrt(4*n-1); t=np.arange(-r,r+1); u=t/(2*math.sqrt(n)); h=H[4*n-t*t]
        y=p/q; a=4*math.pi*math.sqrt(y)
        S=(U(k-2,u)*np.exp(1j*a*u)*h).sum()
        c0=2*(1j)**(k-2)*(k-1)*J(k-1,a)/a
        main=c0*(2*(p+1)*(q+1)-2-2*min(p,q))
        print("twist k",k,"p",p, np.round(S,1), np.round(main,1), " diff/sqrt(n)",np.round(abs(S-main)/math.sqrt(n),2))

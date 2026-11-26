"""Confirms the normalization of the Petersson formula with nebentypus (equation (3)) numerically, q = 11, k = 3."""
import numpy as np, math, cmath
from trace_prime import hurwitz_table, char_table, trace
H=hurwitz_table(3000)
def J(nu,x):
    th=np.linspace(0,np.pi,4001); return np.trapezoid(np.cos(nu*th-x*np.sin(th)),th)/np.pi
def kl(m,n,c,chi,q,conv):
    d=np.arange(c); g=np.gcd(d,c); d=d[g==1]
    db=np.array([pow(int(x),-1,c) for x in d])
    ch=np.array([chi[x%q] for x in d])
    if conv=="chi(d)": w=ch
    else: w=np.conj(ch)
    return np.sum(w*np.exp(2j*np.pi*(m*db+n*d)/c))
def Delta(m,n,chi,q,k,conv,C=1500):
    s=0
    for r in range(1,C+1):
        c=q*r; s+=kl(m,n,c,chi,q,conv)/c*J(k-1,4*math.pi*math.sqrt(m*n)/c)
    return (1 if m==n else 0)+2*math.pi*(1j)**(-k)*s
q,k=11,3
for j in [1,3]:
    chi=char_table(q,j)
    lam=lambda n: trace(n,q,chi,k,H)/n**((k-1)/2)
    for conv in ["chi(d)","chibar(d)"]:
        D11=Delta(1,1,chi,q,k,conv)
        for (m,n) in [(2,11),(3,2)]:
            print(j,conv,(m,n),"Delta(m,n)/Delta(1,1)=",np.round(Delta(m,n,chi,q,k,conv)/D11,4),
                  " lam(m)conj lam(n)=",np.round(lam(m)*np.conj(lam(n)),4),
                  " conj lam(m) lam(n)=",np.round(np.conj(lam(m))*lam(n),4))

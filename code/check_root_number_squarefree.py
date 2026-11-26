"""Checks the root number formula (4) at composite level N = 21."""
import numpy as np, math, cmath
from trace_prime import hurwitz_table
from trace_squarefree import *
H=hurwitz_table(3000); hw=hw_table(H,3000)
for (N,k) in [(21,2),(33,2),(39,2),(15,4),(21,4)]:
    for tab,prim in chars_mod(N):
        if not prim or abs(tab[N-1]-(-1)**k)>1e-9: continue
        if round(trace(1,N,tab,k,H,hw).real)!=1: continue
        a=[0]+[trace(n,N,tab,k,H,hw) for n in range(1,46)]
        f=lambda z,c: sum(c[n]*cmath.exp(2j*math.pi*n*z) for n in range(1,46))
        ac=[x.conjugate() for x in a]
        g=sum(tab[x]*cmath.exp(2j*math.pi*x/N) for x in range(N))
        for th in [1.2,1.8]:
            z=cmath.exp(1j*th)/math.sqrt(N); w=-1/(N*z)
            eta=N**(k/2)*(N*z)**(-k)*f(w,a)/f(z,ac)
            print(N,k,"eta/(chi(-1) tau N^(k/2-1)/a(N)) =",np.round(eta/(tab[N-1]*g*N**(k/2-1)/a[N]),6), " |a(N)|^2/N^(k-1)=",round(abs(a[N])**2/N**(k-1),6))

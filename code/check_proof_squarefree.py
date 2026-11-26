"""Checks the two inputs to Theorem 12: the Eisenstein main term 2N(p+1) and smallness of psi-twisted sums."""
import numpy as np, math
from trace_squarefree import chars_mod, factor_sqfree
from family_prime import primes_upto
H=np.load("H_10000000.npy").astype(np.float64)
for N in [865,899,471]:
    ls=factor_sqfree(N); chars=chars_mod(N)
    for p in [int(x) for x in primes_upto(3*N) if x>N//2][:3]:
        n=p*N; r=math.isqrt(4*n-1); t=np.arange(-r,r+1); h=H[4*n-t*t]
        cop=np.gcd(t,N)==1
        main=2*N*(p+1)
        S0=h[cop].sum()
        # nontrivial psi: max |sum_t H psi(t)| over psi nonprincipal with psi(-1)=1
        R=np.bincount(t%N,weights=h,minlength=N)
        mx=0
        for tab,prim in chars:
            tab=np.array(tab)
            if abs(tab[N-1]-1)>1e-9: continue
            if np.allclose(tab[np.gcd(np.arange(N),N)==1],1): continue
            mx=max(mx,abs(tab@R))
        print(f"N={N}={ls} p={p}: sum_(t,N)=1 H = {S0:.0f}, 2N(p+1) = {main},  rel.diff={(S0-main)/main:+.4f};  max|sum H psi| / sqrt(n) = {mx/math.sqrt(n):.2f}  (vs main/sqrt(n) = {main/math.sqrt(n):.0f})")

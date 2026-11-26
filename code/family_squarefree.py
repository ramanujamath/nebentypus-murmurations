"""Natural-weight family sums at odd squarefree level N, averaged over primitive characters (Theorem 12).
Usage: python family_squarefree.py k N1 N2 ...   (prints the fitted constant C_N; Table 3)."""
import numpy as np, math, cmath, sys
from trace_squarefree import factor_sqfree, chars_mod, Pk, trace, hw_table, mobius
from family_prime import Pk_vec, primes_upto
H=np.load("H_10000000.npy").astype(np.float64)
def J(nu,x):
    th=np.linspace(0,np.pi,2001); return np.trapezoid(np.cos(nu*th-np.outer(np.atleast_1d(x),np.sin(th))),th,axis=1)/np.pi
def family(N,k,pmax_ratio=3.0):
    ps=factor_sqfree(N)
    from trace_prime import hurwitz_table
    Hsmall=hurwitz_table(100); hwsmall=hw_table(Hsmall,100)
    chars=[tab for tab,prim in chars_mod(N) if prim and abs(tab[N-1]-(-1)**k)<1e-9]
    C=np.array(chars)                         # (nchi, N)
    tau=C@np.exp(2j*np.pi*np.arange(N)/N)     # Gauss sums
    D=sum(trace(1,N,list(c),k,Hsmall,hwsmall).real for c in chars)
    divsN=[d for d in range(1,N+1) if N%d==0]
    out=[]
    for p in primes_upto(int(pmax_ratio*N)):
        p=int(p)
        if N%p==0: continue
        n=p*N; r=math.isqrt(4*n-1); t=np.arange(-r,r+1)
        w=Pk_vec(t,n,k)*H[4*n-t*t]
        R=np.bincount(t%N,weights=w,minlength=N)   # sum over t in each class mod N
        main=-0.5*(np.conj(C)@R)                   # chibar(t), zero when gcd(t,N)>1
        # A3 for chibar
        A3=np.zeros(len(chars),complex)
        for d in [d1*d2 for d1 in (1,p) for d2 in divsN]:
            e=n//d; inner=np.zeros(len(chars),complex)
            for tt in divsN:
                M=N//tt; y=next(y for y in range(N) if (y-d)%tt==0 and (y-e)%M==0)
                inner+=np.conj(C[:,y])
            A3+=min(d,e)**(k-1)*inner
        tr=main-0.5*A3
        s=np.sum((1j)**(-k)*tau/math.sqrt(N)*tr)/n**((k-1)/2)
        out.append((p,s.real/D,abs(s.imag)/D))
    return out,D,len(chars)
if __name__=="__main__":
    for N in [int(x) for x in sys.argv[2:]]:
        k=int(sys.argv[1])
        out,D,nc=family(N,k)
        ps=np.array([o[0] for o in out]); v=np.array([o[1] for o in out])
        pred=12/(math.pi*math.sqrt(N))*(1+1/ps)*J(k-1,4*math.pi*np.sqrt(ps/N))
        Cfit=(v@pred)/(pred@pred)
        res=v-Cfit*pred
        ls=factor_sqfree(N)
        print(f"N={N}={ls} k={k} #chi={nc} D={D:.0f}  fitted C_N={Cfit:.4f}  rel.resid={np.sqrt(np.mean(res**2))/np.sqrt(np.mean(pred**2)):.3f}  max imag={max(o[2] for o in out):.1e}")

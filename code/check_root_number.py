"""Numerically determines the Fricke pseudo-eigenvalue in one-dimensional spaces S_k(q, chi) and compares it
with the Atkin-Li formula (Lemma 4 of the paper)."""
from trace_prime import *
H=hurwitz_table(2000)
for (q,k) in [(13,2),(17,2),(7,4),(5,4),(5,5),(7,6),(5,6)]:
  for j in range(1,q-1):
    chi=char_table(q,j)
    if abs(chi[q-1]-(-1)**k)>1e-9: continue
    if round(trace(1,q,chi,k,H).real)!=1: continue
    a=[0]+[trace(n,q,chi,k,H) for n in range(1,31)]
    f=lambda z,c: sum(c[n]*cmath.exp(2j*math.pi*n*z) for n in range(1,31))
    ac=[x.conjugate() for x in a]
    z=cmath.exp(1.3j)/math.sqrt(q); w=-1/(q*z)
    eta=q**(k/2)*(q*z)**(-k)*f(w,a)/f(z,ac)
    g=sum(chi[x]*cmath.exp(2j*math.pi*x/q) for x in range(q))
    print(q,k,j,"ratio eta/(g q^(k/2-1)/a(q)) =",np.round(eta/(q**(k/2-1)*g/a[q]),6))

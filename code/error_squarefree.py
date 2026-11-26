"""Error test for Theorem 12 at levels N = 3l and 5l."""
import numpy as np, math
from family_squarefree import family, J
from trace_squarefree import factor_sqfree
for k in [3,4]:
  for N in [159,303,597,879, 5*37, 5*79, 5*157]:
    out,D,nc=family(N,k)
    ps=np.array([o[0] for o in out]); v=np.array([o[1] for o in out])
    C=math.prod(l*l/(l*l-1) for l in factor_sqfree(N))
    pred=12/(math.pi*math.sqrt(N))*C*(1+1/ps)*J(k-1,4*math.pi*np.sqrt(ps/N))
    e=np.abs(v-pred)/(1+ps/N)
    print(f"k={k} N={N}={factor_sqfree(N)}  max err*N/(1+p/N) = {e.max()*N:.2f}   mean = {e.mean()*N:.2f}")

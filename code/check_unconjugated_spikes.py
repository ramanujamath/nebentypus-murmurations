"""Tests Theorem 2: the unconjugated statistic matches the conjugated one when p = 2q +- 1."""
import numpy as np, math
from family_prime import *
H=np.load("H_10000000.npy").astype(np.float64)
P=set(primes_upto(4000).tolist())
for k in [2,3]:
  n=0
  for q in sorted(x for x in P if 450<=x<=900):
    for p in [2*q+1,2*q-1,q+2]:
        if p in P and n<8:
            n+=1
            D=family_dim(q,k)
            print(k,q,p,"p mod q=",p%q if p%q<q//2 else p%q-q,
                  " unconj",round(math.sqrt(q)*family_sum(q,p,k,H).real/D,4),
                  " conj",round(math.sqrt(q)*family_sum_conj(q,p,k,H).real/D,4))

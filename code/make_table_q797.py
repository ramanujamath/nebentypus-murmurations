"""Produces Table 1 (natural averages at q = 797)."""
import numpy as np, math
from family_prime import *
H=np.load("H_10000000.npy").astype(np.float64)
def J(nu,x):
    th=np.linspace(0,np.pi,4001); return np.trapezoid(np.cos(nu*th-x*np.sin(th)),th)/np.pi
P=primes_upto(3000)
q=797
for k in [2,3,4]:
    D=family_dim(q,k)
    for target in [0.1,0.5,1.0,1.5,2.5]:
        PP=P[P!=q]; p=int(PP[np.argmin(abs(PP/q-target))])
        v=math.sqrt(q)*family_sum_conj(q,p,k,H).real/D
        pr=12/math.pi*J(k-1,4*math.pi*math.sqrt(p/q))
        print(f"{k} & {p} & {p/q:.3f} & {v:.4f} & {pr:.4f} \\\\")

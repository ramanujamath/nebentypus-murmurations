"""Produces Figure 1 (nebentypus_murmurations.png)."""
import numpy as np, math
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from family_prime import *
H=np.load("H_10000000.npy").astype(np.float64)
def J(nu,x):
    th=np.linspace(0,np.pi,2001); return np.trapezoid(np.cos(nu*th-np.outer(x,np.sin(th))),th,axis=1)/np.pi
X=450; c=3.0
P=primes_upto(int(c*2*X)+10); qs=[int(q) for q in P if X<=q<=2*X]
fig,axs=plt.subplots(2,2,figsize=(11,7.5),sharex=True)
for ax,k in zip(axs.flat,[2,3,4,6]):
    xs=[];v2=[];v1=[]
    for q in qs:
        D=family_dim(q,k)
        for p in P:
            p=int(p)
            if p>=c*q: break
            if p==q: continue
            xs.append(p/q)
            v2.append(math.sqrt(q)*family_sum_conj(q,p,k,H).real/D)
            v1.append(math.sqrt(q)*family_sum(q,p,k,H).real/D)
    xs=np.array(xs);v2=np.array(v2);v1=np.array(v1)
    bins=np.linspace(0,c,76);mid=(bins[1:]+bins[:-1])/2;idx=np.digitize(xs,bins)-1
    m2=np.array([v2[idx==b].mean() for b in range(75)]);m1=np.array([v1[idx==b].mean() for b in range(75)])
    y=np.linspace(0.001,c,600); pred=12/math.pi*J(k-1,4*math.pi*np.sqrt(y))
    ax.scatter(xs,v2,s=1,alpha=0.08,color="tab:blue")
    ax.plot(mid,m2,"o",ms=3,color="tab:blue",label=r"$\sqrt{q}\,\langle \varepsilon_f\overline{\lambda_f(p)}\rangle$ (binned)")
    ax.plot(mid,m1,"s",ms=2.5,color="tab:orange",label=r"$\sqrt{q}\,\langle \varepsilon_f\lambda_f(p)\rangle$")
    ax.plot(y,pred,"k-",lw=1.2,label=r"$\frac{12}{\pi}J_{k-1}(4\pi\sqrt{p/q})$")
    ax.set_ylim(-3,3.5); ax.axhline(0,color="gray",lw=0.5); ax.set_title(f"weight k={k}")
    pv=12/math.pi*J(k-1,4*math.pi*np.sqrt(mid))
    print(k,"rms(binned-pred)",np.sqrt(np.nanmean((m2-pv)**2)).round(4),"rms pred",np.sqrt(np.mean(pv**2)).round(3),"rms unconj",np.sqrt(np.mean(m1**2)).round(4))
for ax in axs[1]: ax.set_xlabel("p / q")
axs[0,0].legend(fontsize=8,loc="upper right")
fig.suptitle(r"Murmurations for $S_k(q,\chi)$, q prime in [450,900], all nontrivial $\chi$ with $\chi(-1)=(-1)^k$",fontsize=11)
fig.tight_layout(); fig.savefig("nebentypus_murmurations.png",dpi=140)

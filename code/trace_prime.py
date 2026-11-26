"""Eichler-Selberg / Cohen trace formula on S_k(q, chi), q prime, chi primitive mod q (nontrivial).
Characters are represented by a generator g mod q and an exponent j: chi(g^m) = exp(2 pi i j m/(q-1))."""
import cmath, math
import numpy as np

def hurwitz_table(M):
    """H[n] = Hurwitz class number H(n) for 0<=n<=M (H(0) unused), weights 1/2, 1/3."""
    H = np.zeros(M + 1)
    a = 1
    while 3 * a * a <= M:
        for b in range(-a + 1, a + 1):
            c0 = a if b >= 0 else a + 1
            start = 4 * a * c0 - b * b
            if start > M:
                continue
            H[start::4 * a] += 1.0
        # corrections at c=a: b=0 -> weight 1/2 (n=4a^2), b=a -> weight 1/3 (n=3a^2)
        if 4 * a * a <= M:
            H[4 * a * a] -= 0.5
        H[3 * a * a] -= 2.0 / 3.0
        a += 1
    return H

def primroot(q):
    fac = [p for p in range(2, q) if (q - 1) % p == 0 and all(p % r for r in range(2, int(p**.5) + 1))]
    for g in range(2, q):
        if all(pow(g, (q - 1) // p, q) != 1 for p in fac):
            return g

def char_table(q, j):
    g = primroot(q)
    tab = [0j] * q
    x = 1
    for m in range(q - 1):
        tab[x] = cmath.exp(2j * math.pi * j * m / (q - 1))
        x = x * g % q
    return tab

def Pk(t, n, k):
    # (rho^{k-1}-rhobar^{k-1})/(rho-rhobar) via recurrence u_0=0,u_1=1,u_{m+1}=t u_m - n u_{m-1}
    u0, u1 = 0, 1
    for _ in range(k - 2):
        u0, u1 = u1, t * u1 - n * u0
    return u1

def trace(n, q, chi, k, H):
    tot = 0j
    r = math.isqrt(n)
    if r * r == n:
        tot += n ** (k / 2 - 1) * chi[r % q] * (k - 1) * (q + 1) / 12
    t = 0
    s2 = 0j
    for t in range(-(2 * r + 2), 2 * r + 3):
        D = 4 * n - t * t
        if D <= 0:
            continue
        mu = sum(chi[x] for x in range(q) if (x * x - t * x + n) % q == 0)
        # assumes q does not divide any f with f^2 | D (checked by caller's choice of n)
        s2 += Pk(t, n, k) * H[D] * mu
    tot += -0.5 * s2
    s3 = 0j
    for d in range(1, n + 1):
        if n % d == 0:
            s3 += min(d, n // d) ** (k - 1) * (chi[(n // d) % q] + chi[d % q])
    tot += -0.5 * s3
    return tot

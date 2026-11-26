"""Cohen/Eichler-Selberg trace formula on S_k(Gamma_0(N), chi) for N odd squarefree (any chi mod N)."""
import math, cmath
import numpy as np
from functools import lru_cache

def factor_sqfree(N):
    ps = []; m = N; p = 2
    while p * p <= m:
        if m % p == 0:
            ps.append(p); m //= p
            assert m % p != 0
        p += 1
    if m > 1: ps.append(m)
    return ps

def hw_table(H, M):
    """h_w(-m) for m<=M via Mobius inversion of H(m)=sum_{f^2|m} h_w(-m/f^2)."""
    out = np.zeros(M + 1)
    for m in range(1, M + 1):
        s = 0.0; f = 1
        while f * f <= m:
            if m % (f * f) == 0:
                mu = mobius(f)
                if mu: s += mu * H[m // (f * f)]
            f += 1
        out[m] = s
    return out

@lru_cache(None)
def mobius(n):
    r = 1; p = 2; m = n
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    if m > 1: r = -r
    return r

def Pk(t, n, k):
    u0, u1 = 0, 1
    for _ in range(k - 2):
        u0, u1 = u1, t * u1 - n * u0
    return u1

def trace(n, N, chi, k, H, hw, trivial=False):
    ps = factor_sqfree(N)
    psi = math.prod(l + 1 for l in ps)
    tot = 0j
    r = math.isqrt(n)
    if r * r == n:
        tot += n ** (k / 2 - 1) * chi[r % N] * (k - 1) * psi / 12
    s2 = 0j
    for t in range(-2 * r - 2, 2 * r + 3):
        D = 4 * n - t * t
        if D <= 0: continue
        f = 1
        while f * f <= D:
            if D % (f * f) == 0 and (D // (f * f)) % 4 in (0, 3):
                Nf = math.gcd(N, f)
                mu = 0j
                for x in range(N):
                    ok = True
                    for l in ps:
                        if Nf % l == 0:
                            if (2 * x - t) % l != 0: ok = False; break
                        elif (x * x - t * x + n) % l != 0: ok = False; break
                    if ok: mu += chi[x]
                mu *= math.prod(l + 1 for l in ps if Nf % l == 0)
                s2 += Pk(t, n, k) * hw[D // (f * f)] * mu
            f += 1
    tot += -0.5 * s2
    s3 = 0j
    divsN = [d for d in range(1, N + 1) if N % d == 0]
    for d in range(1, n + 1):
        if n % d: continue
        e = n // d
        inner = 0j
        for tau in divsN:
            M = N // tau
            # y = d mod tau, y = e mod M
            y = next(y for y in range(N) if (y - d) % tau == 0 and (y - e) % M == 0)
            inner += chi[y]
        s3 += min(d, e) ** (k - 1) * inner
    tot += -0.5 * s3
    if trivial and k == 2:
        tot += sum(t for t in range(1, n + 1) if n % t == 0 and math.gcd(N, n // t) == 1)
    return tot

def chars_mod(N):
    """all Dirichlet characters mod N (N odd squarefree) as (table, is_primitive)."""
    ps = factor_sqfree(N)
    from trace_prime import primroot
    gens = [(l, primroot(l)) for l in ps]
    # discrete logs
    logs = []
    for l, g in gens:
        lg = {}; x = 1
        for m in range(l - 1): lg[x] = m; x = x * g % l
        logs.append(lg)
    import itertools
    out = []
    for js in itertools.product(*[range(l - 1) for l in ps]):
        tab = [0j] * N
        for x in range(N):
            v = 1 + 0j
            for (l, g), lg, j in zip(gens, logs, js):
                if x % l == 0: v = 0; break
                v *= cmath.exp(2j * math.pi * j * lg[x % l] / (l - 1))
            tab[x] = v
        out.append((tab, all(j != 0 for j in js)))
    return out

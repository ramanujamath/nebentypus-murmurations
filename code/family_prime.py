"""Root-number-weighted Hecke sums for the family
   F_k(q) = { f newform in S_k(q, chi) : chi mod q nontrivial, chi(-1) = (-1)^k },  q prime.
Every such space is new (chi primitive).  Uses
   eps_f = i^{-k} tau(chi)/sqrt(q) * conj(lambda_f(q))            (verified numerically)
   sum_f a_f(p) conj(a_f(q)) = chi(p) Tr(T_{pq} | S_k(q, chibar))
   Tr(T_{pq} | S_k(q,chibar)) = -1/2 sum_t P_k(t,pq) chibar(t) H(4pq-t^2) - 1 - chibar(p) min(p,q)^{k-1}
and sums over chi in closed form via Gauss sums."""
import numpy as np, math, cmath

def Pk_vec(t, n, k):
    t = t.astype(np.float64); n = float(n)
    u0, u1 = np.zeros_like(t), np.ones_like(t)
    for _ in range(k - 2):
        u0, u1 = u1, t * u1 - n * u0
    return u1

def Gfun(ybar, q, k):
    """sum over chi (nontrivial, parity (-1)^k) of tau(chi)/sqrt(q) * chi(y), given ybar = y^{-1} mod q."""
    e = np.exp(2j * np.pi * ybar / q)
    val = (q - 1) / 2 * (e + (-1) ** k * np.conj(e)) + (1 if k % 2 == 0 else 0)
    return val / math.sqrt(q)

def family_sum(q, p, k, H):
    """sum over chi and f in S_k(q,chi) of eps_f * lambda_f(p)."""
    n = p * q
    r = math.isqrt(4 * n - 1)
    t = np.arange(-r, r + 1)
    t = t[t % q != 0]
    pbar = pow(p, -1, q)
    main = -0.5 * np.sum(Pk_vec(t, n, k) * H[4 * n - t * t] * Gfun((t * pbar) % q, q, k))
    tot = main - Gfun(pbar, q, k) - min(p, q) ** (k - 1) * Gfun(1, q, k)
    return (1j) ** (-k) * tot / (n ** ((k - 1) / 2))

def family_dim(q, k):
    nchi = (q - 1) // 2 - (1 if k % 2 == 0 else 0)
    def C(x):
        x %= q
        return (q - 1) / 2 * ((x == 1) + (-1) ** k * (x == q - 1)) - (1 if k % 2 == 0 else 0)
    Hs = {0: 0.5, 1: 1 / 3, -1: 1 / 3}
    s = 0
    for t in (-1, 0, 1):
        P = Pk_vec(np.array([t]), 1, k)[0]
        s += P * Hs[t] * sum(C(x) for x in range(q) if (x * x - t * x + 1) % q == 0)
    return (k - 1) * (q + 1) / 12 * nchi - 0.5 * s - nchi

def primes_upto(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def family_sum_conj(q, p, k, H):
    """sum over chi and f in S_k(q,chi) of eps_f * conj(lambda_f(p))  (= eps_f * lambda_{fbar}(p))."""
    n = p * q
    r = math.isqrt(4 * n - 1)
    t = np.arange(-r, r + 1)
    t = t[t % q != 0]
    main = -0.5 * np.sum(Pk_vec(t, n, k) * H[4 * n - t * t] * Gfun(t % q, q, k))
    tot = main - Gfun(1, q, k) - min(p, q) ** (k - 1) * Gfun(p % q, q, k)
    return (1j) ** (-k) * tot / (n ** ((k - 1) / 2))

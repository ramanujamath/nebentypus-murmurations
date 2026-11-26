"""Hurwitz class numbers H(m), m <= M, by enumerating reduced binary quadratic forms.
Usage: python hurwitz_sieve.py 1e7   (writes H_10000000.npy, about 25 s; needed by most scripts)."""
import numpy as np, sys, time
def hurwitz(M, dtype=np.float32):
    H = np.zeros(M + 1, dtype=dtype)
    a = 1
    while 3 * a * a <= M:
        s = 4 * a
        for b in range(-a + 1, a + 1):
            c0 = a if b >= 0 else a + 1
            st = 4 * a * c0 - b * b
            if st <= M:
                H[st::s] += 1
        if 4 * a * a <= M: H[4 * a * a] -= 0.5
        H[3 * a * a] -= 2/3
        a += 1
    return H
if __name__ == "__main__":
    M = int(float(sys.argv[1])); t = time.time()
    H = hurwitz(M); print("time", time.time() - t)
    np.save(f"H_{M}.npy", H)
    print(H[[3,4,23,1000003 % (M+1)]])

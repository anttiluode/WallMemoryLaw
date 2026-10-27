"""The calibration law: geometry -> interaction onsets -> window degree -> calibration set.

A written patch is a local scatterer.  In the multiple-scattering (Born / Neumann)
series of the cavity, the k-patch interaction term of a wall state collects only the
paths that scatter off every one of those k patches, so it is zero until the shortest
such path can reach the receiver:

    onset(S) ~ tour(S) / STEP,  tour(S) = min over orders and wall points of
               |source - p1| + |p1 - p2| + ... + |pk - receiver|,  p_i on patch S_i

(straight chords: the cavity is convex, and any extra reflection only adds length).
Inside a listening window that ends at t1, every interaction whose onset is later is
absent, so the echo is a polynomial in the wall bits of degree

    d(t1) = max { k : min_{|S|=k} tour(S)/STEP - margin <= t1 }

and the whole memory is pinned down by the sum_{k<=d} C(K,k) states with at most d
written patches, instead of 2^K.
"""
import itertools
import numpy as np

from .cavity import STEP


# ---------------------------------------------------------------- tours
def tours(cav, source=0, max_points=None):
    """Shortest scattering tour (in steps) for EVERY subset of patches (Held-Karp over
    patches, with each patch represented by up to `max_points` of its wall cells,
    evenly spaced).  Returns array T[mask] in steps (T[0] = direct path)."""
    K = cav.K
    S = np.array(cav.sources[source], float)
    R = np.array(cav.receiver, float)
    pts, owner = [], []
    for k in range(K):
        p = cav.patch_points(k)
        if max_points and len(p) > max_points:
            # order along the wall by elliptic angle, keep evenly spaced cells incl. both ends
            order = np.argsort(np.arctan2((p[:, 0] - 47.5) / 27.0, (p[:, 1] - 47.5) / 36.0))
            p = p[order][np.round(np.linspace(0, len(p) - 1, max_points)).astype(int)]
        pts.append(p)
        owner += [k] * len(p)
    P = np.concatenate(pts)
    owner = np.array(owner)
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    dS = np.linalg.norm(P - S, axis=1)
    dR = np.linalg.norm(P - R, axis=1)
    M = 1 << K
    INF = np.inf
    dp = np.full((M, len(P)), INF)                           # min length source -> ... ending at point
    for i in range(len(P)):
        dp[1 << owner[i], i] = dS[i]
    bits_of = np.array([[(m >> k) & 1 for k in range(K)] for m in range(M)], dtype=bool)
    popc = bits_of.sum(1)
    for size in range(1, K):
        masks = np.flatnonzero(popc == size)
        cur = dp[masks]                                      # (m, P)
        for q in range(K):
            sel = ~bits_of[masks, q]
            if not sel.any():
                continue
            mq = masks[sel]
            cols = np.flatnonzero(owner == q)
            cand = (cur[sel][:, :, None] + D[None, :, cols]).min(1)   # (m, |cols|)
            tgt = mq | (1 << q)
            dp[tgt[:, None], cols[None, :]] = np.minimum(dp[tgt[:, None], cols[None, :]], cand)
    T = np.full(M, INF)
    T[0] = np.linalg.norm(S - R)
    for m in range(1, M):
        T[m] = (dp[m] + dR).min()
    return T / STEP


def min_tour_by_order(T, K):
    popc = np.array([bin(m).count("1") for m in range(1 << K)])
    return {k: float(T[popc == k].min()) for k in range(K + 1)}


def predicted_degree(T, K, t1, margin):
    m = min_tour_by_order(T, K)
    return max(k for k in range(K + 1) if m[k] - margin <= t1)


# ---------------------------------------------------------------- interactions (Moebius)
def codes(K):
    return ((np.arange(1 << K)[:, None] >> np.arange(K)) & 1).astype(np.int8)


def low_order_masks(K, d):
    return [m for m in range(1 << K) if bin(m).count("1") <= d]


def mobius_from_table(Y, K):
    """Full Moebius transform of a 2^K table (rows indexed by mask): interaction
    coefficients c[S] = sum_{T subset S} (-1)^{|S|-|T|} y[T]."""
    c = np.array(Y, dtype=np.float64, copy=True).reshape(1 << K, -1)
    for i in range(K):
        v = c.reshape(1 << (K - i - 1), 2, 1 << i, -1)
        v[:, 1] -= v[:, 0]
    return c


def zeta(c, K):
    """Inverse of mobius_from_table: y[x] = sum_{S subset x} c[S]."""
    y = np.array(c, copy=True).reshape(1 << K, -1)
    for i in range(K):
        v = y.reshape(1 << (K - i - 1), 2, 1 << i, -1)
        v[:, 1] += v[:, 0]
    return y


def coefficients_from_calibration(cal, K, d):
    """cal: dict mask -> response vector, containing every mask with <= d bits.
    Returns coefficient table (2^K, L) with c[S] = 0 for |S| > d (exact Moebius
    inversion restricted to the calibrated subsets)."""
    L = len(next(iter(cal.values())))
    c = np.zeros((1 << K, L))
    for m in low_order_masks(K, d):
        sub = m
        acc = np.zeros(L)
        k = bin(m).count("1")
        while True:                                          # all submasks of m
            acc += (-1) ** (k - bin(sub).count("1")) * cal[sub]
            if sub == 0:
                break
            sub = (sub - 1) & m
        c[m] = acc
    return c


def decode(model_table, obs, chunk=4096):
    """Nearest model response (exhaustive over all 2^K codes). obs: (n, L)."""
    best = np.full(len(obs), np.inf)
    arg = np.zeros(len(obs), dtype=np.int64)
    o2 = (obs ** 2).sum(1)
    for s in range(0, len(model_table), chunk):
        Mt = model_table[s:s + chunk]
        d = o2[:, None] - 2 * obs @ Mt.T + (Mt ** 2).sum(1)[None, :]
        j = d.argmin(1)
        v = d[np.arange(len(obs)), j]
        upd = v < best
        best[upd] = v[upd]
        arg[upd] = s + j[upd]
    return arg

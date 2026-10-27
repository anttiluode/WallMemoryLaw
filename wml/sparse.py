"""The per-subset law and the amplitude model, for K too large for any 2^K table.

active sets   S is *active* in window [t0, t1) if its own scattering tour arrives in
              time:  tour(S) - MARGIN <= t1 for at least one source.  A superset can
              never have a shorter tour, so active sets are enumerated level by level
              (Apriori): a k-set is a candidate only if all its (k-1)-subsets are active.
amplitude     log10(|I_S| / one-bit) ~ A0 + A1 * lead - A2 * |S|,
              lead = t1 - min_source tour(S).  Fit on the K=16 ellipse (discovery):
              A0 = -0.42, A1 = 0.0992 per step, A2 = 0.75 per patch, r = 0.964.
pruned sets   singles always, plus active sets with predicted size >= EPS.
model         y(x) = sum over kept S subset x of c_S, c_S by Moebius inversion from the
              calibration states (the kept sets are closed under taking subsets).
decode        descent over 1- and 2-bit flips from 18 starts (least-squares rounding of
              the single-patch terms, empty wall, 16 random walls).  A miss is a *search*
              failure if the true wall fits the model better than the found one, otherwise
              a *model* failure; both are reported.  On the K=16 ellipse (discovery) it
              matched exhaustive decoding on all 300 held-out walls (with 32 restarts; 16 is
              frozen for the stadium to fit the compute budget, misses are classified).
"""
import itertools

import numpy as np

from .cavity import STEP

MARGIN = 40.0
A0, A1, A2 = -0.42, 0.0992, 0.75
EPS = 1e-3


def _points(cav):
    return [cav.patch_points(k) for k in range(cav.K)]


def subset_tour(cav, subset, source, pts=None):
    """Exact shortest tour (steps) source -> every patch of `subset` (any order) -> receiver."""
    pts = pts or _points(cav)
    S = np.array(cav.sources[source], float)
    R = np.array(cav.receiver, float)
    sub = list(subset)
    k = len(sub)
    P = [pts[q] for q in sub]
    # Held-Karp inside the subset: dp[mask][j] = min cost ending at each point of patch j
    dp = {}
    for j in range(k):
        dp[(1 << j, j)] = np.linalg.norm(P[j] - S, axis=1)
    for mask in range(1, 1 << k):
        for j in range(k):
            if (mask, j) not in dp:
                continue
            cur = dp[(mask, j)]
            for q in range(k):
                if mask >> q & 1:
                    continue
                D = np.linalg.norm(P[j][:, None, :] - P[q][None, :, :], axis=2)
                cand = (cur[:, None] + D).min(0)
                key = (mask | 1 << q, q)
                dp[key] = np.minimum(dp[key], cand) if key in dp else cand
    full = (1 << k) - 1
    best = min((dp[(full, j)] + np.linalg.norm(P[j] - R, axis=1)).min() for j in range(k))
    return best / STEP


def active_sets(cav, t1, margin=MARGIN, max_order=8):
    """{frozenset: min-over-sources tour} for every set (|S| >= 1) whose tour arrives in time."""
    pts = _points(cav)
    K = cav.K
    tour = {}
    level = []
    for q in range(K):
        t = min(subset_tour(cav, (q,), s, pts) for s in range(len(cav.sources)))
        tour[frozenset([q])] = t
        if t - margin <= t1:
            level.append(frozenset([q]))
    active = {s: tour[s] for s in level}
    for k in range(2, max_order + 1):
        prev = set(level)
        cands = set()
        for a in level:
            for q in range(K):
                if q in a:
                    continue
                c = a | {q}
                if all((c - {r}) in prev for r in c):
                    cands.add(c)
        level = []
        for c in cands:
            t = min(subset_tour(cav, tuple(sorted(c)), s, pts) for s in range(len(cav.sources)))
            if t - margin <= t1:
                level.append(c)
                active[c] = t
        if not level:
            break
    return active


def predicted_log_size(tour_steps, order, t1):
    return A0 + A1 * (t1 - tour_steps) - A2 * order


def pruned_sets(active, t1, eps=EPS):
    keep = {}
    for s, t in active.items():
        if len(s) == 1 or predicted_log_size(t, len(s), t1) >= np.log10(eps):
            keep[s] = t
    return keep


def closure(sets):
    """All subsets (incl. empty) of the given sets, as int masks."""
    out = {0}
    for s in sets:
        m = sum(1 << q for q in s)
        sub = m
        while True:
            out.add(sub)
            if sub == 0:
                break
            sub = (sub - 1) & m
    return sorted(out)


class SparseMemory:
    """y(x) = sum_{S in kept, S subset x} c_S."""

    def __init__(self, K, table):
        """table: dict int-mask -> feature vector, closed under subsets."""
        self.K = K
        self.masks = np.array(sorted(table), dtype=np.int64)
        coef = []
        for m in self.masks:
            k = bin(int(m)).count("1")
            acc = 0.0
            sub = int(m)
            while True:
                acc = acc + (-1) ** (k - bin(sub).count("1")) * table[sub]
                if sub == 0:
                    break
                sub = (sub - 1) & int(m)
            coef.append(acc)
        self.coef = np.array(coef)

    def predict(self, X):
        """X: (n,) int masks -> (n, L).  Sparse: only kept sets inside each x contribute."""
        from scipy.sparse import csr_matrix
        X = np.atleast_1d(np.asarray(X, dtype=np.int64))
        r, c = np.nonzero((X[:, None] & self.masks[None, :]) == self.masks[None, :])
        inc = csr_matrix((np.ones(len(r)), (r, c)), shape=(len(X), len(self.masks)))
        return np.asarray(inc @ self.coef)

    def sizes(self):
        one = np.median([np.linalg.norm(c) for m, c in zip(self.masks, self.coef) if bin(int(m)).count("1") == 1])
        return {int(m): float(np.linalg.norm(c) / one) for m, c in zip(self.masks, self.coef)}

    def decode(self, obs, restarts=16, seed=0):
        """Descent over 1- and 2-bit flips from several starts (least-squares rounding of
        the single-patch terms, empty wall, random walls).  Returns (best_mask, best_cost)."""
        rng = np.random.default_rng(seed)
        K = self.K
        flips = (1 << np.arange(K)).astype(np.int64)
        pairs = np.array([(1 << a) | (1 << b) for a in range(K) for b in range(a + 1, K)], dtype=np.int64)
        moves = np.concatenate([flips, pairs])
        y0 = self.predict([0])[0]
        single = np.array([self.predict([1 << q])[0] - y0 for q in range(K)])
        w = np.linalg.lstsq(single.T, obs - y0, rcond=None)[0]
        starts = [int(sum(1 << q for q in range(K) if w[q] > 0.5)), 0]
        starts += [int(rng.integers(0, 1 << K)) for _ in range(restarts)]
        best, best_c = None, np.inf
        for x in starts:
            c = float(((self.predict([x])[0] - obs) ** 2).sum())
            while True:
                nb = x ^ moves
                cs = ((self.predict(nb) - obs) ** 2).sum(1)
                j = int(cs.argmin())
                if cs[j] < c - 1e-15:
                    x, c = int(nb[j]), float(cs[j])
                else:
                    break
            if c < best_c:
                best, best_c = x, c
        return best, best_c

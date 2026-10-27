"""Confirmatory run 2: new geometry, K = 24, the per-subset law and the amplitude model.
See PREREG_STADIUM.md.  Nothing here was tuned on the stadium.

    python confirm_stadium.py plan        # geometry only: window, active and pruned sets
    python confirm_stadium.py calibrate   # every state in the tour-law closure, 3 sources
    python confirm_stadium.py heldout     # 300 random 24-bit walls
    python confirm_stadium.py inactive    # direct inclusion-exclusion on tour-INACTIVE sets
    python confirm_stadium.py decide      # -> results/STADIUM.json
"""
import json
import sys
import time

import numpy as np

from wml import Cavity
from wml.sparse import (active_sets, pruned_sets, closure, subset_tour, predicted_log_size,
                        SparseMemory, MARGIN, EPS, A0, A1, A2)

K, FRACTION, SHAPE = 24, 0.5, "stadium"
HELD, HELD_SEED = 300, 2027
INACTIVE_SEED = 11
NOISE = 0.04
BATCH = 64
OUT = "results/stadium_"


def cav():
    return Cavity(K, FRACTION, SHAPE)


def bits_of(masks):
    m = np.asarray(masks, dtype=np.int64)[:, None]
    return ((m >> np.arange(K)) & 1).astype(float)


def plan():
    c = cav()
    pts = [c.patch_points(k) for k in range(K)]
    single = np.array([[subset_tour(c, (q,), s, pts) for s in range(3)] for q in range(K)]).min(1)
    t0 = int(np.floor(single.min() - MARGIN))
    t1 = int(np.ceil(single.max() + 20))
    act = active_sets(c, t1, max_order=12)
    keep = pruned_sets(act, t1)
    P = {"t0": t0, "t1": t1, "single_tours": single.tolist(),
         "active": {",".join(map(str, sorted(s))): t for s, t in act.items()},
         "pruned": sorted(",".join(map(str, sorted(s))) for s in keep),
         "tour_law_states": len(closure(act)), "pruned_states": len(closure(keep))}
    json.dump(P, open(OUT + "plan.json", "w"))
    print({k: v for k, v in P.items() if k in ("t0", "t1", "tour_law_states", "pruned_states")},
          "active", len(act), "pruned", len(keep))


def load_plan():
    P = json.load(open(OUT + "plan.json"))
    act = {frozenset(map(int, k.split(","))): v for k, v in P["active"].items()}
    keep = {frozenset(map(int, k.split(","))) for k in P["pruned"]}
    return P, act, keep


def simulate(masks, t0, t1, tag):
    c = cav()
    out = []
    t = time.time()
    for s in range(3):
        rows = []
        for i in range(0, len(masks), BATCH):
            rows.append(c.run(bits_of(masks[i:i + BATCH]), source=s, steps=t1)[:, t0:t1, 0])
            if (i // BATCH) % 25 == 0:
                print(f"{tag} source {s}: {i + len(rows[-1])}/{len(masks)}  {time.time() - t:.0f}s", flush=True)
        out.append(np.concatenate(rows))
    return np.concatenate(out, axis=1)                         # (n, 3*(t1-t0)) float64


def calibrate():
    P, act, keep = load_plan()
    masks = closure(act)
    F = simulate(masks, P["t0"], P["t1"], "calibrate")
    np.savez_compressed(OUT + "calibration.npz", masks=np.array(masks, dtype=np.int64), F=F)


def heldout():
    P, _, _ = load_plan()
    rng = np.random.default_rng(HELD_SEED)
    masks = rng.integers(0, 1 << K, HELD, dtype=np.int64)
    F = simulate(masks, P["t0"], P["t1"], "heldout")
    np.savez_compressed(OUT + "heldout.npz", masks=masks, F=F)


def inactive():
    """All tour-inactive pairs and 24 random tour-inactive triples whose 2-subsets are all
    active (the inactive sets closest to the boundary), measured directly in float64."""
    P, act, _ = load_plan()
    import itertools
    pairs = [frozenset(p) for p in itertools.combinations(range(K), 2) if frozenset(p) not in act]
    tri = [frozenset(t) for t in itertools.combinations(range(K), 3)
           if frozenset(t) not in act and all(frozenset(p) in act for p in itertools.combinations(t, 2))]
    rng = np.random.default_rng(INACTIVE_SEED)
    tri = [tri[i] for i in rng.choice(len(tri), min(24, len(tri)), replace=False)] if tri else []
    sets = pairs + tri
    rows = []
    one = None
    for S in sets:
        sub = sorted(S)
        ms, sg = [], []
        for r in range(len(sub) + 1):
            for on in itertools.combinations(sub, r):
                ms.append(sum(1 << q for q in on))
                sg.append((-1) ** (len(sub) - r))
        F = simulate(ms, P["t0"], P["t1"], "inactive")
        I = (np.array(sg)[:, None] * F).sum(0)
        rows.append({"set": sub, "size_abs": float(np.linalg.norm(I))})
    json.dump(rows, open(OUT + "inactive.json", "w"))
    print(len(pairs), "inactive pairs,", len(tri), "inactive triples measured")


def decide():
    P, act, keep = load_plan()
    cal = np.load(OUT + "calibration.npz")
    ho = np.load(OUT + "heldout.npz")
    table = {int(m): f for m, f in zip(cal["masks"], cal["F"])}
    Fh, hm = ho["F"], ho["masks"]
    one_v = np.array([table[1 << q] - table[0] for q in range(K)])
    one = float(np.median(np.linalg.norm(one_v, axis=1)))
    sigma = NOISE * one / np.sqrt(Fh.shape[1])
    noisy = Fh + np.random.default_rng(99).normal(0, sigma, Fh.shape)
    out = {"K": K, "shape": SHAPE, "window": [P["t0"], P["t1"]], "all_states": 1 << K,
           "frozen": {"margin": MARGIN, "eps": EPS, "A0": A0, "A1": A1, "A2": A2},
           "heldout": HELD, "features": int(Fh.shape[1]), "models": {}}
    models = {"tour_law": closure(act), "pruned": closure(keep),
              "tour_law_order2": closure([s for s in act if len(s) <= 2]),
              "additive": closure([s for s in act if len(s) == 1])}
    for name, masks in models.items():
        mem = SparseMemory(K, {m: table[m] for m in masks})
        t = time.time()
        r0 = [mem.decode(o) for o in Fh]
        rn = [mem.decode(o) for o in noisy]
        found0 = np.array([r[0] for r in r0])
        foundn = np.array([r[0] for r in rn])
        true0 = ((mem.predict(hm) - Fh) ** 2).sum(1)
        cost0 = np.array([r[1] for r in r0])
        miss = found0 != hm
        err = np.linalg.norm(Fh - mem.predict(hm), axis=1) / np.linalg.norm(Fh - table[0], axis=1)
        out["models"][name] = {"calibration_states": len(masks), "exact": int((~miss).sum()),
                               "exact_at_noise": int((foundn == hm).sum()),
                               "search_failures": int((miss & (cost0 > true0 + 1e-12)).sum()),
                               "model_failures": int((miss & (cost0 <= true0 + 1e-12)).sum()),
                               "median_rel_model_error": float(np.median(err))}
        print(name, out["models"][name], f"{time.time() - t:.0f}s", flush=True)
    # amplitude model vs measured sizes on every active set of order >= 2
    full = SparseMemory(K, {m: table[m] for m in models["tour_law"]})
    sz = {m: float(np.linalg.norm(c) / one) for m, c in zip(full.masks.tolist(), full.coef)}
    meas, pred, pruned_away = [], [], []
    keepm = {sum(1 << q for q in s) for s in keep}
    for s, tt in act.items():
        if len(s) < 2:
            continue
        m = sum(1 << q for q in s)
        lp = predicted_log_size(tt, len(s), P["t1"])
        if sz[m] > 1e-5:
            meas.append(np.log10(sz[m]))
            pred.append(lp)
        if m not in keepm:
            pruned_away.append(sz[m])
    meas, pred = np.array(meas), np.array(pred)
    out["amplitude"] = {"n": int(len(meas)), "r": float(np.corrcoef(meas, pred)[0, 1]),
                        "median_bias_dex": float(np.median(meas - pred)),
                        "mad_dex": float(np.median(np.abs(meas - pred))),
                        "pruned_away": len(pruned_away), "max_pruned_away_size": float(max(pruned_away)),
                        "pruned_away_above_1e-2": int(sum(v > 1e-2 for v in pruned_away))}
    ina = json.load(open(OUT + "inactive.json"))
    out["inactive"] = {"n": len(ina), "max_rel_size": float(max(r["size_abs"] for r in ina) / one)}
    m = out["models"]
    out["decisions"] = {
        "N1_tour_law_exact": bool(m["tour_law"]["exact"] >= 0.99 * HELD),
        "N2_pruned_exact": bool(m["pruned"]["exact"] >= 0.99 * HELD),
        "N3_inactive_negligible": bool(out["inactive"]["max_rel_size"] <= 1e-5),
        "N4a_amplitude_transfers": bool(out["amplitude"]["r"] >= 0.7),
        "N4b_pruning_safe": bool(out["amplitude"]["pruned_away_above_1e-2"] == 0),
    }
    json.dump(out, open("results/STADIUM.json", "w"), indent=1)
    print(json.dumps(out["decisions"], indent=1), json.dumps(out["amplitude"]), json.dumps(out["inactive"]))


if __name__ == "__main__":
    {"plan": plan, "calibrate": calibrate, "heldout": heldout, "inactive": inactive, "decide": decide}[sys.argv[1]]()

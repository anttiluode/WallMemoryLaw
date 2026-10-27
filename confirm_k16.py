"""Confirmatory run at K=16 (65,536 wall states; never simulated in full).  See PREREG.md.

    python confirm_k16.py calibrate   # all states with <= SAFE_DEGREE written patches, 3 sources
    python confirm_k16.py heldout     # HELD random states with >= SAFE_DEGREE+1 patches
    python confirm_k16.py onsets      # direct inclusion-exclusion on sampled subsets (P2)
    python confirm_k16.py decide      # decisions -> results/K16.json
"""
import itertools
import json
import os
import sys
import time

import numpy as np

from wml import Cavity, tours, min_tour_by_order, predicted_degree, zeta, coefficients_from_calibration, decode

K, FRACTION = 16, 0.5
T0, T1 = 145, 230             # listening window [T0, T1)
MARGIN = 40.0                 # frozen at K=8 discovery
SAFE_DEGREE = 6               # predicted from geometry (PREREG.md)
HELD, HELD_SEED = 300, 2026
NOISE = 0.04                  # Kavitaatio's definition: 4% of median one-bit response RMS
ONSET_SAMPLES, ONSET_SEED, ONSET_STEPS = 10, 7, 600
BATCH = 64


def masks_upto(d):
    return [m for m in range(1 << K) if bin(m).count("1") <= d]


def bits_of(masks):
    m = np.asarray(masks)[:, None]
    return ((m >> np.arange(K)) & 1).astype(float)


def simulate(cav, masks, steps, tag):
    out = []
    t = time.time()
    for s in range(3):
        rows = []
        for i in range(0, len(masks), BATCH):
            rows.append(cav.run(bits_of(masks[i:i + BATCH]), source=s, steps=steps)[:, :, 0].astype(np.float32))
            if (i // BATCH) % 40 == 0:
                print(f"{tag} source {s}: {i + len(rows[-1])}/{len(masks)}  {time.time() - t:.0f}s", flush=True)
        out.append(np.concatenate(rows))
    return np.stack(out)                                     # (3, n, steps)


def features(Y):
    """(3, n, steps) -> (n, 3*(T1-T0))"""
    return np.concatenate([Y[s][:, T0:T1] for s in range(3)], axis=1).astype(np.float64)


def calibrate():
    cav = Cavity(K, FRACTION)
    masks = masks_upto(SAFE_DEGREE)
    Y = simulate(cav, masks, T1, "calibrate")
    np.savez_compressed("results/k16_calibration.npz", masks=np.array(masks), Y=Y)


def heldout():
    cav = Cavity(K, FRACTION)
    rng = np.random.default_rng(HELD_SEED)
    pool = [m for m in range(1 << K) if bin(m).count("1") > SAFE_DEGREE]
    masks = sorted(rng.choice(pool, HELD, replace=False).tolist())
    Y = simulate(cav, masks, T1, "heldout")
    np.savez_compressed("results/k16_heldout.npz", masks=np.array(masks), Y=Y)


def onsets():
    cav = Cavity(K, FRACTION)
    T = tours(cav, source=1)
    rng = np.random.default_rng(ONSET_SEED)
    rows = []
    for k in range(1, SAFE_DEGREE + 1):
        subs = [tuple(sorted(rng.choice(K, k, replace=False).tolist())) for _ in range(ONSET_SAMPLES)]
        for sub in subs:
            subm = [sum(1 << i for i in on) for r in range(k + 1) for on in itertools.combinations(sub, r)]
            sign = [(-1) ** (k - r) for r in range(k + 1) for on in itertools.combinations(sub, r)]
            Y = cav.run(bits_of(subm), source=1, steps=ONSET_STEPS)[:, :, 0]
            I = (np.array(sign)[:, None] * Y).sum(0)
            env = np.abs(I)
            pk = env.max()
            on = int(np.argmax(env > 1e-6 * pk)) if pk > 1e-12 else None
            m = sum(1 << i for i in sub)
            rows.append({"subset": list(sub), "k": k, "onset": on, "tour": float(T[m])})
        print("onsets order", k, "done", flush=True)
    json.dump(rows, open("results/k16_onsets.json", "w"))


def decide():
    from scipy.stats import pearsonr
    cal = np.load("results/k16_calibration.npz")
    ho = np.load("results/k16_heldout.npz")
    Fc, Fh = features(cal["Y"]), features(ho["Y"])
    cmask, hmask = cal["masks"].tolist(), ho["masks"]
    table = {m: Fc[i] for i, m in enumerate(cmask)}
    one = np.array([table[1 << k] - table[0] for k in range(K)])
    sigma = NOISE * np.median(np.linalg.norm(one, axis=1)) / np.sqrt(Fc.shape[1])
    rng = np.random.default_rng(99)
    noisy = Fh + rng.normal(0, sigma, Fh.shape)
    out = {"K": K, "fraction": FRACTION, "window": [T0, T1], "margin": MARGIN, "safe_degree": SAFE_DEGREE,
           "calibration_states": len(cmask), "all_states": 1 << K, "heldout": int(len(hmask)),
           "features": int(Fc.shape[1]), "noise_sigma": float(sigma), "by_order": {}}
    for d in range(1, SAFE_DEGREE + 1):
        sub = {m: v for m, v in table.items() if bin(m).count("1") <= d}
        model = zeta(coefficients_from_calibration(sub, K, d), K)
        err = np.linalg.norm(Fh - model[hmask], axis=1) / np.linalg.norm(Fh - table[0], axis=1)
        p0 = decode(model, Fh)
        pn = decode(model, noisy)
        wrong_bits = [bin(int(a) ^ int(b)).count("1") for a, b in zip(p0, hmask)]
        out["by_order"][d] = {"calibration_states": len(sub), "exact_noise_free": int((p0 == hmask).sum()),
                              "exact_at_noise": int((pn == hmask).sum()),
                              "median_rel_model_error": float(np.median(err)),
                              "mean_wrong_bits_noise_free": float(np.mean(wrong_bits))}
        print(d, out["by_order"][d], flush=True)
    needed = next((d for d in range(1, SAFE_DEGREE + 1)
                   if out["by_order"][d]["exact_noise_free"] >= 0.99 * len(hmask)), None)
    out["needed_order"] = needed
    rows = json.load(open("results/k16_onsets.json"))
    on = np.array([r["onset"] for r in rows], float)
    tu = np.array([r["tour"] for r in rows])
    out["onsets"] = {"n": len(rows), "corr": float(pearsonr(on, tu)[0]), "max_early": float((tu - on).max()),
                     "median_early": float(np.median(tu - on))}
    s = out["by_order"][SAFE_DEGREE]
    out["decisions"] = {
        "P1_safe_degree_is_exact": bool(s["exact_noise_free"] >= 0.99 * len(hmask) and s["median_rel_model_error"] <= 1e-3),
        "P2_onsets_follow_tours": bool(out["onsets"]["corr"] >= 0.95 and out["onsets"]["max_early"] <= MARGIN),
        "P3_needed_order_in_range": bool(needed is not None and SAFE_DEGREE - 2 <= needed <= SAFE_DEGREE),
    }
    json.dump(out, open("results/K16.json", "w"), indent=1)
    print(json.dumps(out["decisions"], indent=1), json.dumps(out["onsets"]))


if __name__ == "__main__":
    {"calibrate": calibrate, "heldout": heldout, "onsets": onsets, "decide": decide}[sys.argv[1]]()

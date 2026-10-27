"""Discovery at K=8, where every one of the 256 wall states can be simulated.

1. All 256 states x 3 sources x 1000 steps at the receiver.
2. Every interaction (Moebius transform of the table) -> onset step (first |I| > 1e-6 of its peak).
3. Tours from geometry (Held-Karp, every wall cell) -> margin = max(tour - onset), rounded up.
4. Window sweep: for windows [145, t1], predicted degree vs the lowest calibration order
   that decodes every held-out state.
Writes results/discovery_k8.json and results/k8_tables.npz.
"""
import json
import time

import numpy as np

from wml import Cavity, tours, min_tour_by_order, predicted_degree, codes, mobius_from_table, \
    coefficients_from_calibration, zeta, decode

K, STEPS, T0 = 8, 1000, 145
WINDOWS = [245, 295, 345, 395, 445, 545, 745]


def onset(trace, rel=1e-6):
    env = np.abs(trace)
    pk = env.max()
    return int(np.argmax(env > rel * pk)) if pk > 1e-12 else None


def main():
    cav = Cavity(K, 1.0)
    C = codes(K)
    t = time.time()
    Y = np.stack([cav.run(C, source=s, steps=STEPS)[:, :, 0] for s in range(3)])   # (3, 256, STEPS)
    print(f"simulated 3 x 256 x {STEPS} in {time.time() - t:.0f}s", flush=True)
    np.savez_compressed("results/k8_tables.npz", Y=Y.astype(np.float32))
    out = {"K": K, "steps": STEPS, "sources": cav.sources, "receiver": cav.receiver}

    # interactions and onsets, per source
    rows = []
    for s in range(3):
        I = mobius_from_table(Y[s], K)                      # (256, STEPS)
        T = tours(cav, source=s)
        for m in range(1, 1 << K):
            k = bin(m).count("1")
            rows.append({"source": s, "mask": m, "k": k, "onset": onset(I[m]), "tour": float(T[m]),
                         "rms_late": float(np.sqrt(np.mean(I[m, 700:] ** 2)))})
        out[f"min_tour_by_order_src{s}"] = min_tour_by_order(T, K)
    on = np.array([r["onset"] for r in rows], float)
    tu = np.array([r["tour"] for r in rows])
    early = tu - on
    out["onset_vs_tour"] = {"n": len(rows), "corr": float(np.corrcoef(on, tu)[0, 1]),
                            "slope": float(np.polyfit(tu, on, 1)[0]),
                            "max_early_steps": float(early.max()), "median_early_steps": float(np.median(early)),
                            "by_order": {k: {"median_onset": float(np.median(on[[r['k'] == k for r in rows]])),
                                             "median_tour": float(np.median(tu[[r['k'] == k for r in rows]])),
                                             "max_early": float(early[[r['k'] == k for r in rows]].max())}
                                         for k in range(1, K + 1)}}
    margin = float(np.ceil(early.max() / 5) * 5)
    out["margin_frozen"] = margin
    print("onset vs tour:", out["onset_vs_tour"]["corr"], "max early", early.max(), "-> margin", margin)

    # window sweep: predicted degree (worst over the 3 sources) vs needed calibration order
    Ts = [tours(cav, source=s) for s in range(3)]
    sweep = []
    for t1 in WINDOWS:
        dhat = max(predicted_degree(T, K, t1, margin) for T in Ts)
        F = np.concatenate([Y[s][:, T0:t1] for s in range(3)], axis=1)   # (256, 3*(t1-T0))
        res = {}
        for d in range(1, K + 1):
            cal = {m: F[m] for m in range(1 << K) if bin(m).count("1") <= d}
            model = zeta(coefficients_from_calibration(cal, K, d), K)
            held = [m for m in range(1 << K) if bin(m).count("1") > d]
            if not held:
                res[d] = {"exact": None}
                continue
            pred = decode(model, F[held])
            err = np.median(np.linalg.norm(F[held] - model[held], axis=1) / np.linalg.norm(F[held] - F[0], axis=1))
            res[d] = {"exact": int((pred == np.array(held)).sum()), "held": len(held), "median_rel_err": float(err)}
        needed = next((d for d in range(1, K + 1) if res[d]["exact"] is None or res[d]["exact"] >= 0.99 * res[d]["held"]), K)
        sweep.append({"t1": t1, "predicted_degree": dhat, "needed_order": needed, "by_order": res,
                      "calibration_states_at_prediction": int(sum(1 for m in range(1 << K) if bin(m).count('1') <= dhat))})
        print(f"window [145,{t1}]: predicted degree {dhat}, needed order {needed}  "
              + "  ".join(f"d{d}:{v['exact']}/{v.get('held', '-')}" for d, v in res.items() if d <= dhat + 1))
    out["window_sweep"] = sweep
    json.dump(out, open("results/discovery_k8.json", "w"), indent=1)
    json.dump(rows, open("results/k8_onsets.json", "w"))


if __name__ == "__main__":
    main()

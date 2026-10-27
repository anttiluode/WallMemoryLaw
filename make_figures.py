"""Figures for the README."""
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wml import Cavity

ORDER_COL = plt.cm.viridis(np.linspace(0.05, 0.95, 8))


def fig_onsets():
    rows = json.load(open("results/k8_onsets.json"))
    k16 = json.load(open("results/k16_onsets.json"))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for k in range(1, 9):
        r = [x for x in rows if x["k"] == k]
        ax[0].scatter([x["tour"] for x in r], [x["onset"] for x in r], s=9, color=ORDER_COL[k - 1], label=f"k={k}")
    for k in range(1, 7):
        r = [x for x in k16 if x["k"] == k]
        ax[1].scatter([x["tour"] for x in r], [x["onset"] for x in r], s=18, color=ORDER_COL[k - 1], label=f"k={k}")
    for a, t in zip(ax, ["K = 8: all 765 interactions (3 sources)", "K = 16: 60 sampled interactions (source 1)"]):
        lo, hi = a.get_xlim()
        xs = np.linspace(lo, hi, 2)
        a.plot(xs, xs, "k-", lw=.8, label="onset = tour")
        a.plot(xs, xs - 40, "k:", lw=.8, label="tour − 40 (margin)")
        a.set_xlabel("shortest scattering tour / wave speed (steps)")
        a.set_ylabel("measured onset, |I| > 1e-6 of peak (steps)")
        a.set_title(t, fontsize=11)
        a.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig("figures/onsets.png", dpi=110)
    plt.close(fig)


def fig_geometry():
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for a, (K, f) in zip(ax, [(8, 1.0), (16, 0.5)]):
        cav = Cavity(K, f)
        img = np.where(cav.inside, 0.15, 0.0)
        rgb = np.dstack([img] * 3)
        for k in range(K):
            m = cav.counts[k] > 0
            rgb[m] = plt.cm.hsv(k / K)[:3]
        a.imshow(rgb, origin="upper")
        for (y, x), c, lab in [(cav.sources[0], "c", "sources"), (cav.receiver, "w", "receiver")]:
            a.plot(x, y, "o", color=c, ms=4)
        for y, x in cav.sources[1:]:
            a.plot(x, y, "o", color="c", ms=3)
        a.set_title(f"K = {K}" + (" (full sectors, Kavitaatio)" if f == 1 else " (half sectors, gaps)"), fontsize=11)
        a.axis("off")
    fig.tight_layout()
    fig.savefig("figures/geometry.png", dpi=110)
    plt.close(fig)


def fig_k16():
    R = json.load(open("results/K16.json"))
    D = json.load(open("results/discovery_k8.json"))
    orders = sorted(R["by_order"], key=int)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ex = [R["by_order"][o]["exact_noise_free"] / R["heldout"] for o in orders]
    exn = [R["by_order"][o]["exact_at_noise"] / R["heldout"] for o in orders]
    st = [R["by_order"][o]["calibration_states"] for o in orders]
    x = np.arange(len(orders))
    ax[0].bar(x - .2, ex, .4, label="noise-free", color="#e0a030")
    ax[0].bar(x + .2, exn, .4, label="at Kavitaatio noise", color="#4c78a8")
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([f"order {o}\n{s:,} states" for o, s in zip(orders, st)], fontsize=8)
    ax[0].axvline(x[-1] + .5, color="k", lw=.6)
    ax[0].set_ylabel(f"held-out exact recovery (of {R['heldout']})")
    ax[0].set_title(f"K = 16, window {R['window']}: predicted degree {R['safe_degree']}, "
                    f"{R['calibration_states']:,} of 65,536 states", fontsize=10)
    ax[0].legend(fontsize=8)
    w = D["window_sweep"]
    ax[1].plot([v["t1"] for v in w], [v["predicted_degree"] for v in w], "o-", color="#e45756", label="predicted degree (geometry)")
    ax[1].plot([v["t1"] for v in w], [v["needed_order"] for v in w], "s-", color="#2f7d32", label="order needed to decode ≥ 99 %")
    ax[1].set_xlabel("window end t₁ (steps)")
    ax[1].set_ylabel("polynomial degree in the wall bits")
    ax[1].set_title("K = 8: the window sets the degree (all 256 states known)", fontsize=10)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/degree.png", dpi=110)
    plt.close(fig)


def fig_sparse():
    """Post hoc: every k16 interaction (orders 2-6) by its own tour vs the window."""
    from wml import coefficients_from_calibration
    import confirm_k16 as C
    cal = np.load("results/k16_calibration.npz")
    F = C.features(cal["Y"])
    table = {m: F[i] for i, m in enumerate(cal["masks"].tolist())}
    coef = coefficients_from_calibration(table, 16, 6)
    one = np.median([np.linalg.norm(coef[1 << k]) for k in range(16)])
    T = np.load("results/k16_tours.npy").min(0)
    pop = np.array([bin(m).count("1") for m in range(1 << 16)])
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    for k in range(2, 7):
        m = np.flatnonzero(pop == k)
        sz = np.linalg.norm(coef[m], axis=1) / one
        ax.scatter(T[m] - 40 - C.T1, np.maximum(sz, 1e-9), s=4, color=ORDER_COL[k - 1], alpha=.6, label=f"order {k}")
    ax.axvline(0, color="k", lw=.8)
    ax.axhline(np.finfo(np.float32).eps * 30, color="k", ls=":", lw=.8)
    ax.text(2, 3e-6, "float32 storage floor", fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("own tour − 40 − window end (steps)   ← predicted active | predicted inactive →")
    ax.set_ylabel("interaction size / one-bit response")
    ax.set_title("K = 16, window [145, 230): every interaction, placed by its own scattering tour (post hoc)", fontsize=10)
    ax.legend(fontsize=8, markerscale=3)
    fig.tight_layout()
    fig.savefig("figures/sparse_law.png", dpi=110)
    plt.close(fig)


def fig_stadium():
    import os
    if not os.path.exists("results/STADIUM.json"):
        return
    from wml.sparse import SparseMemory, predicted_log_size
    import confirm_stadium as S
    R = json.load(open("results/STADIUM.json"))
    P, act, keep = S.load_plan()
    cal = np.load(S.OUT + "calibration.npz")
    table = {int(m): f for m, f in zip(cal["masks"], cal["F"])}
    mem = SparseMemory(S.K, table)
    one = np.median([np.linalg.norm(c) for m, c in zip(mem.masks, mem.coef) if bin(int(m)).count("1") == 1])
    size = {int(m): np.linalg.norm(c) / one for m, c in zip(mem.masks, mem.coef)}
    keepm = {sum(1 << q for q in s) for s in keep}
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for kept, col, lab in ((False, "#9aa0a8", "pruned away"), (True, "#e0a030", "kept (ε = 1e-3)")):
        xs, ys = [], []
        for s_, t in act.items():
            if len(s_) < 2:
                continue
            m = sum(1 << q for q in s_)
            if (m in keepm) != kept:
                continue
            xs.append(predicted_log_size(t, len(s_), P["t1"]))
            ys.append(np.log10(max(size[m], 1e-12)))
        ax[0].scatter(xs, ys, s=4, color=col, alpha=.6, label=lab)
    lo, hi = -12, 2
    ax[0].plot([lo, hi], [lo, hi], "k-", lw=.8, label="measured = predicted")
    ax[0].axvline(-3, color="k", ls=":", lw=.8)
    ax[0].axhline(-2, color="#e45756", ls="--", lw=.8, label="N4b limit (1e-2)")
    ax[0].set_xlim(lo, hi); ax[0].set_ylim(-12.5, 1)
    ax[0].set_xlabel("predicted log10 size (model fitted on the K=16 ellipse)")
    ax[0].set_ylabel("measured log10 size / one-bit")
    ax[0].set_title(f"Stadium, every active set (orders 2-10): r = {R['amplitude']['r']:.2f} above 1e-5", fontsize=10)
    ax[0].legend(fontsize=8, markerscale=3)
    names = ["all walls\n2^24", "degree ceiling\n(order 10)", "per-subset\ntour law", "amplitude-\npruned", "tour law\ncut at pairs", "additive"]
    vals = [1 << 24, 4540386, R["models"]["tour_law"]["calibration_states"], R["models"]["pruned"]["calibration_states"],
            R["models"]["tour_law_order2"]["calibration_states"], R["models"]["additive"]["calibration_states"]]
    ex = [None, None, R["models"]["tour_law"]["exact"], R["models"]["pruned"]["exact"],
          R["models"]["tour_law_order2"]["exact"], R["models"]["additive"]["exact"]]
    cols = ["#9aa0a8", "#9aa0a8", "#2f7d32", "#2f7d32", "#e45756", "#e45756"]
    ax[1].bar(range(6), vals, color=cols)
    ax[1].set_yscale("log")
    for i, (v, e) in enumerate(zip(vals, ex)):
        ax[1].text(i, v * 1.4, f"{v:,}" + (f"\n{e}/300 exact" if e is not None else "\nnot simulated"), ha="center", fontsize=8)
    ax[1].set_xticks(range(6)); ax[1].set_xticklabels(names, fontsize=8)
    ax[1].set_ylim(5, 2e9)
    ax[1].set_ylabel("calibration states (log)")
    ax[1].set_title("What it costs to read a 24-bit wall in the stadium", fontsize=10)
    fig.tight_layout()
    fig.savefig("figures/stadium.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    fig_geometry()
    fig_onsets()
    fig_k16()
    fig_sparse()
    fig_stadium()

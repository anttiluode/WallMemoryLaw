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


if __name__ == "__main__":
    fig_geometry()
    fig_onsets()
    fig_k16()

# Pre-registration — K = 16

Pushed **before** any K = 16 wall state was simulated. Every number below the
"frozen" line was computed from geometry or from the K = 8 discovery run
(`discovery_k8.py`, `results/discovery_k8.json`), and is also hard-coded in
`confirm_k16.py`.

## The law under test

A written wall patch is a local scatterer. In the cavity's multiple-scattering series,
the interaction of a set S of patches (its Möbius / inclusion–exclusion coefficient)
collects only the paths that scatter off every patch in S. So it stays negligible until
the shortest such path, tour(S), can reach the receiver. A listening window that ends
at t₁ therefore sees the echo as a polynomial in the wall bits of degree

    d(t₁) = max { k : min_{|S| = k} tour(S) − margin ≤ t₁ }

and the whole memory is determined by the Σ_{k ≤ d} C(K, k) states with at most d
written patches, instead of 2^K.

## What K = 8 discovery showed (all 256 states, 3 sources, 1000 steps)

- Onset (first |I| > 1e-6 of peak) vs tour/c over all 765 interactions: **r = 0.998**.
  Onsets run up to 37 steps *early* (pulse tails and stencil precursors). This is not
  exact zero: 1e-12 of peak arrives ~54 steps early (`tests/test_gate0.py`).
- Margin frozen at **40 steps** (37 rounded up to 5).
- Window sweep [145, t₁], t₁ ∈ {245, 295, 345, 395, 445, 545}: calibrating to the
  predicted degree decoded every held-out state every time. The order actually
  needed was **1–2 below** the predicted degree (2 at the first five windows, 1 at
  545). So the law is a safe ceiling; the needed order is softer, because an
  interaction that has only just arrived is still small.

## Frozen for K = 16

- 16 patches, each covering the central **50 %** of its angular sector (gaps between
  memory cells), the same cavity, sources and receiver as Kavitaatio.
- Window **[145, 230)**, point receiver, 3 sources, 255 features. Chosen from geometry
  alone: by step 216 every patch's first echo has reached the receiver from some source.
- Predicted (safe) degree **6** (max over sources; per source 5, 6, 6).
  Calibration = all **14,893** states with ≤ 6 written patches (23 % of 65,536).
- Held-out: **300** random states with ≥ 7 written patches (seed 2026); none of them,
  nor any state with > 6 patches, is in the calibration set.
- Decoding: exhaustive nearest neighbour over all 65,536 model-predicted echoes.

## Decisions

| id | claim | pass rule |
|---|---|---|
| P1 | the safe degree is exact | order-6 model decodes ≥ 99 % of the 300 held-out states exactly (noise-free), and median relative model error ≤ 1e-3 |
| P2 | interactions switch on at their tours | 60 sampled subsets (10 per order 1–6, source 1, direct inclusion–exclusion): corr(onset, tour) ≥ 0.95 and no onset earlier than tour − 40 |
| P3 | the needed order is predictable | lowest calibration order that decodes ≥ 99 % lies in **[4, 6]** (safe − 2 … safe) |

Reported, not decided: exact decoding at Kavitaatio's noise level (4 % of the median
one-bit response RMS) for every order; additive (order 1) and pairwise (order 2)
models, both expected to fail.

## What would change the claim

- P1 fails → the geometric degree is not a safe ceiling at K = 16; the law is
  wrong or the margin is K-dependent. The README says so.
- P2 fails → onsets do not follow tours at K = 16 (e.g. gaps change the scattering).
- P3 fails low (needed < 4) → the predictor over-pays by more than two orders, so it
  is safe but too loose to be useful. Fails high → impossible if P1 passes.

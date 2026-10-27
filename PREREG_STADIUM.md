# Pre-registration 2 — new geometry, K = 24

Pushed **before** any stadium state was simulated. Everything frozen here was computed
from geometry (`python confirm_stadium.py plan` → `results/stadium_plan.json`) or
fitted on the K = 16 ellipse data that is already public in this repo.

## Why

The first confirmatory run (PREREG.md) passed P1 and P2 and failed P3: the order-level
degree is a safe ceiling but over-pays (14,893 states where 137 sufficed). Post hoc on
that data, two sharper rules appeared:

1. **Per-subset tour law.** Calibrate a patch set only if its *own* scattering tour
   arrives in the window: tour(S) − 40 ≤ t₁ for some source. A superset never has a
   shorter tour, so active sets are enumerated level by level.
2. **Amplitude model.** Inside the active sets,
   log₁₀(|I_S| / one-bit) ≈ −0.42 + 0.0992·(t₁ − tour(S)) − 0.75·|S|
   (K = 16 ellipse, r = 0.964, 0.33 dex scatter). Keep every single patch, plus every
   active set predicted ≥ ε = 1e-3. On the K = 16 data this needed 59 states.

Both were chosen after seeing data. This run tests them where they were not fitted.

## Frozen

- **Geometry:** a Bunimovich stadium (straight half-length 18, cap radius 22), with no
  foci and chaotic rays. Source and receiver are off every symmetry axis. Same grid,
  wave speed, damping and Robin strength as Kavitaatio. 24 patches, each covering the
  central 50 % of its angular sector.
- **Window rule:** t₀ = min single-patch tour − 40, t₁ = max single-patch tour + 20,
  giving **[155, 283)**. 3 sources, point receiver, 384 features.
- **Margin** 40 steps; **ε** = 1e-3; **amplitude coefficients** −0.42 / 0.0992 / 0.75;
  all as above.
- **Decoder:** descent over 1- and 2-bit flips, 18 starts (least-squares rounding,
  empty, 16 random). It matched exhaustive decoding 300/300 on K = 16. A miss counts as
  a search failure if the true wall fits the model better than the one found, otherwise
  as a model failure. Both are reported, and both count as misses.
- **Held-out:** 300 walls with uniformly random bits (seed 2027), about 12 written
  patches each.
- 2^24 = 16,777,216 walls; the order-level ceiling (degree 10) would need 4,540,386.

| model | calibration states | share of 2^24 |
|---|---|---|
| per-subset tour law | 7,164 | 0.043 % |
| amplitude-pruned (ε = 1e-3) | 848 | 0.005 % |

## Decisions

| id | claim | pass rule | prediction |
|---|---|---|---|
| N1 | the tour law is exact at K = 24 in a new shape | ≥ 99 % of 300 held-out walls decoded exactly (noise-free) | pass |
| N2 | the amplitude-pruned model (848 states) is exact | ≥ 99 % exact (noise-free) | uncertain |
| N3 | tour-inactive sets are empty | all 16 inactive pairs plus 24 random inactive triples whose pairs are all active (the hardest inactive cases), measured directly in float64: every one ≤ 1e-5 of the median one-bit response | pass |
| N4a | the amplitude model transfers across geometry | measured vs predicted log size over active sets above 1e-5: r ≥ 0.7 | uncertain |
| N4b | pruning is safe | no active set pruned away by the model measures > 1e-2 | uncertain |

Reported, not decided: decoding at Kavitaatio's noise level (4 % of the median one-bit
response RMS) for each model; the order-2 truncation of the tour law; the additive
model; bias and absolute deviation of the amplitude model (it is linear in lead time,
so it cannot know about saturation, which K = 16 never reached).

## What would change the claim

- N1 fails → the tour law is not a geometry-independent safe rule. Headline withdrawn.
- N2 fails but N1 passes → the amplitude model does not transfer. The README keeps the
  tour law and drops the 848-state claim.
- N3 fails → "inactive means empty" is false in a chaotic cavity.
- N4b fails → pruning by the frozen model can throw away real structure. N2 can then
  pass only by luck of the held-out draw.

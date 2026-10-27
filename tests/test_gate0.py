"""Gate 0: the identities the law rests on."""
import itertools
import os
import sys

import numpy as np
import pytest

from wml import Cavity, codes, mobius_from_table, zeta, coefficients_from_calibration, decode, tours

RNG = np.random.default_rng(0)


def test_mobius_zeta_inverse():
    K = 6
    Y = RNG.standard_normal((1 << K, 5))
    assert np.abs(zeta(mobius_from_table(Y, K), K) - Y).max() < 1e-12


def test_low_order_calibration_exact_on_low_degree_table():
    K, d = 7, 3
    c = np.zeros((1 << K, 4))
    for m in range(1 << K):
        if bin(m).count("1") <= d:
            c[m] = RNG.standard_normal(4)
    Y = zeta(c, K)
    cal = {m: Y[m] for m in range(1 << K) if bin(m).count("1") <= d}
    model = zeta(coefficients_from_calibration(cal, K, d), K)
    assert np.abs(model - Y).max() < 1e-10                     # degree-d table fully pinned by <=d states


def test_decode_exhaustive():
    T = RNG.standard_normal((64, 9))
    assert (decode(T, T[[3, 17, 63]] + 1e-9) == [3, 17, 63]).all()


@pytest.mark.skipif(not os.path.exists("/home/claude/anttiluode/kavitaatio/experiment.py") and
                    not os.environ.get("KAVITAATIO"), reason="Kavitaatio checkout not available")
def test_matches_kavitaatio():
    sys.path.insert(0, os.environ.get("KAVITAATIO", "/home/claude/anttiluode/kavitaatio"))
    import experiment as E
    cav = Cavity(8, 1.0)
    g = E.geometry()
    bits = RNG.integers(0, 2, (2, 8))
    mine = cav.run(bits, source=1, steps=350)
    for i, b in enumerate(bits):
        tr, _ = E.simulate(b, 1, g)
        assert np.abs(mine[i, 145:345:2, 0] - tr[:, 0]).max() < 1e-12


def test_batch_equals_single_and_linearity_in_wall_free_part():
    cav = Cavity(8, 1.0)
    bits = RNG.integers(0, 2, (3, 8))
    b = cav.run(bits, steps=120)
    for i in range(3):
        assert np.abs(cav.run(bits[i:i + 1], steps=120)[0] - b[i]).max() < 1e-12


def test_causality_no_interaction_before_its_tour():
    """A 2-patch interaction stays below 1e-6 of its peak until 40 steps before its tour.
    Not exactly zero: the Gaussian pulse's tails and the discrete stencil leak
    exponentially small precursors (1e-12 of peak arrives ~54 steps early)."""
    cav = Cavity(8, 1.0)
    T = tours(cav, source=0)
    C4 = codes(8)[[0, 1, 2, 3]]                                # {}, {0}, {1}, {0,1}
    Y = cav.run(C4, source=0, steps=320)[:, :, 0]
    I = Y[3] - Y[1] - Y[2] + Y[0]
    t = int(T[0b11] - 40)
    assert np.abs(I[:t]).max() < 1e-6 * np.abs(I).max()


def test_k16_patches_have_gaps_and_cells():
    cav = Cavity(16, 0.5)
    assert all(len(cav.patch_points(k)) >= 3 for k in range(16))
    full = Cavity(16, 1.0)
    assert cav.counts.sum() < full.counts.sum()

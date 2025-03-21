import numpy as np
import pytest

from curvekit import engine

pytestmark = pytest.mark.skipif(not engine.available(), reason="C++ library not built")


def test_bootstrap_reprices_and_matches_hand_calculation():
    c = engine.call("curve", mode="log_linear")
    assert c["max_error"] < 1e-13 and c["mode"] == "log_linear"
    pillars = np.array(c["rows"]["pillar"])
    assert len(pillars) == 12
    # first two pillars by hand: simple-rate deposits
    assert abs(pillars[0, 2] - 1 / (1 + 0.0310 * 0.25)) < 1e-12
    assert abs(pillars[1, 2] - 1 / (1 + 0.0325 * 0.5)) < 1e-12
    # the forward rate agreement fixes the next one
    assert abs(pillars[2, 2] - pillars[1, 2] / (1 + 0.0352 * 0.5)) < 1e-12
    assert np.all(np.diff(pillars[:, 2]) < 0)


def test_custom_quotes(tmp_path):
    q = tmp_path / "q.txt"
    q.write_text("deposit 0 1 0.05 1\nswap 0 2 0.05 1\n")
    c = engine.call("curve", q, mode="linear_zero")
    pillars = np.array(c["rows"]["pillar"])
    # annual 5% everywhere: discount factors 1 / 1.05^t
    assert np.allclose(pillars[:, 2], 1.05 ** -pillars[:, 0], atol=1e-12)


def test_bond_numbers_are_consistent():
    b = engine.call("bond", maturity=5, coupon=0.03, frequency=1, shifts=[-1, 1])
    assert 0 < b["modified"] < b["macaulay"] < 5
    lo, hi = b["rows"]["shift"]
    assert abs((lo[1] - hi[1]) / 2 - b["dv01"]) < 1e-6
    assert abs(sum(r[1] for r in b["rows"]["krd"]) - b["effective_duration"]) < 1e-4


def test_hull_white_tree_and_calibration():
    h = engine.call("hw", steps=[16, 64])
    coarse, fine = h["rows"]["tree"]
    assert abs(h["parity_gap"]) < 1e-13 and fine[2] < 1e-12
    assert abs(fine[5] / h["payer"] - 1) < 2e-3 and fine[7] > fine[5]
    c = engine.call("calibrate", a=0.1, sigma=0.008)
    assert abs(c["fitted_a"] - 0.1) < 1e-4 and abs(c["fitted_sigma"] - 0.008) < 1e-6


def test_errors_are_reported(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("bond 0 1 0.05 1\n")
    with pytest.raises(Exception):  # noqa: B017
        engine.call("curve", bad)

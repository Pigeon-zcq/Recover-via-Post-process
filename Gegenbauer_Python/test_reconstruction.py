"""Lightweight regression tests for the translated Gegenbauer routines."""

from __future__ import annotations

import unittest

from mpmath import mp

from gegenbauer_reconstruction import (
    composite_gauss_fourier_coefficients,
    endpoint_power,
    gegenbauer_value,
    paper_example_function,
    reconstruct_from_function,
    transform_x_to_y,
    transform_y_to_x,
)


class GegenbauerTranslationTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.dps = 60

    def test_transform_round_trip(self) -> None:
        for q in (2, 3):
            bet = mp.one / q
            for x in (mp.mpf("-1"), mp.mpf("-0.75"), mp.mpf("0.125"), mp.mpf("1")):
                y = transform_x_to_y(x, bet)
                recovered = transform_y_to_x(y, bet)
                self.assertLess(abs(recovered - x), mp.mpf("1e-50"))

    def test_general_interval_round_trip(self) -> None:
        bet = mp.mpf("0.5")
        a = mp.mpf("-0.5")
        b = mp.zero
        for x in (a, mp.mpf("-0.375"), mp.mpf("-0.125"), b):
            y = transform_x_to_y(x, bet, a, b)
            recovered = transform_y_to_x(y, bet, a, b)
            self.assertLess(abs(recovered - x), mp.mpf("1e-50"))

    def test_gegenbauer_known_polynomial(self) -> None:
        x = mp.mpf("0.3")
        lam = 2
        expected = 2 * lam * (lam + 1) * x**2 - lam
        self.assertLess(abs(gegenbauer_value(2, lam, x) - expected), mp.mpf("1e-50"))

    def test_endpoint_power_zero(self) -> None:
        self.assertEqual(endpoint_power(mp.zero, mp.mpf("0.5")), mp.zero)

    def test_fourier_coefficients_have_conjugate_symmetry(self) -> None:
        coeffs = composite_gauss_fourier_coefficients(paper_example_function(1, 2), 4, panels=40)
        for k in range(1, 5):
            self.assertLess(abs(coeffs[-k] - mp.conj(coeffs[k])), mp.mpf("1e-50"))

    def test_small_reconstruction_pipeline(self) -> None:
        result = reconstruct_from_function(
            paper_example_function(1, 2),
            n=8,
            m=1,
            lam=1,
            q=2,
            dps=50,
            panels=80,
            nzn=8,
        )
        self.assertEqual(len(result.reconstructed), 9)
        self.assertTrue(mp.isfinite(result.max_error))
        self.assertLess(result.max_error, mp.mpf("10"))


if __name__ == "__main__":
    unittest.main()

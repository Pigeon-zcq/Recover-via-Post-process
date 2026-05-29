"""Regression tests for the high-precision Fourier-GLOF implementation."""

from __future__ import annotations

import unittest

from mpmath import mp

from GLOF_Python.glof_reconstruction import (
    fourier_coefficients,
    fourier_glof_reconstruct,
    fourier_partial_sum,
    fourier_to_glof_coefficients_stable,
    glof_basis_matrix,
    glof_normalization,
    glof_quadrature,
    laguerre_gauss_nodes,
    laguerre_poly_matrix,
    model_log_singularity,
    stable_b_sigma,
)


class FourierGLOFTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.dps = 80

    def assertSmall(self, value: mp.mpf, tolerance: str = "1e-40") -> None:
        self.assertLess(abs(value), mp.mpf(tolerance))

    def test_laguerre_gauss_integrates_monomials(self) -> None:
        q = 4
        alpha = mp.mpf("0.5")
        y, w = laguerre_gauss_nodes(q, alpha)

        for degree in range(2 * q + 2):
            approx = mp.fsum(w_j * y_j**degree for y_j, w_j in zip(y, w))
            exact = mp.gamma(degree + alpha + 1)
            self.assertSmall(approx - exact)

    def test_glof_basis_is_orthogonal_under_glof_quadrature(self) -> None:
        alpha = mp.mpf("0.5")
        beta = mp.mpf("1.7")
        lambda_value = mp.mpf("0.2")
        m = 3
        q = 8

        quadrature = glof_quadrature(alpha, beta, lambda_value, q)
        phi = glof_basis_matrix(quadrature.x, quadrature.y, alpha, beta, lambda_value, m)
        gamma_l = glof_normalization(alpha, beta, m)

        for ell in range(m + 1):
            for order in range(m + 1):
                inner = mp.fsum(
                    phi[row, ell] * phi[row, order] * quadrature.chi[row]
                    for row in range(q + 1)
                )
                expected = gamma_l[ell] if ell == order else mp.zero
                self.assertSmall(inner - expected)

    def test_fourier_partial_sum_uses_negative_to_positive_ordering(self) -> None:
        values = fourier_partial_sum([0, 0, 1], [mp.mpf("0.25")], nf=1)
        self.assertSmall(values[0] - 1j)

    def test_fourier_coefficients_match_known_zero_mode(self) -> None:
        r = mp.mpf("0.5")
        k = 1
        func = model_log_singularity(r, k)
        coeffs = fourier_coefficients(func, 0, dps=80)
        exact = mp.gamma(k + 1) / (r + 1) ** (k + 1)
        self.assertSmall(coeffs[0] - exact, "1e-35")

    def test_stable_b_coefficients_match_direct_integral(self) -> None:
        alpha = 0
        beta = mp.mpf(1)
        lambda_value = 0
        m = 3
        q = 160
        mode = 1
        f_hat = {-1: 0, 0: 0, 1: 1}

        coefficients, gamma_l, _ = fourier_to_glof_coefficients_stable(
            f_hat, alpha, beta, lambda_value, m, q, nf=1
        )
        sigma = stable_b_sigma(beta, lambda_value)

        for ell in range(m + 1):
            direct = mp.quad(
                lambda z: mp.exp(2 * mp.pi * 1j * mode * mp.exp(-z / (sigma * (beta + 1))))
                * laguerre_poly_matrix([z / sigma], alpha, ell)[0, ell]
                * mp.exp(-z),
                [0, mp.inf],
            )
            direct /= gamma_l[ell] * (beta + 1) ** (alpha + 1) * sigma ** (alpha + 1)
            self.assertSmall(coefficients[ell] - direct, "1e-20")

    def test_reconstruct_constant_zero_mode(self) -> None:
        result = fourier_glof_reconstruct(
            [1],
            r=0,
            z_eval=[mp.mpf("0.2"), mp.mpf("0.8")],
            alpha=0,
            theta=0,
            gamma_fraction=0,
            q=4,
            m=0,
            lambda_value=0,
            beta=0,
            dps=80,
            return_details=True,
        )

        self.assertEqual(result.coefficient_method, "stable_b")
        self.assertEqual(result.params.q, 4)
        self.assertEqual(result.params.q_b, 4)
        self.assertSmall(result.values[0] - 1)
        self.assertSmall(result.values[1] - 1)

    def test_stable_b_uses_separate_qb_and_reports_adaptive_diagnostics(self) -> None:
        result = fourier_glof_reconstruct(
            [1],
            r=0,
            z_eval=[mp.mpf("0.25")],
            alpha=0,
            theta=0,
            gamma_fraction=0,
            q=2,
            q_b=3,
            m=0,
            lambda_value=0,
            beta=0,
            dps=80,
            adaptive_b=True,
            b_tolerance="1e-30",
            q_b_max=6,
            return_details=True,
        )

        self.assertEqual(result.params.q, 2)
        self.assertGreaterEqual(result.params.q_b, 3)
        self.assertIsNotNone(result.b_diagnostics)
        self.assertTrue(result.b_diagnostics.converged)
        self.assertSmall(result.values[0] - 1)


if __name__ == "__main__":
    unittest.main()

"""High-precision Fourier-GLOF reconstruction.

This module implements the algorithm described in ``code.md`` for recovering
endpoint log-singular functions from finitely many Fourier coefficients.
All numerical scalar arithmetic is routed through ``mpmath`` so callers can
raise ``mp.dps`` for difficult parameter regimes.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from mpmath import mp


Number = Union[int, float, str, mp.mpf]
ComplexNumber = Union[int, float, complex, str, mp.mpf, mp.mpc]
RealFunction = Callable[[mp.mpf], Union[mp.mpf, mp.mpc]]

DEFAULT_DPS = 80
CI = mp.mpc(0, 1)


@dataclass(frozen=True)
class GLOFParameters:
    """Derived parameters used by a Fourier-GLOF reconstruction."""

    nf: int
    m: int
    q: int
    q_b: int
    r: mp.mpf
    k: int
    alpha: mp.mpf
    theta: mp.mpf
    gamma_fraction: mp.mpf
    lambda_value: mp.mpf
    beta: mp.mpf
    nu: mp.mpf
    dps: int


@dataclass(frozen=True)
class StableBDiagnostics:
    """Convergence diagnostics for the stable ``B_{ell n}`` quadrature."""

    converged: bool
    q_b: int
    q_b_history: Tuple[int, ...]
    absolute_differences: Tuple[mp.mpf, ...]
    relative_differences: Tuple[mp.mpf, ...]
    error_estimate: Optional[mp.mpf] = None
    relative_error_estimate: Optional[mp.mpf] = None


@dataclass(frozen=True)
class GLOFQuadrature:
    """Laguerre-Gauss and mapped GLOF quadrature data."""

    y: List[mp.mpf]
    omega: List[mp.mpf]
    x: List[mp.mpf]
    chi: List[mp.mpf]


@dataclass(frozen=True)
class ReconstructionResult:
    """Detailed output from ``fourier_glof_reconstruct``."""

    values: List[mp.mpc]
    coefficients: List[mp.mpc]
    params: GLOFParameters
    quadrature: GLOFQuadrature
    gamma_l: List[mp.mpf]
    coefficient_method: str = "stable_b"
    b_diagnostics: Optional[StableBDiagnostics] = None


def set_precision(dps: int = DEFAULT_DPS) -> None:
    """Set the global mpmath decimal precision."""

    mp.dps = int(dps)


def as_mpf(value: Number) -> mp.mpf:
    """Convert numeric inputs to ``mp.mpf`` without preserving float noise."""

    if isinstance(value, mp.mpf):
        return value
    if isinstance(value, int):
        return mp.mpf(value)
    return mp.mpf(str(value))


def as_mpc(value: ComplexNumber) -> mp.mpc:
    """Convert numeric inputs to ``mp.mpc``."""

    if isinstance(value, mp.mpc):
        return value
    if isinstance(value, mp.mpf):
        return mp.mpc(value)
    if isinstance(value, complex):
        return mp.mpc(str(value.real), str(value.imag))
    if isinstance(value, int):
        return mp.mpc(value)
    return mp.mpc(value)


def infer_nf_from_coefficients(f_hat: Union[Sequence[ComplexNumber], Mapping[int, ComplexNumber]]) -> int:
    """Infer ``N_F`` from Fourier coefficients ordered by modes ``-N_F..N_F``."""

    if isinstance(f_hat, Mapping):
        if not f_hat:
            raise ValueError("f_hat mapping must not be empty")
        nf = max(abs(int(mode)) for mode in f_hat.keys())
        expected = set(range(-nf, nf + 1))
        missing = expected.difference(int(mode) for mode in f_hat.keys())
        if missing:
            raise ValueError(f"f_hat mapping is missing Fourier modes: {sorted(missing)[:8]}")
        return nf

    length = len(f_hat)
    if length == 0 or length % 2 == 0:
        raise ValueError("f_hat sequence length must be odd and ordered as -N_F,...,N_F")
    return (length - 1) // 2


def coefficient_items(
    f_hat: Union[Sequence[ComplexNumber], Mapping[int, ComplexNumber]],
    nf: Optional[int] = None,
) -> List[Tuple[int, mp.mpc]]:
    """Return Fourier coefficient pairs ``(n, f_hat_n)`` sorted by mode."""

    if nf is None:
        nf = infer_nf_from_coefficients(f_hat)
    nf = int(nf)

    if isinstance(f_hat, Mapping):
        converted = {int(mode): as_mpc(value) for mode, value in f_hat.items()}
        return [(n, converted[n]) for n in range(-nf, nf + 1)]

    expected_length = 2 * nf + 1
    if len(f_hat) != expected_length:
        raise ValueError(f"expected {expected_length} Fourier coefficients, got {len(f_hat)}")
    return [(n, as_mpc(f_hat[n + nf])) for n in range(-nf, nf + 1)]


def validate_parameters(alpha: mp.mpf, beta: mp.mpf, r: mp.mpf, lambda_value: mp.mpf, m: int, q: int) -> None:
    """Check the basic admissibility conditions from the algorithm document."""

    if alpha <= -1:
        raise ValueError(f"need alpha > -1, got {alpha}")
    if beta <= -1:
        raise ValueError(f"need beta > -1, got {beta}")
    if 2 * r + lambda_value <= -1:
        raise ValueError(f"need 2*r + lambda > -1, got {2 * r + lambda_value}")
    if m < 0:
        raise ValueError(f"need M >= 0, got {m}")
    if q < 0:
        raise ValueError(f"need Q >= 0, got {q}")
    if m > q:
        raise ValueError(f"need M <= Q for quadrature stability, got M={m}, Q={q}")


def validate_b_quadrature(m: int, q_b: int) -> None:
    """Check the basic admissibility of the stable ``B_{ell n}`` quadrature."""

    if q_b < 0:
        raise ValueError(f"need Q_B >= 0, got {q_b}")
    if m > q_b:
        raise ValueError(f"need M <= Q_B for stable B quadrature, got M={m}, Q_B={q_b}")


def laguerre_gauss_nodes(q: int, alpha: Number) -> Tuple[List[mp.mpf], List[mp.mpf]]:
    """Return ``Q+1`` generalized Laguerre-Gauss nodes and weights.

    The returned rule integrates ``p(y) * y**alpha * exp(-y)`` exactly for
    all polynomials ``p`` of degree at most ``2Q+1``.
    """

    q = int(q)
    alpha_m = as_mpf(alpha)
    if q < 0:
        raise ValueError(f"need Q >= 0, got {q}")
    if alpha_m <= -1:
        raise ValueError(f"need alpha > -1, got {alpha_m}")

    nodes, weights = mp.gauss_quadrature(q + 1, "glaguerre", alpha=alpha_m)
    y = [mp.mpf(nodes[i]) for i in range(q + 1)]
    omega = [mp.mpf(weights[i]) for i in range(q + 1)]
    return y, omega


def laguerre_poly_matrix(y_values: Iterable[Number], alpha: Number, m: int) -> mp.matrix:
    """Evaluate ``L_l^(alpha)(y_j)`` for ``l=0..M`` via three-term recurrence."""

    alpha_m = as_mpf(alpha)
    m = int(m)
    if m < 0:
        raise ValueError(f"need M >= 0, got {m}")

    y = [as_mpf(value) for value in y_values]
    n_points = len(y)
    values = mp.matrix(n_points, m + 1)

    for row in range(n_points):
        values[row, 0] = mp.one
        if m >= 1:
            values[row, 1] = alpha_m + 1 - y[row]

    for ell in range(1, m):
        ell_m = mp.mpf(ell)
        for row, y_row in enumerate(y):
            values[row, ell + 1] = (
                (2 * ell_m + alpha_m + 1 - y_row) * values[row, ell]
                - (ell_m + alpha_m) * values[row, ell - 1]
            ) / (ell_m + 1)

    return values


def glof_basis_matrix(
    x_values: Iterable[Number],
    y_values: Iterable[Number],
    alpha: Number,
    beta: Number,
    lambda_value: Number,
    m: int,
) -> mp.matrix:
    """Build ``Phi[j,l] = x_j^nu L_l^(alpha)(y_j)``."""

    x = [as_mpf(value) for value in x_values]
    y = [as_mpf(value) for value in y_values]
    if len(x) != len(y):
        raise ValueError("x_values and y_values must have the same length")

    beta_m = as_mpf(beta)
    lambda_m = as_mpf(lambda_value)
    nu = (beta_m - lambda_m) / 2
    laguerre_values = laguerre_poly_matrix(y, alpha, m)

    basis = mp.matrix(len(x), int(m) + 1)
    for row, x_row in enumerate(x):
        if x_row <= 0:
            raise ValueError("GLOF nodes must be positive")
        x_nu = mp.exp(nu * mp.log(x_row))
        for ell in range(int(m) + 1):
            basis[row, ell] = x_nu * laguerre_values[row, ell]
    return basis


def glof_normalization(alpha: Number, beta: Number, m: int) -> List[mp.mpf]:
    """Return ``gamma_l^(alpha,beta)`` for ``l=0..M`` using log-gamma."""

    alpha_m = as_mpf(alpha)
    beta_m = as_mpf(beta)
    m = int(m)
    if beta_m <= -1:
        raise ValueError(f"need beta > -1, got {beta_m}")

    gamma_l: List[mp.mpf] = []
    log_beta_scale = (alpha_m + 1) * mp.log(beta_m + 1)
    for ell in range(m + 1):
        ell_m = mp.mpf(ell)
        log_gamma_l = mp.loggamma(ell_m + alpha_m + 1) - mp.loggamma(ell_m + 1) - log_beta_scale
        gamma_l.append(mp.exp(log_gamma_l))
    return gamma_l


def glof_quadrature(alpha: Number, beta: Number, lambda_value: Number, q: int) -> GLOFQuadrature:
    """Generate mapped nodes ``x_j`` and GLOF weights ``chi_j``."""

    alpha_m = as_mpf(alpha)
    beta_m = as_mpf(beta)
    lambda_m = as_mpf(lambda_value)
    if beta_m <= -1:
        raise ValueError(f"need beta > -1, got {beta_m}")

    y, omega = laguerre_gauss_nodes(q, alpha_m)
    beta_plus_one = beta_m + 1
    scale = mp.power(beta_plus_one, -(alpha_m + 1))

    x: List[mp.mpf] = []
    chi: List[mp.mpf] = []
    for y_j, omega_j in zip(y, omega):
        x_j = mp.exp(-y_j / beta_plus_one)
        x.append(x_j)
        chi.append(mp.exp((lambda_m - beta_m) * mp.log(x_j)) * scale * omega_j)

    return GLOFQuadrature(y=y, omega=omega, x=x, chi=chi)


def fourier_partial_sum(
    f_hat: Union[Sequence[ComplexNumber], Mapping[int, ComplexNumber]],
    x_values: Iterable[Number],
    nf: Optional[int] = None,
) -> List[mp.mpc]:
    """Evaluate ``sum_{n=-N_F}^{N_F} f_hat_n exp(2*pi*i*n*x)`` at nonuniform nodes."""

    items = coefficient_items(f_hat, nf)
    values: List[mp.mpc] = []
    for x in x_values:
        x_m = as_mpf(x)
        values.append(mp.fsum(coeff * mp.exp(2 * mp.pi * CI * n * x_m) for n, coeff in items))
    return values


def stable_b_sigma(beta: Number, lambda_value: Number) -> mp.mpf:
    """Return the damping ``sigma`` in the stable ``B_{ell n}`` formula."""

    beta_m = as_mpf(beta)
    lambda_m = as_mpf(lambda_value)
    if beta_m <= -1:
        raise ValueError(f"need beta > -1, got {beta_m}")
    sigma = (beta_m + lambda_m + 2) / (2 * (beta_m + 1))
    if sigma <= 0:
        raise ValueError(f"stable B formula needs sigma > 0, got {sigma}")
    return sigma


def stable_b_quadrature(alpha: Number, beta: Number, lambda_value: Number, q: int) -> GLOFQuadrature:
    """Generate quadrature data for the scaled stable ``B_{ell n}`` formula.

    Here ``y`` stores the Laguerre variable ``z`` and ``x`` stores
    ``x_j^(B) = exp(-z_j / (sigma * (beta+1)))``.  The ``chi`` entries are the
    common scaled weights for the unnormalized matrix element ``D_{ell n}``.
    """

    alpha_m = as_mpf(alpha)
    beta_m = as_mpf(beta)
    sigma = stable_b_sigma(beta_m, lambda_value)
    z, omega = laguerre_gauss_nodes(q, alpha_m)

    scale = mp.power(beta_m + 1, -(alpha_m + 1)) * mp.power(sigma, -(alpha_m + 1))
    x_phase: List[mp.mpf] = []
    chi: List[mp.mpf] = []
    for z_j, omega_j in zip(z, omega):
        x_phase.append(mp.exp(-z_j / (sigma * (beta_m + 1))))
        chi.append(scale * omega_j)
    return GLOFQuadrature(y=z, omega=omega, x=x_phase, chi=chi)


def fourier_to_glof_coefficients_stable(
    f_hat: Union[Sequence[ComplexNumber], Mapping[int, ComplexNumber]],
    alpha: Number,
    beta: Number,
    lambda_value: Number,
    m: int,
    q: int,
    nf: Optional[int] = None,
) -> Tuple[List[mp.mpc], List[mp.mpf], GLOFQuadrature]:
    """Compute GLOF coefficients by directly applying stable ``B_{ell n}``.

    This implements the scaled formula from ``code.md``:

    ``B_{ell n} = [gamma_l (beta+1)^(alpha+1) sigma^(alpha+1)]^-1
       sum_j omega_j exp(2*pi*i*n*x_j^B) L_ell^(alpha)(z_j/sigma)``.

    It deliberately does not first evaluate the Fourier partial sum
    ``f_N(x_j)`` at GLOF quadrature nodes.
    """

    items = coefficient_items(f_hat, nf)
    alpha_m = as_mpf(alpha)
    beta_m = as_mpf(beta)
    sigma = stable_b_sigma(beta_m, lambda_value)
    m = int(m)
    q = int(q)

    quadrature = stable_b_quadrature(alpha_m, beta_m, lambda_value, q)
    scaled_nodes = [z_j / sigma for z_j in quadrature.y]
    laguerre_values = laguerre_poly_matrix(scaled_nodes, alpha_m, m)
    gamma_l = glof_normalization(alpha_m, beta_m, m)

    coefficients: List[mp.mpc] = []
    for ell in range(m + 1):
        normalizer = 1 / (gamma_l[ell] * mp.power(beta_m + 1, alpha_m + 1) * mp.power(sigma, alpha_m + 1))
        mode_terms: List[mp.mpc] = []
        for mode, f_hat_n in items:
            integral_sum = mp.fsum(
                quadrature.omega[row]
                * laguerre_values[row, ell]
                * mp.exp(2 * mp.pi * CI * mode * quadrature.x[row])
                for row in range(q + 1)
            )
            mode_terms.append(f_hat_n * normalizer * integral_sum)
        coefficients.append(mp.fsum(mode_terms))

    return coefficients, gamma_l, quadrature


def coefficient_max_norm(coefficients: Sequence[ComplexNumber]) -> mp.mpf:
    """Return the max norm of a coefficient vector."""

    if not coefficients:
        return mp.zero
    return max(abs(as_mpc(value)) for value in coefficients)


def coefficient_difference_norm(
    left: Sequence[ComplexNumber],
    right: Sequence[ComplexNumber],
) -> mp.mpf:
    """Return the max norm of the difference between two coefficient vectors."""

    if len(left) != len(right):
        raise ValueError("coefficient vectors must have the same length")
    if not left:
        return mp.zero
    return max(abs(as_mpc(left[idx]) - as_mpc(right[idx])) for idx in range(len(left)))


def fourier_to_glof_coefficients_stable_adaptive(
    f_hat: Union[Sequence[ComplexNumber], Mapping[int, ComplexNumber]],
    alpha: Number,
    beta: Number,
    lambda_value: Number,
    m: int,
    q_b: int,
    *,
    nf: Optional[int] = None,
    tolerance: Number = "1e-8",
    q_b_max: Optional[int] = None,
    max_steps: int = 4,
    growth_factor: Number = 2,
) -> Tuple[List[mp.mpc], List[mp.mpf], GLOFQuadrature, StableBDiagnostics]:
    """Compute stable ``B_{ell n}`` coefficients with a ``Q_B`` convergence check.

    The check compares the coefficient vector produced by the stable
    Fourier-to-GLOF matrix at successive ``Q_B`` values.  This catches the main
    practical failure mode: the bounded stable phase is safe to evaluate, but
    the oscillatory Laguerre-Gauss integral may still need many more nodes than
    the ordinary GLOF projection rule.
    """

    q_current = int(q_b)
    validate_b_quadrature(int(m), q_current)
    tol = as_mpf(tolerance)
    if tol < 0:
        raise ValueError(f"need tolerance >= 0, got {tol}")
    max_steps_i = int(max_steps)
    if max_steps_i < 0:
        raise ValueError(f"need max_steps >= 0, got {max_steps}")
    growth = as_mpf(growth_factor)
    if growth <= 1:
        raise ValueError(f"need growth_factor > 1, got {growth}")

    if q_b_max is None:
        q_limit = max(q_current, int(mp.ceil(q_current * mp.power(growth, max_steps_i))))
    else:
        q_limit = int(q_b_max)
        if q_limit < q_current:
            raise ValueError(f"need q_b_max >= q_b, got {q_limit} < {q_current}")

    coefficients, gamma_l, quadrature = fourier_to_glof_coefficients_stable(
        f_hat, alpha, beta, lambda_value, m, q_current, nf=nf
    )
    q_history = [q_current]
    absolute_differences: List[mp.mpf] = []
    relative_differences: List[mp.mpf] = []
    converged = False

    for _ in range(max_steps_i):
        if q_current >= q_limit:
            break
        q_next = int(mp.ceil(q_current * growth))
        q_next = max(q_current + 1, min(q_next, q_limit))
        validate_b_quadrature(int(m), q_next)

        next_coefficients, next_gamma_l, next_quadrature = fourier_to_glof_coefficients_stable(
            f_hat, alpha, beta, lambda_value, m, q_next, nf=nf
        )
        diff = coefficient_difference_norm(coefficients, next_coefficients)
        scale = max(mp.one, coefficient_max_norm(next_coefficients))
        rel = diff / scale

        q_history.append(q_next)
        absolute_differences.append(diff)
        relative_differences.append(rel)

        coefficients = next_coefficients
        gamma_l = next_gamma_l
        quadrature = next_quadrature
        q_current = q_next

        if rel <= tol:
            converged = True
            break

    diagnostics = StableBDiagnostics(
        converged=converged,
        q_b=q_current,
        q_b_history=tuple(q_history),
        absolute_differences=tuple(absolute_differences),
        relative_differences=tuple(relative_differences),
        error_estimate=absolute_differences[-1] if absolute_differences else None,
        relative_error_estimate=relative_differences[-1] if relative_differences else None,
    )
    return coefficients, gamma_l, quadrature, diagnostics


def fourier_to_glof_coefficients_via_nodes(
    f_hat: Union[Sequence[ComplexNumber], Mapping[int, ComplexNumber]],
    alpha: Number,
    beta: Number,
    lambda_value: Number,
    m: int,
    q: int,
    nf: Optional[int] = None,
) -> Tuple[List[mp.mpc], List[mp.mpf], GLOFQuadrature]:
    """Legacy path: evaluate ``f_N`` at GLOF nodes and then project."""

    if nf is None:
        nf = infer_nf_from_coefficients(f_hat)
    quadrature = glof_quadrature(alpha, beta, lambda_value, q)
    f_nodes = fourier_partial_sum(f_hat, quadrature.x, nf)
    phi = glof_basis_matrix(quadrature.x, quadrature.y, alpha, beta, lambda_value, m)
    gamma_l = glof_normalization(alpha, beta, m)
    coefficients = glof_coefficients(f_nodes, phi, quadrature.chi, gamma_l)
    return coefficients, gamma_l, quadrature


def glof_coefficients(
    f_values: Sequence[ComplexNumber],
    phi: mp.matrix,
    chi: Sequence[Number],
    gamma_l: Sequence[Number],
) -> List[mp.mpc]:
    """Compute ``c_tilde = G^{-1} Phi^T W f``."""

    f_m = [as_mpc(value) for value in f_values]
    chi_m = [as_mpf(value) for value in chi]
    gamma_m = [as_mpf(value) for value in gamma_l]

    rows, cols = phi.rows, phi.cols
    if len(f_m) != rows or len(chi_m) != rows:
        raise ValueError("f_values, chi, and phi row count must agree")
    if len(gamma_m) != cols:
        raise ValueError("gamma_l length must equal phi column count")

    coefficients: List[mp.mpc] = []
    for ell in range(cols):
        numerator = mp.fsum(f_m[row] * phi[row, ell] * chi_m[row] for row in range(rows))
        coefficients.append(numerator / gamma_m[ell])
    return coefficients


def glof_evaluate(
    coefficients: Sequence[ComplexNumber],
    z_eval: Iterable[Number],
    alpha: Number,
    beta: Number,
    lambda_value: Number,
) -> List[mp.mpc]:
    """Evaluate a GLOF expansion at points in ``[0, 1]``."""

    coeffs = [as_mpc(value) for value in coefficients]
    m = len(coeffs) - 1
    if m < 0:
        return []

    alpha_m = as_mpf(alpha)
    beta_m = as_mpf(beta)
    lambda_m = as_mpf(lambda_value)
    nu = (beta_m - lambda_m) / 2

    z = [as_mpf(value) for value in z_eval]
    output: List[mp.mpc] = []
    positive_z: List[mp.mpf] = []
    positive_index: List[int] = []

    for idx, z_i in enumerate(z):
        if z_i < 0 or z_i > 1:
            raise ValueError(f"evaluation point must be in [0, 1], got {z_i}")
        if z_i == 0:
            if nu > 0:
                output.append(mp.mpc(0))
            else:
                raise ValueError("cannot evaluate at x=0 when nu <= 0")
        else:
            output.append(mp.mpc(0))
            positive_z.append(z_i)
            positive_index.append(idx)

    if not positive_z:
        return output

    y_eval = [-(beta_m + 1) * mp.log(z_i) for z_i in positive_z]
    laguerre_values = laguerre_poly_matrix(y_eval, alpha_m, m)

    for local_idx, z_i in enumerate(positive_z):
        z_nu = mp.exp(nu * mp.log(z_i))
        value = mp.fsum(coeffs[ell] * z_nu * laguerre_values[local_idx, ell] for ell in range(m + 1))
        output[positive_index[local_idx]] = value

    return output


def derive_parameters(
    f_hat: Union[Sequence[ComplexNumber], Mapping[int, ComplexNumber]],
    r: Number,
    k: int = 0,
    alpha: Number = 0,
    theta: Number = "0.05",
    gamma_fraction: Optional[Number] = None,
    gamma: Optional[Number] = None,
    q: Optional[int] = None,
    Q: Optional[int] = None,
    q_b: Optional[int] = None,
    QB: Optional[int] = None,
    m: Optional[int] = None,
    M: Optional[int] = None,
    lambda_value: Optional[Number] = None,
    lambda_: Optional[Number] = None,
    beta: Optional[Number] = None,
    dps: int = DEFAULT_DPS,
) -> GLOFParameters:
    """Derive the parameter set recommended in ``code.md``."""

    nf = infer_nf_from_coefficients(f_hat)
    r_m = as_mpf(r)
    alpha_m = as_mpf(alpha)
    theta_m = as_mpf(theta)

    if gamma_fraction is not None and gamma is not None:
        raise ValueError("use either gamma_fraction or gamma, not both")
    gamma_input = gamma_fraction if gamma_fraction is not None else gamma
    gamma_fraction_m = as_mpf("0.1" if gamma_input is None else gamma_input)

    if lambda_value is not None and lambda_ is not None:
        raise ValueError("use either lambda_value or lambda_, not both")
    lambda_input = lambda_value if lambda_value is not None else lambda_

    if q is not None and Q is not None:
        raise ValueError("use either q or Q, not both")
    q_input = q if q is not None else Q

    if q_b is not None and QB is not None:
        raise ValueError("use either q_b or QB, not both")
    q_b_input = q_b if q_b is not None else QB

    if m is not None and M is not None:
        raise ValueError("use either m or M, not both")
    m_input = m if m is not None else M

    if lambda_input is None:
        lambda_m = mp.mpf(mp.floor(theta_m * nf))
    else:
        lambda_m = as_mpf(lambda_input)

    if beta is None:
        beta_m = lambda_m + 2 * r_m
    else:
        beta_m = as_mpf(beta)

    if m_input is None:
        m_value = int(mp.floor(gamma_fraction_m * nf))
    else:
        m_value = int(m_input)

    if q_input is None:
        q_value = max(2 * m_value, 2 * nf)
    else:
        q_value = int(q_input)

    if q_b_input is None:
        q_b_value = q_value
    else:
        q_b_value = int(q_b_input)

    nu = (beta_m - lambda_m) / 2
    validate_parameters(alpha_m, beta_m, r_m, lambda_m, m_value, q_value)
    validate_b_quadrature(m_value, q_b_value)
    return GLOFParameters(
        nf=nf,
        m=m_value,
        q=q_value,
        q_b=q_b_value,
        r=r_m,
        k=int(k),
        alpha=alpha_m,
        theta=theta_m,
        gamma_fraction=gamma_fraction_m,
        lambda_value=lambda_m,
        beta=beta_m,
        nu=nu,
        dps=int(dps),
    )


def fourier_glof_reconstruct(
    f_hat: Union[Sequence[ComplexNumber], Mapping[int, ComplexNumber]],
    r: Number,
    z_eval: Iterable[Number],
    *,
    k: int = 0,
    alpha: Number = 0,
    theta: Number = "0.05",
    gamma_fraction: Optional[Number] = None,
    gamma: Optional[Number] = None,
    q: Optional[int] = None,
    Q: Optional[int] = None,
    q_b: Optional[int] = None,
    QB: Optional[int] = None,
    m: Optional[int] = None,
    M: Optional[int] = None,
    lambda_value: Optional[Number] = None,
    lambda_: Optional[Number] = None,
    beta: Optional[Number] = None,
    dps: int = DEFAULT_DPS,
    coefficient_method: str = "stable_b",
    adaptive_b: bool = False,
    b_tolerance: Number = "1e-8",
    q_b_max: Optional[int] = None,
    b_max_steps: int = 4,
    b_growth_factor: Number = 2,
    return_details: bool = False,
) -> Union[List[mp.mpc], ReconstructionResult]:
    """Reconstruct values from Fourier coefficients using Fourier-GLOF.

    ``f_hat`` may be a sequence ordered as ``[-N_F, ..., N_F]`` or a mapping
    from integer Fourier modes to coefficients.  By default the recommended
    matching choice ``beta = lambda + 2*r`` is used.

    The default ``coefficient_method="stable_b"`` directly applies the stable
    Fourier-to-GLOF matrix ``B_{ell n}`` from ``code.md``.  Use
    ``coefficient_method="nodes"`` only for comparison with the older path
    that first evaluates the Fourier partial sum at GLOF nodes.

    ``Q``/``q`` controls the ordinary GLOF quadrature used by the node path.
    ``QB``/``q_b`` controls the Laguerre-Gauss rule used in the stable
    ``B_{ell n}`` integrals.  Set ``adaptive_b=True`` to compare successive
    ``Q_B`` values and return a ``StableBDiagnostics`` object.
    """

    set_precision(dps)
    params = derive_parameters(
        f_hat=f_hat,
        r=r,
        k=k,
        alpha=alpha,
        theta=theta,
        gamma_fraction=gamma_fraction,
        gamma=gamma,
        q=q,
        Q=Q,
        q_b=q_b,
        QB=QB,
        m=m,
        M=M,
        lambda_value=lambda_value,
        lambda_=lambda_,
        beta=beta,
        dps=dps,
    )

    method = coefficient_method.lower().replace("-", "_")
    b_diagnostics: Optional[StableBDiagnostics] = None
    if method in {"stable", "stable_b", "b", "direct_b"}:
        if adaptive_b:
            coefficients, gamma_l, quadrature, b_diagnostics = fourier_to_glof_coefficients_stable_adaptive(
                f_hat,
                params.alpha,
                params.beta,
                params.lambda_value,
                params.m,
                params.q_b,
                nf=params.nf,
                tolerance=b_tolerance,
                q_b_max=q_b_max,
                max_steps=b_max_steps,
                growth_factor=b_growth_factor,
            )
            params = replace(params, q_b=b_diagnostics.q_b)
        else:
            coefficients, gamma_l, quadrature = fourier_to_glof_coefficients_stable(
                f_hat,
                params.alpha,
                params.beta,
                params.lambda_value,
                params.m,
                params.q_b,
                params.nf,
            )
        method = "stable_b"
    elif method in {"nodes", "node", "legacy"}:
        coefficients, gamma_l, quadrature = fourier_to_glof_coefficients_via_nodes(
            f_hat,
            params.alpha,
            params.beta,
            params.lambda_value,
            params.m,
            params.q,
            params.nf,
        )
        method = "nodes"
    else:
        raise ValueError("coefficient_method must be 'stable_b' or 'nodes'")

    values = glof_evaluate(coefficients, z_eval, params.alpha, params.beta, params.lambda_value)

    if return_details:
        return ReconstructionResult(
            values=values,
            coefficients=coefficients,
            params=params,
            quadrature=quadrature,
            gamma_l=gamma_l,
            coefficient_method=method,
            b_diagnostics=b_diagnostics,
        )
    return values


def model_log_singularity(r: Number, k: int) -> RealFunction:
    """Return ``f(x) = x**r * (-log(x))**k`` with endpoint limits handled."""

    r_m = as_mpf(r)
    k_i = int(k)
    if r_m <= 0:
        raise ValueError("model helper assumes r > 0")
    if k_i < 0:
        raise ValueError("need k >= 0")

    def func(x: mp.mpf) -> mp.mpf:
        x_m = as_mpf(x)
        if x_m == 0:
            return mp.zero
        if x_m == 1:
            return mp.one if k_i == 0 else mp.zero
        return mp.power(x_m, r_m) * mp.power(-mp.log(x_m), k_i)

    return func


def fourier_coefficients(
    func: RealFunction,
    nf: int,
    *,
    dps: Optional[int] = None,
    intervals: Optional[Sequence[Number]] = None,
) -> List[mp.mpc]:
    """Compute Fourier coefficients ``hat f_n`` for modes ``-N_F..N_F``.

    The integral convention is
    ``hat f_n = integral_0^1 f(x) exp(-2*pi*i*n*x) dx``.
    """

    if dps is not None:
        set_precision(dps)
    nf = int(nf)
    if nf < 0:
        raise ValueError(f"need nf >= 0, got {nf}")

    if intervals is None:
        points = [mp.zero, mp.mpf("1e-12"), mp.mpf("1e-8"), mp.mpf("1e-5"), mp.mpf("1e-3"), mp.mpf("0.05"), mp.one]
    else:
        points = [as_mpf(point) for point in intervals]
    if points[0] != 0 or points[-1] != 1:
        raise ValueError("integration intervals must start at 0 and end at 1")

    coefficients: List[mp.mpc] = []
    for n in range(-nf, nf + 1):
        n_i = int(n)

        def integrand(x: mp.mpf) -> mp.mpc:
            return as_mpc(func(x)) * mp.exp(-2 * mp.pi * CI * n_i * x)

        coefficients.append(mp.quad(integrand, points))
    return coefficients


def linspace(start: Number, stop: Number, count: int) -> List[mp.mpf]:
    """Small mpmath-only replacement for ``numpy.linspace``."""

    count = int(count)
    if count <= 0:
        return []
    start_m = as_mpf(start)
    stop_m = as_mpf(stop)
    if count == 1:
        return [start_m]
    step = (stop_m - start_m) / (count - 1)
    return [start_m + step * i for i in range(count)]


def max_abs_error(reference: Sequence[ComplexNumber], approximation: Sequence[ComplexNumber]) -> mp.mpf:
    """Return max pointwise absolute error."""

    if len(reference) != len(approximation):
        raise ValueError("reference and approximation lengths must agree")
    if not reference:
        return mp.zero
    return max(abs(as_mpc(a) - as_mpc(b)) for a, b in zip(reference, approximation))


def rms_error(reference: Sequence[ComplexNumber], approximation: Sequence[ComplexNumber]) -> mp.mpf:
    """Return root mean square pointwise error."""

    if len(reference) != len(approximation):
        raise ValueError("reference and approximation lengths must agree")
    if not reference:
        return mp.zero
    mean_square = mp.fsum(abs(as_mpc(a) - as_mpc(b)) ** 2 for a, b in zip(reference, approximation)) / len(reference)
    return mp.sqrt(mean_square)

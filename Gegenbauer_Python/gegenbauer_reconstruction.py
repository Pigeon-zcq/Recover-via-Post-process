"""High-precision Gegenbauer reconstruction translated from the Fortran code.

The routines in this module keep the numerical structure of the original
Fortran programs:

* composite three-point Gauss integration for Fourier coefficients,
* Chebyshev-Gauss quadrature for the Gegenbauer coefficients,
* the singularity-removing map between x and y,
* optional Fourier-Galerkin time stepping used by the transport example.

All scalar real arithmetic is routed through ``mpmath``.  The functions accept
plain Python integers for counters and convert numerical data to ``mp.mpf`` or
``mp.mpc`` at the boundary.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from mpmath import mp


RealFunction = Callable[[mp.mpf], mp.mpf]
ComplexCoefficients = Dict[int, mp.mpc]


DEFAULT_DPS = 80
MPSML = mp.mpf("1e-16")
CI = mp.mpc(0, 1)


def as_mpf(value: object) -> mp.mpf:
    """Convert literals through strings when possible to avoid binary floats."""

    if isinstance(value, mp.mpf):
        return value
    if isinstance(value, int):
        return mp.mpf(value)
    return mp.mpf(str(value))


def set_precision(dps: int = DEFAULT_DPS) -> None:
    mp.dps = int(dps)


def _fortran_int(value: mp.mpf) -> int:
    """Fortran INT truncates toward zero."""

    return int(value)


def _periodic_two(value: mp.mpf) -> mp.mpf:
    """Mirror ``tem - dble(INT(tem/2)*2)`` from the Fortran sources."""

    value = as_mpf(value)
    return value - mp.mpf(2) * _fortran_int(value / 2)


def endpoint_power(base: mp.mpf, exponent: mp.mpf, eps: mp.mpf = MPSML) -> mp.mpf:
    base = as_mpf(base)
    exponent = as_mpf(exponent)
    if abs(base) < eps:
        return mp.zero
    return mp.power(base, exponent)


def fs_periodic_singularity(
    z: mp.mpf, p: mp.mpf = mp.mpf(1), q: mp.mpf = mp.mpf(2)
) -> mp.mpf:
    """Fortran ``fs``: periodicized ``(1+x) ** (p/q)`` on period two."""

    z = as_mpf(z)
    p = as_mpf(p)
    q = as_mpf(q)
    tt = _periodic_two(z + 1 + 2)
    return endpoint_power(tt, p / q)


def fv_transport_exact(
    z: mp.mpf, p: mp.mpf = mp.mpf(1), q: mp.mpf = mp.mpf(2)
) -> mp.mpf:
    """Fortran ``fv``: exact value in the variable-coefficient example."""

    z = as_mpf(z)
    p = as_mpf(p)
    q = as_mpf(q)
    tt = _periodic_two(2 * z + 1 + 2)
    return endpoint_power(tt, p / q)


def paper_example_function(p: int, q: int) -> RealFunction:
    """Return ``cos(x) + sin(x) * (1+x) ** (p/q)`` from the 2015 paper."""

    p_m = mp.mpf(p)
    q_m = mp.mpf(q)

    def func(x: mp.mpf) -> mp.mpf:
        x = as_mpf(x)
        return mp.cos(x) + mp.sin(x) * endpoint_power(1 + x, p_m / q_m)

    return func


def transform_y_to_x(
    y: mp.mpf,
    bet: mp.mpf,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
    eps: mp.mpf = MPSML,
) -> mp.mpf:
    """Fortran ``transf`` with flag 1: singular map from y-space to x-space."""

    y = as_mpf(y)
    bet = as_mpf(bet)
    a = as_mpf(a)
    b = as_mpf(b)
    temp = mp.mpf("0.5") * mp.power(b - a, bet)
    if abs(1 + y) < eps:
        return a
    if abs(1 - y) < eps:
        return b
    return mp.power(temp * (1 + y), 1 / bet) + a


def transform_x_to_y(
    x: mp.mpf,
    bet: mp.mpf,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
    eps: mp.mpf = MPSML,
) -> mp.mpf:
    """Fortran ``transf`` with flag 2: inverse map from x-space to y-space."""

    x = as_mpf(x)
    bet = as_mpf(bet)
    a = as_mpf(a)
    b = as_mpf(b)
    temp = mp.mpf("0.5") * mp.power(b - a, bet)
    if abs(x - a) < eps:
        return mp.mpf(-1)
    if abs(b - x) < eps:
        return mp.one
    return -1 + mp.power(x - a, bet) / temp


def transform_values_y_to_x(
    y_values: Iterable[mp.mpf],
    bet: mp.mpf,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
) -> List[mp.mpf]:
    return [transform_y_to_x(y, bet, a, b) for y in y_values]


def transform_values_x_to_y(
    x_values: Iterable[mp.mpf],
    bet: mp.mpf,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
) -> List[mp.mpf]:
    return [transform_x_to_y(x, bet, a, b) for x in x_values]


def rising_over_factorial(lam: int, max_order: int) -> List[mp.mpf]:
    """Return ``(lambda)_j / j!`` for ``0 <= j <= max_order``."""

    lam_m = mp.mpf(lam)
    gt = [mp.one]
    for j in range(1, max_order + 1):
        gt.append(gt[-1] * (lam_m + j - 1) / j)
    return gt


def gegenbauer_value(order: int, lam: int, y: mp.mpf) -> mp.mpf:
    """Evaluate ``C_order^lam(y)`` using the recurrence in ``gegenP``."""

    y = as_mpf(y)
    lam_m = mp.mpf(lam)
    if order == 0:
        return mp.one
    if order == 1:
        return 2 * lam_m * y

    ta = mp.one
    tb = 2 * lam_m * y
    tc = mp.zero
    for j in range(2, order + 1):
        j_m = mp.mpf(j)
        tc = 2 * (j_m - 1 + lam_m) * y * tb
        tc -= (j_m - 2 + 2 * lam_m) * ta
        tc /= j_m
        ta = tb
        tb = tc
    return tc


def gegenbauer_grid_values(
    order: int, nn: int, lam: int, gt: Optional[Sequence[mp.mpf]] = None
) -> Tuple[List[mp.mpf], List[mp.mpf]]:
    """Fortran ``gegen`` on the Chebyshev grid ``cos(pi*i/nn)``."""

    if gt is None or len(gt) <= order:
        gt = rising_over_factorial(lam, order)
    y_values: List[mp.mpf] = []
    c_values: List[mp.mpf] = []
    for i in range(nn + 1):
        theta = mp.pi * i / nn
        y_values.append(mp.cos(theta))
        cn = mp.zero
        for j in range(order + 1):
            cn += gt[j] * gt[order - j] * mp.cos(theta * (order - 2 * j))
        c_values.append(cn)
    return y_values, c_values


def gegenbauer_norm(order: int, lam: int) -> mp.mpf:
    """Normalization ``h_l^lambda`` used by the original ``gcoe``."""

    order_m = mp.mpf(order)
    lam_m = mp.mpf(lam)
    return (
        mp.pi
        * mp.power(mp.mpf("0.5"), 2 * lam_m - 1)
        * mp.gamma(order_m + 2 * lam_m)
        / ((order_m + lam_m) * mp.gamma(lam_m) ** 2 * mp.gamma(order_m + 1))
    )


def composite_gauss_fourier_coefficients(
    func: RealFunction, n: int, panels: Optional[int] = None
) -> ComplexCoefficients:
    """Translate ``getPointValuefromFC``: first ``2N+1`` Fourier coefficients."""

    x1 = mp.sqrt(mp.mpf("0.6"))
    wt1 = mp.mpf(5) / 9
    wt0 = mp.mpf(8) / 9
    panel_count = int(panels if panels is not None else min(100 * n, 4000))
    dx = mp.one / panel_count
    wt1 *= dx
    wt0 *= dx
    coeffs: ComplexCoefficients = {}
    for k in range(0, n + 1):
        acc = mp.mpc(0)
        for j in range(-panel_count + 1, panel_count, 2):
            xn = (mp.mpf(j) - x1) * dx
            xm = mp.mpf(j) * dx
            xp = (mp.mpf(j) + x1) * dx
            phase_n = mp.exp(-CI * k * mp.pi * xn)
            phase_m = mp.exp(-CI * k * mp.pi * xm)
            phase_p = mp.exp(-CI * k * mp.pi * xp)
            acc += wt1 * (func(xn) * phase_n + func(xp) * phase_p)
            acc += wt0 * func(xm) * phase_m
        acc *= mp.mpf("0.5")
        coeffs[k] = mp.mpc(acc)
        coeffs[-k] = mp.conj(acc)
    return coeffs


def quad_fourier_coefficients(func: RealFunction, n: int) -> ComplexCoefficients:
    """High-precision Fourier coefficients via adaptive ``mp.quad``."""

    coeffs: ComplexCoefficients = {}
    intervals = [mp.mpf(-1), mp.mpf("-0.5"), mp.zero, mp.mpf("0.5"), mp.one]
    for k in range(0, n + 1):
        integrand = lambda x, kk=k: func(x) * mp.exp(-CI * kk * mp.pi * x)
        coeff = mp.mpf("0.5") * mp.quad(integrand, intervals)
        coeffs[k] = mp.mpc(coeff)
        coeffs[-k] = mp.conj(coeff)
    return coeffs


def fourier_coefficients(
    func: RealFunction, n: int, method: str = "gauss3", panels: Optional[int] = None
) -> ComplexCoefficients:
    if method == "gauss3":
        return composite_gauss_fourier_coefficients(func, n, panels=panels)
    if method == "quad":
        return quad_fourier_coefficients(func, n)
    raise ValueError(f"unknown Fourier coefficient method: {method}")


def evaluate_fourier_series(
    coeffs: Mapping[int, mp.mpc], x: mp.mpf, n: int
) -> mp.mpf:
    x = as_mpf(x)
    total = mp.mpc(coeffs.get(0, mp.mpc(0)))
    for k in range(1, n + 1):
        phase = mp.exp(CI * k * mp.pi * x)
        total += coeffs[k] * phase + coeffs[-k] / phase
    return mp.re(total)


def gegenbauer_coefficients_from_fourier(
    coeffs: Mapping[int, mp.mpc],
    n: int,
    m: int,
    lam: int,
    q: int,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
    nt: Optional[int] = None,
) -> List[mp.mpf]:
    """Translate the current Fortran ``gcoe`` for Fourier coefficients."""

    if lam <= 0:
        raise ValueError("lambda must be positive for this Gegenbauer basis")
    nt = int(nt if nt is not None else 3 * n)
    bet = mp.one / q
    xi = [mp.zero] * (nt + 1)
    xr = [mp.zero] * (nt + 1)
    ff = [mp.zero] * (nt + 1)
    for i in range(1, nt):
        theta = mp.pi * i / nt
        xi[i] = mp.cos(theta)
        xr[i] = transform_y_to_x(xi[i], bet, a, b)
        ff[i] = evaluate_fourier_series(coeffs, xr[i], n)

    gt = rising_over_factorial(lam, m)
    hg: List[mp.mpf] = []
    for order in range(m + 1):
        _, cnl = gegenbauer_grid_values(order, nt, lam, gt)
        acc = mp.zero
        for i in range(1, nt):
            # Chebyshev-Gauss quadrature converts the Gegenbauer weight
            # (1-y^2)^(lambda-1/2) into the sampled factor (1-y^2)^lambda.
            acc += cnl[i] * ff[i] * mp.power(1 - xi[i], lam) * mp.power(1 + xi[i], lam)
        hg.append(acc * mp.pi / nt / gegenbauer_norm(order, lam))
    return hg


def reconstruct_value(
    x: mp.mpf,
    hg: Sequence[mp.mpf],
    lam: int,
    q: int,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
) -> mp.mpf:
    bet = mp.one / q
    y = transform_x_to_y(x, bet, a, b)
    total = mp.zero
    for order, coeff in enumerate(hg):
        total += coeff * gegenbauer_value(order, lam, y)
    return total


def reconstruct_on_y_grid(
    hg: Sequence[mp.mpf],
    lam: int,
    q: int,
    nzn: int,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
) -> Tuple[List[mp.mpf], List[mp.mpf], List[mp.mpf]]:
    """Reconstruct on the same Chebyshev-in-y grid used by the Fortran code."""

    bet = mp.one / q
    y_values = [mp.zero] * (nzn + 1)
    x_values = [mp.zero] * (nzn + 1)
    values = [mp.zero] * (nzn + 1)
    gt = rising_over_factorial(lam, len(hg) - 1)
    for order, coeff in enumerate(hg):
        y_grid, cnl = gegenbauer_grid_values(order, nzn, lam, gt)
        for i in range(nzn + 1):
            y_values[i] = y_grid[i]
            values[i] += coeff * cnl[i]
    for i, y in enumerate(y_values):
        x_values[i] = transform_y_to_x(y, bet, a, b)
    return x_values, y_values, values


@dataclass
class ReconstructionResult:
    n: int
    m: int
    lam: int
    q: int
    max_error: mp.mpf
    mean_abs_error: mp.mpf
    x_at_max: mp.mpf
    x_values: List[mp.mpf]
    y_values: List[mp.mpf]
    reconstructed: List[mp.mpf]
    exact: List[mp.mpf]
    errors: List[mp.mpf]
    fourier_coefficients: ComplexCoefficients
    gegenbauer_coefficients: List[mp.mpf]


def reconstruct_from_function(
    func: RealFunction,
    n: int,
    m: int,
    lam: int,
    q: int,
    *,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
    nzn: Optional[int] = None,
    dps: int = DEFAULT_DPS,
    fourier_method: str = "gauss3",
    panels: Optional[int] = None,
) -> ReconstructionResult:
    """Full Fourier-coefficient to Gegenbauer reconstruction pipeline."""

    set_precision(dps)
    nzn = int(n if nzn is None else nzn)
    coeffs = fourier_coefficients(func, n, method=fourier_method, panels=panels)
    hg = gegenbauer_coefficients_from_fourier(coeffs, n, m, lam, q, a, b)
    x_values, y_values, reconstructed = reconstruct_on_y_grid(hg, lam, q, nzn, a, b)
    exact = [func(x) for x in x_values]
    errors = [abs(u - ex) for u, ex in zip(reconstructed, exact)]
    max_index = max(range(len(errors)), key=lambda i: errors[i])
    mean_abs_error = mp.fsum(errors) / len(errors)
    return ReconstructionResult(
        n=n,
        m=m,
        lam=lam,
        q=q,
        max_error=errors[max_index],
        mean_abs_error=mean_abs_error,
        x_at_max=x_values[max_index],
        x_values=x_values,
        y_values=y_values,
        reconstructed=reconstructed,
        exact=exact,
        errors=errors,
        fourier_coefficients=coeffs,
        gegenbauer_coefficients=hg,
    )


@dataclass
class FourierSolutionResult:
    n: int
    q: int
    max_error: mp.mpf
    mean_abs_error: mp.mpf
    x_at_max: mp.mpf
    x_values: List[mp.mpf]
    y_values: List[mp.mpf]
    values: List[mp.mpf]
    exact: List[mp.mpf]
    errors: List[mp.mpf]
    fourier_coefficients: ComplexCoefficients


def paper_linear_parameters(example: str, n: int) -> Tuple[int, int, int]:
    """Return ``(p, q, lambda, m)`` for Examples 6.1 and 6.2."""

    key = example.lower().replace("example", "").replace(".", "").strip()
    if key in {"61", "6 1"}:
        return 1, 2, max(1, n // 16), max(1, (3 * n) // 80)
    if key in {"62", "6 2"}:
        return 1, 3, max(1, n // 8), max(1, n // 32)
    raise ValueError("example must be '6.1' or '6.2'")


def run_paper_example(
    example: str,
    n: int,
    *,
    dps: int = DEFAULT_DPS,
    fourier_method: str = "gauss3",
    panels: Optional[int] = None,
    nzn: Optional[int] = None,
) -> ReconstructionResult:
    p, q, lam, m = paper_linear_parameters(example, n)
    return reconstruct_from_function(
        paper_example_function(p, q),
        n=n,
        m=m,
        lam=lam,
        q=q,
        dps=dps,
        fourier_method=fourier_method,
        panels=panels,
        nzn=nzn,
    )


@dataclass
class ParameterSearchResult:
    n: int
    q: int
    best_m: int
    best_lam: int
    max_error: mp.mpf
    mean_abs_error: mp.mpf
    x_at_max: mp.mpf


def find_best_parameters(
    func: RealFunction,
    n: int,
    q: int,
    *,
    max_lam: Optional[int] = None,
    max_m: Optional[int] = None,
    dps: int = DEFAULT_DPS,
    fourier_method: str = "gauss3",
    panels: Optional[int] = None,
    nzn: Optional[int] = None,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
) -> ParameterSearchResult:
    """Translate the optimize program: scan lambda and partial sums in m."""

    set_precision(dps)
    max_lam = int(max_lam if max_lam is not None else n // 2)
    max_m = int(max_m if max_m is not None else n // 2)
    nzn = int(nzn if nzn is not None else n)
    coeffs = fourier_coefficients(func, n, method=fourier_method, panels=panels)
    best: Optional[ParameterSearchResult] = None

    for lam in range(1, max_lam + 1):
        hg = gegenbauer_coefficients_from_fourier(coeffs, n, max_m, lam, q, a, b)
        y_values = [mp.cos(mp.pi * i / nzn) for i in range(nzn + 1)]
        x_values = [transform_y_to_x(y, mp.one / q, a, b) for y in y_values]
        exact = [func(x) for x in x_values]
        reconstructed = [mp.zero] * (nzn + 1)
        gt = rising_over_factorial(lam, max_m)
        for m in range(0, max_m + 1):
            _, cnl = gegenbauer_grid_values(m, nzn, lam, gt)
            for i in range(nzn + 1):
                reconstructed[i] += hg[m] * cnl[i]
            errors = [abs(u - ex) for u, ex in zip(reconstructed, exact)]
            max_index = max(range(len(errors)), key=lambda i: errors[i])
            max_error = errors[max_index]
            if best is None or max_error < best.max_error:
                best = ParameterSearchResult(
                    n=n,
                    q=q,
                    best_m=m,
                    best_lam=lam,
                    max_error=max_error,
                    mean_abs_error=mp.fsum(errors) / len(errors),
                    x_at_max=x_values[max_index],
                )
    if best is None:
        raise ValueError("empty parameter search")
    return best


def collocation_points(
    nfc: int,
    npt: int,
    n: int,
    *,
    points_are_y: bool,
    q: int,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
) -> Tuple[int, Dict[int, mp.mpf], Dict[int, mp.mpf]]:
    """Collocation grid generation used by the Fortran main programs."""

    bet = mp.one / q
    x: Dict[int, mp.mpf] = {}
    y: Dict[int, mp.mpf] = {}
    if nfc == 1:
        nl = -n + 1 if npt == 2 else -n
        if npt == 2:
            raw = {i: mp.mpf(i) / n for i in range(-n + 1, n + 1)}
        else:
            raw = {i: mp.mpf(2 * i) / (2 * n + 1) for i in range(-n, n + 1)}
    else:
        nl = 0
        if npt == 1:
            raw = {i: mp.cos(mp.pi * (2 * i + 1) / (2 * n + 2)) for i in range(n + 1)}
        elif npt == 2:
            raw = {i: mp.cos(mp.pi * i / n) for i in range(n + 1)}
        else:
            raw = {i: mp.cos(mp.pi * (2 * i) / (2 * n + 1)) for i in range(n + 1)}

    for i, value in raw.items():
        if points_are_y:
            y[i] = value
            x[i] = transform_y_to_x(value, bet, a, b)
        else:
            x[i] = value
            y[i] = transform_x_to_y(value, bet, a, b)
    return nl, x, y


def collocation_spectral_coefficients(
    values: Mapping[int, mp.mpf],
    x_nodes: Mapping[int, mp.mpf],
    nfc: int,
    n: int,
    npt: int,
) -> Tuple[ComplexCoefficients, List[mp.mpf]]:
    """Translate ``ccoe`` for Fourier or Chebyshev collocation coefficients."""

    hf: ComplexCoefficients = {}
    hhf = [mp.zero] * (n + 1)
    if nfc == 1:
        if npt == 2:
            for k in range(-n, n + 1):
                acc = mp.mpc(0)
                for ell in range(-n + 1, n + 1):
                    acc += values[ell] * mp.exp(-CI * k * mp.pi * x_nodes[ell])
                hf[k] = acc / (2 * n)
            hf[-n] *= mp.mpf("0.5")
            hf[n] *= mp.mpf("0.5")
        else:
            for k in range(-n, n + 1):
                acc = mp.mpc(0)
                for ell in range(-n, n + 1):
                    acc += values[ell] * mp.exp(-CI * k * mp.pi * x_nodes[ell])
                hf[k] = acc / (2 * n + 1)
        return hf, hhf

    if npt == 1:
        for k in range(n + 1):
            acc = mp.zero
            for ell in range(n + 1):
                acc += values[ell] * mp.cos(mp.pi * (2 * ell + 1) * k / (2 * n + 2))
            hhf[k] = acc * 2 / (n + 1)
        hhf[0] *= mp.mpf("0.5")
        hhf[n] *= mp.mpf("0.5")
    elif npt == 2:
        for k in range(n + 1):
            acc = mp.zero
            for ell in range(1, n):
                acc += values[ell] * mp.cos(mp.pi * ell * k / n)
            acc += values[0] * mp.mpf("0.5")
            acc += values[n] * (mp.one if k % 2 == 0 else -mp.one) * mp.mpf("0.5")
            hhf[k] = acc * 2 / n
        hhf[0] *= mp.mpf("0.5")
        hhf[n] *= mp.mpf("0.5")
    else:
        for k in range(n + 1):
            acc = values[0] * mp.mpf("0.5")
            for ell in range(1, n + 1):
                acc += values[ell] * mp.cos(mp.pi * 2 * ell * k / (2 * n + 1))
            hhf[k] = acc * 4 / (2 * n + 1)
        hhf[0] *= mp.mpf("0.5")
        hhf[n] *= mp.mpf("0.5")
    return hf, hhf


def old_collocation_gegenbauer_coefficients(
    hf: Mapping[int, mp.mpc],
    hhf: Sequence[mp.mpf],
    nfc: int,
    n: int,
    m: int,
    lam: int,
    q: int,
    *,
    npro: int = 2,
    a: mp.mpf = mp.mpf(-1),
    b: mp.mpf = mp.mpf(1),
    nt: Optional[int] = None,
) -> List[mp.mpf]:
    """Translate ``doublescollofindbestparameters.f`` ``gcoe``."""

    nt = int(nt if nt is not None else 3 * n)
    bet = mp.one / q
    xi = [mp.zero] * (nt + 1)
    xr = [mp.zero] * (nt + 1)
    ff = [mp.zero] * (nt + 1)
    for i in range(1, nt):
        theta = mp.pi * i / nt
        xi[i] = mp.cos(theta)
        xr[i] = transform_y_to_x(xi[i], bet, a, b)

    if nfc == 1:
        for k in range(-n, n + 1):
            for i in range(1, nt):
                ff[i] += mp.re(hf[k] * mp.exp(CI * k * mp.pi * xr[i]))
    else:
        for k in range(n + 1):
            for i in range(1, nt):
                ff[i] += hhf[k] * mp.cos(k * mp.acos(xr[i]))

    gt = rising_over_factorial(lam, m)
    hg: List[mp.mpf] = []
    for order in range(m + 1):
        _, cnl = gegenbauer_grid_values(order, nt, lam, gt)
        acc = mp.zero
        for i in range(1, nt):
            if npro == 2:
                factor = mp.sqrt(1 - xi[i]) * mp.power(1 + xi[i], mp.mpf(q) / 2)
            else:
                factor = mp.sqrt(1 - xi[i]) * mp.sqrt(1 + xi[i])
            acc += cnl[i] * ff[i] * factor
        hg.append(acc * mp.pi / nt / gegenbauer_norm(order, lam))
    return hg


def fourier_galerkin_transport_coefficients(
    n: int,
    *,
    p: int = 1,
    q: int = 2,
    cfl: mp.mpf = mp.mpf("0.5"),
    panels: Optional[int] = None,
) -> ComplexCoefficients:
    """Translate ``FourierGalerkinTransport`` for ``u_t - x u_x = 0``."""

    initial = lambda x: fs_periodic_singularity(x, mp.mpf(p), mp.mpf(q))
    ak = composite_gauss_fourier_coefficients(initial, n, panels=panels)
    a1 = mp.mpf(16) * mp.log(10)
    r1 = int(mp.sqrt(2 * n))
    filt = {
        i: mp.exp(-a1 * mp.power(mp.mpf(i) / n, r1))
        for i in range(-n, n + 1)
    }

    indices = list(range(-n, n + 1))
    matrix: Dict[Tuple[int, int], mp.mpf] = {}
    for i in indices:
        for j in indices:
            if i == j:
                matrix[(i, j)] = mp.zero
            else:
                sign = mp.one if (j - i) % 2 == 0 else -mp.one
                matrix[(i, j)] = mp.mpf(j) / (j - i) * sign

    def matvec(vec: Mapping[int, mp.mpc]) -> ComplexCoefficients:
        out: ComplexCoefficients = {}
        for i in indices:
            acc = mp.mpc(0)
            for j in indices:
                acc += matrix[(i, j)] * vec[j]
            out[i] = acc
        return out

    def add_scaled(base: Mapping[int, mp.mpc], scale: mp.mpf, inc: Mapping[int, mp.mpc]) -> ComplexCoefficients:
        return {i: base[i] + scale * inc[i] for i in indices}

    tfinal = mp.log(2)
    dt_initial = cfl * tfinal / (2 * n + 1)
    t = mp.zero
    while t < tfinal:
        dt = min(dt_initial, tfinal - t)
        ak1 = matvec(ak)
        ak2 = matvec(add_scaled(ak, mp.mpf("0.5") * dt, ak1))
        ak3 = matvec(add_scaled(ak, mp.mpf("0.5") * dt, ak2))
        ak4 = matvec(add_scaled(ak, dt, ak3))
        for i in indices:
            ak[i] += dt * (ak1[i] + 2 * ak2[i] + 2 * ak3[i] + ak4[i]) / 6
            ak[i] *= filt[i]
        t += dt
    return ak


def evaluate_transport_fourier_solution(
    n: int,
    *,
    a: mp.mpf = mp.mpf("-0.5"),
    b: mp.mpf = mp.zero,
    q: int = 2,
    dps: int = DEFAULT_DPS,
    panels: Optional[int] = None,
    nzn: Optional[int] = None,
    coeffs: Optional[Mapping[int, mp.mpc]] = None,
) -> FourierSolutionResult:
    """Evaluate the unprocessed Fourier-Galerkin transport solution."""

    set_precision(dps)
    nzn = int(nzn if nzn is not None else n)
    coeff_dict = (
        dict(coeffs)
        if coeffs is not None
        else fourier_galerkin_transport_coefficients(n, q=q, panels=panels)
    )
    bet = mp.one / q
    y_values = [mp.cos(mp.pi * i / nzn) for i in range(nzn + 1)]
    x_values = [transform_y_to_x(y, bet, a, b) for y in y_values]
    values = [evaluate_fourier_series(coeff_dict, x, n) for x in x_values]
    exact_func = lambda x: fv_transport_exact(x, mp.one, mp.mpf(q))
    exact = [exact_func(x) for x in x_values]
    errors = [abs(u - ex) for u, ex in zip(values, exact)]
    max_index = max(range(len(errors)), key=lambda i: errors[i])
    return FourierSolutionResult(
        n=n,
        q=q,
        max_error=errors[max_index],
        mean_abs_error=mp.fsum(errors) / len(errors),
        x_at_max=x_values[max_index],
        x_values=x_values,
        y_values=y_values,
        values=values,
        exact=exact,
        errors=errors,
        fourier_coefficients=coeff_dict,
    )


def reconstruct_variable_transport(
    n: int,
    m: int,
    lam: int,
    *,
    a: mp.mpf = mp.mpf("-0.5"),
    b: mp.mpf = mp.zero,
    q: int = 2,
    dps: int = DEFAULT_DPS,
    panels: Optional[int] = None,
    nzn: Optional[int] = None,
    coeffs: Optional[Mapping[int, mp.mpc]] = None,
) -> ReconstructionResult:
    """Run the variable-coefficient transport post-processing path."""

    set_precision(dps)
    nzn = int(nzn if nzn is not None else n)
    coeff_dict = (
        dict(coeffs)
        if coeffs is not None
        else fourier_galerkin_transport_coefficients(n, q=q, panels=panels)
    )
    hg = gegenbauer_coefficients_from_fourier(coeff_dict, n, m, lam, q, a, b)
    x_values, y_values, reconstructed = reconstruct_on_y_grid(hg, lam, q, nzn, a, b)
    exact_func = lambda x: fv_transport_exact(x, mp.one, mp.mpf(q))
    exact = [exact_func(x) for x in x_values]
    errors = [abs(u - ex) for u, ex in zip(reconstructed, exact)]
    max_index = max(range(len(errors)), key=lambda i: errors[i])
    return ReconstructionResult(
        n=n,
        m=m,
        lam=lam,
        q=q,
        max_error=errors[max_index],
        mean_abs_error=mp.fsum(errors) / len(errors),
        x_at_max=x_values[max_index],
        x_values=x_values,
        y_values=y_values,
        reconstructed=reconstructed,
        exact=exact,
        errors=errors,
        fourier_coefficients=coeff_dict,
        gegenbauer_coefficients=hg,
    )


def find_best_transport_parameters(
    n: int,
    *,
    max_lam: Optional[int] = None,
    max_m: Optional[int] = None,
    a: mp.mpf = mp.mpf("-0.5"),
    b: mp.mpf = mp.zero,
    q: int = 2,
    dps: int = DEFAULT_DPS,
    panels: Optional[int] = None,
    nzn: Optional[int] = None,
    coeffs: Optional[Mapping[int, mp.mpc]] = None,
) -> ParameterSearchResult:
    """Scan the transport post-processing parameters after one FG solve."""

    set_precision(dps)
    max_lam = int(max_lam if max_lam is not None else n // 2)
    max_m = int(max_m if max_m is not None else n // 2)
    nzn = int(nzn if nzn is not None else n)
    coeff_dict = (
        dict(coeffs)
        if coeffs is not None
        else fourier_galerkin_transport_coefficients(n, q=q, panels=panels)
    )
    y_values = [mp.cos(mp.pi * i / nzn) for i in range(nzn + 1)]
    x_values = [transform_y_to_x(y, mp.one / q, a, b) for y in y_values]
    exact = [fv_transport_exact(x, mp.one, mp.mpf(q)) for x in x_values]
    best: Optional[ParameterSearchResult] = None

    for lam in range(1, max_lam + 1):
        hg = gegenbauer_coefficients_from_fourier(coeff_dict, n, max_m, lam, q, a, b)
        reconstructed = [mp.zero] * (nzn + 1)
        gt = rising_over_factorial(lam, max_m)
        for m in range(0, max_m + 1):
            _, cnl = gegenbauer_grid_values(m, nzn, lam, gt)
            for i in range(nzn + 1):
                reconstructed[i] += hg[m] * cnl[i]
            errors = [abs(u - ex) for u, ex in zip(reconstructed, exact)]
            max_index = max(range(len(errors)), key=lambda i: errors[i])
            max_error = errors[max_index]
            if best is None or max_error < best.max_error:
                best = ParameterSearchResult(
                    n=n,
                    q=q,
                    best_m=m,
                    best_lam=lam,
                    max_error=max_error,
                    mean_abs_error=mp.fsum(errors) / len(errors),
                    x_at_max=x_values[max_index],
                )
    if best is None:
        raise ValueError("empty parameter search")
    return best


def mp_format(value: mp.mpf, digits: int = 12) -> str:
    return mp.nstr(value, n=digits, min_fixed=-5, max_fixed=5)

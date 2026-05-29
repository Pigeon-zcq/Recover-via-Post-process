"""Command-line demo for high-precision Fourier-GLOF reconstruction."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from mpmath import mp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from GLOF_Python.glof_reconstruction import (  # noqa: E402
    fourier_coefficients,
    fourier_glof_reconstruct,
    linspace,
    max_abs_error,
    model_log_singularity,
    rms_error,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a Fourier-GLOF reconstruction demo.")
    parser.add_argument("--NF", type=int, default=12, help="Fourier half-bandwidth.")
    parser.add_argument("--M", type=int, default=None, help="GLOF highest order. Defaults to floor(gamma*NF).")
    parser.add_argument("--Q", type=int, default=None, help="Quadrature highest index. Defaults to max(2M, 2NF).")
    parser.add_argument("--QB", type=int, default=None, help="Stable B quadrature highest index. Defaults to Q.")
    parser.add_argument("--r", default="0.5", help="Power singularity exponent.")
    parser.add_argument("--k", type=int, default=1, help="Log power.")
    parser.add_argument("--alpha", default="0.5", help="Laguerre weight parameter.")
    parser.add_argument("--theta", default="0.05", help="lambda = floor(theta*NF).")
    parser.add_argument("--gamma", dest="gamma_fraction", default="0.1", help="M = floor(gamma*NF).")
    parser.add_argument("--dps", type=int, default=80, help="mpmath decimal precision.")
    parser.add_argument("--points", type=int, default=80, help="Number of evaluation points.")
    parser.add_argument("--z-min", default="0.005", help="Left evaluation endpoint.")
    parser.add_argument("--z-max", default="0.995", help="Right evaluation endpoint.")
    parser.add_argument("--adaptive-B", action="store_true", help="Adaptively increase QB until stable B coefficients converge.")
    parser.add_argument("--B-tol", default="1e-8", help="Relative tolerance for --adaptive-B.")
    parser.add_argument("--QB-max", type=int, default=None, help="Maximum QB for --adaptive-B.")
    parser.add_argument(
        "--coefficient-method",
        choices=("stable_b", "nodes"),
        default="stable_b",
        help="Use stable direct B_{ell n} coefficients or the legacy node-value path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mp.dps = args.dps

    f_exact = model_log_singularity(args.r, args.k)
    z_eval = linspace(args.z_min, args.z_max, args.points)

    print("=== Fourier-GLOF Reconstruction ===")
    print(f"f(x) = x^{args.r} * (-log x)^{args.k}")
    print(f"Computing {2 * args.NF + 1} Fourier coefficients with mp.dps={args.dps} ...")

    f_hat = fourier_coefficients(f_exact, args.NF, dps=args.dps)
    result = fourier_glof_reconstruct(
        f_hat,
        args.r,
        z_eval,
        k=args.k,
        alpha=args.alpha,
        theta=args.theta,
        gamma_fraction=args.gamma_fraction,
        q=args.Q,
        QB=args.QB,
        m=args.M,
        dps=args.dps,
        coefficient_method=args.coefficient_method,
        adaptive_b=args.adaptive_B,
        b_tolerance=args.B_tol,
        q_b_max=args.QB_max,
        return_details=True,
    )

    f_true = [f_exact(z) for z in z_eval]
    e_inf = max_abs_error(f_true, result.values)
    e_rms = rms_error(f_true, result.values)
    params = result.params

    print(f"NF={params.nf}, M={params.m}, Q={params.q}, QB={params.q_b}")
    print(f"coefficient method={result.coefficient_method}")
    if result.b_diagnostics is not None:
        diag = result.b_diagnostics
        rel = "None" if diag.relative_error_estimate is None else mp.nstr(diag.relative_error_estimate, 12)
        print(f"B convergence={diag.converged}, QB history={diag.q_b_history}, rel estimate={rel}")
    print(f"alpha={mp.nstr(params.alpha, 12)}, lambda={mp.nstr(params.lambda_value, 12)}, beta={mp.nstr(params.beta, 12)}")
    print(f"max error: {mp.nstr(e_inf, 12)}")
    print(f"rms error: {mp.nstr(e_rms, 12)}")
    print("first coefficients:")
    for ell, coeff in enumerate(result.coefficients[: min(6, len(result.coefficients))]):
        print(f"  c[{ell}] = {mp.nstr(coeff, 12)}")


if __name__ == "__main__":
    main()

"""Command-line runners for the translated Gegenbauer reconstruction code."""

from __future__ import annotations

import argparse
from typing import Iterable, List

from mpmath import mp

try:
    from .gegenbauer_reconstruction import (
        DEFAULT_DPS,
        find_best_parameters,
        fs_periodic_singularity,
        mp_format,
        paper_example_function,
        paper_linear_parameters,
        reconstruct_from_function,
        reconstruct_variable_transport,
        run_paper_example,
        set_precision,
    )
except ImportError:  # pragma: no cover - enables direct script execution.
    from gegenbauer_reconstruction import (
        DEFAULT_DPS,
        find_best_parameters,
        fs_periodic_singularity,
        mp_format,
        paper_example_function,
        paper_linear_parameters,
        reconstruct_from_function,
        reconstruct_variable_transport,
        run_paper_example,
        set_precision,
    )


def _print_result_rows(results: Iterable[object]) -> None:
    print("N        lambda  m        Linf error        L1 mean          x@max")
    print("-" * 72)
    for result in results:
        print(
            f"{result.n:<8d} {result.lam:<7d} {result.m:<8d} "
            f"{mp_format(result.max_error, 10):<17s} "
            f"{mp_format(result.mean_abs_error, 10):<16s} "
            f"{mp_format(result.x_at_max, 10)}"
        )


def run_paper(args: argparse.Namespace) -> None:
    set_precision(args.dps)
    results = [
        run_paper_example(
            args.example,
            n,
            dps=args.dps,
            fourier_method=args.method,
            panels=args.panels,
            nzn=args.nzn,
        )
        for n in args.N
    ]
    _print_result_rows(results)


def run_scan(args: argparse.Namespace) -> None:
    set_precision(args.dps)
    p, q, _, _ = paper_linear_parameters(args.example, args.N)
    result = find_best_parameters(
        paper_example_function(p, q),
        args.N,
        q,
        max_lam=args.max_lambda,
        max_m=args.max_m,
        dps=args.dps,
        fourier_method=args.method,
        panels=args.panels,
        nzn=args.nzn,
    )
    print("N        best lambda  best m    Linf error        L1 mean          x@max")
    print("-" * 78)
    print(
        f"{result.n:<8d} {result.best_lam:<12d} {result.best_m:<9d} "
        f"{mp_format(result.max_error, 10):<17s} "
        f"{mp_format(result.mean_abs_error, 10):<16s} "
        f"{mp_format(result.x_at_max, 10)}"
    )


def run_periodic(args: argparse.Namespace) -> None:
    set_precision(args.dps)
    func = lambda x: fs_periodic_singularity(x, mp.one, mp.mpf(args.q))
    result = reconstruct_from_function(
        func,
        n=args.N,
        m=args.m,
        lam=args.lam,
        q=args.q,
        dps=args.dps,
        fourier_method=args.method,
        panels=args.panels,
        nzn=args.nzn,
    )
    _print_result_rows([result])


def run_transport(args: argparse.Namespace) -> None:
    set_precision(args.dps)
    result = reconstruct_variable_transport(
        args.N,
        args.m,
        args.lam,
        a=mp.mpf(str(args.a)),
        b=mp.mpf(str(args.b)),
        q=args.q,
        dps=args.dps,
        panels=args.panels,
        nzn=args.nzn,
    )
    _print_result_rows([result])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="High-precision Python translation of the Gegenbauer Fortran code."
    )
    parser.set_defaults(func=None)
    subparsers = parser.add_subparsers(dest="command")

    paper = subparsers.add_parser("paper", help="Run Chen-Shu 2015 Example 6.1 or 6.2.")
    paper.add_argument("--example", choices=["6.1", "6.2"], default="6.1")
    paper.add_argument("--N", type=int, nargs="+", default=[40])
    paper.add_argument("--dps", type=int, default=DEFAULT_DPS)
    paper.add_argument("--method", choices=["gauss3", "quad"], default="gauss3")
    paper.add_argument("--panels", type=int, default=None)
    paper.add_argument("--nzn", type=int, default=None)
    paper.set_defaults(func=run_paper)

    scan = subparsers.add_parser("scan", help="Search the best lambda and m.")
    scan.add_argument("--example", choices=["6.1", "6.2"], default="6.1")
    scan.add_argument("--N", type=int, default=40)
    scan.add_argument("--max-lambda", type=int, default=None)
    scan.add_argument("--max-m", type=int, default=None)
    scan.add_argument("--dps", type=int, default=DEFAULT_DPS)
    scan.add_argument("--method", choices=["gauss3", "quad"], default="gauss3")
    scan.add_argument("--panels", type=int, default=None)
    scan.add_argument("--nzn", type=int, default=None)
    scan.set_defaults(func=run_scan)

    periodic = subparsers.add_parser(
        "periodic", help="Run the current Fortran fs periodic singularity path."
    )
    periodic.add_argument("--N", type=int, required=True)
    periodic.add_argument("--m", type=int, required=True)
    periodic.add_argument("--lambda", dest="lam", type=int, required=True)
    periodic.add_argument("--q", type=int, default=2)
    periodic.add_argument("--dps", type=int, default=DEFAULT_DPS)
    periodic.add_argument("--method", choices=["gauss3", "quad"], default="gauss3")
    periodic.add_argument("--panels", type=int, default=None)
    periodic.add_argument("--nzn", type=int, default=None)
    periodic.set_defaults(func=run_periodic)

    transport = subparsers.add_parser(
        "transport", help="Run the Fourier-Galerkin variable-coefficient example."
    )
    transport.add_argument("--N", type=int, required=True)
    transport.add_argument("--m", type=int, required=True)
    transport.add_argument("--lambda", dest="lam", type=int, required=True)
    transport.add_argument("--a", default="-0.5")
    transport.add_argument("--b", default="0")
    transport.add_argument("--q", type=int, default=2)
    transport.add_argument("--dps", type=int, default=DEFAULT_DPS)
    transport.add_argument("--panels", type=int, default=None)
    transport.add_argument("--nzn", type=int, default=None)
    transport.set_defaults(func=run_transport)

    return parser


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.func is None:
        parser.print_help()
        return 2
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Compare Examples 6.1 and 6.2 against the published linear-choice tables."""

from __future__ import annotations

import argparse
from typing import Dict, List

from mpmath import mp

try:
    from .gegenbauer_reconstruction import DEFAULT_DPS, mp_format, run_paper_example
except ImportError:  # pragma: no cover
    from gegenbauer_reconstruction import DEFAULT_DPS, mp_format, run_paper_example


REFERENCE_LINEAR: Dict[str, Dict[int, mp.mpf]] = {
    "6.1": {
        40: mp.mpf("5.91e-1"),
        80: mp.mpf("2.95e-1"),
        160: mp.mpf("1.56e-2"),
        320: mp.mpf("1.33e-5"),
        640: mp.mpf("6.46e-12"),
    },
    "6.2": {
        40: mp.mpf("7.13e-1"),
        80: mp.mpf("3.64e-1"),
        160: mp.mpf("2.91e-1"),
        320: mp.mpf("1.17e-2"),
        640: mp.mpf("1.83e-6"),
    },
}


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example", choices=["6.1", "6.2"], nargs="+", default=["6.1", "6.2"])
    parser.add_argument("--N", type=int, nargs="+", default=[40])
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--method", choices=["gauss3", "quad"], default="gauss3")
    parser.add_argument("--panels", type=int, default=None)
    parser.add_argument("--nzn", type=int, default=None)
    args = parser.parse_args(argv)

    print("Example  N        lambda  m        computed       reference      ratio")
    print("-" * 78)
    for example in args.example:
        for n in args.N:
            result = run_paper_example(
                example,
                n,
                dps=args.dps,
                fourier_method=args.method,
                panels=args.panels,
                nzn=args.nzn,
            )
            ref = REFERENCE_LINEAR[example].get(n)
            if ref is None:
                ref_text = "n/a"
                ratio_text = "n/a"
            else:
                ref_text = mp_format(ref, 8)
                ratio_text = mp_format(result.max_error / ref, 8)
            print(
                f"{example:<8s} {n:<8d} {result.lam:<7d} {result.m:<8d} "
                f"{mp_format(result.max_error, 8):<14s} {ref_text:<14s} {ratio_text}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

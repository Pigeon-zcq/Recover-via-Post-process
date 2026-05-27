"""Reproduce Chen-Shu 2015 Examples 6.1 and 6.2 tables and Fig. 1.

The script recomputes the Gegenbauer reconstructions from the first 2N+1
Fourier coefficients and writes:

* Table 1 and Table 2 comparison files in CSV and Markdown.
* Pointwise error CSV files used for plotting.
* Fig. 1a, Fig. 1b, and a combined Fig. 1 PNG.

Full N=40,80,160,320,640 runs are intentionally expensive with mpmath.  Use
``--quick`` or pass a short ``--N`` list while developing.

Note: in the OCR Markdown table, the final two optimal-choice columns are
labelled as lambda and m, but the numbers match the Fortran program output
order ``best pair(m,lambda)``.  The constants below store them as
``optimal_lam`` and ``optimal_m`` after swapping that order.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from mpmath import mp
from PIL import Image, ImageDraw, ImageFont

try:
    from .gegenbauer_reconstruction import (
        DEFAULT_DPS,
        find_best_parameters,
        fourier_coefficients,
        gegenbauer_coefficients_from_fourier,
        mp_format,
        paper_example_function,
        reconstruct_on_y_grid,
        set_precision,
    )
except ImportError:  # pragma: no cover - allows direct script execution.
    from gegenbauer_reconstruction import (
        DEFAULT_DPS,
        find_best_parameters,
        fourier_coefficients,
        gegenbauer_coefficients_from_fourier,
        mp_format,
        paper_example_function,
        reconstruct_on_y_grid,
        set_precision,
    )


ROOT = Path(__file__).resolve().parent
DEFAULT_PAPER_MD = (
    ROOT
    / "Chen和Shu_2015_Recovering_Exponential_Accuracy_in_Fourier_Spectral_Methods_Involving_Piecewise_Smooth_Functions_wit.md"
)
DEFAULT_OUTPUT_DIR = ROOT / "reproductions" / "chen_shu_2015"
FULL_N_VALUES = [40, 80, 160, 320, 640]
QUICK_N_VALUES = [40, 80]


@dataclass(frozen=True)
class PaperRow:
    n: int
    linear_error: str
    order: Optional[str]
    linear_lam: int
    linear_m: int
    optimal_error: str
    optimal_lam: int
    optimal_m: int


@dataclass(frozen=True)
class PaperExample:
    key: str
    table_number: int
    figure_label: str
    p: int
    q: int
    title: str
    rows: Sequence[PaperRow]


EXAMPLES: Mapping[str, PaperExample] = {
    "6.1": PaperExample(
        key="6.1",
        table_number=1,
        figure_label="Fig. 1a",
        p=1,
        q=2,
        title="Example 6.1: cos(x) + sin(x) sqrt(1+x)",
        rows=[
            PaperRow(40, "5.91e-1", None, 2, 1, "1.21e-3", 18, 8),
            PaperRow(80, "2.95e-1", "1.00", 5, 3, "5.25e-5", 40, 12),
            PaperRow(160, "1.56e-2", "4.25", 10, 6, "5.57e-10", 80, 20),
            PaperRow(320, "1.33e-5", "10.19", 20, 12, "2.75e-11", 51, 22),
            PaperRow(640, "6.46e-12", "20.98", 40, 24, "1.11e-12", 44, 23),
        ],
    ),
    "6.2": PaperExample(
        key="6.2",
        table_number=2,
        figure_label="Fig. 1b",
        p=1,
        q=3,
        title="Example 6.2: cos(x) + sin(x) (1+x)^(1/3)",
        rows=[
            PaperRow(40, "7.13e-1", None, 5, 1, "1.83e-1", 12, 5),
            PaperRow(80, "3.64e-1", "0.97", 10, 2, "1.98e-2", 33, 7),
            PaperRow(160, "2.91e-1", "0.33", 20, 5, "1.24e-3", 76, 12),
            PaperRow(320, "1.17e-2", "4.64", 40, 10, "5.50e-6", 97, 17),
            PaperRow(640, "1.83e-6", "12.64", 80, 20, "3.15e-7", 89, 21),
        ],
    ),
}


@dataclass
class RunResult:
    n: int
    lam: int
    m: int
    max_error: mp.mpf
    mean_abs_error: mp.mpf
    x_at_max: mp.mpf
    x_values: List[mp.mpf]
    errors: List[mp.mpf]


def sci(value: mp.mpf | str | None, digits: int = 8) -> str:
    if value is None:
        return ""
    return mp.nstr(mp.mpf(value), n=digits, min_fixed=-4, max_fixed=4)


def ratio_text(computed: mp.mpf, reference: str) -> str:
    return sci(computed / mp.mpf(reference), 8)


def compute_reconstruction(
    *,
    example: PaperExample,
    n: int,
    lam: int,
    m: int,
    coeffs: Mapping[int, mp.mpc],
    nzn: Optional[int],
) -> RunResult:
    func = paper_example_function(example.p, example.q)
    grid_size = n if nzn is None else nzn
    hg = gegenbauer_coefficients_from_fourier(coeffs, n, m, lam, example.q)
    x_values, _, reconstructed = reconstruct_on_y_grid(hg, lam, example.q, grid_size)
    exact = [func(x) for x in x_values]
    errors = [abs(got - want) for got, want in zip(reconstructed, exact)]
    max_index = max(range(len(errors)), key=lambda i: errors[i])
    return RunResult(
        n=n,
        lam=lam,
        m=m,
        max_error=errors[max_index],
        mean_abs_error=mp.fsum(errors) / len(errors),
        x_at_max=x_values[max_index],
        x_values=x_values,
        errors=errors,
    )


def computed_order(previous: Optional[mp.mpf], current: mp.mpf) -> str:
    if previous is None:
        return ""
    return sci(mp.log(previous / current, 2), 5)


def selected_rows(example: PaperExample, n_values: Sequence[int]) -> List[PaperRow]:
    wanted = set(n_values)
    return [row for row in example.rows if row.n in wanted]


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_markdown_table(path: Path, title: str, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    headers = list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(f"# {title}\n\n")
        handle.write("| " + " | ".join(headers) + " |\n")
        handle.write("| " + " | ".join(["---"] * len(headers)) + " |\n")
        for row in rows:
            handle.write("| " + " | ".join(str(row[header]) for header in headers) + " |\n")


def write_pointwise_csv(
    path: Path, example: PaperExample, results: Sequence[RunResult]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["example", "N", "lambda", "m", "x", "error", "log10_error"],
        )
        writer.writeheader()
        for result in results:
            for x, error in sorted(zip(result.x_values, result.errors), key=lambda item: item[0]):
                log_error = mp.log10(max(error, mp.mpf("1e-300")))
                writer.writerow(
                    {
                        "example": example.key,
                        "N": result.n,
                        "lambda": result.lam,
                        "m": result.m,
                        "x": sci(x, 18),
                        "error": sci(error, 18),
                        "log10_error": sci(log_error, 18),
                    }
                )


def load_font(size: int = 14) -> ImageFont.ImageFont:
    for name in ("Arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, size=size)
        except OSError:
            pass
    return ImageFont.load_default()


def nice_log_bounds(errors: Iterable[mp.mpf]) -> Tuple[int, int]:
    positive_logs = [mp.log10(error) for error in errors if error > 0]
    if not positive_logs:
        return -16, 0
    bottom = int(mp.floor(min(positive_logs)))
    top = int(mp.ceil(max(positive_logs)))
    if bottom == top:
        bottom -= 1
        top += 1
    return bottom, top


def draw_error_plot(
    *,
    title: str,
    results: Sequence[RunResult],
    width: int = 980,
    height: int = 680,
) -> Image.Image:
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font = load_font(14)
    title_font = load_font(18)
    small_font = load_font(12)
    left, right, top, bottom = 95, 230, 58, 82
    plot_left = left
    plot_right = width - right
    plot_top = top
    plot_bottom = height - bottom
    colors = ["#c13f3f", "#2c6bb1", "#258a51", "#8a4fb2", "#d28a20", "#3b8c91"]

    all_errors = [error for result in results for error in result.errors]
    log_min, log_max = nice_log_bounds(all_errors)
    y_floor = mp.power(10, log_min)

    def px(x: mp.mpf) -> int:
        return int(plot_left + (x + 1) * (plot_right - plot_left) / 2)

    def py(error: mp.mpf) -> int:
        value = mp.log10(max(error, y_floor))
        return int(plot_bottom - (value - log_min) * (plot_bottom - plot_top) / (log_max - log_min))

    draw.rectangle([plot_left, plot_top, plot_right, plot_bottom], outline="black", width=1)
    draw.text((plot_left, 20), title, fill="black", font=title_font)
    draw.text((plot_left, height - 36), "x", fill="black", font=font)
    draw.text((18, plot_top + 12), "log10 pointwise error", fill="black", font=font)

    for tick in [-1, -0.5, 0, 0.5, 1]:
        x_pos = px(mp.mpf(str(tick)))
        draw.line([x_pos, plot_bottom, x_pos, plot_bottom + 6], fill="black")
        draw.text((x_pos - 13, plot_bottom + 10), str(tick), fill="black", font=small_font)

    for exponent in range(log_min, log_max + 1):
        y_pos = py(mp.power(10, exponent))
        draw.line([plot_left - 6, y_pos, plot_left, y_pos], fill="black")
        draw.line([plot_left, y_pos, plot_right, y_pos], fill="#e8e8e8")
        draw.text((plot_left - 74, y_pos - 8), f"1e{exponent}", fill="black", font=small_font)

    for index, result in enumerate(results):
        color = colors[index % len(colors)]
        ordered = sorted(zip(result.x_values, result.errors), key=lambda item: item[0])
        points = [(px(x), py(error)) for x, error in ordered]
        if len(points) > 1:
            draw.line(points, fill=color, width=2)
        for x_pos, y_pos in points[:: max(1, len(points) // 28)]:
            draw.ellipse([x_pos - 2, y_pos - 2, x_pos + 2, y_pos + 2], fill=color)

        legend_y = plot_top + 6 + index * 24
        legend_x = plot_right + 18
        draw.line([legend_x, legend_y + 8, legend_x + 30, legend_y + 8], fill=color, width=3)
        draw.text(
            (legend_x + 38, legend_y),
            f"N={result.n}, lambda={result.lam}, m={result.m}",
            fill="black",
            font=small_font,
        )

    return image


def save_figures(output_dir: Path, plot_results: Mapping[str, Sequence[RunResult]]) -> None:
    images: List[Tuple[str, Image.Image]] = []
    for example_key, results in plot_results.items():
        example = EXAMPLES[example_key]
        image = draw_error_plot(title=f"{example.figure_label} {example.title}", results=results)
        panel = "a" if example_key == "6.1" else "b"
        path = output_dir / f"figure_1{panel}_example_{example_key.replace('.', '_')}.png"
        image.save(path)
        images.append((example_key, image))

    if len(images) == 2:
        gap = 28
        combined = Image.new(
            "RGB",
            (images[0][1].width * 2 + gap, images[0][1].height),
            "white",
        )
        combined.paste(images[0][1], (0, 0))
        combined.paste(images[1][1], (images[0][1].width + gap, 0))
        combined.save(output_dir / "figure_1_examples_6_1_6_2.png")


def run_example(
    *,
    example: PaperExample,
    n_values: Sequence[int],
    dps: int,
    method: str,
    panels: Optional[int],
    nzn: Optional[int],
    output_dir: Path,
    search_optimal: bool,
) -> Tuple[List[Mapping[str, object]], List[RunResult]]:
    set_precision(dps)
    func = paper_example_function(example.p, example.q)
    table_rows: List[Mapping[str, object]] = []
    plot_results: List[RunResult] = []
    previous_linear_error: Optional[mp.mpf] = None

    for paper_row in selected_rows(example, n_values):
        print(
            f"Computing Example {example.key}, N={paper_row.n}, "
            f"linear(lambda={paper_row.linear_lam}, m={paper_row.linear_m}), "
            f"paper-opt(lambda={paper_row.optimal_lam}, m={paper_row.optimal_m})"
        )
        coeffs = fourier_coefficients(func, paper_row.n, method=method, panels=panels)
        linear = compute_reconstruction(
            example=example,
            n=paper_row.n,
            lam=paper_row.linear_lam,
            m=paper_row.linear_m,
            coeffs=coeffs,
            nzn=nzn,
        )
        optimal_pair = compute_reconstruction(
            example=example,
            n=paper_row.n,
            lam=paper_row.optimal_lam,
            m=paper_row.optimal_m,
            coeffs=coeffs,
            nzn=nzn,
        )

        row: Dict[str, object] = {
            "N": paper_row.n,
            "linear_lambda": paper_row.linear_lam,
            "linear_m": paper_row.linear_m,
            "computed_linear_Linf": sci(linear.max_error),
            "paper_linear_Linf": sci(paper_row.linear_error),
            "linear_ratio": ratio_text(linear.max_error, paper_row.linear_error),
            "computed_order": computed_order(previous_linear_error, linear.max_error),
            "paper_order": paper_row.order or "",
            "computed_linear_L1_mean": sci(linear.mean_abs_error),
            "linear_x_at_max": sci(linear.x_at_max),
            "paper_opt_lambda": paper_row.optimal_lam,
            "paper_opt_m": paper_row.optimal_m,
            "computed_Linf_at_paper_opt_pair": sci(optimal_pair.max_error),
            "paper_opt_Linf": sci(paper_row.optimal_error),
            "paper_opt_pair_ratio": ratio_text(optimal_pair.max_error, paper_row.optimal_error),
        }

        if search_optimal:
            found = find_best_parameters(
                func,
                paper_row.n,
                example.q,
                max_lam=paper_row.n // 2,
                max_m=paper_row.n // 2,
                dps=dps,
                fourier_method=method,
                panels=panels,
                nzn=nzn,
            )
            row.update(
                {
                    "searched_opt_lambda": found.best_lam,
                    "searched_opt_m": found.best_m,
                    "searched_opt_Linf": sci(found.max_error),
                }
            )

        table_rows.append(row)
        plot_results.append(linear)
        previous_linear_error = linear.max_error

    table_stem = f"table_{example.table_number}_example_{example.key.replace('.', '_')}"
    write_csv(output_dir / f"{table_stem}.csv", table_rows)
    write_markdown_table(
        output_dir / f"{table_stem}.md",
        f"Table {example.table_number} reproduction for Example {example.key}",
        table_rows,
    )
    write_pointwise_csv(
        output_dir / f"pointwise_errors_example_{example.key.replace('.', '_')}.csv",
        example,
        plot_results,
    )
    return table_rows, plot_results


def write_summary(
    output_dir: Path,
    metadata: Mapping[str, object],
    all_rows: Mapping[str, Sequence[Mapping[str, object]]],
) -> None:
    payload = {
        "metadata": metadata,
        "tables": {
            key: [dict(row) for row in rows]
            for key, rows in all_rows.items()
        },
    }
    with (output_dir / "summary.json").open("w") as handle:
        json.dump(payload, handle, indent=2)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-md", type=Path, default=DEFAULT_PAPER_MD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--example", choices=["6.1", "6.2"], nargs="+", default=["6.1", "6.2"])
    parser.add_argument("--N", type=int, nargs="+", default=None)
    parser.add_argument("--quick", action="store_true", help="Use N=40,80 unless --N is given.")
    parser.add_argument("--full", action="store_true", help="Use N=40,80,160,320,640 unless --N is given.")
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--method", choices=["gauss3", "quad"], default="gauss3")
    parser.add_argument("--panels", type=int, default=None)
    parser.add_argument("--nzn", type=int, default=None)
    parser.add_argument(
        "--search-optimal",
        action="store_true",
        help="Also scan 1 <= lambda,m <= N/2. This is very slow for large N.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not args.paper_md.exists():
        raise FileNotFoundError(f"Paper Markdown file not found: {args.paper_md}")

    if args.N is not None:
        n_values = args.N
    elif args.full:
        n_values = FULL_N_VALUES
    else:
        n_values = QUICK_N_VALUES
    if args.quick and args.N is None and not args.full:
        n_values = QUICK_N_VALUES

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_rows: Dict[str, Sequence[Mapping[str, object]]] = {}
    all_plot_results: Dict[str, Sequence[RunResult]] = {}

    for example_key in args.example:
        rows, plot_results = run_example(
            example=EXAMPLES[example_key],
            n_values=n_values,
            dps=args.dps,
            method=args.method,
            panels=args.panels,
            nzn=args.nzn,
            output_dir=args.output_dir,
            search_optimal=args.search_optimal,
        )
        all_rows[example_key] = rows
        all_plot_results[example_key] = plot_results

    save_figures(args.output_dir, all_plot_results)
    write_summary(
        args.output_dir,
        {
            "paper_md": str(args.paper_md),
            "N_values": n_values,
            "dps": args.dps,
            "method": args.method,
            "panels": args.panels,
            "nzn": args.nzn,
            "search_optimal": args.search_optimal,
        },
        all_rows,
    )
    print(f"Wrote reproduction artifacts to: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Runnable example for conceptual fluorobenzene gas-phase nitration modeling."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from reaction import ArrheniusParameters, ReactionProperties
from reactor import ReactorParameters, ReactorResult, simulate_fixed_bed, summarize_result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Simulate isothermal and non-isothermal 1-D fixed-bed nitration models."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Directory for generated CSV files and SVG figures.",
    )
    return parser


def _scale(value: float, minimum: float, maximum: float, start: float, end: float) -> float:
    if maximum == minimum:
        return (start + end) / 2.0
    return start + (value - minimum) * (end - start) / (maximum - minimum)


def write_svg_plot(
    path: Path,
    title: str,
    x_label: str,
    y_label: str,
    series: list[tuple[str, list[float], list[float], str]],
    target_y: float | None = None,
) -> None:
    """Write a small dependency-free SVG line plot."""

    width, height = 820, 520
    left, right, top, bottom = 80, 30, 45, 70
    xs = [x for _, x_values, _, _ in series for x in x_values]
    ys = [y for _, _, y_values, _ in series for y in y_values]
    if target_y is not None:
        ys.append(target_y)
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    padding = (max_y - min_y) * 0.08 or 1.0
    min_y -= padding
    max_y += padding

    def points(x_values: list[float], y_values: list[float]) -> str:
        return " ".join(
            f"{_scale(x, min_x, max_x, left, width - right):.1f},"
            f"{_scale(y, min_y, max_y, height - bottom, top):.1f}"
            for x, y in zip(x_values, y_values)
        )

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="25" text-anchor="middle" font-size="20">{title}</text>',
        f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black"/>',
        f'<text x="{width / 2}" y="{height-20}" text-anchor="middle" font-size="14">{x_label}</text>',
        f'<text x="20" y="{height / 2}" text-anchor="middle" font-size="14" transform="rotate(-90 20 {height / 2})">{y_label}</text>',
        f'<text x="{left}" y="{height-bottom+20}" text-anchor="middle" font-size="12">{min_x:.2f}</text>',
        f'<text x="{width-right}" y="{height-bottom+20}" text-anchor="middle" font-size="12">{max_x:.2f}</text>',
        f'<text x="{left-10}" y="{height-bottom}" text-anchor="end" font-size="12">{min_y:.2f}</text>',
        f'<text x="{left-10}" y="{top+5}" text-anchor="end" font-size="12">{max_y:.2f}</text>',
    ]
    if target_y is not None:
        y = _scale(target_y, min_y, max_y, height - bottom, top)
        lines.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" stroke="black" stroke-dasharray="6,6"/>'
        )
    legend_y = 58
    for label, x_values, y_values, color in series:
        lines.append(
            f'<polyline points="{points(x_values, y_values)}" fill="none" stroke="{color}" stroke-width="2.5"/>'
        )
        lines.append(
            f'<text x="{width-210}" y="{legend_y}" font-size="13" fill="{color}">{label}</text>'
        )
        legend_y += 20
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_profiles_csv(path: Path, isothermal: ReactorResult, non_isothermal: ReactorResult) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "z_m",
                "residence_time_s",
                "conversion_isothermal",
                "temperature_isothermal_K",
                "conversion_non_isothermal",
                "temperature_non_isothermal_K",
            ]
        )
        for row in zip(
            isothermal.axial_position,
            isothermal.residence_time,
            isothermal.conversion,
            isothermal.temperature,
            non_isothermal.conversion,
            non_isothermal.temperature,
        ):
            writer.writerow([f"{value:.8g}" for value in row])


def write_outputs(isothermal: ReactorResult, non_isothermal: ReactorResult, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_profiles_csv(output_dir / "profiles.csv", isothermal, non_isothermal)
    write_svg_plot(
        output_dir / "conversion_time.svg",
        "Conversion-time comparison",
        "Residence time (s)",
        "Fluorobenzene conversion",
        [
            ("Isothermal", isothermal.residence_time, isothermal.conversion, "#1f77b4"),
            ("Non-isothermal", non_isothermal.residence_time, non_isothermal.conversion, "#d62728"),
        ],
        target_y=0.98,
    )
    write_svg_plot(
        output_dir / "temperature_profile.svg",
        "Fixed-bed temperature distribution",
        "Bed axial position (m)",
        "Temperature (K)",
        [
            ("Isothermal", isothermal.axial_position, isothermal.temperature, "#1f77b4"),
            ("Non-isothermal", non_isothermal.axial_position, non_isothermal.temperature, "#d62728"),
        ],
    )


def main() -> None:
    args = build_parser().parse_args()

    kinetics = ArrheniusParameters()
    properties = ReactionProperties(target_conversion=0.98, total_yield=0.92)
    reactor = ReactorParameters()

    isothermal = simulate_fixed_bed(reactor, kinetics, properties, mode="isothermal")
    non_isothermal = simulate_fixed_bed(reactor, kinetics, properties, mode="non_isothermal")

    for result in (isothermal, non_isothermal):
        summary = summarize_result(result, properties)
        print(
            f"{summary['mode']}: conversion={summary['outlet_conversion']:.3f} "
            f"(target {summary['target_conversion']:.2f}), "
            f"yield={summary['total_yield']:.3f} (target {summary['target_yield']:.2f}), "
            f"outlet T={summary['outlet_temperature_K']:.1f} K"
        )

    write_outputs(isothermal, non_isothermal, args.output_dir)
    print(f"Outputs saved to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()

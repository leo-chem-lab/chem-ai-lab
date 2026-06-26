"""Run the gas-phase nitration reactor simulation and save profile plots."""

from pathlib import Path
from xml.sax.saxutils import escape

from reactor import run_reactor


def _write_svg(x_values, y_values, ylabel: str, title: str, path: Path) -> None:
    width, height = 760, 470
    left, right, top, bottom = 78, 24, 42, 64
    x_min, x_max = min(x_values), max(x_values)
    y_min, y_max = min(y_values), max(y_values)
    pad = (y_max - y_min) * 0.08 or 1.0
    y_min -= pad
    y_max += pad

    def sx(x):
        return left + (x - x_min) / (x_max - x_min) * (width - left - right)

    def sy(y):
        return height - bottom - (y - y_min) / (y_max - y_min) * (height - top - bottom)

    points = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(x_values, y_values))
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="white"/>
  <text x="{width/2}" y="26" text-anchor="middle" font-family="Arial" font-size="19">{escape(title)}</text>
  <line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#222"/>
  <line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#222"/>
  <polyline points="{points}" fill="none" stroke="#1f77b4" stroke-width="3"/>
  <text x="{width/2}" y="{height-18}" text-anchor="middle" font-family="Arial" font-size="15">Reactor length (m)</text>
  <text x="18" y="{height/2}" transform="rotate(-90 18 {height/2})" text-anchor="middle" font-family="Arial" font-size="15">{escape(ylabel)}</text>
  <text x="{left}" y="{height-bottom+24}" text-anchor="middle" font-family="Arial" font-size="12">{x_min:.1f}</text>
  <text x="{width-right}" y="{height-bottom+24}" text-anchor="middle" font-family="Arial" font-size="12">{x_max:.1f}</text>
  <text x="{left-8}" y="{sy(y_min+pad):.1f}" text-anchor="end" font-family="Arial" font-size="12">{y_min+pad:.2f}</text>
  <text x="{left-8}" y="{sy(y_max-pad):.1f}" text-anchor="end" font-family="Arial" font-size="12">{y_max-pad:.2f}</text>
</svg>
"""
    path.write_text(svg)


def plot_profiles(results, output_dir: Path = Path("outputs")) -> None:
    """Create conversion, temperature, selectivity, and pressure profile plots."""

    output_dir.mkdir(exist_ok=True)
    z = results["z"]
    profiles = [
        ("conversion_vs_length.svg", [v * 100.0 for v in results["conversion"]], "Benzene conversion (%)", "Conversion vs reactor length"),
        ("temperature_profile.svg", results["temperature"], "Temperature (K)", "Temperature profile"),
        ("selectivity_profile.svg", [v * 100.0 for v in results["selectivity"]], "Selectivity to nitrobenzene (%)", "Selectivity profile"),
        ("pressure_profile.svg", [v / 100_000.0 for v in results["pressure"]], "Pressure (bar)", "Pressure profile"),
    ]
    for filename, y_values, ylabel, title in profiles:
        _write_svg(z, y_values, ylabel, title, output_dir / filename)


def main() -> None:
    results = run_reactor()
    plot_profiles(results)
    print("Industrial gas-phase nitration simulation complete.")
    print(f"Outlet conversion: {results['conversion'][-1] * 100:.2f}%")
    print(f"Outlet temperature: {results['temperature'][-1]:.2f} K")
    print(f"Outlet selectivity: {results['selectivity'][-1] * 100:.2f}%")
    print("Plots written to outputs/.")


if __name__ == "__main__":
    main()

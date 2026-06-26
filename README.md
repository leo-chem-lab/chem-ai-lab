# chem-ai-lab

A runnable Python example for a conceptual fluorobenzene gas-phase nitration
model in a one-dimensional fixed-bed reactor.

> The default kinetic and heat-transfer values are placeholders for model
> structure demonstration. They are not validated design or operating data.

## Features

- 1-D fixed-bed plug-flow style reactor simulation.
- Target assumptions for 98% conversion and 92% total yield.
- Arrhenius kinetic expression with replaceable placeholder parameters.
- Isothermal and simplified non-isothermal heat-balance comparison.
- Dependency-free SVG/CSV output for conversion vs. residence time and temperature vs. bed length.
- Modular code split into `reaction.py`, `reactor.py`, and `main.py`.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt  # optional; no third-party packages are required
python main.py --output-dir outputs
```

Generated files:

- `outputs/profiles.csv`
- `outputs/conversion_time.svg`
- `outputs/temperature_profile.svg`

## Module overview

- `reaction.py` defines Arrhenius kinetics, reaction properties, and rates.
- `reactor.py` defines fixed-bed parameters, profile simulation, and summaries.
- `main.py` runs both thermal cases and writes the example plots.

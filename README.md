# chem-ai-lab

A compact packed-bed reactor simulation for industrial gas-phase nitration. The
model includes a desired nitration pathway, NO2 side decomposition, heat release,
non-isothermal mass/energy coupling, heat transfer to a coolant, pressure drop,
and selectivity tracking.

## Run

```bash
python main.py
```

The script writes plots to `outputs/`:

- `conversion_vs_length.svg`
- `temperature_profile.svg`
- `selectivity_profile.svg`
- `pressure_profile.svg`

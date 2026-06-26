# chem-ai-labwork

A small lab for chemistry-oriented simulation utilities.

## Nitration reactor simulation

This repository includes a deterministic semi-batch reactor model for aromatic nitration studies. The model tracks substrate, nitrating acid, nitro product, water, and reactor temperature while accounting for:

- Arrhenius temperature-dependent nitration kinetics.
- Semi-batch acid dosing with configurable feed concentration and temperature.
- Heat generation from nitration and heat removal to a jacket.
- Safety diagnostics for temperature, residual acid inventory, and estimated heat-release rate.

Run the built-in example:

```bash
python -m chem_ai_lab.nitration
```

Run tests:

```bash
python -m pytest
```

"""Reaction kinetics for a conceptual fluorobenzene gas-phase nitration model.

This module intentionally uses placeholder kinetic parameters. The equations are
for reactor-model demonstration only and are not a validated process recipe.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ArrheniusParameters:
    """Arrhenius kinetic parameters for pseudo-first-order nitration."""

    pre_exponential: float = 4.2e4  # 1/s, placeholder
    activation_energy: float = 55_000.0  # J/mol, placeholder
    gas_constant: float = 8.314462618  # J/(mol K)


@dataclass(frozen=True)
class ReactionProperties:
    """Thermochemical and selectivity assumptions used by the reactor model."""

    heat_of_reaction: float = -95_000.0  # J/mol fluorobenzene converted, placeholder
    target_conversion: float = 0.98
    total_yield: float = 0.92

    @property
    def selectivity(self) -> float:
        """Overall product selectivity implied by conversion and total yield."""

        return self.total_yield / self.target_conversion


def arrhenius_rate_constant(
    temperature: float, parameters: ArrheniusParameters | None = None
) -> float:
    """Return the Arrhenius rate constant at ``temperature`` in kelvin."""

    if temperature <= 0.0:
        raise ValueError("temperature must be greater than zero kelvin")
    params = parameters or ArrheniusParameters()
    return params.pre_exponential * math.exp(
        -params.activation_energy / (params.gas_constant * temperature)
    )


def reaction_rate(
    fluorobenzene_concentration: float,
    temperature: float,
    parameters: ArrheniusParameters | None = None,
) -> float:
    """Return pseudo-first-order fluorobenzene consumption rate in mol/(m^3 s)."""

    if fluorobenzene_concentration < 0.0:
        raise ValueError("fluorobenzene_concentration cannot be negative")
    return arrhenius_rate_constant(temperature, parameters) * fluorobenzene_concentration


def product_rate(
    fluorobenzene_concentration: float,
    temperature: float,
    kinetic_parameters: ArrheniusParameters | None = None,
    reaction_properties: ReactionProperties | None = None,
) -> float:
    """Return desired nitration-product formation rate after yield/selectivity loss."""

    props = reaction_properties or ReactionProperties()
    return props.selectivity * reaction_rate(
        fluorobenzene_concentration, temperature, kinetic_parameters
    )

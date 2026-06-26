"""Reaction kinetics for an industrial gas-phase aromatic nitration model.

The model tracks benzene (A), nitrogen dioxide (B), desired nitrobenzene
product (P), side-product lump (S), temperature, and pressure.  It is a
compact teaching model rather than a plant-design package, but the equations
are arranged so rate expressions and thermodynamic data can be changed without
rewriting the reactor balance.
"""

from dataclasses import dataclass
from math import exp
from typing import Dict

R_GAS = 8.314462618  # J mol-1 K-1


@dataclass(frozen=True)
class ReactionParameters:
    """Kinetic and thermal parameters for the gas-phase nitration network."""

    # Main nitration: A + B -> P
    k0_main: float = 2.5e-5  # m3 mol-1 s-1
    ea_main: float = 68_000.0  # J mol-1
    delta_h_main: float = -118_000.0  # J mol-1, exothermic

    # NO2 side decomposition / oxidation lump: 2 B -> S
    k0_side: float = 1.8e-5  # m3 mol-1 s-1
    ea_side: float = 82_000.0  # J mol-1
    delta_h_side: float = -56_000.0  # J mol-1 of side-product extent

    reference_temperature: float = 600.0  # K, used to scale Arrhenius values


def arrhenius_rate_constant(k0: float, ea: float, temperature: float) -> float:
    """Return an Arrhenius rate constant with temperature clipping for stability."""

    safe_temperature = min(max(temperature, 250.0), 1_200.0)
    return k0 * exp(-ea / R_GAS * (1.0 / safe_temperature - 1.0 / 600.0))


def reaction_rates(concentrations: Dict[str, float], temperature: float, params: ReactionParameters) -> Dict[str, float]:
    """Calculate main and side volumetric reaction rates.

    Concentrations are expected in mol m-3. Negative concentrations can occur as
    tiny numerical artifacts in stiff integrations, so they are clipped to zero.
    """

    c_a = max(concentrations.get("benzene", 0.0), 0.0)
    c_b = max(concentrations.get("no2", 0.0), 0.0)

    k_main = arrhenius_rate_constant(params.k0_main, params.ea_main, temperature)
    k_side = arrhenius_rate_constant(params.k0_side, params.ea_side, temperature)

    return {
        "main": k_main * c_a * c_b,
        "side": k_side * c_b * c_b,
    }


def selectivity(rates: Dict[str, float]) -> float:
    """Return instantaneous selectivity to desired product P over side product S."""

    side = max(rates.get("side", 0.0), 0.0)
    main = max(rates.get("main", 0.0), 0.0)
    return main / (main + side) if main + side > 0.0 else 1.0

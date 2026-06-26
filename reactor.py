"""One-dimensional fixed-bed reactor models for fluorobenzene nitration.

The implementation compares an isothermal plug-flow approximation with a
simplified non-isothermal energy balance. Parameters are illustrative and should
be replaced with experimentally validated data before any engineering use.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from reaction import ArrheniusParameters, ReactionProperties, reaction_rate


Mode = Literal["isothermal", "non_isothermal"]


@dataclass(frozen=True)
class ReactorParameters:
    """Geometry, feed, and heat-transfer settings for a 1-D fixed bed."""

    length: float = 1.0  # m
    n_nodes: int = 200
    superficial_velocity: float = 0.25  # m/s
    inlet_concentration: float = 12.0  # mol/m^3 fluorobenzene
    inlet_temperature: float = 620.0  # K
    coolant_temperature: float = 600.0  # K
    gas_density: float = 1.2  # kg/m^3
    heat_capacity: float = 1_150.0  # J/(kg K)
    bed_void_fraction: float = 0.45
    heat_transfer_coefficient: float = 170.0  # W/(m^3 K), volumetric placeholder
    max_time: float = 4.0  # s, used for conversion-time output

    def __post_init__(self) -> None:
        if self.length <= 0.0:
            raise ValueError("length must be positive")
        if self.n_nodes < 3:
            raise ValueError("n_nodes must be at least 3")
        if self.superficial_velocity <= 0.0:
            raise ValueError("superficial_velocity must be positive")


@dataclass(frozen=True)
class ReactorResult:
    """Computed reactor profiles."""

    axial_position: list[float]
    residence_time: list[float]
    concentration: list[float]
    conversion: list[float]
    temperature: list[float]
    mode: Mode

    @property
    def outlet_conversion(self) -> float:
        return self.conversion[-1]

    @property
    def outlet_temperature(self) -> float:
        return self.temperature[-1]


def simulate_fixed_bed(
    reactor_parameters: ReactorParameters | None = None,
    kinetic_parameters: ArrheniusParameters | None = None,
    reaction_properties: ReactionProperties | None = None,
    mode: Mode = "isothermal",
) -> ReactorResult:
    """Simulate a 1-D fixed-bed reactor by explicit axial marching."""

    params = reactor_parameters or ReactorParameters()
    kinetics = kinetic_parameters or ArrheniusParameters()
    props = reaction_properties or ReactionProperties()

    dz = params.length / (params.n_nodes - 1)
    z = [i * dz for i in range(params.n_nodes)]
    tau = [position / params.superficial_velocity for position in z]
    concentration = [params.inlet_concentration] + [0.0] * (params.n_nodes - 1)
    temperature = [params.inlet_temperature] + [0.0] * (params.n_nodes - 1)

    for i in range(1, params.n_nodes):
        previous_c = max(concentration[i - 1], 0.0)
        previous_t = temperature[i - 1]
        rate = reaction_rate(previous_c, previous_t, kinetics)

        dcdz = -rate / params.superficial_velocity
        concentration[i] = max(previous_c + dcdz * dz, 0.0)

        if mode == "isothermal":
            temperature[i] = params.inlet_temperature
        elif mode == "non_isothermal":
            reaction_heat = (-props.heat_of_reaction * rate) / (
                params.gas_density * params.heat_capacity * params.superficial_velocity
            )
            heat_removal = (
                params.heat_transfer_coefficient
                * (previous_t - params.coolant_temperature)
                / (params.gas_density * params.heat_capacity * params.superficial_velocity)
            )
            temperature[i] = previous_t + (reaction_heat - heat_removal) * dz
        else:
            raise ValueError("mode must be 'isothermal' or 'non_isothermal'")

    conversion = [1.0 - c / params.inlet_concentration for c in concentration]
    return ReactorResult(z, tau, concentration, conversion, temperature, mode)


def summarize_result(result: ReactorResult, reaction_properties: ReactionProperties) -> dict[str, float | str]:
    """Return key scalar metrics for display or testing."""

    return {
        "mode": result.mode,
        "outlet_conversion": result.outlet_conversion,
        "target_conversion": reaction_properties.target_conversion,
        "total_yield": result.outlet_conversion * reaction_properties.selectivity,
        "target_yield": reaction_properties.total_yield,
        "outlet_temperature_K": result.outlet_temperature,
    }

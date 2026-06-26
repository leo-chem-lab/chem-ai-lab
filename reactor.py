"""Packed-bed plug-flow reactor balances for gas-phase nitration."""

from dataclasses import dataclass, field
from typing import Dict, List

from reaction import R_GAS, ReactionParameters, reaction_rates, selectivity


@dataclass(frozen=True)
class ReactorParameters:
    """Geometry, transport, and feed settings for the reactor model."""

    length: float = 8.0  # m
    diameter: float = 0.45  # m
    catalyst_bulk_density: float = 850.0  # kg catalyst m-3 bed
    particle_diameter: float = 0.004  # m
    bed_void_fraction: float = 0.42
    gas_viscosity: float = 3.2e-5  # Pa s
    heat_transfer_coefficient: float = 320.0  # W m-2 K-1
    coolant_temperature: float = 560.0  # K
    heat_capacity: Dict[str, float] = field(
        default_factory=lambda: {
            "benzene": 145.0,
            "no2": 48.0,
            "nitrobenzene": 170.0,
            "side_product": 70.0,
        }
    )  # J mol-1 K-1
    inlet_temperature: float = 585.0  # K
    inlet_pressure: float = 650_000.0  # Pa
    inlet_flows: Dict[str, float] = field(
        default_factory=lambda: {
            "benzene": 8.0,
            "no2": 9.2,
            "nitrobenzene": 0.0,
            "side_product": 0.0,
        }
    )  # mol s-1

    @property
    def area(self) -> float:
        return 3.141592653589793 * self.diameter**2 / 4.0

    @property
    def heat_transfer_area_per_length(self) -> float:
        return 3.141592653589793 * self.diameter


SPECIES = ("benzene", "no2", "nitrobenzene", "side_product")


def _state_to_flows(y: List[float]) -> Dict[str, float]:
    return {name: max(value, 0.0) for name, value in zip(SPECIES, y[:4])}


def reactor_odes(z: float, y: List[float], reactor: ReactorParameters, kinetics: ReactionParameters) -> List[float]:
    """Mass, energy, and pressure balances along reactor length z."""

    flows = _state_to_flows(y)
    temperature = min(max(y[4], 250.0), 1_200.0)
    pressure = max(y[5], 25_000.0)
    total_flow = max(sum(flows.values()), 1e-12)

    volumetric_flow = total_flow * R_GAS * temperature / pressure
    concentrations = {name: flow / volumetric_flow for name, flow in flows.items()}
    rates = reaction_rates(concentrations, temperature, kinetics)

    catalyst_factor = reactor.catalyst_bulk_density * reactor.area
    main_extent = rates["main"] * catalyst_factor
    side_extent = rates["side"] * catalyst_factor

    heat_released = -(kinetics.delta_h_main * main_extent + kinetics.delta_h_side * side_extent)
    heat_removed = reactor.heat_transfer_coefficient * reactor.heat_transfer_area_per_length * (
        temperature - reactor.coolant_temperature
    )
    mixture_heat_capacity = sum(flows[name] * reactor.heat_capacity[name] for name in SPECIES)
    d_temperature = (heat_released - heat_removed) / max(mixture_heat_capacity, 1e-9)

    superficial_velocity = volumetric_flow / reactor.area
    molecular_weight_avg = 0.078  # kg mol-1, representative aromatic/NO2 gas mixture
    gas_density = pressure * molecular_weight_avg / (R_GAS * temperature)
    eps = reactor.bed_void_fraction
    dp = reactor.particle_diameter
    ergun_laminar = 150.0 * (1.0 - eps) ** 2 * reactor.gas_viscosity * superficial_velocity / (eps**3 * dp**2)
    ergun_inertial = 1.75 * (1.0 - eps) * gas_density * superficial_velocity**2 / (eps**3 * dp)

    return [
        -main_extent,
        -main_extent - 2.0 * side_extent,
        main_extent,
        side_extent,
        d_temperature,
        -(ergun_laminar + ergun_inertial),
    ]


def _rk4_step(z: float, y: List[float], dz: float, reactor: ReactorParameters, kinetics: ReactionParameters) -> List[float]:
    def add_scaled(base, deriv, scale):
        return [b + scale * d for b, d in zip(base, deriv)]

    k1 = reactor_odes(z, y, reactor, kinetics)
    k2 = reactor_odes(z + dz / 2.0, add_scaled(y, k1, dz / 2.0), reactor, kinetics)
    k3 = reactor_odes(z + dz / 2.0, add_scaled(y, k2, dz / 2.0), reactor, kinetics)
    k4 = reactor_odes(z + dz, add_scaled(y, k3, dz), reactor, kinetics)
    stepped = [value + dz / 6.0 * (a + 2.0 * b + 2.0 * c + d) for value, a, b, c, d in zip(y, k1, k2, k3, k4)]
    return [max(stepped[i], 0.0) if i < 4 else stepped[i] for i in range(len(stepped))]


def _instant_selectivity(y: List[float], kinetics: ReactionParameters) -> float:
    flows = _state_to_flows(y)
    temperature = min(max(y[4], 250.0), 1_200.0)
    pressure = max(y[5], 25_000.0)
    total_flow = max(sum(flows.values()), 1e-12)
    volumetric_flow = total_flow * R_GAS * temperature / pressure
    concentrations = {name: flow / volumetric_flow for name, flow in flows.items()}
    return selectivity(reaction_rates(concentrations, temperature, kinetics))


def run_reactor(reactor: ReactorParameters | None = None, kinetics: ReactionParameters | None = None, n_steps: int = 800):
    """Integrate the non-isothermal packed-bed model and return profiles."""

    reactor = reactor or ReactorParameters()
    kinetics = kinetics or ReactionParameters()
    dz = reactor.length / n_steps
    z_values = [0.0]
    y_values = [[*(reactor.inlet_flows[name] for name in SPECIES), reactor.inlet_temperature, reactor.inlet_pressure]]

    for step in range(n_steps):
        z = step * dz
        y_values.append(_rk4_step(z, y_values[-1], dz, reactor, kinetics))
        z_values.append((step + 1) * dz)

    inlet_benzene = reactor.inlet_flows["benzene"]
    flows = {name: [row[i] for row in y_values] for i, name in enumerate(SPECIES)}
    return {
        "z": z_values,
        "flows": flows,
        "temperature": [row[4] for row in y_values],
        "pressure": [row[5] for row in y_values],
        "conversion": [min(max((inlet_benzene - row[0]) / inlet_benzene, 0.0), 1.0) for row in y_values],
        "selectivity": [_instant_selectivity(row, kinetics) for row in y_values],
    }

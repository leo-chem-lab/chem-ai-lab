"""Semi-batch nitration reactor simulation.

The model is intentionally compact and deterministic so it can be used for
screening studies and agent-driven optimization experiments. It is not a
replacement for validated process-safety calculations.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp
from typing import Iterable

GAS_CONSTANT_J_PER_MOL_K = 8.314462618


@dataclass(frozen=True)
class NitrationParameters:
    """Configurable process and kinetic parameters for nitration simulation."""

    reactor_volume_l: float = 1.0
    acid_feed_rate_mol_min: float = 0.018
    acid_feed_duration_min: float = 45.0
    acid_feed_temperature_k: float = 293.15
    pre_exponential_l_mol_min: float = 1.8e7
    activation_energy_j_mol: float = 54_000.0
    reaction_enthalpy_j_mol: float = -115_000.0
    heat_capacity_j_k: float = 4_200.0
    heat_transfer_j_min_k: float = 95.0
    jacket_temperature_k: float = 298.15
    max_temperature_k: float = 323.15
    max_heat_release_j_min: float = 2_000.0

    def validate(self) -> None:
        """Validate parameters before numerical integration."""

        positive_fields = {
            "reactor_volume_l": self.reactor_volume_l,
            "acid_feed_duration_min": self.acid_feed_duration_min,
            "pre_exponential_l_mol_min": self.pre_exponential_l_mol_min,
            "activation_energy_j_mol": self.activation_energy_j_mol,
            "heat_capacity_j_k": self.heat_capacity_j_k,
        }
        for field_name, value in positive_fields.items():
            if value <= 0:
                raise ValueError(f"{field_name} must be positive")
        if self.acid_feed_rate_mol_min < 0:
            raise ValueError("acid_feed_rate_mol_min cannot be negative")
        if self.heat_transfer_j_min_k < 0:
            raise ValueError("heat_transfer_j_min_k cannot be negative")


@dataclass(frozen=True)
class NitrationState:
    """Instantaneous reactor state."""

    time_min: float
    substrate_mol: float
    acid_mol: float
    product_mol: float
    water_mol: float
    temperature_k: float
    reaction_rate_mol_min: float
    heat_release_j_min: float

    @property
    def conversion(self) -> float:
        initial_substrate = self.substrate_mol + self.product_mol
        if initial_substrate <= 0:
            return 0.0
        return self.product_mol / initial_substrate


@dataclass(frozen=True)
class SimulationResult:
    """Full trajectory and safety summaries for a nitration run."""

    states: tuple[NitrationState, ...]
    warnings: tuple[str, ...]

    @property
    def final_state(self) -> NitrationState:
        return self.states[-1]

    @property
    def peak_temperature_k(self) -> float:
        return max(state.temperature_k for state in self.states)

    @property
    def peak_heat_release_j_min(self) -> float:
        return max(state.heat_release_j_min for state in self.states)


def _rate_constant(temperature_k: float, params: NitrationParameters) -> float:
    return params.pre_exponential_l_mol_min * exp(
        -params.activation_energy_j_mol / (GAS_CONSTANT_J_PER_MOL_K * temperature_k)
    )


def _acid_feed(time_min: float, params: NitrationParameters) -> float:
    if time_min <= params.acid_feed_duration_min:
        return params.acid_feed_rate_mol_min
    return 0.0


def _derivatives(state: NitrationState, params: NitrationParameters) -> tuple[float, ...]:
    feed = _acid_feed(state.time_min, params)
    acid_concentration = max(state.acid_mol, 0.0) / params.reactor_volume_l
    substrate_concentration = max(state.substrate_mol, 0.0) / params.reactor_volume_l
    rate = _rate_constant(state.temperature_k, params) * acid_concentration * substrate_concentration
    rate = min(rate * params.reactor_volume_l, max(state.substrate_mol, 0.0), max(state.acid_mol + feed, 0.0))
    heat_release = -params.reaction_enthalpy_j_mol * rate
    heat_removal = params.heat_transfer_j_min_k * (state.temperature_k - params.jacket_temperature_k)
    feed_sensible_heat = feed * 75.0 * (params.acid_feed_temperature_k - state.temperature_k)
    d_temperature = (heat_release - heat_removal + feed_sensible_heat) / params.heat_capacity_j_k
    return (-rate, feed - rate, rate, rate, d_temperature, rate, heat_release)


def simulate_nitration(
    initial_substrate_mol: float = 1.0,
    initial_acid_mol: float = 0.02,
    initial_temperature_k: float = 298.15,
    duration_min: float = 120.0,
    step_min: float = 0.25,
    params: NitrationParameters | None = None,
) -> SimulationResult:
    """Simulate a semi-batch nitration reactor with fourth-order Runge-Kutta.

    Returns the full state trajectory plus warnings when process limits are
    exceeded. Inputs are mol, kelvin, minutes, and joules.
    """

    params = params or NitrationParameters()
    params.validate()
    if initial_substrate_mol <= 0:
        raise ValueError("initial_substrate_mol must be positive")
    if initial_acid_mol < 0:
        raise ValueError("initial_acid_mol cannot be negative")
    if duration_min <= 0 or step_min <= 0:
        raise ValueError("duration_min and step_min must be positive")

    current = NitrationState(
        time_min=0.0,
        substrate_mol=initial_substrate_mol,
        acid_mol=initial_acid_mol,
        product_mol=0.0,
        water_mol=0.0,
        temperature_k=initial_temperature_k,
        reaction_rate_mol_min=0.0,
        heat_release_j_min=0.0,
    )
    states = [current]

    while current.time_min < duration_min:
        dt = min(step_min, duration_min - current.time_min)
        current = _rk4_step(current, params, dt)
        states.append(current)

    warnings = _build_warnings(states, params)
    return SimulationResult(tuple(states), tuple(warnings))


def _rk4_step(state: NitrationState, params: NitrationParameters, dt: float) -> NitrationState:
    def shifted(base: NitrationState, delta: Iterable[float], scale: float) -> NitrationState:
        ds, da, dp, dw, dtemp, rate, heat = delta
        return NitrationState(
            time_min=base.time_min + dt * scale,
            substrate_mol=max(base.substrate_mol + ds * dt * scale, 0.0),
            acid_mol=max(base.acid_mol + da * dt * scale, 0.0),
            product_mol=max(base.product_mol + dp * dt * scale, 0.0),
            water_mol=max(base.water_mol + dw * dt * scale, 0.0),
            temperature_k=base.temperature_k + dtemp * dt * scale,
            reaction_rate_mol_min=rate,
            heat_release_j_min=heat,
        )

    k1 = _derivatives(state, params)
    k2 = _derivatives(shifted(state, k1, 0.5), params)
    k3 = _derivatives(shifted(state, k2, 0.5), params)
    k4 = _derivatives(shifted(state, k3, 1.0), params)
    combined = tuple((a + 2 * b + 2 * c + d) / 6 for a, b, c, d in zip(k1, k2, k3, k4))
    ds, da, dp, dw, dtemp, rate, heat = combined
    return NitrationState(
        time_min=state.time_min + dt,
        substrate_mol=max(state.substrate_mol + ds * dt, 0.0),
        acid_mol=max(state.acid_mol + da * dt, 0.0),
        product_mol=max(state.product_mol + dp * dt, 0.0),
        water_mol=max(state.water_mol + dw * dt, 0.0),
        temperature_k=state.temperature_k + dtemp * dt,
        reaction_rate_mol_min=max(rate, 0.0),
        heat_release_j_min=max(heat, 0.0),
    )


def _build_warnings(states: list[NitrationState], params: NitrationParameters) -> list[str]:
    warnings: list[str] = []
    peak_temperature = max(state.temperature_k for state in states)
    peak_heat_release = max(state.heat_release_j_min for state in states)
    final_acid = states[-1].acid_mol
    if peak_temperature > params.max_temperature_k:
        warnings.append(
            f"Peak temperature {peak_temperature:.2f} K exceeds limit {params.max_temperature_k:.2f} K"
        )
    if peak_heat_release > params.max_heat_release_j_min:
        warnings.append(
            f"Peak heat release {peak_heat_release:.0f} J/min exceeds limit {params.max_heat_release_j_min:.0f} J/min"
        )
    if final_acid > 0.05:
        warnings.append(f"Residual acid inventory remains high at {final_acid:.3f} mol")
    return warnings


def _format_state(state: NitrationState) -> str:
    return (
        f"t={state.time_min:6.1f} min | T={state.temperature_k:6.2f} K | "
        f"conversion={state.conversion:5.1%} | acid={state.acid_mol:6.3f} mol | "
        f"q={state.heat_release_j_min:7.1f} J/min"
    )


def main() -> None:
    result = simulate_nitration()
    for state in result.states[:: max(1, len(result.states) // 10)]:
        print(_format_state(state))
    print(_format_state(result.final_state))
    if result.warnings:
        print("Warnings:")
        for warning in result.warnings:
            print(f"- {warning}")


if __name__ == "__main__":
    main()

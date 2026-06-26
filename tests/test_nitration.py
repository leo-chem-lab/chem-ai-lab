import pytest

from chem_ai_lab.nitration import NitrationParameters, simulate_nitration


def test_simulation_converts_substrate_and_tracks_heat():
    result = simulate_nitration(duration_min=60, step_min=0.5)

    assert result.final_state.product_mol > 0.1
    assert result.final_state.substrate_mol < 1.0
    assert result.peak_heat_release_j_min > 0
    assert result.final_state.conversion == pytest.approx(
        result.final_state.product_mol
        / (result.final_state.product_mol + result.final_state.substrate_mol)
    )


def test_safety_warnings_include_temperature_excursion():
    params = NitrationParameters(
        acid_feed_rate_mol_min=0.08,
        heat_transfer_j_min_k=5.0,
        max_temperature_k=300.0,
    )

    result = simulate_nitration(duration_min=30, step_min=0.25, params=params)

    assert any("Peak temperature" in warning for warning in result.warnings)


def test_validation_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="initial_substrate_mol"):
        simulate_nitration(initial_substrate_mol=0)

    with pytest.raises(ValueError, match="reactor_volume_l"):
        simulate_nitration(params=NitrationParameters(reactor_volume_l=0))

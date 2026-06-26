"""Chem AI Lab simulation utilities."""

__all__ = [
    "NitrationParameters",
    "NitrationState",
    "SimulationResult",
    "simulate_nitration",
]


def __getattr__(name: str):
    if name in __all__:
        from . import nitration

        return getattr(nitration, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

"""Simple reactor utilities for running example reactions."""

from reaction import Reaction


class Reactor:
    """A minimal reactor that formats reaction results."""

    def __init__(self, temperature_celsius: float = 25.0) -> None:
        self.temperature_celsius = temperature_celsius

    def run(self, reaction: Reaction) -> str:
        """Run a reaction and return a short status message."""
        return (
            f"Ran {reaction.name} at {self.temperature_celsius:.1f} °C: "
            f"{reaction.equation()}"
        )

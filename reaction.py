"""Reaction domain models for the Chem AI Lab example project."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Reaction:
    """A minimal representation of a chemical reaction."""

    name: str
    reactants: tuple[str, ...]
    products: tuple[str, ...]

    def equation(self) -> str:
        """Return a human-readable reaction equation."""
        return f"{' + '.join(self.reactants)} -> {' + '.join(self.products)}"

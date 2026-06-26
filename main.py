"""Minimal runnable example for the Chem AI Lab project."""

from reaction import Reaction
from reactor import Reactor


def build_example_reaction() -> Reaction:
    """Create the demonstration reaction used by the CLI entry point."""
    return Reaction(
        name="Water formation",
        reactants=("2 H2", "O2"),
        products=("2 H2O",),
    )


def main() -> None:
    """Run the minimal working example."""
    reactor = Reactor(temperature_celsius=25.0)
    reaction = build_example_reaction()
    print(reactor.run(reaction))


if __name__ == "__main__":
    main()

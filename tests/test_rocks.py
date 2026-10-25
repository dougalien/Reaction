"""Rocks and glasses as reactants (stoichiometry only; no thermodynamic data needed)."""
import pytest

from reaction.app import rocks
from reaction.app.assistant import interpret_rules
from reaction.app.scenario import Scenario
from reaction.thermo.eosdb import ATOMIC_WEIGHTS


def test_a_gram_of_basalt_glass_weighs_a_gram():
    el, oxygen = rocks.elements_per_gram(rocks.ROCKS["basalt_glass"])
    mass = sum(n * ATOMIC_WEIGHTS["Fe" if e == "Fe3" else e] for e, n in el.items())
    mass += sum(oxygen.values()) * ATOMIC_WEIGHTS["O"]
    # FeO is total iron: the ferric part carries the extra oxygen of Fe2O3 on top of the gram
    assert mass == pytest.approx(1.0 + 0.5 * el["Fe3"] * ATOMIC_WEIGHTS["O"], rel=1e-9)
    assert el["Fe3"] / (el["Fe"] + el["Fe3"]) == pytest.approx(0.10)
    assert el["Si"] == pytest.approx(0.505 / 60.0843 / 1.0 * 100 / sum(rocks.ROCKS["basalt_glass"].oxides.values()), rel=1e-6)


def test_names():
    for name in ("basalt glass", "Basalt glass", "MORB", "basaltic glass", "basalt_glass"):
        assert rocks.find(name).key == "basalt_glass", name
    assert rocks.find("albite") is None


def test_basalt_glass_from_the_chat():
    got = interpret_rules("react 100 g basalt glass with seawater at 350 C and 500 bar", Scenario())
    card, problems = Scenario().changed(got.scenario_changes)
    assert not problems and card.reactants == [{"name": "basalt glass", "amount": 100.0, "unit": "g"}]
    assert card.t_start_c == card.t_end_c == 350.0 and got.run

"""Ideal solid solutions: a(endmember) = x^sites, formed in place of the pure endmembers."""
import math

import pytest

from reaction.app import config
from reaction.app.assistant import interpret_rules
from reaction.app.engine import Engine
from reaction.app.narrate import story
from reaction.app.scenario import Scenario
from reaction.equilibrium import System, equilibrate
from reaction.solidsolution import by_name, table


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


def test_table():
    names = {s.name for s in table()}
    assert {"plagioclase", "carbonate", "chlorite", "olivine"} <= names
    assert by_name()["chlorite"].sites == 5 and by_name()["garnet"].endmembers == ("grossular", "andradite")


def test_endmembers_exist(engine):
    missing = [e for s in table() for e in s.endmembers if e not in engine.db.species]
    assert not missing


def test_mass_action_and_mass_balance(engine):
    """Ca-Mn carbonate from a supersaturated fluid: SI of each endmember = sites * log x, and
    the elements add up."""
    system = System(engine.db, ["Na+", "Cl-", "Ca+2", "Mn+2", "HCO3-"], solid_solutions=[by_name()["carbonate"]])
    bulk = {"H2O": 55.508, "H+": 0.0, "Na+": 0.1, "Cl-": 0.122, "Ca+2": 0.01, "Mn+2": 0.001, "HCO3-": 0.02}
    st = equilibrate(system, 25.0, bulk)
    assert not st.minerals and set(st.solid_solutions) == {"carbonate"}
    members = st.solid_solutions["carbonate"]
    total = sum(members.values())
    for name, moles in members.items():
        assert st.speciation.saturation[name] == pytest.approx(math.log10(moles / total), abs=1e-6)
    mn_left = st.water_kg * st.speciation.component_totals["Mn+2"]
    assert mn_left + members["rhodochrosite"] == pytest.approx(0.001, rel=1e-6)
    assert members["rhodochrosite"] / total > 0.1          # Mn goes into calcite strongly
    assert st.solids["carbonate"] == pytest.approx(total)


def test_without_the_switch_minerals_stay_pure(engine):
    system = System(engine.db, ["Na+", "Cl-", "Ca+2", "Mn+2", "HCO3-"])
    bulk = {"H2O": 55.508, "H+": 0.0, "Na+": 0.1, "Cl-": 0.122, "Ca+2": 0.01, "Mn+2": 0.001, "HCO3-": 0.02}
    st = equilibrate(system, 25.0, bulk)
    assert not st.solid_solutions and st.minerals


def test_greenschist_assemblage(engine):
    card, problems = Scenario().changed({"t_start_c": 350, "t_end_c": 350, "p_start_bar": 500, "p_end_bar": 500,
                                         "steps": 8, "reactants": "500 g basalt", "solid_solutions": True})
    assert not problems and "with solid solutions" in card.describe()
    result = engine.run(card)
    run = result["modes"]["static"]
    end = {k: v[-1] for k, v in run["minerals"].items() if v[-1]}
    assert {"chlorite", "epidote-clinozoisite", "actinolite", "albite"} <= set(end)
    assert not {"clinochlore", "epidote", "tremolite"} & set(end)          # taken over by the mixtures
    x = run["composition"]["chlorite:daphnite"][-1]
    assert 0.05 < x < 0.95
    text = " ".join(story(result))
    assert "clinochlore" in text and "% at the end" in text
    assert any("mix ideally" in n for n in result["notes"])


def test_flow_through_keeps_mass(engine):
    card = Scenario().changed({"t_start_c": 350, "t_end_c": 350, "p_start_bar": 500, "p_end_bar": 500,
                               "steps": 5, "reactants": "200 g basalt", "solid_solutions": True,
                               "mode": "flow_through"})[0]
    run = engine.run(card)["modes"]["flow_through"]
    chlorite = [v or 0.0 for v in run["minerals"].get("chlorite", [])]
    assert chlorite and chlorite[-1] > 0 and all(b >= a * (1 - 1e-6) for a, b in zip(chlorite, chlorite[1:]))


def test_chat():
    ch = interpret_rules("react seawater with basalt at 350 C with solid solutions", Scenario()).scenario_changes
    assert ch["solid_solutions"] is True
    assert interpret_rules("no solid solutions", Scenario()).scenario_changes["solid_solutions"] is False

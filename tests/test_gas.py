"""The free gas phase: Peng-Robinson fugacity coefficients, gas solubility, degassing, and methane held
apart from sulfate.  Needs the thermodynamic data (with OBIGT organic_aq.csv and organic_gas.csv)."""
import pytest

from reaction.app import config
from reaction.app.engine import Engine
from reaction.app.narrate import story
from reaction.app.scenario import Scenario
from reaction.equilibrium import System, equilibrate
from reaction.thermo.gas import log_phi


def test_peng_robinson_pure_gases():
    # CO2 at 50 C and 100 bar: phi 0.63 from the Span-Wagner equation (NIST); CH4 0.88 at 100 C, 200 bar
    assert 10 ** log_phi(["CO2"], [1.0], 323.15, 100.0)[0] == pytest.approx(0.63, abs=0.03)
    assert 10 ** log_phi(["CH4"], [1.0], 373.15, 200.0)[0] == pytest.approx(0.88, abs=0.03)
    assert 10 ** log_phi(["H2"], [1.0], 298.15, 1.0)[0] == pytest.approx(1.0, abs=1e-3)
    assert log_phi(["Xe"], [1.0], 300.0, 100.0)[0] == 0.0             # unknown gases stay ideal


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "organic_gas.csv").exists():
        pytest.skip("thermodynamic data with methane not found; set REACTION_DATA_DIR")
    return Engine()


@pytest.mark.parametrize("t_c, p_bar, duan", [(25, 50, 1.13), (50, 100, 0.98), (50, 200, 1.21),
                                               (100, 100, 0.72)])
def test_co2_solubility_in_water_against_duan_and_sun(engine, t_c, p_bar, duan):
    """CO2-saturated water (Duan and Sun 2003, mol/kg) within 10% to 200 bar."""
    system = System(engine.db, ["HCO3-", "Na+", "Cl-"], gas_phase=True)
    bulk = {"H2O": 55.508 - 5, "H+": 5.0, "HCO3-": 5.0, "Na+": 1e-6, "Cl-": 1e-6}
    st = equilibrate(system, t_c, bulk, pressure_bar=p_bar)
    assert st.speciation.species["CO2"].molality == pytest.approx(duan, rel=0.10)
    y_h2o = st.gas["H2O(g)"] / sum(st.gas.values())
    assert 0.0 < y_h2o < 0.03


def test_no_gas_phase_when_undersaturated_or_switched_off(engine):
    bulk = {"H2O": 55.508 - 0.1, "H+": 0.1, "HCO3-": 0.1, "Na+": 1e-6, "Cl-": 1e-6}
    on = System(engine.db, ["HCO3-", "Na+", "Cl-"], gas_phase=True)
    assert equilibrate(on, 50, bulk, pressure_bar=100).gas == {}
    off = System(engine.db, ["HCO3-", "Na+", "Cl-"])
    rich = dict(bulk, **{"H2O": 50.5, "H+": 5.0, "HCO3-": 5.0})
    assert equilibrate(off, 50, rich, pressure_bar=100).gas == {}


def test_methane_held_apart_from_sulfate(engine):
    base = {"H2O": 55.508, "H+": 0.0, "HCO3-": 2e-3, "SO4-2": 28e-3, "Na+": 0.526, "Cl-": 0.47,
            "O2(g)": 0.0}
    apart = System(engine.db, ["HCO3-", "SO4-2", "Na+", "Cl-", "O2(g)"], gas_phase=True, hold_apart=("CH4",))
    st = equilibrate(apart, 100, dict(base, CH4=0.01), pressure_bar=50)
    assert st.speciation.species["CH4"].molality == pytest.approx(0.01, rel=1e-3)
    assert st.speciation.component_totals["SO4-2"] == pytest.approx(28e-3, rel=1e-3)
    joined = System(engine.db, ["HCO3-", "SO4-2", "Na+", "Cl-", "O2(g)"], gas_phase=True)
    b = dict(base)
    b["HCO3-"] += 0.01; b["H+"] += 0.01; b["H2O"] += 0.01; b["O2(g)"] -= 0.02    # the same CH4
    st = equilibrate(joined, 100, b, pressure_bar=50)
    sp = st.speciation.species                       # the methane reduces as much sulfate to sulfide
    assert sp["CH4"].molality < 1e-6
    assert sp["H2S"].molality + sp["HS-"].molality == pytest.approx(0.01, rel=0.05)


def test_decompression_degasses_and_the_story_says_so(engine):
    card = Scenario().changed({"fluid": "custom", "analysis": "Na 0.5\nCl 0.5\nHCO3 1.0", "units": "molal",
                               "ph": 5.0, "fluid_t_c": 150, "t_start_c": 150, "t_end_c": 150,
                               "p_start_bar": 300, "p_end_bar": 20, "steps": 8, "mode": "both"})[0]
    result = engine.run(card)
    static, flow = result["modes"]["static"]["gas"], result["modes"]["flow_through"]["gas"]
    assert all(v is not None and v > 0 for v in static["CO2(g)"])
    # at 20 bar and 150 C the gas is about a quarter steam (vapour pressure 4.8 bar)
    y = static["H2O(g)"][-1] / (static["H2O(g)"][-1] + static["CO2(g)"][-1])
    assert y == pytest.approx(0.24, abs=0.04)
    # flow-through vents the gas as it forms: the cumulative CO2 is about what the static run holds
    assert flow["CO2(g)"][-1] == pytest.approx(static["CO2(g)"][-1], rel=0.05)
    assert result["catalog"]["gas"]["CO2(g)"]["label"] == "CO₂ gas"
    text = " ".join(story(result))
    assert "A gas phase forms from the start" in text and "is vented" in text


def test_gases_as_reactants_and_methane_card_switch(engine):
    base = {"fluid": "seawater", "t_start_c": 25, "t_end_c": 150, "p_start_bar": 100, "p_end_bar": 100,
            "steps": 6, "mode": "static", "reactants": [{"name": "methane gas", "amount": 0.2, "unit": "mol"}]}
    apart = engine.run(Scenario().changed(base)[0])
    react = engine.run(Scenario().changed({**base, "methane_sulfate": "react"})[0])
    so4_apart = apart["modes"]["static"]["species"]["SO4-2"][-1]
    so4_react = react["modes"]["static"]["species"]["SO4-2"][-1]
    assert so4_apart > 1e-3 and so4_react < 1e-6
    assert apart["modes"]["static"]["gas"]["CH4(g)"][-1] > 0.05
    co2 = engine.run(Scenario().changed({**base, "t_end_c": 25, "reactants": [
        {"name": "CO2(g)", "amount": 2, "unit": "mol"}]})[0])
    assert co2["modes"]["static"]["gas"]["CO2(g)"][-1] > 0.5


def test_card_fields_for_the_gas_phase():
    card, problems = Scenario().changed({"gas_phase": "no", "methane_sulfate": "equilibrate"})
    assert not problems and card.gas_phase is False and card.methane_sulfate == "react"
    _, problems = Scenario().changed({"methane_sulfate": "sometimes"})
    assert problems

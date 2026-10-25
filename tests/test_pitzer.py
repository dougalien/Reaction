"""Pitzer activity coefficients, brines and evaporation."""
import math

import numpy as np
import pytest

from reaction.app import config
from reaction.app.assistant import interpret_rules
from reaction.app.engine import Engine
from reaction.app.narrate import story
from reaction.app.scenario import Scenario
from reaction.equilibrium import System, equilibrate
from reaction.pitzer import Pitzer, at_temperature

A_PHI_25 = 0.3915


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


@pytest.mark.parametrize("m, gamma, phi", [(0.1, 0.778, 0.932), (1.0, 0.657, 0.936), (3.0, 0.714, 1.045),
                                           (6.0, 0.986, 1.271)])
def test_sodium_chloride_against_robinson_and_stokes(m, gamma, phi):
    p = Pitzer(["Na+", "Cl-"], np.array([1.0, -1.0]), 25.0, A_PHI_25)
    mm = np.array([m, m])
    ln_g, _ = p.compute(mm)
    assert math.exp(ln_g.mean()) == pytest.approx(gamma, rel=0.01)
    assert p.osmotic(mm) == pytest.approx(phi, rel=0.005)


@pytest.mark.parametrize("m, gamma, phi", [(0.1, 0.518, 0.854), (1.0, 0.500, 1.046), (3.0, 1.483, 1.769)])
def test_calcium_chloride_mixing_charges(m, gamma, phi):
    p = Pitzer(["Ca+2", "Cl-"], np.array([2.0, -1.0]), 25.0, A_PHI_25)
    mm = np.array([m, 2 * m])
    ln_g, _ = p.compute(mm)
    assert math.exp((ln_g[0] + 2 * ln_g[1]) / 3) == pytest.approx(gamma, rel=0.02)
    assert p.osmotic(mm) == pytest.approx(phi, rel=0.01)


def test_gibbs_duhem():
    """d ln a_w and the solutes' d ln gamma must agree: sum m_i d ln(m_i gamma_i) = -55.51 d ln a_w."""
    names = ["Na+", "Mg+2", "Cl-", "SO4-2", "CO2"]
    p = Pitzer(names, np.array([1.0, 2.0, -1.0, -2.0, 0.0]), 60.0, 0.41)
    m = np.array([2.0, 0.5, 2.4, 0.3, 0.05])
    dm = m * 1e-6
    g0, w0 = p.compute(m)
    g1, w1 = p.compute(m + dm)
    lhs = np.sum(m * (np.log(m + dm) - np.log(m) + g1 - g0))
    assert lhs == pytest.approx(-55.508 * (w1 - w0), rel=1e-3)


def test_temperature_function():
    # NaCl beta0 from its PHREEQC terms: 0.0753 at 25 C, about 0.10 at 100 C (Pitzer and others, 1984)
    a = (7.534e-2, 9598.4, 35.48, -5.8731e-2, 1.798e-5, -5e5)
    assert at_temperature(a, 298.15) == pytest.approx(0.07534)
    assert 0.09 < at_temperature(a, 373.15) < 0.11


def test_seawater_water_activity(engine):
    sp, _ = engine.starting_fluid(Scenario().changed({"activity_model": "pitzer"})[0])
    assert 10 ** sp.log_a_water == pytest.approx(0.9814, abs=0.0005)


@pytest.mark.parametrize("t_c, molal", [(25.0, 6.15), (100.0, 6.71)])
def test_halite_solubility(engine, t_c, molal):
    system = System(engine.db, ["Na+", "Cl-"], activity_model="pitzer")
    st = equilibrate(system, t_c, {"H2O": 55.508, "H+": 0.0, "Na+": 10.0, "Cl-": 10.0},
                     minerals=["halite"], eligible=["halite"])
    assert st.speciation.totals["Na+"] == pytest.approx(molal, rel=0.02)


def test_seawater_evaporation_sequence(engine):
    card, problems = Scenario().changed({"fluid": "seawater", "temperature_c": 25, "pressure_bar": 1, "steps": 31,
                                         "evaporate": 40, "activity_model": "pitzer"})
    assert not problems
    result = engine.run(card)
    run = result["modes"]["static"]
    factor = result["axes"]["evaporation"]["values"]
    assert factor[0] == 1 and factor[-1] == pytest.approx(40)

    def first(name):
        series = run["minerals"][name]
        return factor[next(i for i, v in enumerate(series) if v and v > 1e-6)]

    assert 3.0 < first("gypsum") < 6.0                       # about 3.8 in seawater
    assert 9.5 < first("halite") < 12.5                      # about 10.6
    assert run["properties"]["log_aw"][-1] < math.log10(0.8)
    text = " ".join(story(result))
    assert "times concentrated" in text and "Pitzer" in " ".join(result["notes"])


def test_brine_without_pitzer_is_flagged(engine):
    card = Scenario().changed({"fluid": "seawater", "temperature_c": 25, "pressure_bar": 1, "steps": 4,
                               "evaporate": 4})[0]
    notes = " ".join(engine.run(card)["notes"])
    assert "choose Pitzer" in notes


def test_card_and_chat():
    card, problems = Scenario().changed({"evaporate": "20", "activity_model": "pitzer"})
    assert not problems and card.evaporating and "evaporated to 20 times" in card.describe()
    assert Scenario().changed({"evaporate": 0.5})[0].problems()
    ch = interpret_rules("evaporate seawater to 50x with pitzer at 25 C", Scenario()).scenario_changes
    assert ch["evaporate"] == 50 and ch["activity_model"] == "pitzer" and ch["fluid"] == "seawater"
    assert ch["temperature_c"] == 25 and "t_end_c" not in ch

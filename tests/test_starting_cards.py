"""The ready-made cards and their fluids on the real engine.  Needs the thermodynamic data."""
import pytest

from reaction.app import config
from reaction.app.engine import Engine
from reaction.app.fluids import PRESETS
from reaction.app.scenario import STARTING_CARDS, Scenario, starting_card

ATMOSPHERE_LOG_FCO2 = -3.38                 # 420 ppm


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


@pytest.fixture(scope="module")
def runs(engine):
    # the column card (stream) runs for a minute; test_clays.py runs it shortened
    return {key: engine.run(starting_card(key)) for key, c in STARTING_CARDS.items()
            if c.changes.get("kind") != "transport"}


def test_fresh_water_and_brine_sit_at_the_atmospheres_co2(engine):
    fresh, _ = engine.starting_fluid(Scenario(fluid="freshwater"))
    assert fresh.log_fugacity["CO2(g)"] == pytest.approx(ATMOSPHERE_LOG_FCO2, abs=0.02)
    brine, notes = engine.starting_fluid(Scenario(fluid="gsl_brine", activity_model="pitzer"))
    assert brine.log_fugacity["CO2(g)"] == pytest.approx(ATMOSPHERE_LOG_FCO2, abs=0.05)
    assert -0.3 < brine.saturation["halite"] < 0                       # close to halite saturation


def test_hydrothermal_fluid_is_the_end_of_the_seawater_heating_run(engine, runs):
    """The preset was copied from the default card's last point; if the engine changes, so must it."""
    run = runs["seawater"]["modes"]["static"]
    entries = dict(line.split() for line in PRESETS["hydrothermal"].analysis.splitlines())
    names = {"Na": "Na+", "K": "K+", "Ca": "Ca+2", "Mg": "Mg+2", "Cl": "Cl-", "S": "SO4-2", "C": "HCO3-"}
    for element, species in names.items():
        assert float(entries[species]) / 1000 == pytest.approx(run["elements"][element][-1], rel=2e-3), element
    assert float(entries["O2"]) / 1000 == pytest.approx(run["species"]["O2"][-1], rel=2e-3)
    assert PRESETS["hydrothermal"].ph == pytest.approx(run["properties"]["pH"][-1], abs=0.005)
    hot, _ = engine.starting_fluid(Scenario(fluid="hydrothermal"))
    for mineral in ("anhydrite", "brucite"):                          # saturated, as it left them
        assert abs(hot.saturation[mineral]) < 0.01


def test_every_starting_card_runs_and_gives_its_minerals(runs):
    expect = {"seawater": {"anhydrite", "brucite"}, "freshwater": {"calcite", "anhydrite"},
              "brine": {"halite", "anhydrite", "epsomite"}, "hydrothermal": {"anhydrite"}}
    for key, result in runs.items():
        run = result["modes"]["static"]
        assert all(v is not None for v in run["properties"]["pH"]), key
        formed = {m for m, v in run["minerals"].items() if max(v) > 0}
        assert expect[key] <= formed, (key, formed)
    mixed = runs["hydrothermal"]["modes"]["static"]["minerals"]["anhydrite"]
    assert max(mixed) > 1e-3 and mixed[-1] < max(mixed) / 10           # forms in the hot mix, dissolves as it cools

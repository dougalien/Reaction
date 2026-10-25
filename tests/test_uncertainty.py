"""Uncertainty bands: the run repeated with every log K shifted at random by its uncertainty."""
import numpy as np
import pytest

from reaction.app import config
from reaction.app.assistant import interpret_rules
from reaction.app.engine import Engine, _bands
from reaction.app.scenario import Scenario
from reaction.equilibrium import System
from reaction.thermo import validity


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


def test_sigma_follows_the_reliability_check(engine):
    db = engine.db
    cool = validity.logk_sigma(db, 25.0, 1.0)
    assert cool == pytest.approx(validity.BASE_SIGMA, abs=0.02)
    near_critical = validity.logk_sigma(db, 400.0, 300.0)       # low density: at least the caution width
    assert near_critical >= validity.CAUTION


def test_perturb_shifts_log_k_and_undoes(engine):
    system = System(engine.db, ["Na+", "Cl-", "Ca+2", "HCO3-"])
    m0 = system.model(100.0, 100.0)
    net = system.net
    system.perturb(np.ones(len(net.dep_names)), np.full(len(net.phase_names), -1.0), lambda t, p: 0.2)
    m1 = system.model(100.0, 100.0)
    assert np.allclose(m1.logk_dep - m0.logk_dep, 0.2) and np.allclose(m1.logk_ph - m0.logk_ph, -0.2)
    system.perturb()
    assert np.allclose(system.model(100.0, 100.0).logk_dep, m0.logk_dep)


def test_bands_take_in_the_line():
    run = {"properties": {"pH": [7.0, None, 6.0]}, "minerals": {"calcite": [0.0, 1.0, 2.0]}, "warnings": []}
    samples = [{"properties": {"pH": [7.0 + d, 6.5, 6.0 - d]}, "minerals": {}} for d in (0.1, 0.2, 0.3, -0.1)]
    bands = _bands(run, samples)
    lo, hi = bands["properties"]["pH"]
    assert lo[1] is None and lo[0] < 7.0 < hi[0] and hi[0] <= 7.3
    assert bands["minerals"]["calcite"] == [[0.0, 0.0, 0.0], [0.0, 1.0, 2.0]]   # missing in samples = none


def test_seawater_heating_band(engine):
    card, problems = Scenario().changed({"t_end_c": 250, "steps": 6, "uncertainty_runs": 6})
    assert not problems and "uncertainty range from 6 runs" in card.describe()
    result = engine.run(card)
    run = result["modes"]["static"]
    assert run["band_runs"] == 6
    ph, (lo, hi) = run["properties"]["pH"], run["bands"]["properties"]["pH"]
    assert all(a <= v <= b for a, v, b in zip(lo, ph, hi))
    assert 0.05 < hi[-1] - lo[-1] < 1.5                  # about 0.1 log unit on every log K
    assert any("Shaded ranges" in n for n in result["notes"])
    again = engine.run(card)["modes"]["static"]["bands"]["properties"]["pH"]
    assert again == [lo, hi]                             # the same card gives the same range


def test_card_and_chat():
    assert Scenario().changed({"uncertainty_runs": 3})[0].problems()
    assert Scenario().changed({"uncertainty_runs": "off"})[0].uncertainty_runs == 0
    ch = interpret_rules("heat seawater to 300 C with uncertainty bands", Scenario()).scenario_changes
    assert ch["uncertainty_runs"] == 20
    assert interpret_rules("error bars from 50 runs", Scenario()).scenario_changes["uncertainty_runs"] == 50
    assert interpret_rules("no uncertainty", Scenario()).scenario_changes["uncertainty_runs"] == 0


def test_kinetics_chat_words():
    # these patterns once held literal backspaces instead of word boundaries
    ch = interpret_rules("react with 1 g calcite instantly", Scenario()).scenario_changes
    assert ch["kinetics"] is False
    assert interpret_rules("use measured rates for 10 g albite", Scenario()).scenario_changes["kinetics"] is True

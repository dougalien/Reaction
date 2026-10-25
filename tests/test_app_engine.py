"""The heating run on the real engine.  Needs the thermodynamic data (see app/config.py)."""
import math

import pytest

from reaction.app import config
from reaction.app.engine import Engine
from reaction.app.narrate import story
from reaction.app.scenario import Scenario


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


@pytest.fixture(scope="module")
def result(engine):
    card = Scenario(t_start_c=25.0, t_end_c=250.0, steps=10, p_start_bar=500.0, p_end_bar=500.0, mode="both")
    return engine.run(card)


def test_both_modes_are_returned_on_the_same_grid(result):
    assert set(result["modes"]) == {"static", "flow_through"}
    assert len(result["x"]["values"]) == 10
    for run in result["modes"].values():
        assert all(v is not None for v in run["properties"]["pH"])
        assert len(run["species"]["Mg+2"]) == 10


def test_anhydrite_appears_on_heating_near_150c(result):
    temps = result["x"]["values"]
    amounts = result["modes"]["static"]["minerals"]["anhydrite"]
    first = next(i for i, v in enumerate(amounts) if v > 1e-6)
    assert 100 <= temps[first] <= 175                 # Bischoff and Seyfried (1978): about 150 C
    sat = result["modes"]["static"]["saturation"]["anhydrite"]
    assert sat[first - 1] < 0 and abs(sat[first]) < 1e-6   # undersaturated before, at saturation after


def test_ph_falls_on_heating_and_blocked_carbonates_stay_out(result):
    ph = result["modes"]["static"]["properties"]["pH"]
    assert ph[0] > 7.5 and ph[-1] < ph[0] - 2
    minerals = result["modes"]["static"]["minerals"]
    assert "dolomite" not in minerals and "magnesite" not in minerals
    assert max(result["modes"]["static"]["saturation"]["dolomite"]) > 0      # supersaturated yet blocked


def test_mass_is_conserved_in_the_static_run(result):
    """Sulfur dissolved plus sulfur in anhydrite equals the sulfur the seawater started with."""
    run = result["modes"]["static"]
    start = result["start_elements"]["S"]
    for i in range(10):
        dissolved = run["elements"]["S"][i] * run["properties"]["water_kg"][i]
        solid = run["minerals"].get("anhydrite", [0.0] * 10)[i] + run["minerals"].get("kieserite", [0.0] * 10)[i]
        assert dissolved + solid == pytest.approx(start, rel=2e-4)


def test_flow_through_keeps_what_static_redissolves(result):
    st, ft = result["modes"]["static"]["minerals"], result["modes"]["flow_through"]["minerals"]
    # calcite forms on warming and static heating redissolves it as the fluid turns acid;
    # in flow-through it was left behind and stays
    assert ft["calcite"][-1] > 0
    assert st.get("calcite", [0.0])[-1] < ft["calcite"][-1]


def test_story_and_catalog_are_consistent(result):
    text = " ".join(story(result))
    assert "anhydrite" in text
    for mode_run in result["modes"].values():
        for name in mode_run["species"]:
            assert name in result["catalog"]["species"]
        for table in ("minerals", "saturation"):
            for name in mode_run[table]:
                assert name in result["catalog"]["minerals"]
    assert math.isfinite(result["elapsed_s"])


def test_a_bad_card_is_refused(engine):
    with pytest.raises(ValueError, match="boils"):
        engine.run(Scenario(p_start_bar=100.0, p_end_bar=100.0))


def test_flow_through_runs_above_the_critical_temperature(engine):
    """Leaving minerals behind must not need a saturation pressure, which does not exist above 374 C."""
    card = Scenario(t_start_c=300.0, t_end_c=450.0, steps=4, p_start_bar=1000.0, p_end_bar=1000.0,
                    mode="flow_through")
    run = engine.run(card)["modes"]["flow_through"]
    assert all(v is not None for v in run["properties"]["pH"])
    assert run["minerals"]["anhydrite"][-1] > 0


def test_a_failed_step_leaves_nothing_behind_halfway(engine, monkeypatch):
    """When a step is retried in two halves, the halfway point is only a starting guess: minerals are
    left behind at the steps of the card and nowhere else."""
    from reaction.app import engine as engine_module
    from reaction.speciation import SpeciationError

    card = Scenario(t_start_c=25.0, t_end_c=350.0, steps=8, p_start_bar=500.0, p_end_bar=500.0,
                    mode="flow_through")
    target = card.points()[5][0]
    real_solve, real_row = engine_module.equilibrate, engine_module._composition_row
    failed, removals = [], []

    def flaky(system, t_c, *args, **kwargs):
        if abs(t_c - target) < 1e-9 and not failed:
            failed.append(t_c)
            raise SpeciationError("simulated convergence failure")
        return real_solve(system, t_c, *args, **kwargs)

    def counted(*args):
        removals.append(args[1])
        return real_row(*args)

    monkeypatch.setattr(engine_module, "_composition_row", counted)
    plain = engine.run(card)["modes"]["flow_through"]
    plain_removals = len(removals)
    removals.clear()
    monkeypatch.setattr(engine_module, "equilibrate", flaky)
    retried = engine.run(card)["modes"]["flow_through"]
    assert failed
    assert len(removals) == plain_removals
    assert set(retried["minerals"]) == set(plain["minerals"])
    for name, amounts in plain["minerals"].items():
        assert retried["minerals"][name] == pytest.approx(amounts, rel=1e-6, abs=1e-12)

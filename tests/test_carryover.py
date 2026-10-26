"""Going on from the fluid at the end of the last run: "take the result of the heating run and cool it to 25"."""
import pytest

from reaction.app import carryover, config
from reaction.app.assistant import interpret_rules
from reaction.app.engine import Engine
from reaction.app.scenario import Scenario
from reaction.app.session import Session

from test_app import FakeEngine, ScriptedModel, fake_result, wait

DOUGLAS = "take the result of the seawater heating fluid and cool it back down to 25"


@pytest.fixture
def session(tmp_path):
    engine = FakeEngine()
    s = Session(tmp_path, engine, ScriptedModel({}, broken=True))
    s.engine_stub = engine
    return s


def run(s: Session, text: str, calls: int):
    s.handle_chat(text)
    assert wait(lambda: len(s.engine_stub.calls) >= calls and s.run_state["state"] == "idle")
    return s.engine_stub.calls[calls - 1] if len(s.engine_stub.calls) >= calls else None


def heated(s: Session):
    card = run(s, "heat seawater to 350 C at 500 bar, static", 1)
    assert (card["fluid"], card["t_start_c"], card["t_end_c"]) == ("seawater", 2.0, 350.0)


def test_rules_read_a_cooling_of_the_last_runs_fluid():
    got = interpret_rules(DOUGLAS, Scenario())
    assert got.from_result and got.direction == "cool" and got.run
    assert got.scenario_changes == {"t_end_c": 25.0}          # "seawater" names the old run, not a new fluid
    assert interpret_rules("cool it back down to 25.", Scenario()).scenario_changes == {"t_end_c": 25.0}
    assert interpret_rules("cool to 500 bar", Scenario()).scenario_changes == {"pressure_bar": 500.0}
    assert not interpret_rules("heat seawater to 300 C", Scenario()).from_result


def test_the_fluid_at_the_end_of_the_run_becomes_the_starting_fluid(session):
    heated(session)
    card = run(session, DOUGLAS, 2)
    assert card["fluid"] == "custom"
    assert (card["t_start_c"], card["t_end_c"], card["p_start_bar"], card["p_end_bar"]) == (350.0, 25.0, 500.0, 500.0)
    assert (card["fluid_t_c"], card["fluid_p_bar"], card["units"], card["charge_balance"]) == (350.0, 500.0, "molal", "none")
    end = fake_result(Scenario(**{**card, "fluid": "seawater"}))["modes"]["static"]
    assert card["ph"] == end["properties"]["pH"][-1]
    assert f"SO4-2 {end['elements']['S'][-1]:.6g}" in card["analysis"].splitlines()
    assert "Took the fluid at the end of the last run (350 °C, 500 bar" in " ".join(c["text"] for c in session.chat[-3:])


def test_cool_never_heats_a_fresh_fluid(session):
    heated(session)
    card = run(session, "cool to 25 C", 2)                     # no "result" words: only the last run is hotter
    assert (card["fluid"], card["t_start_c"], card["t_end_c"]) == ("custom", 350.0, 25.0)


def test_cooling_with_nothing_hotter_is_refused_not_run(session):
    session.handle_chat("static, cool to 25 C")
    assert session.engine_stub.calls == []
    assert "would heat the fluid from 2 to 25 °C" in session.chat[-1]["text"]


def test_asking_for_the_last_fluid_without_a_run_is_refused(session):
    session.handle_chat("static. " + DOUGLAS)
    assert session.engine_stub.calls == []
    assert "no run yet" in session.chat[-1]["text"]


def test_naming_a_fresh_fluid_afterwards_starts_it_fresh(session):
    heated(session)
    run(session, DOUGLAS, 2)
    card = run(session, "heat seawater to 300 C", 3)
    assert (card["fluid"], card["t_start_c"], card["t_end_c"]) == ("seawater", 2.0, 300.0)


def test_the_language_model_cannot_put_seawater_back(tmp_path):
    engine = FakeEngine()
    s = Session(tmp_path, engine, ScriptedModel({"scenario_changes": {"fluid": "seawater", "t_end_c": 25},
                                                 "run": True}))
    s.engine_stub = engine
    heated(s)
    card = run(s, DOUGLAS, 2)
    assert (card["fluid"], card["t_start_c"], card["t_end_c"]) == ("custom", 350.0, 25.0)


@pytest.mark.skipif(not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists(),
                    reason="thermodynamic data not found; set REACTION_DATA_DIR")
def test_the_carried_fluid_is_the_same_fluid_on_the_real_engine():
    engine = Engine()
    hot = engine.run(Scenario(t_start_c=2.0, t_end_c=350.0, steps=8, mode="static"))
    changes, _ = carryover.fluid_from_result(hot)
    card, problems = Scenario().changed({**changes, "t_end_c": 25.0, "mode": "static"})
    assert not problems and not card.problems()
    cooled = engine.run(card)
    a, b = hot["modes"]["static"], cooled["modes"]["static"]
    assert b["properties"]["pH"][0] == pytest.approx(a["properties"]["pH"][-1], abs=1e-4)
    assert b["properties"]["log_fO2"][0] == pytest.approx(a["properties"]["log_fO2"][-1], abs=1e-4)
    for el, values in a["elements"].items():
        assert b["elements"][el][0] == pytest.approx(values[-1], rel=1e-5)
    assert cooled["axes"]["temperature_c"]["values"][-1] == 25.0

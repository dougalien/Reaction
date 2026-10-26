"""A reactant named in the chat without an amount: "react seawater with olivine (Fo88) at W/R=1 at 300C"."""
import pytest

from reaction.app import config, rocks
from reaction.app.assistant import interpret_rules
from reaction.app.engine import Engine
from reaction.app.scenario import Scenario
from reaction.app.session import Session

from test_app import FakeEngine, ScriptedModel, wait

DOUGLAS = "react seawater with olivine ((Fo88) at W/R=1 at 300C and 350 bar"


@pytest.fixture
def session(tmp_path):
    engine = FakeEngine()
    s = Session(tmp_path, engine, ScriptedModel({}, broken=True))
    s.engine_stub = engine
    return s


def test_a_water_rock_ratio_sets_the_amount_and_the_temperature_is_held():
    got = interpret_rules(DOUGLAS, Scenario())
    assert got.scenario_changes == {"fluid": "seawater", "reactants": "1000 g olivine fo88",
                                    "temperature_c": 300.0, "pressure_bar": 350.0}
    assert not got.question
    got = interpret_rules("react seawater with basalt glass at a water to rock ratio of 10 at 350 C", Scenario())
    assert got.scenario_changes["reactants"] == "100 g basalt glass"
    got = interpret_rules("react seawater with forsterite, water:rock 2:1, at 250 C", Scenario())
    assert got.scenario_changes["reactants"] == "500 g forsterite"


def test_a_reactant_without_an_amount_is_asked_about_not_dropped(session):
    session.handle_chat("static. react seawater with olivine at 300 C")
    assert session.engine_stub.calls == []
    assert "How much olivine?" in session.chat[-1]["text"]


def test_other_with_phrases_are_not_reactants():
    for text in ("heat seawater to 300 C with isotopes", "react the fluid with 50 meq of exchange sites at 25 C",
                 "mix seawater with vent fluid at 350 C"):
        got = interpret_rules(text, Scenario())
        assert "reactants" not in got.scenario_changes and not got.question, text


def test_olivine_by_its_forsterite_content():
    for name in ("olivine fo88", "olivine (Fo88)", "Fo88", "fo88 olivine"):
        assert rocks.find(name).key == "olivine_fo88"
    assert rocks.find("olivine").label == "Olivine (Fo90)"
    assert rocks.find("forsterite") is None                  # the mineral, from the database
    el, _ = rocks.elements_per_gram(rocks.get("olivine_fo88"))
    assert el["Mg"] / el["Fe"] == pytest.approx(88 / 12)
    assert el["Si"] / (el["Mg"] + el["Fe"]) == pytest.approx(0.5)


def test_the_reply_says_what_the_card_holds_not_what_the_model_said(tmp_path):
    engine = FakeEngine()
    s = Session(tmp_path, engine, ScriptedModel({"reply": "Reacting seawater with olivine (Fo88).",
                                                 "scenario_changes": {"t_end_c": 300}, "run": True}))
    s.handle_chat("static. react seawater with olivine at 300 C")
    assert "Reacting seawater with olivine" not in s.chat[-1]["text"]
    assert "How much olivine?" in s.chat[-1]["text"] and engine.calls == []


def test_an_unchanged_card_is_not_run_again(session):
    session.handle_chat("heat seawater to 300 C at 350 bar, static")
    assert wait(lambda: len(session.engine_stub.calls) == 1 and session.run_state["state"] == "idle")
    session.handle_chat("run the olivine reaction")
    assert "Nothing on the card changed" in session.chat[-1]["text"]
    session.handle_chat("run it again")
    assert wait(lambda: len(session.engine_stub.calls) == 2 and session.run_state["state"] == "idle")


@pytest.mark.skipif(not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists(),
                    reason="thermodynamic data not found; set REACTION_DATA_DIR")
def test_seawater_serpentinizes_fo88_on_the_real_engine():
    card, problems = Scenario().changed(interpret_rules(DOUGLAS, Scenario()).scenario_changes)
    card.steps = 8
    assert not problems and not card.problems()
    run = Engine().run(card)["modes"]["static"]
    formed = {k for k, v in run["minerals"].items() if v[-1] > 0}
    assert {"brucite", "magnetite"} <= formed
    assert any(k in formed for k in ("antigorite", "chrysotile", "lizardite"))

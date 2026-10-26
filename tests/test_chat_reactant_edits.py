"""Changing reactants in the chat, as people say it (Douglas's transcript of 2026-10-26 and others)."""
import pytest

from reaction.app import config, rocks
from reaction.app.assistant import interpret, interpret_rules
from reaction.app.engine import Engine
from reaction.app.scenario import Scenario
from reaction.app.session import Session

from test_app import FakeEngine, ScriptedModel, wait

HAS_DATA = (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists()
needs_data = pytest.mark.skipif(not HAS_DATA, reason="thermodynamic data not found; set REACTION_DATA_DIR")
CLAIMS = ScriptedModel({"reply": "Reacting hydrothermal fluid with forsterite (Fo88) in a static run. Running now.",
                        "run": True})


def edit(text, reactants=None, **card):
    sc = Scenario(**card)
    if reactants:
        sc, _ = sc.changed({"reactants": reactants})
    got = interpret_rules(text, sc)
    new, problems = sc.changed(got.scenario_changes)
    assert not problems
    return [(r["amount"], r["name"]) for r in new.reactants], got


def test_douglas_transcript_reads_each_message():
    have, got = edit("React the hydrothermal fluid with forsterite(fo88) at a W/R=1 at 350C and 500 bar. Static.")
    assert have == [(1000.0, "forsterite fo88")] and got.scenario_changes["fluid"] == "hydrothermal"
    assert got.scenario_changes["temperature_c"] == 350.0 and got.run
    assert rocks.find("forsterite fo88").key == "olivine_fo88"
    assert edit("react Forsterite then", "1000 g forsterite fo88")[0] == [(1000.0, "forsterite")]
    assert edit("React olivine then", "1000 g forsterite fo88")[0] == [(1000.0, "olivine")]
    assert edit("Do not react Forsterite Fo88 then. React Forsterite instead.",
                "1000 g forsterite fo88")[0] == [(1000.0, "forsterite")]


@pytest.mark.parametrize("text, before, after", [
    ("now use basalt glass instead", "100 g olivine fo88", [(100.0, "basalt glass")]),
    ("switch to quartz", "100 g olivine", [(100.0, "quartz")]),
    ("replace calcite with albite", "100 g basalt glass; 10 g calcite", [(100.0, "basalt glass"), (10.0, "albite")]),
    ("add 10 g calcite too", "100 g basalt glass", [(100.0, "basalt glass"), (10.0, "calcite")]),
    ("remove the calcite", "100 g basalt glass; 10 g calcite", [(100.0, "basalt glass")]),
    ("no reactants", "100 g basalt glass", []),
    ("same but at a water/rock ratio of 5", "10 g quartz", [(200.0, "quartz")]),
    ("react it with 10 g quartz at 250 C", "100 g basalt glass; 10 g calcite", [(10.0, "quartz")]),
    ("react seawater with 100 g olivine (Fo88) at 300 C", None, [(100.0, "olivine fo88")]),
    ("react 1 kg basalt glass with seawater at 300 C", None, [(1000.0, "basalt glass")]),
    ("dissolve 5 g calcite in rainwater at 25 C", None, [(5.0, "calcite")]),
    ("react seawater with 20 g anorthite and 10 g diopside at 300 C", None, [(20.0, "anorthite"), (10.0, "diopside")]),
])
def test_everyday_reactant_phrasings(text, before, after):
    assert edit(text, before)[0] == after


@pytest.mark.parametrize("text", ["heat seawater to 300 C with isotopes", "use a linear axis", "hide the calcite",
                                  "mix seawater with vent fluid at 350 C", "plot pH and anhydrite",
                                  "react the fluid with 50 meq of exchange sites at 25 C"])
def test_other_messages_leave_the_reactants_alone(text):
    have, got = edit(text, "100 g basalt glass; 10 g calcite")
    assert have == [(100.0, "basalt glass"), (10.0, "calcite")] and "reactants" not in got.scenario_changes


def test_the_model_never_speaks_for_the_run():
    got = interpret("React olivine then", Scenario().changed({"reactants": "1000 g forsterite fo88"})[0], [], [], CLAIMS)
    assert "Running now" not in got.reply and "Reacting" not in got.reply
    assert got.scenario_changes["reactants"] == "1000 g olivine"
    got = interpret("what does W/R mean?", Scenario(), [], [], CLAIMS)
    assert got.reply == "" and not got.run


def test_session_follows_the_transcript(tmp_path):
    engine = FakeEngine()
    s = Session(tmp_path, engine, CLAIMS)
    for i, (text, name) in enumerate([
            ("React the hydrothermal fluid with forsterite(fo88) at a W/R=1 at 350C and 500 bar. Static.", "forsterite fo88"),
            ("react Forsterite then", "forsterite"), ("React olivine then", "olivine")], start=1):
        s.handle_chat(text)
        assert wait(lambda: len(engine.calls) == i and s.run_state["state"] == "idle"), s.chat[-1]["text"]
        assert engine.calls[-1]["reactants"][0]["name"] == name
        assert s.chat[-2]["text"].count("Running:") == 1 and "Running now" not in s.chat[-2]["text"]


@needs_data
def test_an_unknown_name_is_refused_before_running_and_the_old_reactant_stays(tmp_path):
    s = Session(tmp_path, Engine(), None)
    s.handle_chat("static, react seawater with 10 g albite at 200 C")
    assert wait(lambda: s.current and s.run_state["state"] == "idle", 120)
    runs = len(s.list_runs())
    s.handle_chat("react it with unobtanium instead")
    assert "not a mineral, rock or species" in s.chat[-1]["text"] and "Running now" not in s.chat[-1]["text"]
    assert s.card.reactants[0]["name"] == "albite" and len(s.list_runs()) == runs
    s.handle_chat("react plagioclase instead")
    assert "mineral group" in s.chat[-1]["text"] and s.card.reactants[0]["name"] == "albite"


@needs_data
def test_hydrothermal_fluid_reacts_with_forsterite_fo88():
    sc, _ = Scenario().changed(interpret_rules(
        "React the hydrothermal fluid with forsterite(fo88) at a W/R=1 at 350C and 500 bar. Static.",
        Scenario()).scenario_changes)
    sc.steps = 6
    assert not sc.problems()
    run = Engine().run(sc)["modes"]["static"]
    formed = {k for k, v in run["minerals"].items() if v[-1] > 0}
    assert formed and run["properties"]["pH"][-1] is not None
    serp, _ = Scenario().changed({"reactants": "50 g serpentine", "temperature_c": 300.0})
    assert Engine()._reactant_names(serp) == {"serpentine": "chrysotile"}

"""A table of everyday chat messages and what each must do to the card, including negations.

Each row: the card to start from, the message, the card fields that must result, and whether a run starts.
The session uses the stub engine and no language model, so this tests the reading of the words alone.
"""
import pytest

from reaction.app.scenario import Scenario, starting_card
from reaction.app.session import Session

from test_app import FakeEngine, ScriptedModel, wait

OLIVINE = [{"name": "olivine", "amount": 1000.0, "unit": "g"}]


def hydrothermal_card():
    """The Hydrothermal fluid starting card (it mixes in 9 kg of 2 C seawater), with olivine on it."""
    card = starting_card("hydrothermal")
    card.reactants = [dict(r) for r in OLIVINE]
    return card


def plain(**kw):
    return Scenario(mode="static", **kw)


ROWS = [
    # Douglas, 2026-10-26: a starting card's mixing must not ride along, and "not ... mix" must remove it
    (starting_card("hydrothermal"), "react the hydrothermal fluid with olivine at 300C and 500 bar at a W/R=1",
     {"fluid": "hydrothermal", "mix_with": "", "reactants": OLIVINE, "t_start_c": 300.0, "t_end_c": 300.0,
      "p_start_bar": 500.0}, True),
    (hydrothermal_card(), "We are not supposed to mix with seawater. Just the hydrothermal fluid recting with olivine",
     {"fluid": "hydrothermal", "mix_with": "", "reactants": OLIVINE}, True),
    (hydrothermal_card(), "don't mix with seawater", {"mix_with": "", "reactants": OLIVINE}, False),
    (hydrothermal_card(), "no mixing", {"mix_with": "", "reactants": OLIVINE}, False),
    (hydrothermal_card(), "without any seawater mixing, react it at 300 C", {"mix_with": "", "reactants": OLIVINE}, True),
    (hydrothermal_card(), "stop mixing in seawater", {"mix_with": ""}, False),
    (plain(), "mix seawater with vent fluid at 350 C", {"fluid": "seawater", "mix_with": "vent_fluid", "mix_t_c": 350.0}, True),
    (plain(), "mix 9 kg of seawater into the hydrothermal fluid", {"mix_with": "seawater"}, True),
    # negated fluids and reactants
    (plain(), "react seawater, not the hydrothermal fluid, with 10 g calcite",
     {"fluid": "seawater", "reactants": [{"name": "calcite", "amount": 10.0, "unit": "g"}]}, True),
    (plain(reactants=OLIVINE), "we should not react olivine", {"reactants": []}, False),
    (plain(reactants=OLIVINE), "do not add calcite", {"reactants": OLIVINE}, False),
    (plain(reactants=OLIVINE), "no calcite please", {"reactants": OLIVINE}, False),
    (plain(reactants=OLIVINE), "without olivine", {"reactants": []}, False),
    (plain(reactants=OLIVINE), "no reactants", {"reactants": []}, False),
    (plain(reactants=OLIVINE), "drop the olivine", {"reactants": []}, False),
    (plain(reactants=OLIVINE), "don't run yet, but use basalt glass instead",
     {"reactants": [{"name": "basalt glass", "amount": 1000.0, "unit": "g"}]}, False),
    # a named fluid starts a new experiment: earlier processes go, the reply says so
    (hydrothermal_card(), "heat seawater to 300 C", {"fluid": "seawater", "mix_with": "", "reactants": []}, True),
    (plain(evaporate=5.0), "react rainwater with 10 g calcite at 25 C",
     {"fluid": "rainwater", "evaporate": 1.0, "t_start_c": 25.0}, True),
    (hydrothermal_card(), "cool it to 25 C", {"mix_with": "seawater", "reactants": OLIVINE}, True),   # no fluid named
    # reactant edits
    (plain(reactants=OLIVINE), "react forsterite instead",
     {"reactants": [{"name": "forsterite", "amount": 1000.0, "unit": "g"}]}, True),
    (plain(reactants=OLIVINE), "also add 5 g calcite",
     {"reactants": OLIVINE + [{"name": "calcite", "amount": 5.0, "unit": "g"}]}, True),
    (plain(reactants=OLIVINE), "W/R 10", {"reactants": [{"name": "olivine", "amount": 100.0, "unit": "g"}]}, True),
    (plain(), "react seawater with olivine (Fo88) at W/R=1 at 300C and 350 bar",
     {"reactants": [{"name": "olivine fo88", "amount": 1000.0, "unit": "g"}], "t_start_c": 300.0, "p_end_bar": 350.0},
     True),
    (plain(), "react seawater with olivine at 300 C", {"reactants": []}, False),       # asks how much
    # paths
    (plain(), "heat seawater from 2 to 350 C at 500 bar", {"t_start_c": 2.0, "t_end_c": 350.0, "p_end_bar": 500.0}, True),
    (plain(t_start_c=350.0, t_end_c=350.0), "decompress from 1000 to 200 bar at 300 C",
     {"t_start_c": 300.0, "p_start_bar": 1000.0, "p_end_bar": 200.0}, True),
    (plain(), "don't run, just set the end temperature to 200 C", {"t_end_c": 200.0}, False),
    (plain(), "evaporate seawater 10 times at 25 C", {"evaporate": 10.0, "t_start_c": 25.0, "t_end_c": 25.0}, True),
    (plain(evaporate=10.0), "no evaporation", {"evaporate": 1.0}, False),
    (plain(), "heat seawater to 300 C without isotopes", {"isotopes": False}, True),
    (plain(isotopes=True), "turn off isotopes", {"isotopes": False}, False),
]


@pytest.mark.parametrize("card, text, expect, runs", ROWS, ids=[r[1][:60] for r in ROWS])
def test_phrasing(tmp_path, card, text, expect, runs):
    engine = FakeEngine()
    s = Session(tmp_path, engine, ScriptedModel(broken=True))
    s.card, s.mode_confirmed = card, True
    s.handle_chat(text)
    if runs:
        assert wait(lambda: engine.calls and s.run_state["state"] == "idle"), s.chat[-1]["text"]
    else:
        assert engine.calls == [], s.chat[-1]["text"]
    reply = s.chat[-1]["text"] if not runs else next(c["text"] for c in reversed(s.chat) if "Running:" in c["text"])
    got = s.card.to_dict()
    for key, value in expect.items():
        assert got[key] == value, f"{key}: {got[key]!r} != {value!r}; reply: {reply}"
    assert "Running now" not in reply
    if runs:
        assert s.card.describe() in reply                  # the reply names exactly what runs


def test_cleared_processes_are_named_in_the_reply(tmp_path):
    s = Session(tmp_path, FakeEngine(), ScriptedModel(broken=True))
    s.card, s.mode_confirmed = hydrothermal_card(), True
    s.handle_chat("heat seawater to 300 C")
    reply = next(c["text"] for c in reversed(s.chat) if "Running:" in c["text"])
    assert "Cleared from the last card: 1000 g olivine; mixing with 9 kg of seawater per kg." in reply


def test_the_model_cannot_add_fields_it_was_not_offered(tmp_path):
    engine = FakeEngine()
    model = ScriptedModel({"reply": "Mixing with seawater. Running now.", "scenario_changes": {"mix_with": "seawater"},
                           "run": True})
    s = Session(tmp_path, engine, model)
    s.card, s.mode_confirmed = starting_card("hydrothermal"), True
    s.handle_chat("react the hydrothermal fluid with olivine at 300C and 500 bar at a W/R=1")
    assert wait(lambda: engine.calls and s.run_state["state"] == "idle")
    assert engine.calls[-1]["mix_with"] == "" and engine.calls[-1]["reactants"] == OLIVINE
    reply = next(c["text"] for c in reversed(s.chat) if "Running:" in c["text"])
    assert "Running now" not in reply and "Cleared from the last card: mixing with 9 kg of seawater per kg." in reply

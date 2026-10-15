"""The scenario card, chat interpretation, shared view, session flow and local server.

These tests use a stub engine and a scripted language model, so they need no thermodynamic data and
no Ollama.  ``test_app_engine.py`` runs the real engine when the data is present.
"""
import http.client
import json
import threading
import time

import pytest

from reaction.app import view as viewmod
from reaction.app.assistant import Interpretation, Ollama, interpret, interpret_rules
from reaction.app.engine import _pack
from reaction.app.names import build_index, pretty_mineral, pretty_species, resolve
from reaction.app.narrate import story
from reaction.app.scenario import Scenario
from reaction.app.server import make_server
from reaction.app.session import Session


# ---- helpers ------------------------------------------------------------------------------------
def fake_result(card: Scenario, run_id: str = "20260101-000000-test") -> dict:
    n = 6
    temps = [card.t_start_c + i * (card.t_end_c - card.t_start_c) / (n - 1) for i in range(n)]
    rows = [{"pH": 8.0 - 0.5 * i, "ionic_strength": 0.6, "water_kg": 1.0, "log_aw": -0.02,
             "species": {"Ca+2": 0.010 - 0.001 * i, "CaCl+": 0.001 * i, "Mg+2": 0.05, "SO4-2": 0.028 - 0.004 * i},
             "elements": {"Ca": 0.0102, "Mg": 0.053 - 0.004 * i, "S": 0.028 - 0.004 * i},
             "minerals": {"anhydrite": 0.0 if i < 2 else 0.002 * i, "dolomite": 0.0},
             "saturation": {"anhydrite": -1.0 + 0.5 * i, "dolomite": 0.8}, "warnings": []}
            for i in range(n)]
    modes = ["static", "flow_through"] if card.mode == "both" else [card.mode]
    runs = {m: _pack(rows) for m in modes}
    for run in runs.values():
        run.pop("warnings")
    catalog = {
        "species": {k: {"label": pretty_species(k, z), "charge": z}
                    for k, z in (("Ca+2", 2), ("CaCl+", 1), ("Mg+2", 2), ("SO4-2", -2))},
        "elements": {k: {"label": k} for k in ("Ca", "Mg", "S")},
        "minerals": {"anhydrite": {"label": "Anhydrite", "blocked": False},
                     "dolomite": {"label": "Dolomite", "blocked": True}},
    }
    return {"id": run_id, "created": "2026-01-01T00:00:00", "scenario": card.to_dict(),
            "summary": card.describe(), "database": "stub",
            "fluid": {"key": "seawater", "label": "Seawater", "description": "", "measured_t_c": 25.0,
                      "source": "", "totals_molal": {}},
            "x": {"key": "temperature_c", "label": "Temperature", "unit": "°C", "values": temps},
            "start_elements": {"Ca": 0.0102, "Mg": 0.053, "S": 0.028}, "catalog": catalog,
            "modes": runs, "notes": ["stub"], "elapsed_s": 0.0}


class FakeEngine:
    _db = object()

    def __init__(self):
        self.run_lock = threading.Lock()
        self.calls = []

    def status(self):
        return {"ready": True, "message": "stub"}

    def run(self, card, progress=None):
        self.calls.append(card.to_dict())
        if progress:
            progress(0.5, "halfway")
            progress(1.0, "done")
        return fake_result(card, f"20260101-00000{len(self.calls)}-test")


class ScriptedModel(Ollama):
    """Answers like a small model that sometimes invents values."""

    def __init__(self, answer: dict | None = None, broken: bool = False):
        super().__init__("http://unused", "stub")
        self.answer, self.broken, self.seen = answer or {}, broken, []

    def available(self, timeout=1.5):
        return not self.broken

    def chat(self, messages, timeout=120.0):
        self.seen.append(messages)
        if self.broken:
            raise OSError("down")
        base = {"reply": "ok", "scenario_changes": {}, "show": [], "hide": [], "axis": "unchanged",
                "run": False, "question": ""}
        return json.dumps({**base, **self.answer})


def wait(pred, timeout=10.0):
    end = time.time() + timeout
    while time.time() < end:
        if pred():
            return True
        time.sleep(0.02)
    return False


# ---- names --------------------------------------------------------------------------------------
def test_pretty_names_and_tolerant_matching():
    assert pretty_species("Ca+2", 2) == "Ca²⁺"
    assert pretty_species("HCO3-", -1) == "HCO₃⁻"
    assert pretty_species("CaSO4", 0) == "CaSO₄(aq)"
    assert pretty_mineral("dolomite,ordered") == "Dolomite (ordered)"
    names = {"Ca+2": pretty_species("Ca+2", 2), "CaCl+": "CaCl⁺", "SO4-2": pretty_species("SO4-2", -2)}
    index = build_index(names)
    for typed in ("Ca++", "ca2+", "Ca²⁺", "Ca+2"):
        assert resolve(typed, index) == "Ca+2"
    assert resolve("cacl+", index) == "CaCl+"
    assert resolve("SO4--", index) == "SO4-2"
    assert resolve("unobtainium", index) is None


# ---- scenario card ------------------------------------------------------------------------------
def test_default_card_is_valid_and_described():
    card = Scenario()
    assert card.problems() == []
    assert "350" in card.describe() and "500 bar" in card.describe()


def test_card_edits_are_validated_field_by_field():
    card, problems = Scenario().changed({"t_end_c": "300", "mode": "Flow-through", "steps": 20.4,
                                         "nonsense": 1, "kind": "mixing"})
    assert (card.t_end_c, card.mode, card.steps) == (300.0, "flow_through", 20)
    assert len(problems) == 2                                    # unknown field, kind not editable
    bad, problems = Scenario().changed({"t_end_c": "hot"})
    assert problems and bad.t_end_c == 350.0                     # unchanged


def test_card_problems_cover_boiling_and_ranges():
    low_p = Scenario(pressure_bar=100.0)                          # water boils at 350 C below ~165 bar
    assert any("boils" in p for p in low_p.problems())
    assert Scenario(t_start_c=300, t_end_c=100).problems()
    assert Scenario(mode="sideways").problems()
    assert Scenario(fluid="mud").problems()
    assert Scenario(pressure_bar=100.0, t_end_c=250).problems() == []


def test_blocked_minerals_accept_names_with_commas():
    card, _ = Scenario().changed({"blocked_minerals": "dolomite,ordered, magnesite"})
    assert card.blocked_minerals == ["dolomite,ordered", "magnesite"]


# ---- reading the chat ---------------------------------------------------------------------------
def test_rules_read_temperature_pressure_mode():
    got = interpret_rules("heat seawater to 300 C at 400 bar, static", Scenario())
    assert got.scenario_changes == {"t_end_c": 300.0, "pressure_bar": 400.0, "mode": "static", "fluid": "seawater"}
    assert got.run and got.mode_mentioned
    got = interpret_rules("run it from 25 to 200 °C with 20 steps, both", Scenario())
    assert got.scenario_changes["t_start_c"] == 25.0 and got.scenario_changes["t_end_c"] == 200.0
    assert got.scenario_changes["steps"] == 20 and got.scenario_changes["mode"] == "both"
    got = interpret_rules("use 30 MPa", Scenario())
    assert got.scenario_changes["pressure_bar"] == 300.0


def test_rules_read_plot_requests():
    got = interpret_rules("plot CaCl+ instead of Ca++", Scenario())
    assert got.show == ["CaCl+"] and got.hide == ["Ca++"] and not got.run
    got = interpret_rules("show anhydrite and brucite on a linear axis", Scenario())
    assert got.show == ["anhydrite", "brucite"] and got.axis == "linear"
    assert interpret_rules("hide sulfate", Scenario()).hide == ["sulfate"]


def test_model_inventions_are_dropped():
    # the user said nothing about mode or pressure; the model made them up
    model = ScriptedModel({"scenario_changes": {"mode": "flow_through", "pressure_bar": 123, "t_end_c": 300},
                           "run": True})
    got = interpret("React seawater at 300 C", Scenario(), [], [], model)
    assert got.scenario_changes == {"t_end_c": 300.0, "fluid": "seawater"}
    assert got.source == "model" and not got.mode_mentioned


def test_model_failure_falls_back_to_rules():
    got = interpret("heat seawater to 250 C, static", Scenario(), [], [], ScriptedModel(broken=True))
    assert got.source == "rules" and got.scenario_changes["t_end_c"] == 250.0
    assert "did not answer" in got.reply


# ---- the shared view ----------------------------------------------------------------------------
def test_view_default_and_unit_groups():
    result = fake_result(Scenario(mode="both"))
    view = viewmod.default_view(result)
    assert [i["key"] for i in view["items"]] == ["Mg+2", "Ca+2", "SO4-2"]
    assert view["modes"] == ["static", "flow_through"]
    # adding something with another unit replaces the species curves
    view = viewmod.apply_patch(view, {"add": [{"group": "minerals", "key": "anhydrite"}]}, result)
    assert view["items"] == [{"group": "minerals", "key": "anhydrite"}]
    # toggling removes it again
    view = viewmod.apply_patch(view, {"toggle": {"group": "minerals", "key": "anhydrite"}}, result)
    assert view["items"] == []


def test_view_drops_items_the_result_cannot_plot():
    result = fake_result(Scenario())
    view = viewmod.sanitize({"items": [{"group": "species", "key": "Ca+2"}, {"group": "species", "key": "Xx+9"}],
                             "modes": ["static", "flow_through"], "unit": "furlongs"}, result)
    assert view["items"] == [{"group": "species", "key": "Ca+2"}]
    assert view["modes"] == ["static"] and view["unit"] == "mmol"


def test_resolve_item_finds_species_minerals_elements_properties():
    result = fake_result(Scenario())
    r = lambda text, **kw: viewmod.resolve_item(text, result, **kw)
    assert r("CaCl+") == {"group": "species", "key": "CaCl+"}
    assert r("Ca++") == {"group": "species", "key": "Ca+2"}
    assert r("anhydrite") == {"group": "minerals", "key": "anhydrite"}
    assert r("saturation index of anhydrite") == {"group": "saturation", "key": "anhydrite"}
    assert r("SI_anhydrite") == {"group": "saturation", "key": "anhydrite"}
    assert r("dolomite") == {"group": "saturation", "key": "dolomite"}      # never present, so its SI
    assert r("magnesium") == {"group": "elements", "key": "Mg"}
    assert r("pH") == {"group": "property", "key": "pH"}
    assert r("unobtainium") is None


# ---- narration ----------------------------------------------------------------------------------
def test_story_quotes_the_tables():
    result = fake_result(Scenario(mode="both"))
    text = " ".join(story(result))
    assert "anhydrite from about" in text
    assert "pH goes from 8.0" in text and "5.5" in text
    assert "dolomite" in text.lower() and "blocked" in text            # supersaturated but held out
    assert "Where the two modes differ" in text or "almost the same" in text


# ---- session flow -------------------------------------------------------------------------------
@pytest.fixture
def session(tmp_path):
    engine = FakeEngine()
    s = Session(tmp_path, engine, ScriptedModel({}, broken=True))
    s.engine_stub = engine
    return s


def chat_and_wait(s: Session, text: str):
    n = len(s.chat)
    s.handle_chat(text)
    assert len(s.chat) >= n + 2


def test_it_asks_for_the_mode_before_running_then_runs(session):
    chat_and_wait(session, "heat seawater to 300 C")
    assert "Static, flow-through, or both?" in session.chat[-1]["text"]
    assert session.engine_stub.calls == []
    chat_and_wait(session, "both")
    assert wait(lambda: session.engine_stub.calls and session.run_state["state"] == "idle")
    assert session.engine_stub.calls[0]["mode"] == "both"
    assert session.engine_stub.calls[0]["t_end_c"] == 300.0
    assert session.current and session.view["items"]


def test_chat_changes_the_plot_and_axis(session):
    session.edit_card({"mode": "static"})
    assert session.start_run()
    assert wait(lambda: session.run_state["state"] == "idle" and session.current)
    chat_and_wait(session, "plot CaCl+ instead of Ca++")
    keys = [i["key"] for i in session.view["items"]]
    assert "CaCl+" in keys and "Ca+2" not in keys
    chat_and_wait(session, "show anhydrite on a linear axis")
    assert session.view["items"] == [{"group": "minerals", "key": "anhydrite"}]
    assert session.view["log"] is False
    chat_and_wait(session, "plot unobtainium")
    assert "could not find" in session.chat[-1]["text"]


def test_a_bad_card_is_reported_not_run(session):
    session.edit_card({"mode": "static", "pressure_bar": 50})
    chat_and_wait(session, "run it")
    assert session.engine_stub.calls == []
    assert "boils" in session.chat[-1]["text"]


def test_second_run_request_during_a_run_is_queued(tmp_path):
    gate = threading.Event()

    class SlowEngine(FakeEngine):
        def run(self, card, progress=None):
            gate.wait(5)
            return super().run(card, progress)

    engine = SlowEngine()
    s = Session(tmp_path, engine, ScriptedModel(broken=True))
    s.edit_card({"mode": "static"})
    assert s.start_run()
    assert wait(lambda: s.run_state["state"] == "running")
    chat_and_wait(s, "run again at 200 C")
    assert "still going" in s.chat[-1]["text"]
    gate.set()
    assert wait(lambda: len(engine.calls) == 2 and s.run_state["state"] == "idle")
    assert engine.calls[1]["t_end_c"] == 200.0


def test_session_survives_a_restart(tmp_path):
    engine = FakeEngine()
    s = Session(tmp_path, engine, ScriptedModel(broken=True))
    s.edit_card({"t_end_c": 222, "mode": "static"})
    s.start_run()
    assert wait(lambda: s.run_state["state"] == "idle" and s.current)
    again = Session(tmp_path, FakeEngine(), ScriptedModel(broken=True))
    assert again.card.t_end_c == 222.0 and again.current == s.current
    assert again.view["items"]
    assert len(again.list_runs()) == 1


# ---- the web server -----------------------------------------------------------------------------
@pytest.fixture
def server(tmp_path):
    httpd = make_server(0, workspace=tmp_path, engine=FakeEngine(), model=ScriptedModel(broken=True))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield httpd
    httpd.shutdown()
    httpd.server_close()


def request(httpd, method, path, body=None, headers=None):
    port = httpd.server_address[1]
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    h = {"Host": f"127.0.0.1:{port}", **(headers or {})}
    if body is not None:
        h.setdefault("Content-Type", "application/json")
    conn.request(method, path, body=json.dumps(body) if body is not None else None, headers=h)
    r = conn.getresponse()
    data = r.read()
    conn.close()
    return r.status, data


def test_server_serves_pages_and_state(server):
    for path in ("/control", "/results", "/static/plot.js", "/static/app.css"):
        status, data = request(server, "GET", path)
        assert status == 200 and data, path
    status, data = request(server, "GET", "/api/state")
    state = json.loads(data)
    assert state["card"]["fluid"] == "seawater" and state["llm"]["ready"] is False


def test_server_edit_run_and_fetch_result(server):
    assert request(server, "POST", "/api/card", {"changes": {"mode": "both", "t_end_c": 240}})[0] == 200
    status, data = request(server, "POST", "/api/run", {})
    assert json.loads(data)["started"] is True
    session = server.session
    assert wait(lambda: session.run_state["state"] == "idle" and session.current)
    status, data = request(server, "GET", f"/api/run/{session.current}")
    result = json.loads(data)
    assert status == 200 and result["story"] and set(result["modes"]) == {"static", "flow_through"}
    assert request(server, "GET", "/api/run/../../etc/passwd")[0] == 404
    status, data = request(server, "POST", "/api/view", {"patch": {"toggle": {"group": "species", "key": "CaCl+"}}})
    assert status == 200
    assert {"group": "species", "key": "CaCl+"} in session.view["items"]


def test_server_rejects_foreign_hosts_and_wrong_content_types(server):
    port = server.server_address[1]
    assert request(server, "GET", "/api/state", headers={"Host": "evil.example"})[0] == 403
    assert request(server, "GET", "/api/state", headers={"Origin": "http://evil.example"})[0] == 403
    assert request(server, "POST", "/api/run", {}, headers={"Origin": "http://evil.example"})[0] == 403
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
    conn.request("POST", "/api/run", body="{}", headers={"Host": f"127.0.0.1:{port}", "Content-Type": "text/plain"})
    assert conn.getresponse().status == 415
    assert request(server, "GET", "/static/../server.py")[0] == 404


def test_server_never_shares_a_busy_port(server):
    port = server.server_address[1]
    other = make_server(port, workspace=server.session.workspace / "other", engine=FakeEngine(),
                        model=ScriptedModel(broken=True))
    try:
        assert other.server_address[1] != port           # fell through to a free port, did not share
    finally:
        other.server_close()


def test_events_stream_pushes_state(server):
    port = server.server_address[1]
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
    conn.request("GET", "/events", headers={"Host": f"127.0.0.1:{port}"})
    r = conn.getresponse()
    assert r.status == 200 and r.getheader("Content-Type") == "text/event-stream"
    first = r.fp.readline().decode()
    assert first.startswith("data: ") and json.loads(first[6:])["card"]["fluid"] == "seawater"
    request(server, "POST", "/api/card", {"changes": {"title": "pushed"}})
    r.fp.readline()                                       # blank line after the event
    second = r.fp.readline().decode()
    assert json.loads(second[6:])["card"]["title"] == "pushed"
    conn.close()


def test_missing_thermodynamic_data_is_reported_clearly(tmp_path):
    from reaction.app.engine import Engine, EngineUnavailable
    engine = Engine(data_dir=tmp_path / "nowhere")
    assert engine.status()["ready"] is False
    with pytest.raises(EngineUnavailable, match="Thermodynamic data not found"):
        engine.run(Scenario())
    s = Session(tmp_path / "ws", engine, ScriptedModel(broken=True))
    s.edit_card({"mode": "static"})
    assert s.start_run()
    assert wait(lambda: s.run_state["state"] == "error")
    assert "not found" in s.run_state["error"] and "failed" in s.chat[-1]["text"]

"""The scenario card, chat interpretation, shared view, session flow and local server.

These tests use a stub engine and a scripted language model, so they need no thermodynamic data and
no Ollama.  ``test_app_engine.py`` runs the real engine when the data is present.
"""
import http.client
import json
import re
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
    points = [(card.t_start_c + (card.t_end_c - card.t_start_c) * i / (n - 1),
               card.p_start_bar + (card.p_end_bar - card.p_start_bar) * i / (n - 1)) for i in range(n)]
    axes = {"temperature_c": {"key": "temperature_c", "label": "Temperature", "unit": "°C",
                              "values": [p[0] for p in points]},
            "pressure_bar": {"key": "pressure_bar", "label": "Pressure", "unit": "bar",
                             "values": [p[1] for p in points]}}
    varies = card.varies
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
            "title": card.display_title, "summary": card.describe(), "database": "stub",
            "fluid": {"key": "seawater", "label": "Seawater", "description": "", "measured_t_c": 25.0,
                      "source": "", "totals_molal": {}},
            "axes": axes, "varies": varies, "x": axes[varies[0] if varies else "temperature_c"],
            "reliability": {"levels": ["good"] * n, "deviation": [0.0] * n, "density": [0.9] * n},
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
    low_p = Scenario(p_start_bar=100.0, p_end_bar=100.0)          # water boils at 350 C below ~165 bar
    assert any("boils" in p for p in low_p.problems())
    assert Scenario(t_end_c=700).problems() and Scenario(p_end_bar=9000).problems()
    assert Scenario(mode="sideways").problems()
    assert Scenario(fluid="mud").problems()
    assert Scenario(p_start_bar=100.0, p_end_bar=100.0, t_end_c=250).problems() == []


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
    session.edit_card({"mode": "static", "pressure_bar": 50})          # sets both ends
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


def test_second_results_window_and_comparison(tmp_path):
    s = Session(tmp_path, FakeEngine(), ScriptedModel(broken=True))
    s.edit_card({"mode": "static"})
    s.start_run()
    assert wait(lambda: s.run_state["state"] == "idle" and s.current)
    first = s.current
    s.edit_card({"t_end_c": 200})
    s.start_run()
    assert wait(lambda: s.run_state["state"] == "idle" and s.current != first)
    second = s.current
    # window 2 follows window 1 until it is changed, then keeps its own run and plot
    assert s.window("2")["current"] == second
    assert s.set_current(first, "2")
    s.set_view({"log": False}, "2")
    assert s.current == second and s.window("2")["current"] == first
    assert s.view["log"] is True and s.window("2")["view"]["log"] is False
    # a comparison names another saved run; itself or a missing run is dropped
    s.set_view({"compare": first})
    assert s.view["compare"] == first
    s.set_view({"compare": second})
    assert s.view["compare"] is None
    s.set_view({"compare": "no-such-run"})
    assert s.view["compare"] is None
    s.set_view({"unit": "mg"})
    assert s.view["unit"] == "mg"
    again = Session(tmp_path, FakeEngine(), ScriptedModel(broken=True))
    assert again.window("2")["current"] == first and again.snapshot()["windows"]["2"]["current"] == first


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


def test_pages_carry_the_dougalien_systems_brand(server):
    for path in ("/control", "/results"):
        page = request(server, "GET", path)[1].decode()
        assert "Dougalien Systems" in page and "/static/dougalien_mark.svg" in page, path
    status, data = request(server, "GET", "/static/dougalien_mark.svg")
    assert status == 200 and data.startswith(b"<svg")


def test_sources_page_lists_the_data_with_links(server):
    for path in ("/control", "/results"):
        assert 'href="/sources"' in request(server, "GET", path)[1].decode(), path
    status, data = request(server, "GET", "/sources")
    page = data.decode()
    assert status == 200
    for name in ("CHNOSZ", "EQ3/6", "B-dot", "Pitzer", "Horita", "Driesner", "Allison", "Bradbury"):
        assert name in page, name
    links = re.findall(r'href="(https?://[^"]+)"', page)
    assert len(links) > 30 and all(u.startswith("https://") for u in links)


def test_clear_chat_empties_the_chat_but_keeps_the_card(server):
    request(server, "POST", "/api/card", {"changes": {"t_end_c": 250}})
    server.session.say("user", "heat seawater to 250 C")
    assert request(server, "POST", "/api/clear-chat", {})[0] == 200
    state = json.loads(request(server, "GET", "/api/state")[1])
    assert [m["role"] for m in state["chat"]] == ["assistant"]
    assert state["card"]["t_end_c"] == 250


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


# ---- pressure as a path variable ----------------------------------------------------------------
def test_pressure_can_vary_alone_or_together_with_temperature():
    assert Scenario().varies == ["temperature_c"]
    dec = Scenario(t_start_c=300, t_end_c=300, p_start_bar=1000, p_end_bar=100)
    assert dec.varies == ["pressure_bar"] and dec.problems() == []
    both = Scenario(t_start_c=400, t_end_c=150, p_start_bar=1000, p_end_bar=200)
    assert both.varies == ["temperature_c", "pressure_bar"] and both.problems() == []
    assert "decompressed from 1000 to 100 bar" in dec.describe()
    assert "cooled" in Scenario(t_start_c=300, t_end_c=100).describe()
    assert "taken from" in both.describe()
    pts = dec.points()
    assert pts[0] == (300.0, 1000.0) and pts[-1] == (300.0, 100.0) and len(pts) == dec.steps
    assert [round(p) for _, p in Scenario(t_start_c=0, t_end_c=100, p_start_bar=100, p_end_bar=500,
                                          steps=5).points()] == [100, 200, 300, 400, 500]


def test_a_path_where_nothing_changes_is_refused():
    assert any("nothing changes" in p for p in Scenario(t_start_c=100, t_end_c=100).problems())


def test_boiling_is_reported_where_it_happens_on_the_path():
    # the line from (100 C, 200 bar) to (300 C, 50 bar) crosses the boiling curve before it ends
    card = Scenario(t_start_c=100, t_end_c=300, p_start_bar=200, p_end_bar=50)
    problems = card.problems()
    assert len(problems) == 1 and "boils" in problems[0]
    import re
    t_report = float(re.search(r"at (\d+) C water boils", problems[0]).group(1))
    assert 250 < t_report < 300
    # the same pressures are fine if the temperature stays low
    assert Scenario(t_start_c=100, t_end_c=200, p_start_bar=200, p_end_bar=50).problems() == []


def test_old_cards_with_a_single_pressure_still_load():
    card = Scenario.from_dict({"title": "Seawater heating", "kind": "heating", "pressure_bar": 300.0,
                               "t_end_c": 250})
    assert (card.p_start_bar, card.p_end_bar) == (300.0, 300.0)
    assert card.kind == "path" and card.title == "" and card.display_title == "Seawater run"


def test_pressure_and_temperature_shorthands_set_both_ends():
    card, problems = Scenario().changed({"pressure_bar": 250, "temperature_c": 200})
    assert (card.p_start_bar, card.p_end_bar, card.t_start_c, card.t_end_c) == (250, 250, 200, 200)
    assert not problems and any("nothing changes" in p for p in card.problems())


def test_rules_read_pressure_paths():
    r = lambda text: interpret_rules(text, Scenario()).scenario_changes
    assert r("decompress from 1000 to 200 bar at 300 C") == {
        "p_start_bar": 1000.0, "p_end_bar": 200.0, "temperature_c": 300.0}
    assert r("compress to 3000 bar at 25 C") == {"p_end_bar": 3000.0, "temperature_c": 25.0}
    assert r("isothermal at 400 C from 2 kbar to 500 bar") == {
        "temperature_c": 400.0, "p_start_bar": 2000.0, "p_end_bar": 500.0}
    assert r("from 25 to 300 C while the pressure goes 500 to 300 bar") == {
        "t_start_c": 25.0, "t_end_c": 300.0, "p_start_bar": 500.0, "p_end_bar": 300.0}
    assert r("heat to 300 C at 400 bar") == {"t_end_c": 300.0, "pressure_bar": 400.0}
    assert interpret_rules("decompress to 200 bar", Scenario()).run


def test_rules_read_the_horizontal_axis():
    got = interpret_rules("plot Ca++ against pressure", Scenario())
    assert got.x_axis == "pressure_bar" and got.show == ["Ca++"]
    assert interpret_rules("show it versus temperature", Scenario()).x_axis == "temperature_c"


def test_the_guard_covers_pressure_and_temperature_inventions():
    model = ScriptedModel({"scenario_changes": {"p_start_bar": 777, "p_end_bar": 400, "t_end_c": 300}})
    got = interpret("heat seawater to 300 C at 400 bar", Scenario(), [], [], model)
    assert "p_start_bar" not in got.scenario_changes                 # 777 was never said
    assert got.scenario_changes["pressure_bar"] == 400.0 and got.scenario_changes["t_end_c"] == 300.0


def test_view_x_axis_follows_what_varies():
    both = fake_result(Scenario(t_start_c=300, t_end_c=100, p_start_bar=1000, p_end_bar=200))
    v = viewmod.default_view(both)
    assert v["x"] == "temperature_c"
    assert viewmod.apply_patch(v, {"x": "pressure_bar"}, both)["x"] == "pressure_bar"
    only_t = fake_result(Scenario())
    v = viewmod.apply_patch(viewmod.default_view(only_t), {"x": "pressure_bar"}, only_t)
    assert v["x"] == "temperature_c"                                  # pressure is constant here
    only_p = fake_result(Scenario(t_start_c=300, t_end_c=300, p_start_bar=1000, p_end_bar=100))
    assert viewmod.default_view(only_p)["x"] == "pressure_bar"


def test_story_describes_a_pressure_path_in_bar(session=None):
    result = fake_result(Scenario(t_start_c=300, t_end_c=300, p_start_bar=1000, p_end_bar=100))
    text = " ".join(story(result))
    assert "pH goes from 8.0 at 1000 bar to 5.5 at 100 bar" in text
    assert "anhydrite at about" in text and "bar" in text


def test_chat_runs_a_decompression_and_keeps_the_pressure_axis(session):
    chat_and_wait(session, "decompress from 1000 to 200 bar at 300 C, static")
    assert wait(lambda: session.engine_stub.calls and session.run_state["state"] == "idle")
    call = session.engine_stub.calls[0]
    assert (call["t_start_c"], call["t_end_c"], call["p_start_bar"], call["p_end_bar"]) == (
        300.0, 300.0, 1000.0, 200.0)
    assert session.view["x"] == "pressure_bar"
    chat_and_wait(session, "plot against temperature")
    assert "does not change" in session.chat[-1]["text"] and session.view["x"] == "pressure_bar"


def test_an_axis_request_made_with_a_new_run_applies_to_that_run(session):
    chat_and_wait(session, "decompress from 1000 to 200 bar at 300 C, static, plot against pressure")
    assert wait(lambda: session.engine_stub.calls and session.run_state["state"] == "idle")
    assert session.view["x"] == "pressure_bar"
    assert not any("does not change" in m["text"] for m in session.chat)


def test_an_axis_the_new_run_does_not_vary_is_explained_after_the_run(session):
    chat_and_wait(session, "heat to 300 C at 400 bar, static, plot against pressure")
    assert wait(lambda: session.engine_stub.calls and session.run_state["state"] == "idle")
    assert wait(lambda: any("does not change" in m["text"] for m in session.chat))
    assert session.view["x"] == "temperature_c"


def test_diagram_requests_in_the_chat_open_and_close_it(session):
    session.handle_chat("Eh-pH diagram for iron at 300 C")
    assert session.diagram["kind"] == "eh_ph" and session.diagram["element"] == "Fe"
    assert session.diagram["t_c"] == 300
    assert "Eh-pH diagram for Fe at 300 C" in session.chat[-1]["text"]
    assert session.snapshot()["diagram"] == session.diagram
    session.handle_chat("close the diagram")
    assert session.diagram is None


# ---- starting cards -----------------------------------------------------------------------------
def test_starting_cards_load_from_the_chat_and_the_reset_button(session):
    snap = session.snapshot()
    assert [c["key"] for c in snap["starting_cards"]] == ["seawater", "freshwater", "brine", "hydrothermal", "stream"]
    chat_and_wait(session, "load the brine card")
    assert session.card.fluid == "gsl_brine" and session.card.activity_model == "pitzer"
    assert "Brine card loaded" in session.chat[-1]["text"]
    chat_and_wait(session, "hydrothermal scenario")
    assert session.card.fluid == "hydrothermal" and session.card.mix_with == "seawater"
    session.reset_card("freshwater")
    assert session.card.fluid == "freshwater" and session.card.t_start_c == 25.0
    session.reset_card()
    assert session.card.to_dict() == Scenario().to_dict()


def test_reset_endpoint_takes_a_starting_card(server):
    assert request(server, "POST", "/api/reset", {"start": "hydrothermal"})[0] == 200
    assert server.session.card.fluid == "hydrothermal"
    assert request(server, "POST", "/api/reset", {"start": "nonsense"})[0] == 400
    assert request(server, "POST", "/api/reset", {})[0] == 200
    assert server.session.card.fluid == "seawater"


# ---- examples ------------------------------------------------------------------------------------
def test_an_example_loads_its_card_and_runs(session):
    from reaction.app.scenario import EXAMPLES, Scenario as Card
    snap = session.snapshot()
    assert [e["key"] for e in snap["examples"]] == list(EXAMPLES) and len(EXAMPLES) == 9
    known = {f for f in Card().to_dict()}
    for e in EXAMPLES.values():                      # every field the note outlines is on the card
        assert {field for field, _ in e.setup} <= known, e.key
        assert not Card(**e.changes).problems(), e.key
    assert session.load_example("serpentinization")
    assert session.card.solid_solutions and session.card.reactants[0]["name"] == "olivine fo90"
    assert session.snapshot()["example"] == "serpentinization"
    assert wait(lambda: session.run_state["state"] == "idle" and session.current)
    session.reset_card("brine")
    assert session.snapshot()["example"] is None


def test_foram_example_opens_on_its_mass_in_micrograms(session, monkeypatch):
    import dataclasses
    from reaction.app import scenario
    from reaction.app.scenario import reactant_text
    assert scenario.EXAMPLES["forams"].plot == ("calcite",) and scenario.EXAMPLES["forams"].unit == "ug"
    # the stub engine's result has anhydrite, not calcite: open on that instead
    monkeypatch.setitem(scenario.EXAMPLES, "forams", dataclasses.replace(scenario.EXAMPLES["forams"], plot=("anhydrite",)))
    assert session.load_example("forams")
    assert reactant_text(session.card.reactants[0]) == "20 µg forams"
    assert wait(lambda: session.run_state["state"] == "idle" and session.current)
    assert session.view["unit"] == "ug" and session.view["items"] == [{"group": "minerals", "key": "anhydrite"}]
    assert not session.view["log"]


def test_example_endpoint(server):
    assert request(server, "POST", "/api/example", {"key": "lake"})[0] == 200
    assert server.session.card.isotopes and server.session.card.evaporate == 10.0
    assert request(server, "POST", "/api/example", {"key": "nonsense"})[0] == 400

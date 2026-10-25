"""Predominance and activity diagrams (needs the thermodynamic data), and reading diagram requests."""
import math

import pytest

from reaction.app import config, diagrams
from reaction.app.engine import Engine
from reaction.app.scenario import Scenario
from reaction.diagram import Axis, DiagramError, predominance


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


@pytest.fixture(scope="module")
def db(engine):
    return engine.db


def test_iron_eh_ph_at_25c_has_the_textbook_fields(db):
    d = predominance(db, "Fe", Axis("pH", 0, 14), Axis("Eh", -1, 1.2), t_c=25)
    assert d.at(2, 1.0) == "Fe+3"
    assert d.at(4, 0.3) == "Fe+2"
    assert d.at(10, 0.3) == "hematite"
    assert d.at(10, -0.6) == "magnetite"
    assert d.at(4, -0.85) == "iron"
    # water: fO2 = 1 at 1.23 V and fH2 = 1 at 0 V at pH 0, both falling 59 mV per pH unit
    up, low = d.water["upper"], d.water["lower"]
    assert up[0][1] == pytest.approx(1.229, abs=0.005)
    assert low[0][1] == pytest.approx(0.0, abs=0.005)
    assert (up[0][1] - up[1][1]) / 14 == pytest.approx(0.0592, abs=0.0005)
    assert d.reliability.level == "good" and d.band == 0.15


def test_fe3_fe2_boundary_sits_at_the_standard_potential(db):
    # with equal activities the Fe+3/Fe+2 line is the standard potential, 0.77 V
    d = predominance(db, "Fe", Axis("pH", 0, 1), Axis("Eh", 0.6, 0.95), t_c=25, include="aqueous",
                     nx=3, ny=351)
    col = d.grid[:, 0]
    k = int((col == d.names.index("Fe+3")).argmax())
    assert d.ys[k] == pytest.approx(0.77, abs=0.01)


def test_diagram_holds_at_high_temperature_and_pressure(db):
    d = predominance(db, "Fe", Axis("pH", 2, 10), Axis("log_fO2", -50, -20), t_c=350, p_bar=500,
                     others={"S": -3})
    assert {"pyrite", "magnetite"} <= set(d.present())
    assert d.p_bar == 500


def test_k_feldspar_activity_diagram_and_quartz_limit(db):
    d = predominance(db, "Al", Axis("log_a", -6, -1, "SiO2"), Axis("log_ratio", 0, 12, "K+"), t_c=25,
                     include="minerals")
    assert {"kaolinite", "muscovite", "K-feldspar"} <= set(d.present())
    assert d.at(-4.5, 1) == "kaolinite"
    assert d.at(-3, 8) == "K-feldspar"
    (lim,) = d.limits
    assert lim["mineral"] == "quartz" and lim["value"] == pytest.approx(-4.0, abs=0.1)


def test_uncertain_band_follows_boundaries(db):
    d = predominance(db, "Fe", Axis("pH", 0, 14), Axis("Eh", -1, 1.2), t_c=25)
    assert d.margin.min() < d.band                       # cells near a boundary
    assert (d.margin > d.band).mean() > 0.8              # most of the plane is clear


def test_bad_requests_are_refused(db):
    with pytest.raises(DiagramError):
        predominance(db, "Fe", Axis("pH", 0, 14), Axis("pH", 0, 14))
    with pytest.raises(DiagramError):
        predominance(db, "Fe", Axis("Eh", 0, 1), Axis("log_fO2", -60, 0))


def test_requests_in_words():
    r = diagrams.parse_request("Eh-pH diagram for Fe at 300 C and 500 bar")
    assert r == {"kind": "eh_ph", "element": "Fe", "t_c": 300.0, "p_bar": 500.0}
    r = diagrams.parse_request("draw an fO2-pH diagram of copper with sulfur")
    assert r["kind"] == "fo2_ph" and r["element"] == "Cu" and r["with"] == ["S"]
    r = diagrams.parse_request("activity diagram for Na at 200 C")
    assert r["kind"] == "activity" and r["cation"] == "Na+" and r["t_c"] == 200
    assert diagrams.parse_request("pourbaix diagram of iron at 25 C")["element"] == "Fe"
    assert diagrams.parse_request("close the diagram") == {"close": True}
    assert diagrams.parse_request("heat seawater to 300 C") is None


def test_merge_checks_values():
    spec = diagrams.merge(None, {"kind": "eh_ph", "element": "Fe"})
    spec = diagrams.merge(spec, {"t_c": "abc", "element": "Xx", "with": "S, C", "kind": "bogus"})
    assert spec["kind"] == "eh_ph" and spec["t_c"] is None and spec["element"] == "Fe"
    assert spec["with"] == ["S", "C"]
    act = diagrams.merge(spec, {"kind": "activity"})
    assert act["element"] == "Al" and act["include"] == "minerals"


def test_diagram_follows_the_run(engine):
    card = Scenario().changed({"fluid": "seawater", "t_end_c": 350, "p_start_bar": 500, "p_end_bar": 500,
                               "steps": 6, "mode": "static",
                               "reactants": [{"name": "basalt glass", "amount": 300, "unit": "g"}]})[0]
    result = engine.run(card)
    out = diagrams.compute(engine.db, diagrams.merge(None, {"kind": "activity", "cation": "K+"}), result)
    assert out["t_c"] == 350 and out["p_bar"] == pytest.approx(500)
    pts = out["path"]["static"]
    assert len(pts) == 6 and all(math.isfinite(v) for p in pts for v in p)
    assert any("Held at the run's values" in n for n in out["notes"])
    fe = diagrams.compute(engine.db, diagrams.merge(None, {"kind": "fo2_ph", "element": "Fe", "with": ["S"]}),
                          result)
    assert {"Pyrite", "Pyrrhotite"} & {f["label"] for f in fe["fields"]}
    assert len(fe["grid"]) == fe["ny"] and len(fe["uncertain"][0]) == fe["nx"]

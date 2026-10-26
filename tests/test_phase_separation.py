"""Brine and vapour of a saline fluid from the NaCl-H2O correlations (Driesner and Heinrich 2007)."""
import math

import pytest

from reaction.app.engine import _phase_split, phase_notes
from reaction.phasesep import grid, phase_state, salinity_from_chloride

SEAWATER = salinity_from_chloride(0.5459)          # about 3.1 wt% NaCl


def test_table_shape_and_halite_saturation():
    t, p, cols = grid()
    assert len(t) == 200 and len(p) == 25
    # halite-saturated brine at 25 C: 26.4 wt% NaCl
    assert phase_state(25.0, 10.0, 0.30).x_liq == pytest.approx(0.264, abs=0.002)


def test_seawater_boils_near_400_c_at_250_bar():
    assert phase_state(350.0, 250.0, SEAWATER).kind == "one phase"
    st = phase_state(401.0, 250.0, SEAWATER)
    assert st.separated
    assert st.x_liq == pytest.approx(0.1936, abs=1e-3) and st.x_vap == pytest.approx(0.00099, rel=0.01)
    # lever rule
    assert st.vapour_mass == pytest.approx((st.x_liq - SEAWATER) / (st.x_liq - st.x_vap))
    assert phase_state(480.0, 250.0, SEAWATER).kind == "vapour + halite"
    assert phase_state(450.0, 450.0, SEAWATER).kind == "one phase"         # above the critical curve
    assert phase_state(400.0, 600.0, SEAWATER).kind == "one phase"
    assert phase_state(400.0, 6000.0, SEAWATER).kind == "outside"


# Bischoff and Rosenbauer (1988), table 1: the two-phase boundary of seawater (3.2 wt% NaCl), C and bar
BISCHOFF_ROSENBAUER = [(195.6, 14.1), (251.9, 40.7), (299.6, 83.4), (309.5, 95.8), (319.6, 109.4), (350.1, 160.8),
                       (370.2, 203.4), (390.2, 253.1), (400.2, 279.7), (405.2, 292.8), (409.8, 304.8),
                       (415.2, 319.3), (429.4, 357.2), (449.8, 409.7)]


def test_seawater_boundary_against_bischoff_and_rosenbauer():
    for t, p in BISCHOFF_ROSENBAUER:
        assert phase_state(t, p - 3.5, 0.032).separated, (t, p)
        assert not phase_state(t, p + 3.5, 0.032).separated, (t, p)
    # their critical point of seawater, 407 C and 298.5 bar: the correlations put it at 408 C, 301 bar
    from reaction import nacl_h2o as dh
    from scipy.optimize import brentq
    tc = brentq(lambda t: dh.x_crit(t) - dh.mole_fraction(0.032), 380.0, 450.0)
    assert tc == pytest.approx(407.0, abs=1.5) and dh.p_crit(tc) == pytest.approx(298.5, abs=3.0)


def test_split_along_a_path_fills_unsolved_points():
    points = [(350.0, 250.0), (390.0, 250.0), (400.0, 250.0), (420.0, 250.0), (480.0, 250.0)]
    run = {"phases": {"salinity": [3.1, None, None, 3.1, 3.1]}}
    _phase_split(run, points)
    assert run["phase_state"] == ["one phase", "liquid + vapour", "liquid + vapour", "liquid + vapour",
                                  "vapour + halite"]
    share = run["phases"]["vapour_share"]
    assert share[0] is None and 0 < share[1] < share[3] < 100
    notes = " ".join(phase_notes(run, points))
    assert "brine and vapour from about 390 °C" in notes and "halite" in notes


def test_salinity_from_chloride():
    assert salinity_from_chloride(1.0) == pytest.approx(58.443 / 1058.443)
    assert math.isclose(salinity_from_chloride(0.0), 0.0)


def _follow(phase):
    from reaction.app.engine import Engine
    from reaction.app.scenario import Scenario
    # seawater is one phase at 400 C and 305 bar and splits at 410 C and 300 bar (Bischoff and Rosenbauer: from
    # 408 C at 300 bar)
    sc, problems = Scenario().changed({"t_start_c": 400, "t_end_c": 410, "p_start_bar": 305, "p_end_bar": 300,
                                       "steps": 2, "follow_phase": phase})
    assert not problems
    return Engine().run(sc)


def test_following_the_brine_keeps_the_brine():
    result = _follow("brine")
    run = result["modes"]["static"]
    ph = run["phases"]
    assert run["phase_state"] == ["one phase", "liquid + vapour"]
    # what is left is the brine of the lever rule: its salinity, and the vapour's mass gone
    assert ph["salinity"][1] == pytest.approx(ph["brine_salinity"][1], rel=0.01)
    water = run["properties"]["water_kg"][1]
    assert water / (1 - ph["salinity"][1] / 100) == pytest.approx(1 - ph["vapour_share"][1] / 100, rel=0.05)
    assert any("chemistry shown is that of the brine" in n for n in result["notes"])


def test_following_the_vapour_leaves_salt_and_minerals_behind():
    run = _follow("vapour")["modes"]["static"]
    ph = run["phases"]
    assert ph["salinity"][1] == pytest.approx(ph["vapour_salinity"][1], rel=0.05)
    assert ph["split_off_kg"][1] == pytest.approx(1 - ph["vapour_share"][1] / 100, rel=0.05)
    # the anhydrite formed before the split stays with the brine: only what the vapour itself makes is left
    assert run["minerals"]["anhydrite"][1] < 0.2 * run["minerals"]["anhydrite"][0]


def test_follow_phase_on_the_card_and_in_chat():
    from reaction.app.assistant import interpret_rules
    from reaction.app.scenario import Scenario
    sc, problems = Scenario().changed({"follow_phase": "vapor"})
    assert not problems and sc.follow_phase == "vapour"
    assert "following the vapour" in sc.describe()
    got = interpret_rules("heat seawater from 380 to 440 C at 250 bar and follow the brine", Scenario())
    assert got.scenario_changes["follow_phase"] == "brine" and got.scenario_changes["fluid"] == "seawater"
    assert interpret_rules("let the vapor boil off", Scenario()).scenario_changes == {"follow_phase": "brine"}
    assert interpret_rules("leave the brine behind", Scenario()).scenario_changes == {"follow_phase": "vapour"}
    assert interpret_rules("follow the whole fluid", Scenario()).scenario_changes == {"follow_phase": "whole"}
    assert "follow_phase" not in interpret_rules("plot brine salinity", Scenario()).scenario_changes


def test_correlations_reproduce_the_table():
    # the lookup table was computed from the same correlations (xThermal): the brine agrees to 0.5%, the vapour to 2%
    from reaction import nacl_h2o as dh
    t, p, cols = grid()
    for tc, pb in ((401.0, 250.0), (451.0, 390.0), (601.0, 470.0), (351.0, 150.0)):
        i, j = list(p).index(pb), list(t).index(tc)
        assert dh.mass_fraction(dh.x_liquid(tc, pb)) == pytest.approx(cols["x_liq"][i, j], rel=0.005)
        assert dh.mass_fraction(dh.x_vapour(tc, pb)) == pytest.approx(cols["x_vap"][i, j], rel=0.02)

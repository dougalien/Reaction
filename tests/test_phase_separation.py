"""Brine and vapour of a saline fluid from the NaCl-H2O table (Driesner and Heinrich 2007)."""
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
    assert phase_state(400.0, 600.0, SEAWATER).kind == "outside"


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

"""One-dimensional reactive transport (transport.py) against PHREEQC 3.7.3's TRANSPORT, and the column run in the app."""
import json

import pytest

from reaction.app.engine import Engine
from reaction.app.narrate import story
from reaction.app.scenario import Scenario
from reaction.transport import Column, mixing_steps

# PHREEQC 3.7.3 (llnl.dat, phreeqc.dat's exchange constants), example 11 on 20 cells: 0.6 mmol CaCl2 displaces
# 1 mmol Na, 0.2 mmol K and 1.2 mmol Br (the tracer) from a column with 1.1 mmol of exchange sites per cell;
# cells 0.004 m, dispersivity 0.002 m, no diffusion, constant-concentration inlet, flux outlet.
# Outlet after each 4th shift: Na, K, Ca, Br, Cl in mmol/kg.
PHREEQC_EX11 = {
    4: [1.0, 0.2, 0.0, 1.2, 0.0], 8: [1.0, 0.2, 0.0, 1.19999, 1e-05], 12: [1.0, 0.2, 0.0, 1.17744, 0.02256],
    16: [0.99998, 0.20002, 0.0, 0.92643, 0.27357], 20: [0.99836, 0.20164, 0.0, 0.47693, 0.72307],
    24: [0.96657, 0.23343, 0.0, 0.16786, 1.03214], 28: [0.71226, 0.48774, 0.0, 0.04488, 1.15512],
    32: [0.2832, 0.9168, 0.0, 0.0099, 1.1901], 36: [0.06433, 0.8355, 0.15008, 0.0019, 1.1981],
    40: [0.00821, 0.10414, 0.54383, 0.00033, 1.19967], 44: [0.00113, 0.01631, 0.59128, 5e-05, 1.19995],
    48: [0.00016, 0.00272, 0.59856, 1e-05, 1.19999], 52: [2e-05, 0.00047, 0.59976, 0.0, 1.2],
    56: [0.0, 8e-05, 0.59996, 0.0, 1.2], 60: [0.0, 1e-05, 0.59999, 0.0, 1.2]}


@pytest.fixture(scope="module")
def engine():
    return Engine()


def test_mixing_steps():
    assert mixing_steps(0.0) == (0, 0.0)
    assert mixing_steps(0.3) == (1, 0.3)
    assert mixing_steps(0.5) == (2, 0.25)
    k, f = mixing_steps(1.0)
    assert k == 3 and f == pytest.approx(1 / 3)


def test_ion_exchange_column_against_phreeqc(engine):
    pore = Scenario(fluid="custom", analysis="Na 1\nK 0.2\nBr 1.2\nCa 1e-9\nCl 1e-9", units="mmol/kg", ph=7.0,
                    redox="none", charge_balance="none", t_start_c=25, t_end_c=25, p_start_bar=1, p_end_bar=1,
                    exchange_meq=1.1)
    inflow = Scenario(fluid="custom", analysis="Ca 0.6\nCl 1.2\nNa 1e-9\nK 1e-9\nBr 1e-9", units="mmol/kg", ph=7.0,
                      redox="none", charge_balance="none")
    sp0, _, _, _, system, bulk0 = engine._prepare(pore, lambda f, m: None)
    entering, _ = engine.starting_fluid(inflow)
    column = Column(system, [(25.0, 1.0)] * 20, bulk0, dict(entering.component_totals), dispersion=0.002 / 0.004,
                    guess=sp0)
    res = column.run(60)
    assert not res.failures and res.pore_volumes[-1] == 3.0
    for shift, expected in PHREEQC_EX11.items():
        c = res.outlet[shift].speciation.component_totals
        got = [1e3 * c[k] for k in ("Na+", "K+", "Ca+2", "Br-", "Cl-")]
        assert got == pytest.approx(expected, abs=0.006), shift
    # the tracer is conserved: what is left in the column and what has left it add up to what was there
    left = sum(1e3 * st.moles_in_solution("Br-") for st in res.profile)
    gone = 1e3 * (res.out["Br-"] + res.back["Br-"]) - 1e3 * 60 * entering.component_totals["Br-"]
    assert left + gone == pytest.approx(sum(1e3 * b["Br-"] for b in [bulk0] * 20), rel=1e-6)
    assert 0 < res.back["Br-"] < 0.05 * res.out["Br-"]          # a little leaves back through the inlet


def test_seawater_basalt_column(engine):
    sc, problems = Scenario().changed({"kind": "transport", "cells": 6, "pore_volumes": 1, "t_start_c": 25,
                                       "t_end_c": 250, "pressure_bar": 400, "reactants": "100 g basalt glass",
                                       "dispersivity_m": 0.05})
    assert not problems and not sc.problems()
    r = json.loads(json.dumps(engine.run(sc)))
    assert list(r["modes"]) == ["outlet", "start", "profile"]
    out, end = r["modes"]["outlet"], r["modes"]["profile"]
    assert out["varies"] == ["pore_volumes"] and len(out["properties"]["pH"]) == 7
    assert end["varies"] == ["distance", "temperature_c"] and len(end["properties"]["pH"]) == 6
    assert end["axes"]["temperature_c"]["values"] == [25, 70, 115, 160, 205, 250]
    assert len(end["reliability"]["levels"]) == 6
    # heated seawater leaves anhydrite in the hot part of the column, not at the cold inlet
    anh = end["minerals"]["anhydrite"]
    assert (anh[0] or 0) < 1e-6 and anh[-1] > 1e-3
    assert any("TRANSPORT" in n for n in r["notes"])
    text = " ".join(story(r))
    assert "At the outlet the pH goes from" in text and "anhydrite (" in text
    assert r["title"] == "Seawater column"


def test_card_and_chat():
    from reaction.app.assistant import interpret_rules
    sc, problems = Scenario().changed({"kind": "column", "pore_volumes": 3, "cells": 30, "column_m": 2,
                                       "dispersivity_m": 0.1})
    assert not problems and sc.transport and sc.shifts == 90 and sc.mixing_fraction == pytest.approx(1.5)
    assert sc.cell_points()[0] == (sc.t_start_c, 500.0) and sc.cell_points()[-1] == (sc.t_end_c, 500.0)
    assert "flowing through a 2 m column of 30 cells" in sc.describe()
    sc2, _ = sc.changed({"kinetics": True, "reactants": "1 g calcite"})
    assert any("not available in a column" in p for p in sc2.problems())
    sc3, _ = sc.changed({"cells": 100, "pore_volumes": 100})
    assert any("equilibrium calculations" in p for p in sc3.problems())
    got = interpret_rules("run rainwater through a 50 cm column filled with seawater, 10 cells, 2.5 pore volumes, "
                          "dispersivity 1 cm", Scenario()).scenario_changes
    assert got == {"fluid": "rainwater", "kind": "transport", "cells": 10, "pore_volumes": 2.5, "column_m": 0.5,
                   "dispersivity_m": 0.01, "pore_fluid": "seawater"}
    assert interpret_rules("no column", sc).scenario_changes == {"kind": "path"}
    assert "mode" not in interpret_rules("flow through a column of 10 g calcite", Scenario()).scenario_changes
    got = interpret_rules("show anhydrite and albite along the column", sc)
    assert got.scenario_changes == {} and got.show == ["anhydrite", "albite"] and got.column_view == "profile"
    assert interpret_rules("plot Ca at the outlet", sc).column_view == "outlet"


def test_views_keep_one_set_of_axes():
    from reaction.app.view import default_view, sanitize
    run = {"species": {"Ca+2": [1.0]}, "properties": {}}
    result = {"catalog": {}, "x": {"key": "pore_volumes"}, "varies": ["pore_volumes"],
              "modes": {"outlet": {**run, "varies": ["pore_volumes"]},
                        "start": {**run, "varies": ["distance", "temperature_c"]},
                        "profile": {**run, "varies": ["distance", "temperature_c"]}}}
    assert default_view(result)["modes"] == ["outlet"]
    view = sanitize({"modes": ["profile", "outlet", "start"], "x": "pore_volumes", "items": []}, result)
    assert view["modes"] == ["profile", "start"] and view["x"] == "distance"
    view = sanitize({"modes": ["profile"], "x": "temperature_c", "items": []}, result)
    assert view["x"] == "temperature_c"

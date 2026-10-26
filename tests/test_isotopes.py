"""Stable isotopes and 87Sr/86Sr as a layer on the reaction path (isotopes.py)."""
import math

import pytest

from reaction import isotopes
from reaction.app.engine import Engine
from reaction.app.narrate import story
from reaction.app.scenario import Scenario


@pytest.fixture(scope="module")
def engine():
    return Engine()


def _run(engine, **card):
    sc, problems = Scenario().changed({"isotopes": True, **card})
    assert not problems
    return engine.run(sc)


def test_fractionation_factors():
    f = isotopes.table()
    # Vho et al. (2020): quartz-water 4.11, -3.33; quartz-calcite 0.39, 0.15
    assert f[("O", "quartz")].delta(25.0) == pytest.approx(35.07, abs=0.02)
    assert f[("O", "calcite")].delta(25.0) == pytest.approx(30.18, abs=0.02)
    # Ohmoto and Rye (1979): calcite-CO2 about 10.1 permil at 25 C; Ohmoto and Lasaga (1982) sulfate-H2S
    assert f[("C", "calcite")].delta(25.0) == pytest.approx(10.1, abs=0.1)
    assert f[("S", "sulfate")].delta(300.0) == pytest.approx(6.463e6 / 573.15 ** 2 + 0.56)
    # liquid-vapour 18O, Horita and Wesolowski (1994): 9.3 permil at 25 C, none above the critical point
    assert isotopes.liquid_vapour_o(25.0) == pytest.approx(9.30, abs=0.02)
    assert isotopes.liquid_vapour_o(400.0) == 0.0


def test_seawater_basalt_at_350_c(engine):
    r = _run(engine, temperature_c=350, pressure_bar=500, steps=6,
             reactants=[{"name": "basalt glass", "amount": 100, "unit": "g"}])
    run = r["modes"]["static"]
    iso = run["isotopes"]
    # the fluid's 87Sr/86Sr is the Sr-weighted mix of seawater (7.8 mg/kg) and 100 g of basalt (115 ppm)
    sw, rock = 7.8e-3 / 87.62, 115e-6 * 100 / 87.62
    assert iso["87Sr/86Sr:fluid"][-1] == pytest.approx((sw * 0.70918 + rock * 0.7027) / (sw + rock), abs=1e-6)
    # each mineral sits at its fractionation from the water
    d = isotopes.table()[("O", "quartz")].delta(350.0)
    assert iso["d18O:quartz"][-1] - iso["d18O:water"][-1] == pytest.approx(d, abs=0.002)
    gap = isotopes.table()[("S", "sulfate")].delta(350.0)
    assert iso["d34S:sulfate"][-1] - iso["d34S:sulfide"][-1] == pytest.approx(gap, abs=0.002)
    # basalt (5.7 permil) raises the water's d18O; sulfate reduction leaves the sulfate heavier
    assert 0 < iso["d18O:water"][-1] < 1.0
    assert iso["d34S:sulfate"][-1] > 21.5
    # sulfur mass balance: the basalt glass brings no sulfur, so the system keeps 21 permil overall
    db = engine.db
    kg = run["properties"]["water_kg"][-1]
    total = inventory = 0.0
    for name, line in run["species"].items():
        s = db.species[name].elements.get("S", 0.0) if name in db.species else 0.0
        if s and line[-1]:
            o = db.species[name].elements.get("O", 0.0)
            key = "d34S:sulfate" if o >= 4 * s else "d34S:sulfide" if o == 0 else None
            if key:
                n = line[-1] * kg * s
                total += n
                inventory += n * iso[key][-1]
    for name, line in run["minerals"].items():
        s = db.species[name].elements.get("S", 0.0)
        if s and line[-1]:
            total += line[-1] * s
            inventory += line[-1] * s * iso[f"d34S:{name}"][-1]
    assert inventory / total == pytest.approx(21.0, abs=0.05)
    assert any("Vho" in n for n in r["notes"])


def test_evaporation_follows_rayleigh(engine):
    r = _run(engine, temperature_c=25, pressure_bar=1, evaporate=10, steps=21)
    water = r["modes"]["static"]["isotopes"]["d18O:water"]
    # Rayleigh: d = -1000 ln(f) * 1000 ln(alpha l-v) / 1000 to first order, f = 0.1
    expected = 1000.0 * (0.1 ** (-(isotopes.liquid_vapour_o(25.0) / 1000.0)) - 1.0)
    assert water[-1] == pytest.approx(expected, rel=0.03)


def test_card_and_chat():
    from reaction.app.assistant import interpret_rules
    sc, problems = Scenario().changed({"isotope_inputs": {"Fluid": {"δ18O": -5, "87sr/86sr": 0.7035}}})
    assert not problems and sc.isotope_inputs == {"fluid": {"d18O": -5.0, "87Sr/86Sr": 0.7035}}
    _, problems = Scenario().changed({"isotope_inputs": {"fluid": {"87Sr/86Sr": 7.035}}})
    assert problems
    card, _ = Scenario().changed({"reactants": [{"name": "basalt glass", "amount": 100, "unit": "g"}]})
    got = interpret_rules("set the fluid d18O to -5 and basalt 87Sr/86Sr 0.7030", card).scenario_changes
    assert got == {"isotopes": True, "isotope_inputs": {"fluid": {"d18O": -5.0}, "basalt glass": {"87Sr/86Sr": 0.703}}}
    assert interpret_rules("no isotopes", card).scenario_changes == {"isotopes": False}
    assert interpret_rules("plot d18O quartz", card).scenario_changes == {}
    assert interpret_rules("Sr+2 0.08 mmol", card).scenario_changes.get("isotopes") is None


def test_given_starting_values_are_used(engine):
    r = _run(engine, t_start_c=25, t_end_c=100, steps=3, isotope_inputs={"fluid": {"d18O": -5.0}})
    water = r["modes"]["static"]["isotopes"]["d18O:water"]
    assert water[0] == pytest.approx(-5.0, abs=0.01)
    assert all(math.isfinite(v) for v in water if v is not None)
    assert not any("the fluid: d18O" in n for n in r["notes"])


def test_hydrogen_isotopes(engine):
    # Horita and Wesolowski (1994) liquid-vapour D/H: about 76 permil at 25 C, 28 at 100 C, slightly inverse
    # (the vapour heavier) near 300 C, none above the critical point
    assert isotopes.liquid_vapour_h(25.0) == pytest.approx(75.8, abs=0.2)
    assert isotopes.liquid_vapour_h(100.0) == pytest.approx(27.9, abs=0.2)
    assert -5.0 < isotopes.liquid_vapour_h(300.0) < 0.0 and isotopes.liquid_vapour_h(400.0) == 0.0
    # evaporating to a tenth of the water at 25 C follows Rayleigh for dD as for d18O
    r = _run(engine, temperature_c=25, pressure_bar=1, evaporate=10, steps=21,
             isotope_inputs={"fluid": {"dD": -54.0, "d18O": -8.0}})
    iso = r["modes"]["static"]["isotopes"]
    f = 0.1
    # residual water: R / R0 = f ** (1 / alpha - 1), alpha = R(liquid) / R(vapour)
    expected = (1000.0 - 54.0) * f ** (math.exp(-isotopes.liquid_vapour_h(25.0) / 1000.0) - 1.0) - 1000.0
    assert iso["dD:water"][-1] == pytest.approx(expected, rel=0.015)
    assert not any(k.startswith("dD:") and k != "dD:water" for k in iso)
    assert r["catalog"]["isotopes"]["dD:water"]["label"] == "δD water"
    # basalt alteration barely moves the water's dD: the minerals hold about 1% of the hydrogen
    r = _run(engine, temperature_c=350, pressure_bar=500, steps=3,
             reactants=[{"name": "basalt glass", "amount": 100, "unit": "g"}])
    assert abs(r["modes"]["static"]["isotopes"]["dD:water"][-1]) < 1.0


def test_mineral_water_hydrogen():
    h = isotopes.mineral_water_h
    # Sheppard and Gilg (1996) kaolinite: -2.2e6/T^2 - 7.7, about -32 permil at 25 C
    assert h("kaolinite", {}, 25.0)[0] == pytest.approx(-32.45, abs=0.05)
    # Suzuoki and Epstein (1976) eq. 1 for muscovite and the general equation (which their table 3 checks at 650 C)
    t = 650.0 + 273.15
    assert h("muscovite", {}, 650.0)[0] == pytest.approx(-22.1e6 / t ** 2 + 19.1)
    phlogopite = {"K": 1, "Mg": 3, "Al": 1, "Si": 3, "O": 12, "H": 2}
    assert h("phlogopite", phlogopite, 650.0)[0] == pytest.approx(-22.4e6 / t ** 2 + 28.2 + 2 * 0.25 - 4 * 0.75)
    annite = {"K": 1, "Fe": 3, "Al": 1, "Si": 3, "O": 12, "H": 2}
    assert h("annite", annite, 650.0)[0] < h("phlogopite", phlogopite, 650.0)[0] - 45   # Fe takes up less D
    # below its 400 C range the general equation is held at its 400 C value
    assert h("phlogopite", phlogopite, 200.0)[0] == h("phlogopite", phlogopite, 400.0)[0]
    # Graham and others (1984): tremolite flat at -21.7 to 650 C, hornblende -23.1, actinolite -29
    assert h("tremolite", {}, 500.0)[0] == -21.7 and h("pargasite", {}, 500.0)[0] == -23.1
    assert h("tremolite", {}, 800.0)[0] == pytest.approx(-31.0e6 / (800 + 273.15) ** 2 + 14.9)
    # molecular water and oxyanion hydrogen take the water's dD
    assert h("laumontite", {"Ca": 1, "Al": 2, "Si": 4, "O": 16, "H": 8}, 200.0)[0] is None
    assert h("gypsum", {"Ca": 1, "S": 1, "O": 6, "H": 4}, 25.0)[0] is None


def test_mineral_dd_in_a_run(engine):
    r = _run(engine, t_start_c=2, t_end_c=350, steps=8, reactants=[{"name": "basalt glass", "amount": 100, "unit": "g"}])
    iso = r["modes"]["static"]["isotopes"]
    minerals = [k for k in iso if k.startswith("dD:") and k != "dD:water"]
    assert minerals
    for k in minerals:                    # every hydrous silicate is lighter than the water it grew from
        for d, w in zip(iso[k], iso["dD:water"]):
            if d is not None:
                assert d < w - 15
    assert any("mineral-water D/H" in n for n in r["notes"])
    assert "lighter in δD than the water" in " ".join(story(r))

"""The heating run on the real engine.  Needs the thermodynamic data (see app/config.py)."""
import math

import pytest

from reaction.app import config
from reaction.app.engine import Engine
from reaction.app.narrate import story
from reaction.app.scenario import Scenario


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


@pytest.fixture(scope="module")
def result(engine):
    card = Scenario(t_start_c=25.0, t_end_c=250.0, steps=10, p_start_bar=500.0, p_end_bar=500.0, mode="both")
    return engine.run(card)


def test_both_modes_are_returned_on_the_same_grid(result):
    assert set(result["modes"]) == {"static", "flow_through"}
    assert len(result["x"]["values"]) == 10
    for run in result["modes"].values():
        assert all(v is not None for v in run["properties"]["pH"])
        assert len(run["species"]["Mg+2"]) == 10


def test_anhydrite_appears_on_heating_near_150c(result):
    temps = result["x"]["values"]
    amounts = result["modes"]["static"]["minerals"]["anhydrite"]
    first = next(i for i, v in enumerate(amounts) if v > 1e-6)
    assert 100 <= temps[first] <= 175                 # Bischoff and Seyfried (1978): about 150 C
    sat = result["modes"]["static"]["saturation"]["anhydrite"]
    assert sat[first - 1] < 0 and abs(sat[first]) < 1e-6   # undersaturated before, at saturation after


def test_ph_falls_on_heating_and_blocked_carbonates_stay_out(result):
    ph = result["modes"]["static"]["properties"]["pH"]
    assert ph[0] > 7.5 and ph[-1] < ph[0] - 2
    minerals = result["modes"]["static"]["minerals"]
    assert "dolomite" not in minerals and "magnesite" not in minerals
    assert max(result["modes"]["static"]["saturation"]["dolomite"]) > 0      # supersaturated yet blocked


def test_mass_is_conserved_in_the_static_run(result):
    """Sulfur dissolved plus sulfur in anhydrite equals the sulfur the seawater started with."""
    run = result["modes"]["static"]
    start = result["start_elements"]["S"]
    for i in range(10):
        dissolved = run["elements"]["S"][i] * run["properties"]["water_kg"][i]
        solid = run["minerals"].get("anhydrite", [0.0] * 10)[i] + run["minerals"].get("kieserite", [0.0] * 10)[i]
        assert dissolved + solid == pytest.approx(start, rel=2e-4)


def test_flow_through_keeps_what_static_redissolves(result):
    st, ft = result["modes"]["static"]["minerals"], result["modes"]["flow_through"]["minerals"]
    # calcite forms on warming and static heating redissolves it as the fluid turns acid;
    # in flow-through it was left behind and stays
    assert ft["calcite"][-1] > 0
    assert st.get("calcite", [0.0])[-1] < ft["calcite"][-1]


def test_story_and_catalog_are_consistent(result):
    text = " ".join(story(result))
    assert "anhydrite" in text
    for mode_run in result["modes"].values():
        for name in mode_run["species"]:
            assert name in result["catalog"]["species"]
        for table in ("minerals", "saturation"):
            for name in mode_run[table]:
                assert name in result["catalog"]["minerals"]
    assert math.isfinite(result["elapsed_s"])


def test_a_bad_card_is_refused(engine):
    with pytest.raises(ValueError, match="boils"):
        engine.run(Scenario(p_start_bar=100.0, p_end_bar=100.0))


def test_flow_through_runs_above_the_critical_temperature(engine):
    """Leaving minerals behind must not need a saturation pressure, which does not exist above 374 C."""
    card = Scenario(t_start_c=300.0, t_end_c=450.0, steps=4, p_start_bar=1000.0, p_end_bar=1000.0,
                    mode="flow_through")
    run = engine.run(card)["modes"]["flow_through"]
    assert all(v is not None for v in run["properties"]["pH"])
    assert run["minerals"]["anhydrite"][-1] > 0


def test_a_failed_step_does_not_change_the_flow_through_result(engine, monkeypatch):
    """A step that fails to converge is retried in smaller steps; minerals are left behind where they
    would start to redissolve, not at the steps, so the retry changes nothing."""
    from reaction import path as path_module
    from reaction.speciation import SpeciationError

    card = Scenario(t_start_c=25.0, t_end_c=350.0, steps=8, p_start_bar=500.0, p_end_bar=500.0,
                    mode="flow_through")
    plain = engine.run(card)["modes"]["flow_through"]
    target = card.points()[5][0]
    real, failed = path_module.equilibrate, []

    def flaky(system, t_c, *args, **kwargs):
        if abs(t_c - target) < 1e-9 and not failed:
            failed.append(t_c)
            raise SpeciationError("simulated convergence failure")
        return real(system, t_c, *args, **kwargs)

    monkeypatch.setattr(path_module, "equilibrate", flaky)
    retried = engine.run(card)["modes"]["flow_through"]
    assert failed
    assert set(retried["minerals"]) == set(plain["minerals"])
    for name, amounts in plain["minerals"].items():
        assert retried["minerals"][name] == pytest.approx(amounts, rel=1e-6, abs=1e-12)
    assert retried["properties"]["pH"] == pytest.approx(plain["properties"]["pH"], abs=1e-6)


def test_flow_through_leaves_calcite_behind_where_it_would_redissolve(engine):
    """Calcite forms on warming and redissolves in the static run; in flow-through the amount left
    behind stops growing where the static run starts to lose it, and never falls."""
    card = Scenario(t_start_c=25.0, t_end_c=250.0, steps=10, p_start_bar=500.0, p_end_bar=500.0,
                    mode="flow_through")
    ft = engine.run(card)["modes"]["flow_through"]["minerals"]["calcite"]
    assert all(b >= a * (1 - 1e-9) for a, b in zip(ft, ft[1:]))
    assert ft[-1] > 0


@pytest.mark.parametrize("key", ["rainwater", "vent_fluid", "formation_brine"])
def test_every_preset_fluid_runs(engine, key):
    card = Scenario(fluid=key, t_start_c=25.0, t_end_c=300.0, steps=6, p_start_bar=500.0,
                    p_end_bar=500.0, mode="static")
    res = engine.run(card)
    assert all(v is not None for v in res["modes"]["static"]["properties"]["pH"])
    assert res["fluid"]["key"] == key and res["fluid"]["analysis"]


def test_vent_fluid_precipitates_silica_on_heating(engine):
    res = engine.run(Scenario(fluid="vent_fluid", t_start_c=25.0, t_end_c=300.0, steps=6,
                              p_start_bar=500.0, p_end_bar=500.0, mode="static"))
    assert res["modes"]["static"]["minerals"]["quartz"][0] > 1e-3          # 17 mmol/kg silica
    assert any("Charge balance closed on Cl" in n for n in res["notes"])


def test_charge_balance_closes_on_the_chosen_ion(engine):
    card = Scenario(fluid="custom", analysis="Na 100 mmol/kg\nCl 80 mmol/kg\nCa 5 mmol/kg", ph=7.0)
    sp, notes = engine.starting_fluid(card)
    assert abs(sp.charge_imbalance) < 1e-10
    assert sp.totals["Cl-"] > 0.08 / (1 - 6.4e-3) and "Cl⁻" in notes[0]
    card.charge_balance = "Na+"
    sp, notes = engine.starting_fluid(card)
    assert abs(sp.charge_imbalance) < 1e-10 and sp.totals["Na+"] < 0.1
    card.charge_balance = "none"
    sp, notes = engine.starting_fluid(card)
    assert sp.charge_imbalance > 0.01 and "out of charge balance" in notes[0]


def test_reacting_albite_makes_clay_and_quartz_and_a_water_rock_axis(engine):
    card = Scenario(t_start_c=25.0, t_end_c=25.0, p_start_bar=1.0, p_end_bar=1.0, steps=6, mode="static",
                    reactants=[{"name": "Albite", "amount": 10.0, "unit": "g"}])
    res = engine.run(card)
    assert res["varies"][:2] == ["reacted", "water_rock"] and res["x"]["key"] == "reacted"
    assert res["axes"]["reacted"]["values"][-1] == pytest.approx(10.0)
    assert res["axes"]["water_rock"]["values"][0] is None and res["axes"]["water_rock"]["values"][-1] == 100.0
    run = res["modes"]["static"]
    assert all(v is not None for v in run["properties"]["pH"])
    assert run["minerals"]["quartz"][-1] > 0.01
    # all the sodium of 10 g of albite stays in solution
    na_added = 10.0 / 262.22
    na = [run["elements"]["Na"][i] * run["properties"]["water_kg"][i] for i in (0, -1)]
    assert na[1] - na[0] == pytest.approx(na_added, rel=2e-3)


def test_titration_axis_is_in_millimoles(engine):
    card = Scenario(t_start_c=25.0, t_end_c=25.0, p_start_bar=1.0, p_end_bar=1.0, steps=6, mode="static",
                    reactants=[{"name": "HCl", "amount": 0.02, "unit": "mol"}])
    res = engine.run(card)
    assert res["varies"] == ["reacted"] and res["axes"]["reacted"]["unit"] == "mmol per kg water"
    ph = res["modes"]["static"]["properties"]["pH"]
    assert ph[-1] < 2.5 < ph[0]


def test_unknown_or_redox_reactants_are_refused_clearly(engine):
    base = dict(t_start_c=25.0, t_end_c=25.0, p_start_bar=1.0, p_end_bar=1.0, steps=4)
    with pytest.raises(ValueError, match="not a mineral"):
        engine.run(Scenario(**base, reactants=[{"name": "unobtainium", "amount": 1.0, "unit": "g"}]))
    with pytest.raises(ValueError, match="oxidation state"):
        engine.run(Scenario(**base, reactants=[{"name": "magnetite", "amount": 1.0, "unit": "g"}]))


def test_solved_state_carries_the_water_activity_it_was_solved_with(engine):
    """The active set ranks minerals by saturation computed from the returned activities; with a stale
    water activity, two polymorphs of a hydrous-reaction mineral (K-feldspar, microcline) each looked
    supersaturated when the other was present, and the search cycled."""
    from reaction.constants import LN10
    from reaction.equilibrium import System, equilibrate
    sp, _ = engine.starting_fluid(Scenario())
    basis = [b for b in sp.basis_names if b not in ("H2O", "H+")]
    state = equilibrate(System(engine.db, basis), 25.0, dict(sp.component_totals), pressure_bar=1.0)
    assert state._u[0] == pytest.approx(state.speciation.log_a_water * LN10, abs=1e-12)
    assert state.speciation.log_a_water < -0.005


def test_mixing_vent_fluid_with_cold_seawater(engine):
    card = Scenario(fluid="vent_fluid", t_start_c=350.0, p_start_bar=250.0, p_end_bar=250.0, steps=11,
                    mode="both", mix_with="seawater", mix_ratio=10.0, mix_t_c=2.0)
    res = engine.run(card)
    assert res["varies"] == ["mixed", "mix_fraction", "temperature_c"]
    temps = res["axes"]["temperature_c"]["values"]
    assert temps[0] == 350.0 and temps[-1] == pytest.approx((350.0 + 10 * 2.0) / 11, rel=1e-5)
    assert res["axes"]["mix_fraction"]["values"][-1] == pytest.approx(10 / 11, rel=1e-5)
    st, ft = res["modes"]["static"], res["modes"]["flow_through"]
    # chloride is in no mineral: what is dissolved at the end is the vent fluid's plus ten seawaters'
    vent, _ = engine.starting_fluid(Scenario(fluid="vent_fluid"))
    sea, _ = engine.starting_fluid(Scenario(fluid="seawater"))
    cl_end = st["elements"]["Cl"][-1] * st["properties"]["water_kg"][-1]
    assert cl_end == pytest.approx(vent.totals["Cl-"] + 10 * sea.totals["Cl-"], rel=1e-4)
    # anhydrite forms where hot vent fluid meets seawater sulfate; cooling redissolves it unless it
    # was left behind
    assert max(v for v in st["minerals"]["anhydrite"] if v is not None) > 1e-3
    assert st["minerals"]["anhydrite"][-1] < 1e-9 < ft["minerals"]["anhydrite"][-1]


def test_mixing_without_its_own_temperature_follows_the_path(engine):
    card = Scenario(t_start_c=25.0, t_end_c=25.0, p_start_bar=1.0, p_end_bar=1.0, steps=5,
                    mix_with="rainwater", mix_ratio=4.0)
    res = engine.run(card)
    assert res["varies"] == ["mixed", "mix_fraction"]
    assert set(res["axes"]["temperature_c"]["values"]) == {25.0}
    na = res["modes"]["static"]["species"]["Na+"]
    assert na[-1] < na[0] / 4

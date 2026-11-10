"""Kinetic reactants: rate laws, the time-stepping and the card (most need the thermodynamic data)."""
import math

import pytest

from reaction import kinetics
from reaction.app import config
from reaction.app.assistant import interpret_rules
from reaction.app.engine import Engine
from reaction.app.narrate import story
from reaction.app.scenario import Scenario
from reaction.kinetics import Mechanism, RateLaw


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


def kinetic_card(**changes):
    base = {"fluid": "rainwater", "temperature_c": 25, "pressure_bar": 1, "steps": 6, "kinetics": True,
            "duration": 1, "duration_unit": "days"}
    card, problems = Scenario().changed({**base, **changes})
    assert not problems
    return card


class _Spec:                     # just what a rate law reads from a speciation
    def __init__(self, t_c, ph, species=None, fugacity=None):
        self.temperature_c, self.ph = t_c, ph
        self.species = species or {}
        self.log_fugacity = fugacity or {}


def test_palandri_kharaka_table():
    laws = kinetics.palandri_kharaka()
    assert len(laws) >= 60
    albite = laws["albite"]
    assert [m.name for m in albite.mechanisms] == ["acid", "neutral", "base"]
    # by hand at pH 5, 25 C: 10^-10.16 aH^0.317 + 10^-12.56 + 10^-15.6 aH^-0.471
    want = 10 ** -10.16 * 1e-5 ** 0.317 + 10 ** -12.56 + 10 ** -15.6 * 1e-5 ** -0.471
    assert albite.far_from_equilibrium(_Spec(25.0, 5.0)) == pytest.approx(want, rel=1e-9)
    assert all(m.e_kj > 0 for m in laws["kyanite"].mechanisms)
    assert laws["albite,low"] == albite and "quartz,beta" in laws


def test_arrhenius_and_affinity():
    m = Mechanism("neutral", -12.0, 60.0)
    assert m.k(298.15) == pytest.approx(1e-12)
    assert m.k(373.15) / m.k(298.15) == pytest.approx(math.exp(-60.0 / kinetics.R_KJ * (1 / 373.15 - 1 / 298.15)))
    law = RateLaw((m,))
    spec = _Spec(25.0, 7.0)
    assert law.rate(spec, 2.0, saturation=-10.0) == pytest.approx(2e-12)
    assert law.rate(spec, 2.0, saturation=0.0) == 0.0
    assert law.rate(spec, 2.0, saturation=math.log10(2.0)) == pytest.approx(-2e-12)   # grows
    assert law.rate(spec, 0.0, saturation=-10.0) == 0.0


def test_calcite_carbonate_mechanism_needs_co2():
    calcite = kinetics.rate_law("calcite")
    without = calcite.far_from_equilibrium(_Spec(25.0, 8.0))
    with_co2 = calcite.far_from_equilibrium(_Spec(25.0, 8.0, fugacity={"CO2(g)": -3.5}))
    assert with_co2 - without == pytest.approx(10 ** -3.48 * 10 ** -3.5)


def test_basalt_glass_constant():
    law = kinetics.basalt_glass(0.0084)
    (m,) = law.mechanisms
    # 10^-5.6 exp(-25.5 kJ / RT) at 25 C
    assert m.k(298.15) == pytest.approx(10 ** -5.6 * math.exp(-25.5 / (kinetics.R_KJ * 298.15)))
    assert law.per_unit == pytest.approx(1 / 0.0084) and not law.affinity


def test_calcite_dissolves_to_saturation(engine):
    result = engine.run(kinetic_card(reactants="1 g calcite"))
    run = result["modes"]["static"]
    left = run["kinetic"]["left:calcite"]
    assert left[0] == 100.0 and 98.0 < left[-1] < 100.0
    assert abs(run["saturation"]["calcite"][-1]) < 0.01
    rate = run["kinetic"]["rate:calcite"]
    assert rate[0] > 1e-8 and abs(rate[-1]) < 1e-3 * rate[0]
    assert run["minerals"]["calcite"][-1] == pytest.approx(left[-1] / 100 / 100.0869, rel=1e-3)
    assert result["x"]["key"] == "time" and result["x"]["unit"] == "days"


def test_answer_does_not_depend_on_the_output_steps(engine):
    ends = []
    for steps in (3, 21):
        run = engine.run(kinetic_card(reactants="10 g albite", duration=10, duration_unit="years",
                                      steps=steps))["modes"]["static"]
        ends.append(100.0 - run["kinetic"]["left:albite"][-1])
    assert ends[0] == pytest.approx(ends[1], rel=0.05)


def test_fast_mineral_dissolves_completely_and_stops(engine):
    run = engine.run(kinetic_card(reactants="0.05 g anhydrite"))["modes"]["static"]
    left = run["kinetic"]["left:anhydrite"]
    assert left[-1] == pytest.approx(0.0, abs=1e-6) and min(left) >= 0.0
    assert run["kinetic"]["rate:anhydrite"][-1] == 0.0


def test_supersaturated_mineral_grows_on_its_seed(engine):
    # surface seawater is supersaturated with calcite
    run = engine.run(kinetic_card(fluid="seawater", reactants="1 g calcite", duration=1,
                                  duration_unit="years"))["modes"]["static"]
    assert run["kinetic"]["left:calcite"][-1] > 100.0
    assert run["kinetic"]["rate:calcite"][1] < 0.0


def test_kinetic_quartz_is_not_replaced_by_its_polymorph(engine):
    result = engine.run(kinetic_card(reactants="10 g albite, 5 g quartz", duration=100, duration_unit="years",
                                     blocked_minerals=["lawsonite"]))
    run = result["modes"]["static"]
    assert "quartz,beta" not in run["minerals"] and "albite,low" not in run["minerals"]
    text = " ".join(story(result))
    assert "At the measured rates" in text and "albite dissolved" in text


def test_reactant_without_rate_is_added_steadily(engine):
    result = engine.run(kinetic_card(reactants="10 g albite, 0.001 mol NaCl"))
    assert any("No measured rate for NaCl" in n for n in result["notes"])


def test_card_and_chat():
    card, problems = Scenario().changed({"kinetics": "yes", "duration": "6", "duration_unit": "months"})
    assert problems == ["duration_unit: 'months' is not a time unit (use seconds, minutes, hours, days, years)"]
    assert card.kinetics is True
    assert "kinetics needs a reactant" in " ".join(card.problems())
    out = interpret_rules("react rainwater with 10 g albite, 1 g calcite and 5 g quartz for 100 years at 25 C",
                          Scenario())
    ch = out.scenario_changes
    assert ch["kinetics"] is True and ch["duration"] == 100 and ch["duration_unit"] == "years"
    assert ch["reactants"] == "10 g albite; 1 g calcite; 5 g quartz"
    ch = interpret_rules("react seawater with 10 g basalt glass over 6 months", Scenario()).scenario_changes
    assert ch["duration"] == pytest.approx(0.5) and ch["duration_unit"] == "years"
    ch = interpret_rules("set the surface area to 0.5 m2/g", Scenario()).scenario_changes
    assert ch["surface_area_cm2_per_g"] == 5000


def test_seawater_carbonate_laws_match_their_fits():
    from reaction.kinetics import SEAWATER_CARBONATE
    cal, arag = SEAWATER_CARBONATE["calcite"], SEAWATER_CARBONATE["aragonite_dong"]
    assert SEAWATER_CARBONATE["aragonite"].regimes == cal.regimes        # aragonite: calcite's fits, its own Omega
    # Naviaux et al. 2019 at 5 C: log k -13.07 mol/cm2/s, n 0.34 above Omega 0.75; -10.01, n 4.81 below
    assert math.log10(cal.rate_per_m2(5.0, math.log10(0.9))) == pytest.approx(-9.07 + 0.34 * math.log10(0.1), abs=1e-3)
    assert math.log10(cal.rate_per_m2(5.0, math.log10(0.3))) == pytest.approx(-6.01 + 4.81 * math.log10(0.7), abs=1e-3)
    # Dong et al. 2019: 0.013 g/g/day (1 - Omega)^1.37 on 2.46 m2/g
    g_per_g_day = arag.rate_per_m2(5.0, math.log10(0.5)) * 2.46 * 100.087 * 86400
    assert g_per_g_day == pytest.approx(0.013 * 0.5 ** 1.37, rel=1e-3)
    # smooth where the fits meet (Omega 0.79), no jump at 0.75
    near = [cal.rate_per_m2(5.0, math.log10(o)) for o in (0.749, 0.751, 0.789, 0.791)]
    assert all(a / b == pytest.approx(1.0, abs=0.06) for a, b in ((near[0], near[1]), (near[2], near[3])))  # the old jump was 2.3x
    assert cal.rate_per_m2(5.0, 0.1) == 0.0 and arag.rate_per_m2(5.0, 0.0) == 0.0     # dissolution only
    assert cal.rate_per_m2(21.0, -0.5) > 1.5 * cal.rate_per_m2(5.0, -0.5)                # Ea 25 kJ/mol


def test_reactant_surface_area_from_text():
    from reaction.app.scenario import parse_reactants, reactant_ssa_m2_g, reactant_text
    r = parse_reactants("20 ug forams at 3 m2/g; 10 g calcite; 5 ug pteropods")
    assert r[0] == {"name": "forams", "amount": 2e-5, "unit": "g", "ssa_m2_g": 3.0}
    assert reactant_text(r[0]) == "20 µg forams at 3 m2/g"
    assert [reactant_ssa_m2_g(x) for x in r] == [3.0, None, 2.0]
    assert reactant_ssa_m2_g(r[2], shells=False) is None
    assert parse_reactants("1 g calcite @ 500 cm2/g")[0]["ssa_m2_g"] == pytest.approx(0.05)
    assert parse_reactants([{"name": "calcite", "amount": 1, "unit": "g", "ssa_m2_g": 2}])[0]["ssa_m2_g"] == 2.0


def test_surface_follows_the_mass_left():
    from reaction.path import Reactant
    grain = Reactant("g", {}, 1.0, area=1.0)
    shell = Reactant("s", {}, 1.0, area=1.0, area_exponent=0.0)
    assert grain.surface(0.5) / 0.5 == pytest.approx(2 ** (1 / 3))        # per gram, a shrinking grain
    assert shell.surface(0.5) / 0.5 == pytest.approx(2.0)                 # a shell keeps its whole surface
    assert shell.surface(1.0) == 0.0


def test_takahashi_be_sinking_speed():
    from reaction.app.engine import sinking_speed
    assert sinking_speed(10.0) * 86400 == pytest.approx(10 ** 2.06 * 10 ** 0.64)       # about 500 m/day
    assert 450 < sinking_speed(10.0) * 86400 < 550

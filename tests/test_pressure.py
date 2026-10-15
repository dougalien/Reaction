"""Pressure as a free variable: checks on the thermodynamics and on runs along T-P paths.

The equilibrium constants come from Gibbs energies at the requested T and P, never from a database
tabulated at one pressure.  These tests check that against things we know independently:
the ionization constant of water (Marshall and Franck, 1981), the signs and sizes of reaction
volumes, mass conservation, and that the reliability indicator flags the right regions.
Needs the thermodynamic data (see app/config.py).
"""
import math

import pytest

from reaction.app import config
from reaction.app.engine import Engine
from reaction.app.scenario import Scenario
from reaction.thermo import validity

R_J = 8.314462618


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


@pytest.fixture(scope="module")
def db(engine):
    return engine.db


# ---- the reference formula itself -----------------------------------------------------------------
def test_marshall_franck_reproduces_the_textbook_value_at_25c():
    log_kw, rho = validity.marshall_franck_log_kw(25.0, 1.0132)
    assert log_kw == pytest.approx(-13.995, abs=0.002)
    assert rho == pytest.approx(0.997, abs=0.001)


# ---- the engine against it: the model's accuracy surface ---------------------------------------------
@pytest.mark.parametrize("t_c, p_bar, limit", [
    (25.0, 1.0132, 0.002), (25.0, 500.0, 0.02), (25.0, 1000.0, 0.04), (100.0, 500.0, 0.02),
    (200.0, 500.0, 0.06), (250.0, 1000.0, 0.12), (300.0, 500.0, 0.06), (300.0, 1000.0, 0.15),
    (350.0, 1000.0, 0.15), (400.0, 1000.0, 0.10), (500.0, 2000.0, 0.35),
])
def test_water_ionization_agrees_with_the_measured_constant(db, t_c, p_bar, limit):
    rel = validity.water_ionization_check(db, t_c, p_bar)
    assert abs(rel.deviation) < limit, rel


def test_the_reliability_flag_marks_the_known_bad_regions(db):
    good = validity.water_ionization_check(db, 250.0, 500.0)
    poor = validity.water_ionization_check(db, 400.0, 500.0)           # near-critical, density 0.58
    assert good.level == "good" and poor.level == "poor"
    assert validity.water_ionization_check(db, 500.0, 300.0).level == "poor"   # below 0.35 g/cm3
    assert validity.water_ionization_check(db, 25.0, 1.0132).level == "good"


def test_the_deviation_varies_smoothly_with_pressure_where_the_model_is_sound(db):
    """No jumps: neighbouring pressures differ by well under 0.05 log units (a smooth surface)."""
    for t_c in (25.0, 100.0, 200.0, 300.0):
        devs = [validity.water_ionization_check(db, t_c, p).deviation for p in range(300, 2001, 100)]
        assert max(abs(a - b) for a, b in zip(devs, devs[1:])) < 0.05


# ---- reaction volumes: the pressure dependence of logK ----------------------------------------------------
def _reaction_volume_cm3(db, species, t_c, p_bar, dp=20.0):
    """Delta V of the reaction as written (species coefficient -1), from d(logK)/dP."""
    hi = db.species[species].logk(t_c, p_bar + dp)
    lo = db.species[species].logk(t_c, p_bar - dp)
    dlogk_dp = (hi - lo) / (2.0 * dp)                                  # per bar
    # d(lnK)/dP = -dV/(RT)  ->  dV = -R T ln10 dlogK/dP   (J/bar/mol; 1 J/bar = 10 cm3)
    return -R_J * (t_c + 273.15) * math.log(10.0) * dlogk_dp * 10.0


def test_water_ionization_volume_at_25c_matches_experiment(db):
    # the database writes OH- = H2O - H+, the reverse of ionization, so its volume change is +22 cm3/mol
    # (ionization itself is -22.1 cm3/mol: Marshall and Franck 1981)
    dv = _reaction_volume_cm3(db, "OH-", 25.0, 100.0)
    assert dv == pytest.approx(22.1, abs=1.5)


def test_mineral_solubility_rises_with_pressure_when_the_volume_change_is_negative(db):
    for mineral in ("anhydrite", "calcite"):
        for t_c in (25.0, 150.0):
            # dissolution Ca2+ + anion -> aqueous ions loses volume (electrostriction): logK increases
            assert db.species[mineral].logk(t_c, 1000.0) > db.species[mineral].logk(t_c, 100.0), (mineral, t_c)
    # calcite + H+ = Ca++ + HCO3-: V(Ca++) + V(HCO3-) - V(calcite) = -18 + 24 - 37 = -31 cm3/mol at
    # 1 bar, smaller in size at 500 bar because the solution is compressed
    dv = _reaction_volume_cm3(db, "calcite", 25.0, 500.0)
    assert -35.0 < dv < -20.0


def test_logk_at_saturation_pressure_is_continuous_with_logk_at_other_pressures(db):
    sp = db.species["anhydrite"]
    base = sp.logk(150.0, db.reference_pressure(150.0))
    assert sp.logk(150.0, db.reference_pressure(150.0) + 0.01) == pytest.approx(base, abs=1e-4)


# ---- runs along paths -----------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def decompression(engine):
    return engine.run(Scenario(t_start_c=300.0, t_end_c=300.0, p_start_bar=1000.0, p_end_bar=150.0,
                               steps=8, mode="static"))


def test_a_pressure_only_path_is_described_in_pressure(decompression):
    r = decompression
    assert r["varies"] == ["pressure_bar"] and r["x"]["key"] == "pressure_bar"
    assert r["axes"]["temperature_c"]["values"] == [300.0] * 8
    assert r["axes"]["pressure_bar"]["values"][0] == 1000.0 and r["axes"]["pressure_bar"]["values"][-1] == 150.0
    assert "decompressed" in r["summary"]
    assert all(v is not None for v in r["modes"]["static"]["properties"]["pH"])


def test_sulfur_is_conserved_along_a_decompression(decompression):
    run, start = decompression["modes"]["static"], decompression["start_elements"]["S"]
    for i in range(8):
        dissolved = run["elements"]["S"][i] * run["properties"]["water_kg"][i]
        solid = sum(run["minerals"].get(m, [0.0] * 8)[i] for m in ("anhydrite", "kieserite"))
        assert dissolved + solid == pytest.approx(start, rel=2e-4)


def test_compression_dissolves_calcite_at_25c(engine):
    r = engine.run(Scenario(t_start_c=25.0, t_end_c=25.0, p_start_bar=1.0, p_end_bar=3000.0,
                            steps=8, mode="static"))
    calcite = r["modes"]["static"]["minerals"]["calcite"]
    assert calcite[0] > 0 and calcite[-1] < calcite[0] * 0.2          # calcite solubility rises with pressure
    assert all(a >= b - 1e-9 for a, b in zip(calcite, calcite[1:]))   # and never comes back


def test_a_combined_path_has_both_axes_and_matches_single_variable_runs_at_the_ends(engine):
    both = engine.run(Scenario(t_start_c=300.0, t_end_c=150.0, p_start_bar=1000.0, p_end_bar=300.0,
                               steps=6, mode="static"))
    assert both["varies"] == ["temperature_c", "pressure_bar"]
    point = engine.run(Scenario(t_start_c=300.0, t_end_c=300.0, p_start_bar=1000.0, p_end_bar=999.0,
                                steps=2, mode="static"))
    a, b = both["modes"]["static"], point["modes"]["static"]
    assert a["properties"]["pH"][0] == pytest.approx(b["properties"]["pH"][0], abs=1e-4)
    assert a["species"]["Mg+2"][0] == pytest.approx(b["species"]["Mg+2"][0], rel=1e-4)


def test_a_supercritical_path_is_flagged_not_hidden(engine):
    r = engine.run(Scenario(t_start_c=380.0, t_end_c=500.0, p_start_bar=800.0, p_end_bar=800.0,
                            steps=6, mode="static"))
    levels = r["reliability"]["levels"]
    assert levels[0] == "good" and levels[-1] == "poor"
    text = " ".join(r["notes"])
    assert "qualitative only" in text
    assert any(math.isfinite(d) for d in r["reliability"]["deviation"] if d is not None)


def test_advice_is_given_before_a_run(engine):
    risky = Scenario(t_start_c=380.0, t_end_c=500.0, p_start_bar=800.0, p_end_bar=800.0)
    safe = Scenario(t_start_c=25.0, t_end_c=250.0)
    assert engine.advisories(risky) and not engine.advisories(safe)
    assert engine.advisories(Scenario(t_start_c=100, t_end_c=100)) == []      # invalid card: no advice

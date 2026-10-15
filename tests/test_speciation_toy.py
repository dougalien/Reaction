"""Speciation checks that need no external data: analytic and hand-computed cases."""
import math

import pytest

from reaction import activity
from reaction.speciation import (Couple, LogActivity, Molal, PH, SpeciationError, speciate)
from toydb import make_toy_db


@pytest.fixture(scope="module")
def db():
    return make_toy_db()


def test_carbonate_system_matches_analytic_solution(db):
    """Fixed pH, ideal solution: carbonate partitioning follows K1 and K2 exactly."""
    ph, total = 8.0, 1.0e-3
    res = speciate(db, 25.0, {"H+": PH(ph), "HCO3-": Molal(total)}, activity_model="ideal")
    ah = 10 ** -ph
    ratio_co2 = ah / 10 ** -6.352
    ratio_co3 = 10 ** -10.329 / ah
    hco3 = total / (1 + ratio_co2 + ratio_co3)
    assert res.molality("HCO3-") == pytest.approx(hco3, rel=1e-9)
    assert res.molality("CO2(aq)") == pytest.approx(hco3 * ratio_co2, rel=1e-9)
    assert res.molality("CO3--") == pytest.approx(hco3 * ratio_co3, rel=1e-9)
    assert res.molality("OH-") == pytest.approx(10 ** -13.995 / ah, rel=1e-9)


def test_mass_balance_closes(db):
    res = speciate(db, 25.0, {"H+": PH(7.4), "HCO3-": Molal(2e-3), "Na+": Molal(0.1),
                              "Cl-": Molal(0.1)})
    carbon = sum(res.molality(n) for n in ("HCO3-", "CO3--", "CO2(aq)"))
    assert carbon == pytest.approx(2e-3, rel=1e-9)


def test_calcite_saturation_index(db):
    res = speciate(db, 25.0, {"H+": PH(8.0), "HCO3-": Molal(2e-3), "Ca++": Molal(1e-3)},
                   activity_model="ideal")
    expected = math.log10(1e-3 * res.molality("HCO3-") / 1e-8) - 1.85
    assert res.saturation["Calcite"] == pytest.approx(expected, abs=1e-9)


def test_nbs_scale_pins_chloride_to_bates_guggenheim(db):
    res = speciate(db, 25.0, {"H+": PH(7.0), "Na+": Molal(0.3), "Cl-": Molal(0.3)})
    expected = activity.log_gamma_chloride_bates_guggenheim(res.ionic_strength, 0.5114)
    assert res.species["Cl-"].log_gamma == pytest.approx(expected, abs=1e-12)


def test_internal_scale_uses_plain_bdot(db):
    res = speciate(db, 25.0, {"H+": PH(7.0), "Na+": Molal(0.3), "Cl-": Molal(0.3)},
                   ph_scale="internal")
    expected = activity.log_gamma_bdot(1.0, 4.0, res.ionic_strength, 0.5114, 0.3288, 0.041)
    assert res.species["Na+"].log_gamma == pytest.approx(expected, abs=1e-12)


def test_log_activity_constraint_fixes_activity(db):
    res = speciate(db, 25.0, {"H+": PH(7.0), "Na+": LogActivity(-1.5), "Cl-": Molal(0.05)})
    assert res.species["Na+"].log_activity == pytest.approx(-1.5, abs=1e-10)


def test_missing_ph_constraint_is_an_error(db):
    with pytest.raises(SpeciationError):
        speciate(db, 25.0, {"Na+": Molal(0.1), "Cl-": Molal(0.1)})


def test_unknown_basis_species_is_an_error(db):
    with pytest.raises(SpeciationError):
        speciate(db, 25.0, {"H+": PH(7.0), "Unobtainium+": Molal(1.0)})


def test_water_activity_limits():
    assert activity.log_activity_water_bdot(0.0, 0.0, 0.5114, 0.3288, 0.041) == 0.0
    # ideal limit: log aw -> -sum(m) / (Omega ln10) as the solution gets dilute
    dilute = activity.log_activity_water_bdot(1e-9, 1e-9, 0.5114, 0.3288, 0.041)
    assert dilute == pytest.approx(-1e-9 / (55.508 * math.log(10)), rel=1e-3)
    # at 1e-4 m the Debye-Hueckel term lowers the deviation from ideal by under 1%
    near = activity.log_activity_water_bdot(1e-4, 1e-4, 0.5114, 0.3288, 0.041)
    ideal = -1e-4 / (55.508 * math.log(10))
    assert 0.98 < near / ideal < 1.0


# --- constraints and model options added for the EQ3NR test library ------------

from reaction.speciation import Equilibrium, MgPerKg, MgPerL  # noqa: E402


def test_charge_balance_on_chloride(db):
    res = speciate(db, 25.0, {"H+": PH(7.0), "Na+": Molal(0.1), "Cl-": Molal(0.5)},
                   balance_on="Cl-")
    assert abs(res.charge_imbalance) < 1e-12
    # the input 0.5 m is only a starting value: chloride is whatever balances Na+ and H+
    assert res.totals["Cl-"] == pytest.approx(0.1 + res.molality("H+") - res.molality("OH-"),
                                              rel=1e-9)


def test_charge_balance_on_hydrogen_computes_ph(db):
    """Balancing on H+ turns the pH input into a starting guess (as in EQ3NR)."""
    res = speciate(db, 25.0, {"H+": PH(3.0), "Na+": Molal(1e-3), "HCO3-": Molal(1e-3)},
                   balance_on="H+", activity_model="ideal")
    assert abs(res.charge_imbalance) < 1e-12
    assert 8.0 < res.ph < 9.0                      # a sodium bicarbonate solution


def test_mineral_equilibrium_fixes_the_basis_species(db):
    res = speciate(db, 25.0, {"H+": PH(8.0), "HCO3-": Molal(2e-3),
                              "Ca++": Equilibrium("Calcite")})
    assert res.saturation["Calcite"] == pytest.approx(0.0, abs=1e-10)


def test_gas_equilibrium_sets_log_fugacity(db):
    """CO2(aq)-HCO3- equilibrium needs a gas in the database; build one on the fly."""
    from toydb import _sp
    from reaction.database.model import Kind
    g = make_toy_db()
    g.species["CO2(g)"] = _sp("CO2(g)", Kind.GAS, 0, {"C": 1, "O": 2},
                              {"CO2(g)": -1, "H2O": -1, "H+": 1, "HCO3-": 1}, -7.8136)
    res = speciate(g, 25.0, {"H+": PH(5.0), "HCO3-": Equilibrium("CO2(g)", -3.5)},
                   activity_model="ideal")
    assert res.log_fugacity["CO2(g)"] == pytest.approx(-3.5, abs=1e-10)
    assert res.molality("HCO3-") == pytest.approx(10 ** (-3.5 - 7.8136 + 5.0), rel=1e-9)


def test_inactive_phase_is_reported_but_cannot_constrain(db):
    res = speciate(db, 25.0, {"H+": PH(8.0), "HCO3-": Molal(2e-3), "Ca++": Molal(1e-3)},
                   suppress_phases=["Calcite"])
    assert "Calcite" in res.saturation
    with pytest.raises(SpeciationError):
        speciate(db, 25.0, {"H+": PH(8.0), "HCO3-": Molal(2e-3),
                            "Ca++": Equilibrium("Calcite")}, suppress_phases=["Calcite"])


def test_mass_unit_conversions(db):
    mw = db.molecular_weight("Na+")
    base = speciate(db, 25.0, {"H+": PH(7.0), "Na+": Molal(1.0 / mw)}, activity_model="ideal")
    kg = speciate(db, 25.0, {"H+": PH(7.0), "Na+": MgPerKg(1000.0)}, activity_model="ideal")
    assert kg.totals["Na+"] == pytest.approx(base.totals["Na+"], rel=1e-12)
    lit = speciate(db, 25.0, {"H+": PH(7.0), "Na+": MgPerL(1100.0)}, density=1.1,
                   activity_model="ideal")
    assert lit.totals["Na+"] == pytest.approx(base.totals["Na+"], rel=1e-12)
    tds = speciate(db, 25.0, {"H+": PH(7.0), "Na+": MgPerKg(1000.0)}, tds_mg_per_kg=100000.0,
                   activity_model="ideal")
    assert tds.totals["Na+"] == pytest.approx(base.totals["Na+"] / 0.9, rel=1e-12)


def test_basis_switch_does_not_change_the_chemistry(db):
    """Replacing HCO3- by CO3-- as the carbon basis species leaves every molality unchanged."""
    plain = speciate(db, 25.0, {"H+": PH(8.0), "HCO3-": Molal(1e-3)})
    switched = speciate(db, 25.0, {"H+": PH(8.0), "CO3--": Molal(1e-3)},
                        switches={"HCO3-": "CO3--"})
    for name, state in plain.species.items():
        assert switched.molality(name) == pytest.approx(state.molality, rel=1e-8), name


def test_decoupled_component_is_independent(db):
    """An auxiliary species given its own total no longer exchanges with its parents."""
    res = speciate(db, 25.0, {"H+": PH(8.0), "HCO3-": Molal(1e-3), "CO2(aq)": Molal(5e-4)},
                   activity_model="ideal")
    assert res.molality("CO2(aq)") == pytest.approx(5e-4, rel=1e-12)
    assert res.molality("HCO3-") + res.molality("CO3--") == pytest.approx(1e-3, rel=1e-9)


def test_suppressing_an_auxiliary_species_removes_it(db):
    res = speciate(db, 25.0, {"H+": PH(8.0), "HCO3-": Molal(1e-3)}, suppress=["CO2(aq)"])
    assert "CO2(aq)" not in res.species and "CO3--" in res.species


def test_davies_model_wiring(db):
    res = speciate(db, 25.0, {"H+": PH(7.0), "Na+": Molal(0.1), "Cl-": Molal(0.1),
                              "HCO3-": Molal(1e-3)},
                   activity_model="davies", ph_scale="internal")
    expected = activity.log_gamma_davies(1.0, res.ionic_strength, 0.5114)
    assert res.species["Na+"].log_gamma == pytest.approx(expected, abs=1e-12)
    # unlike B-dot, Davies gives every neutral species a unit activity coefficient
    assert res.species["CO2(aq)"].log_gamma == 0.0
    bdot = speciate(db, 25.0, {"H+": PH(7.0), "Na+": Molal(0.1), "Cl-": Molal(0.1),
                               "HCO3-": Molal(1e-3)}, ph_scale="internal")
    assert bdot.species["CO2(aq)"].log_gamma != 0.0


def test_redox_needs_oxygen_in_the_database(db):
    """The toy database has no O2(g), so a redox constraint cannot be honoured."""
    with pytest.raises(SpeciationError):
        speciate(db, 25.0, {"H+": PH(7.0), "Na+": Molal(0.1)}, redox=Couple("CO2(aq)"))

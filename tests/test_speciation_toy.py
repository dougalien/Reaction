"""Speciation checks that need no external data: analytic and hand-computed cases."""
import math

import pytest

from reaction import activity
from reaction.speciation import (LogActivity, Molal, PH, SpeciationError, speciate)
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

"""Validation: seawater (major ions) against the EQ3NR reference output.

EQ3/6 test case swmaj (data0.cmp, B-dot, NBS pH scale, 25 C), after
Nordstrom et al. (1979).  The comparison covers every species in EQ3NR's
printed table, the solid saturation indices and the gas fugacities.
"""
import math

import pytest

from reaction.speciation import Eh, LogFO2, MgPerKg, Molal, PH, speciate
from reaction.validation import parse_eq3nr_output

BASIS = {
    "Na+": MgPerKg(10768), "K+": MgPerKg(399.1), "Ca++": MgPerKg(412.3),
    "Mg++": MgPerKg(1291.8), "H+": PH(8.22), "HCO3-": Molal(2.022e-3),
    "Cl-": MgPerKg(19353), "SO4--": MgPerKg(2712),
}


@pytest.fixture(scope="module")
def reference(eq36_tests):
    return parse_eq3nr_output(eq36_tests / "3tlib_cmp" / "Outputs" / "swmaj.3o")


@pytest.fixture(scope="module")
def result(data0_cmp, reference):
    # Fix log fO2 at EQ3NR's value: EQ3/6 converts Eh to pe with a legacy
    # Faraday constant (see test_eh_to_pe), which would otherwise shift every
    # redox-sensitive species by up to ~0.01 log units.
    return speciate(data0_cmp, 25.0, BASIS, redox=LogFO2(reference.log_fo2))


def test_bulk_properties(result, reference):
    assert result.ionic_strength == pytest.approx(reference.ionic_strength, abs=6e-6)
    assert result.sum_molality == pytest.approx(reference.sum_molality, abs=6e-5)
    assert result.log_a_water == pytest.approx(reference.log_a_water, abs=2e-7)
    assert result.osmotic_coefficient == pytest.approx(reference.osmotic_coefficient, abs=6e-6)
    assert result.pmh == pytest.approx(reference.pmh, abs=6e-5)


def test_every_species_in_the_printed_table(result, reference):
    checked = 0
    for name, (_, log_m, log_gamma, log_a) in reference.species.items():
        if log_m < -12:        # the printed table keeps 5 significant figures
            continue
        got = result.species[name]
        assert math.log10(got.molality) == pytest.approx(log_m, abs=2e-4), name
        assert got.log_gamma == pytest.approx(log_gamma, abs=2e-4), name
        assert got.log_activity == pytest.approx(log_a, abs=2e-4), name
        checked += 1
    assert checked >= 30


def test_solid_saturation_indices(result, reference):
    checked = 0
    for name, log_qk in reference.solids.items():
        if name in result.saturation:
            assert result.saturation[name] == pytest.approx(log_qk, abs=5e-4), name
            checked += 1
    assert checked >= 35
    assert result.saturation["Calcite"] == pytest.approx(0.65034, abs=5e-4)


def test_gas_fugacities(result, reference):
    for name in ("CO2(g)", "H2O(g)", "HCl(g)"):
        assert result.log_fugacity[name] == pytest.approx(reference.gases[name], abs=5e-4), name


def test_eh_to_pe(data0_cmp, reference):
    """Reaction uses CODATA R and F; EQ3/6 uses a legacy F, so pe differs by ~6e-5 relative."""
    res = speciate(data0_cmp, 25.0, BASIS, redox=Eh(0.5))
    assert res.pe == pytest.approx(reference.pe, rel=1e-4)
    assert res.pe != pytest.approx(reference.pe, rel=1e-6)

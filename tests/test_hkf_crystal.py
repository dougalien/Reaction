"""HKF aqueous species and Maier-Kelley crystalline phases (published parameters, no files)."""
import math

import pytest
from scipy.integrate import quad

from reaction.thermo import hkf
from reaction.thermo.crystal import CrystalSpecies

CAL = 4.184


def na_plus() -> hkf.HKFSpecies:
    """Na+, Shock and Helgeson (1988); parameters in SUPCRT's scaled units, unscaled here."""
    return hkf.HKFSpecies("Na+", 1.0, -62591 * CAL, 13.96 * CAL,
                          1.839e-1 * CAL, -2.285e2 * CAL, 3.256 * CAL, -2.726e4 * CAL,
                          18.18 * CAL, -2.981e4 * CAL, 0.3306e5 * CAL)


def test_g_function_is_zero_where_defined_so_and_small_in_liquid_water():
    assert hkf.g_function(298.15, 1.0, 0.997) == pytest.approx(0.0, abs=1e-10)
    assert hkf.g_function(573.15, 7000.0, 0.72) == 0.0        # above 6 kbar
    assert hkf.g_function(298.15, 1.0, 1.05) == 0.0           # density of 1 g/cm3 or more
    # at 300 C in dense liquid the correction is a few thousandths of an angstrom, and its sign
    # is fixed by the Shock et al. (1992) constants (checked through Ca2+ logK against SUPCRT92)
    g_300 = hkf.g_function(573.15, 100.0, 0.72)
    assert -0.01 < g_300 < 0.0
    # it grows in magnitude as water expands toward the critical region
    assert abs(hkf.g_function(673.15, 250.0, 0.30)) > 10.0 * abs(g_300)


def test_born_coefficient_responds_to_the_g_function():
    ca = hkf.HKFSpecies("Ca+2", 2.0, -132120 * CAL, -13.5 * CAL, -0.1947e-1 * CAL, -7.252e2 * CAL,
                        5.2966 * CAL, -2.4792e4 * CAL, 9.0 * CAL, -2.522e4 * CAL, 1.2366e5 * CAL)
    assert ca.omega_pt(298.15, 1.0, 0.997) == pytest.approx(ca.omega, rel=1e-9)
    assert ca.omega_pt(673.15, 250.0, 0.30) != pytest.approx(ca.omega, rel=1e-3)
    h = hkf.HKFSpecies("H+", 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    assert h.omega_pt(673.15, 250.0, 0.30) == 0.0              # H+ is exempt


def test_hkf_reference_state_identities():
    sp = na_plus()
    assert sp.gibbs(298.15, 1.0) == pytest.approx(-62591 * CAL, rel=1e-12)
    h = 0.01
    slope = -(sp.gibbs(298.15 + h, 1.0) - sp.gibbs(298.15 - h, 1.0)) / (2 * h)
    assert slope == pytest.approx(13.96 * CAL, rel=2e-4)       # -dG/dT = S


def test_hkf_volume_and_heat_capacity_reproduce_tabulated_values():
    sp = na_plus()
    assert sp.volume_298() == pytest.approx(-1.11, abs=0.15)         # cm3/mol
    assert sp.heat_capacity_298() / CAL == pytest.approx(9.06, abs=0.1)   # cal/mol/K


def test_hkf_pressure_dependence_is_the_volume():
    sp = na_plus()
    h = 0.5
    v = (sp.gibbs(298.15, 1.0 + h) - sp.gibbs(298.15, 1.0 - h + 1e-9)) / (2 * h) / 0.1
    assert v == pytest.approx(sp.volume_298(), abs=0.01)


def test_aqueous_species_are_undefined_in_steam():
    with pytest.raises(hkf.WaterStateError):
        na_plus().gibbs(450.0, 2.0)                          # Psat(450 K) = 9.3 bar


def test_metastable_liquid_within_one_percent_of_saturation_is_accepted():
    na_plus().gibbs(373.15, 1.0132)                          # Psat(100 C) = 1.0142 bar


def calcite() -> CrystalSpecies:
    return CrystalSpecies("calcite", "cr", -269880 * CAL, 22.15 * CAL, 36.934,
                          24.98 * CAL, 0.00524 * CAL, -620000 * CAL, 0.0, 0.0, 1200.0)


def test_crystal_heat_capacity_integrals_match_quadrature():
    c = calcite()
    t = 700.0
    int_cp_dt, int_cp_dlnt = c.integrals(t)
    assert int_cp_dt == pytest.approx(quad(c.heat_capacity, 298.15, t)[0], rel=1e-9)
    assert int_cp_dlnt == pytest.approx(quad(lambda x: c.heat_capacity(x) / x, 298.15, t)[0], rel=1e-9)


def test_crystal_gibbs_energy_reference_state_and_pressure_term():
    c = calcite()
    assert c.gibbs(298.15, 1.0) == pytest.approx(-269880 * CAL, rel=1e-12)
    dg = c.gibbs(298.15, 1001.0) - c.gibbs(298.15, 1.0)
    assert dg == pytest.approx(36.934 * 1000.0 * 0.023901488 * CAL, rel=1e-9)


def test_constant_heat_capacity_fallback():
    sp = CrystalSpecies("x", "cr", 0.0, 10.0, 0.0, 0, 0, 0, 0, 0, 0.0, cp0=50.0)
    assert sp.gibbs(398.15, 1.0) == pytest.approx(
        -10.0 * 100.0 + 50.0 * 100.0 - 398.15 * 50.0 * math.log(398.15 / 298.15), rel=1e-12)


def test_quartz_alpha_beta_transition_volume_terms_are_continuous():
    alpha = CrystalSpecies("quartz", "cr", -204646 * CAL, 9.88 * CAL, 22.688, 11.22 * CAL,
                           0.0082 * CAL, -270000 * CAL, 0, 0, 848.0)
    beta = CrystalSpecies("quartz", "cr2", -204461.901 * CAL, 8.9986 * CAL, 23.06,
                          14.41 * CAL, 0.00194 * CAL, 0, 0, 0, 2000.0)
    # alpha is the stable form below 848 K, beta above; the transition lies near 848 K
    assert alpha.gibbs(600.0, 1.0) < beta.gibbs(600.0, 1.0)
    assert beta.gibbs(1000.0, 1.0) < alpha.gibbs(1000.0, 1.0)
    diff = alpha.gibbs(848.0, 1.0) - beta.gibbs(848.0, 1.0)
    assert abs(diff) < 60.0                                  # J/mol, small next to ~2e5 total

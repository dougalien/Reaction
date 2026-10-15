"""Dielectric constant, Born functions and Debye-Hueckel parameters (no external data)."""
import pytest

from reaction.thermo import electrostatics as el
from reaction.thermo import water


def test_dielectric_constant_at_25c():
    # Johnson and Norton (1991) / SUPCRT92: 78.245
    assert el.epsilon(298.15, 1.0) == pytest.approx(78.245, abs=0.01)


def test_debye_huckel_parameters_match_tabulated_values():
    # EQ3/6 data0 header values (from SUPCRT92) to their printed precision
    for t_c, p, a, b in [(0.01, 1.0132, 0.4939, 0.3253), (25.0, 1.0132, 0.5114, 0.3288),
                         (100.0, 1.0132, 0.5995, 0.3421), (200.0, 15.5365, 0.7994, 0.3639),
                         (300.0, 85.8378, 1.2180, 0.3925)]:
        got_a, got_b = el.debye_huckel_a_b(t_c + 273.15, p)
        assert got_a == pytest.approx(a, abs=3e-4), t_c
        assert got_b == pytest.approx(b, abs=2e-4), t_c


def test_born_y_at_reference_state():
    y, q, x = el.born_functions(298.15, 1.0)
    assert y == pytest.approx(-5.798e-5, rel=2e-3)          # 1/K, SUPCRT92
    assert x == pytest.approx(-3.0e-7, rel=0.1)
    assert q > 0                                           # eps grows with pressure


def test_trends_with_temperature_and_pressure():
    assert el.epsilon(373.15, 500.0) < el.epsilon(298.15, 500.0)
    assert el.epsilon(298.15, 1000.0) > el.epsilon(298.15, 1.0)
    a_cold, _ = el.debye_huckel_a_b(298.15, 500.0)
    a_hot, _ = el.debye_huckel_a_b(623.15, 500.0)
    assert a_hot > 2.0 * a_cold


def test_supercritical_water_is_a_poor_solvent():
    assert el.epsilon(673.15, 250.0) < 10.0                 # low density, low permittivity
    assert el.epsilon(673.15, 1000.0) > 15.0

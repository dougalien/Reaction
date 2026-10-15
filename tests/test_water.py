"""IAPWS-95: coefficients, derivatives, verification values and saturation consistency."""
import pytest

from reaction.thermo import water as w

# Wagner and Pruss (2002), Table 6.6: T (K), rho (kg/m3) -> p (MPa), s (kJ/kg/K), cv (kJ/kg/K)
TABLE_6_6 = [
    (300.0, 0.9965560e3, 0.992418352e-1, 0.393062643, 0.413018112e1),
    (300.0, 0.1005308e4, 0.200022515e2, 0.387405401, 0.406798347e1),
    (300.0, 0.1188202e4, 0.700004704e3, 0.132609616, 0.346135580e1),
    (500.0, 0.4350000, 0.999679423e-1, 0.794488271e1, 0.150817541e1),
    (500.0, 0.4532000e1, 0.999938125, 0.682502725e1, 0.166991025e1),
    (500.0, 0.8380250e3, 0.100003858e2, 0.256690919e1, 0.322106219e1),
    (500.0, 0.1084564e4, 0.700000405e3, 0.203237509e1, 0.307437693e1),
    (647.0, 0.3580000e3, 0.220384756e2, 0.432092307e1, 0.618315728e1),
    (900.0, 0.2410000, 0.100062559, None, None),
    (900.0, 0.5261500e2, 0.200000690e2, None, None),
    (900.0, 0.8707690e3, 0.700000006e3, None, None),
]


def test_coefficient_tables_are_complete():
    assert len(w._N) == 56
    assert len(w._D) == 54 and len(w._T) == 54
    assert len(w._C) == 44


@pytest.mark.parametrize("t,rho,p,s,cv", TABLE_6_6)
def test_pressure_matches_published_values(t, rho, p, s, cv):
    assert w.pressure_mpa(rho, t) == pytest.approx(p, rel=1e-8)


@pytest.mark.parametrize("t,rho,p,s,cv", [row for row in TABLE_6_6 if row[3] is not None])
def test_entropy_and_cv_match_published_values(t, rho, p, s, cv):
    delta, tau = rho / w.RHO_CRIT, w.T_CRIT / t
    i, r = w.phi_ideal(delta, tau), w.phi_residual(delta, tau)
    entropy = (tau * (i[3] + r[3]) - i[0] - r[0]) * w.R_SPECIFIC
    heat = -tau ** 2 * (i[4] + r[4]) * w.R_SPECIFIC
    assert entropy == pytest.approx(s, rel=1e-8)
    assert heat == pytest.approx(cv, rel=1e-8)


@pytest.mark.parametrize("delta,tau", [(0.5, 1.3), (1.0 + 1e-3, 1.03), (0.97, 1.01), (2.1, 0.8),
                                       (2.9, 0.6), (0.01, 2.0)])
def test_analytic_derivatives_agree_with_finite_differences(delta, tau):
    f = w.phi_residual(delta, tau)
    h = 1e-6
    dd = (w.phi_residual(delta + h, tau)[0] - w.phi_residual(delta - h, tau)[0]) / (2 * h)
    dt = (w.phi_residual(delta, tau + h)[0] - w.phi_residual(delta, tau - h)[0]) / (2 * h)
    assert f[1] == pytest.approx(dd, rel=1e-6, abs=1e-9)
    assert f[3] == pytest.approx(dt, rel=1e-6, abs=1e-9)
    ddd = (w.phi_residual(delta + h, tau)[1] - w.phi_residual(delta - h, tau)[1]) / (2 * h)
    dtt = (w.phi_residual(delta, tau + h)[3] - w.phi_residual(delta, tau - h)[3]) / (2 * h)
    ddt = (w.phi_residual(delta + h, tau)[3] - w.phi_residual(delta - h, tau)[3]) / (2 * h)
    assert f[2] == pytest.approx(ddd, rel=1e-6, abs=1e-8)
    assert f[4] == pytest.approx(dtt, rel=1e-6, abs=1e-8)
    assert f[5] == pytest.approx(ddt, rel=1e-6, abs=1e-8)


def _reduced_gibbs(rho, t):
    delta, tau = rho / w.RHO_CRIT, w.T_CRIT / t
    i, r = w.phi_ideal(delta, tau), w.phi_residual(delta, tau)
    return 1.0 + i[0] + r[0] + delta * r[1]


def _saturation_from_eos(t):
    """Phase equilibrium from the EOS alone: equal pressure and equal Gibbs energy."""
    from scipy.optimize import fsolve
    x0 = [w.saturated_liquid_density(t), w.saturated_vapor_density(t)]

    def residual(x):
        rl, rv = x
        return [(w.pressure_mpa(rl, t) - w.pressure_mpa(rv, t)) / w.pressure_mpa(rv, t) * 1e3,
                (_reduced_gibbs(rl, t) - _reduced_gibbs(rv, t)) * 1e3]

    rl, rv = fsolve(residual, x0, xtol=1e-13)
    return w.pressure_mpa(rv, t), rl, rv


def test_saturated_liquid_density_matches_published_value():
    """Wagner and Pruss (2002), Table 13.2: rho_liquid(300 K) = 996.513027 kg/m3."""
    _, rl, _ = _saturation_from_eos(300.0)
    assert rl == pytest.approx(996.513027, rel=1e-8)


@pytest.mark.parametrize("t", [280.0, 373.124, 450.0, 550.0, 620.0, 640.0])
def test_auxiliary_equations_agree_with_the_equation_of_state(t):
    ps, rl, rv = _saturation_from_eos(t)
    assert w.saturation_pressure_mpa(t) == pytest.approx(ps, rel=2e-4)
    assert w.saturated_liquid_density(t) == pytest.approx(rl, rel=2e-3)
    assert w.saturated_vapor_density(t) == pytest.approx(rv, rel=2e-3)


def test_known_saturation_and_density_values():
    assert w.saturation_pressure_mpa(373.124) == pytest.approx(0.101325, rel=1e-4)
    assert w.density(298.15, 1.0) == pytest.approx(997.047, abs=0.01)
    # above the boiling point (373.124 K) at 1 atm the stable phase is steam ...
    assert w.density(373.15, 1.01325) == pytest.approx(0.5976, abs=1e-3)
    # ... but the metastable liquid branch is available, as the EQ3/6 grid needs at 100 C
    assert w.density(373.15, 1.0132, liquid=True) == pytest.approx(958.35, abs=0.1)
    # supercritical and compressed states exist and are ordered sensibly
    assert w.density(673.15, 250.0) < w.density(673.15, 500.0) < w.density(673.15, 1000.0)


def test_density_solver_inverts_the_equation_of_state():
    for t, p in [(300.0, 1.0), (450.0, 500.0), (600.0, 300.0), (650.0, 250.0), (800.0, 1000.0),
                 (500.0, 5.0), (400.0, 0.5)]:
        rho = w.density(t, p)
        assert w.pressure_mpa(rho, t) * 10.0 == pytest.approx(p, rel=1e-9)


def test_vapour_branch_below_saturation():
    rho = w.density(450.0, 5.0)                     # steam: Psat(450 K) = 9.3 bar
    assert 1.0 < rho < 10.0


def test_expansivity_and_compressibility_of_liquid_water():
    alpha, beta = w.expansivity_and_compressibility(298.15, 1.0)
    assert alpha == pytest.approx(2.57e-4, rel=0.01)       # 1/K
    assert beta == pytest.approx(4.52e-5, rel=0.01)        # 1/bar


def test_supcrt_reference_state():
    cal = 4.184
    assert w.gibbs(298.15, 1.0) == pytest.approx(-56688.0 * cal, rel=1e-12)
    # the slope at the reference state is minus the third-law entropy, 16.712 cal/(mol K)
    h = 0.01
    slope = -(w.gibbs(298.15 + h, 1.0) - w.gibbs(298.15 - h, 1.0)) / (2 * h)
    assert slope == pytest.approx(16.712 * cal, rel=1e-6)


def test_heat_capacity_of_liquid_water():
    assert w.heat_capacity_p(298.15, 1.0) == pytest.approx(75.37, rel=2e-3)   # J/(mol K)

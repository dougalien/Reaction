"""The steady pore-water column (porewater.py) against its own mass balance and against Hales and Emerson
(1997, GCA 61, 501-514).  The bottom water is the Ceara Rise water profile (GLODAPv2.2023, station 198)."""
import numpy as np
import pytest

from reaction.app.scenario import example_card
from reaction.porewater import Sediment, d_sediment_m2_s, solve_column


def _bottom(z_m):
    t, s, ta, dic = example_card("ceara").profile_at(z_m)
    return s, t, 1.01325 + z_m * 1033.0 * 9.81 / 1e5, ta, dic


# Hales and Emerson (1997): Station G, 4685 m, CaCO3 35-50% (Fig. 2B), O2 flux 11-14 umol/cm2/yr (Table 2),
# k best fit 125%/day, n 4.5 (Eqn. 4), beta 0.64; the metabolic CO2 scale depth about 0.5 cm (Section 5.1);
# porosity 0.65 from their formation factors 2.5-3 and b 2.3 (Eqn. 1)
STATION_G = Sediment(porosity=0.65, caco3_wt_pct=40.0, k_pct_day=125.0, n=4.5, o2_flux_umol_cm2_yr=13.0,
                     beta=0.64, resp_depth_cm=0.5)


def test_column_converges_with_the_grid_and_closes_its_mass_balance():
    s, t, p, ta, dic = _bottom(4685.0)
    coarse = solve_column(STATION_G, s, t, p, ta, dic, depth_cm=20.0, nodes=80)
    fine = solve_column(STATION_G, s, t, p, ta, dic, depth_cm=20.0, nodes=320, first_cm=1e-4)
    assert coarse.flux_umol_cm2_yr == pytest.approx(fine.flux_umol_cm2_yr, rel=2e-3)
    # the DIC leaving through the seafloor is the calcite dissolved plus the CO2 respired
    d = d_sediment_m2_s(t, 0.65)
    grad = (fine.dic[1] - fine.dic[0]) / ((fine.z_cm[1] - fine.z_cm[0]) / 100.0)
    out = d * grad * 0.65 * 1025.0 * 365.25 * 86400.0 * 1e6 / 1e4          # umol/cm2/yr
    respired = 0.64 * 13.0 * (1.0 - np.exp(-20.0 / 0.5))
    assert out == pytest.approx(fine.flux_umol_cm2_yr + respired, rel=2e-3)
    assert fine.omega[0] < 1.0 and np.all(np.diff(fine.omega) > 0)          # less corrosive with depth


def test_column_dissolution_falls_in_hales_and_emerson_table_2():
    """Their Table 2, beta 0.64: Station A (4000 m) 6-9, Station G (4685 m) 8-23 umol CaCO3/cm2/yr.  Our
    bottom water is GLODAP's (Omega 1.12 at 4000 m, 0.86 below 4619 m) where theirs was 1.02 and 0.78."""
    s, t, p, ta, dic = _bottom(4000.0)
    a = Sediment(porosity=0.65, caco3_wt_pct=65.0, k_pct_day=125.0, n=4.5, o2_flux_umol_cm2_yr=24.0, beta=0.64,
                 resp_depth_cm=0.5)
    assert 6.0 <= solve_column(a, s, t, p, ta, dic).flux_umol_cm2_yr <= 9.0
    s, t, p, ta, dic = _bottom(4685.0)
    assert 8.0 <= solve_column(STATION_G, s, t, p, ta, dic).flux_umol_cm2_yr <= 23.0


def test_no_respiration_and_supersaturated_bottom_water_dissolves_nothing():
    s, t, p, ta, dic = _bottom(3000.0)
    sed = Sediment(porosity=0.65, caco3_wt_pct=65.0, k_pct_day=125.0, n=4.5, o2_flux_umol_cm2_yr=0.0, beta=0.64,
                   resp_depth_cm=0.5)
    col = solve_column(sed, s, t, p, ta, dic)
    assert col.flux_umol_cm2_yr < 1e-6 and np.allclose(col.dic, dic, rtol=1e-9)


def test_biodiffusion_reaches_its_two_limits():
    """No mixing: M = exp(-int lam dz / w).  Fast mixing over L with lam only inside it: below the mixed
    layer M = 1 / (1 + lam L / w) (a well-mixed reactor fed by the rain).  No reaction: M = 1."""
    from reaction.porewater import mixed_mass_left
    z = np.linspace(0.0, 30.0, 6001)
    assert np.allclose(mixed_mass_left(z, np.zeros_like(z), 0.002, 0.2, 10.0), 1.0)
    lam = np.full_like(z, 1e-4)
    m = mixed_mass_left(z, lam, 0.002, 0.0, 10.0)
    assert m[2000] == pytest.approx(np.exp(-1e-4 * z[2000] / 0.002), rel=2e-3)
    lam = np.where(z < 10.0, 1e-4, 0.0)
    assert mixed_mass_left(z, lam, 0.002, 1e3, 10.0)[-1] == pytest.approx(1.0 / (1.0 + 1e-4 * 10.0 / 0.002),
                                                                           rel=2e-3)

"""The 2026 activity-model upgrade: b-gamma (HKF 1981), the IAPWS R11-24 water ionization constant,
Harvie's J integral, and the Pitzer model with pitzer.dat's own constants checked against PHREEQC."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from reaction.pitzer import _j
from reaction.thermo import bgamma as bg
from reaction.thermo import validity


# ---- IAPWS R11-24 (2024): its own verification values (Table 3 of the release) --------------------------
@pytest.mark.parametrize("t_k, rho, pkw", [
    (300.0, 1.0, 13.906672), (600.0, 0.07, 20.161651), (600.0, 0.7, 11.147093),
    (800.0, 0.2, 14.487671), (800.0, 1.2, 6.4058649), (1270.0, 0.0, 35.081557),
])
def test_iapws_r11_24_reproduces_its_check_values(t_k, rho, pkw):
    assert validity.iapws_pkw(t_k, rho) == pytest.approx(pkw, abs=1e-6)


def test_iapws_log_kw_at_25c_is_the_textbook_value():
    log_kw, rho = validity.iapws_log_kw(25.0, 1.01325)
    assert log_kw == pytest.approx(-13.995, abs=0.005)


# ---- the J integral of unsymmetrical mixing -------------------------------------------------------------
def _j_by_quadrature(x):
    """Pitzer (1975, eq. 47): J(x) = x/4 - 1 + (1/x) int_0^inf [1 - exp(-(x/y) e^-y)] y^2 dy."""
    f = lambda y: (1.0 - math.exp(-(x / y) * math.exp(-y))) * y * y
    return x / 4.0 - 1.0 + quad(f, 0.0, np.inf, limit=400)[0] / x


@pytest.mark.parametrize("x", [0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 50.0])
def test_harvie_chebyshev_j_matches_the_integral(x):
    j, dj = _j(np.array([x]))
    assert j[0] == pytest.approx(_j_by_quadrature(x), abs=2e-5)
    h = 1e-4 * x
    slope = (_j_by_quadrature(x + h) - _j_by_quadrature(x - h)) / (2 * h)
    assert dj[0] == pytest.approx(slope, rel=1e-3, abs=1e-5)


# ---- b-gamma --------------------------------------------------------------------------------------------
def test_the_fmm_spline_is_exact_for_a_cubic():
    x = np.array([0.0, 1.0, 2.5, 4.0, 5.0, 7.0])
    f = lambda u: 0.3 * u ** 3 - 2.0 * u ** 2 + u - 4.0
    s = bg.FmmSpline(x, f(x))
    for u in (0.5, 3.3, 6.1, 7.5):
        assert s(u) == pytest.approx(f(u), abs=1e-10)


@pytest.mark.parametrize("p_bar", [500, 1000, 2000, 5000])
def test_b_gamma_returns_the_tabulated_values_at_the_nodes(p_bar):
    temps, values = bg._ISOBARS[p_bar]
    for t, v in zip(temps, values):
        if t > 25:                 # at 25 C the 500 bar node sits close to the Psat node; still a node
            assert bg.bgamma(float(t), float(p_bar)) == pytest.approx(v / 100.0, abs=1e-9), (t, p_bar)


def test_b_gamma_is_continuous_and_flags_extrapolation():
    grid = [bg.bgamma(t, 2000.0) for t in np.arange(25.0, 600.0, 5.0)]
    assert max(abs(a - b) for a, b in zip(grid, grid[1:])) < 0.01
    assert bg.extrapolated(200.0, 1000.0) is None
    assert bg.extrapolated(300.0, 10.0)                     # past the end of the saturation row
    assert bg.extrapolated(300.0, 8000.0)                   # above 5 kbar


# ---- Pitzer with pitzer.dat constants against PHREEQC 3.7.3 ----------------------------------------------
# PHREEQC 3.7.3 (IPhreeqc via phreeqpython 1.6.2), pitzer.dat of phreeqc3 commit be755855 with its
# MEAN_GAMMAS and GAS_BINARY_PARAMETERS blocks removed, -MacInnes false, pH 7, mol/kgw; run 2026-11-11.
PHREEQC_SI = [
    (25, 1, {"Na": 4.0, "Cl": 4.0, "Ca": 0.01, "S(6)": 0.01, "C(4)": 0.002},
     {"halite": -0.589, "anhydrite": -1.6554, "calcite": -0.2226, "CO2(g)": -2.1652}),
    (25, 500, {"Na": 4.0, "Cl": 4.0, "Ca": 0.01, "S(6)": 0.01, "C(4)": 0.002},
     {"halite": -0.6271, "anhydrite": -1.7905, "calcite": -0.3266, "CO2(g)": -2.0385}),
    (100, 1, {"Na": 2.0, "Cl": 2.0, "Ca": 0.02, "S(6)": 0.02, "C(4)": 0.002},
     {"halite": -1.3781, "anhydrite": -0.2864, "calcite": 0.5829, "CO2(g)": -1.6708}),
    (100, 800, {"Na": 2.0, "Cl": 2.0, "Ca": 0.02, "S(6)": 0.02, "C(4)": 0.002},
     {"halite": -1.4163, "anhydrite": -0.5763, "calcite": 0.4721, "CO2(g)": -1.5142}),
    (150, 300, {"Na": 5.0, "Cl": 5.0, "Ca": 0.01, "S(6)": 0.01, "C(4)": 0.001},
     {"halite": -0.4367, "anhydrite": -0.6691, "calcite": 0.4653, "CO2(g)": -1.6695}),
    (25, 1, {"Na": 0.5, "Cl": 0.6, "Mg": 0.05, "Ca": 0.01, "S(6)": 0.03, "C(4)": 0.002},
     {"halite": -2.4609, "anhydrite": -0.9824, "calcite": -0.539, "CO2(g)": -2.1123}),
]


@pytest.fixture(scope="module")
def db():
    from reaction.app import config
    from reaction.app.engine import Engine
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine().db


@pytest.mark.parametrize("t_c, p_atm, comp, ref", PHREEQC_SI)
def test_pitzer_saturation_matches_phreeqc(db, t_c, p_atm, comp, ref):
    from reaction.speciation import PH, Molal, speciate
    names = {"Na": "Na+", "Cl": "Cl-", "Ca": "Ca+2", "S(6)": "SO4-2", "C(4)": "HCO3-", "Mg": "Mg+2"}
    basis = {names[k]: Molal(v) for k, v in comp.items()}
    basis["H+"] = PH(7.0)
    sp = speciate(db, float(t_c), basis, pressure_bar=p_atm * 1.01325, activity_model="pitzer",
                  ph_scale="internal")
    for name, si in ref.items():
        mine = sp.log_fugacity["CO2(g)"] if name == "CO2(g)" else sp.saturation[name]
        assert mine == pytest.approx(si, abs=0.005), name


def test_a_hot_or_deep_b_dot_card_is_pointed_to_b_gamma():
    from reaction.app.engine import bdot_range_note
    from reaction.app.scenario import Scenario
    assert bdot_range_note(Scenario(t_start_c=25.0, t_end_c=250.0)) is None
    assert "b-gamma" in bdot_range_note(Scenario(t_start_c=25.0, t_end_c=350.0))
    assert "b-gamma" in bdot_range_note(Scenario(t_start_c=25.0, t_end_c=200.0, p_start_bar=2000.0,
                                                 p_end_bar=2000.0))
    assert bdot_range_note(Scenario(t_start_c=25.0, t_end_c=350.0, activity_model="bgamma")) is None

"""The marine carbonate system against PyCO2SYS 1.8 (default constants), to 0.1%."""
import pytest

from reaction import marine

# S, T C, P bar (absolute), TA, DIC umol/kg -> pH total, CO3 umol/kg, Omega calcite, Omega aragonite (PyCO2SYS)
CASES = [((35, 25, 1.01325, 2300, 2000), (8.0459, 213.41, 5.1373, 3.3862)),
         ((35, 2, 401.0, 2300, 2000), (8.2544, 190.53, 2.0574, 1.3599)),
         ((34.9, 1.5, 556.0, 2350, 2250), (7.7376, 75.52, 0.6080, 0.4094)),
         ((20, 30, 200.0, 2000, 1900), (7.7066, 89.75, 1.9207, 1.2368))]


@pytest.mark.parametrize("given,expected", CASES)
def test_matches_pyco2sys(given, expected):
    s, t, p, ta, dic = given
    st = marine.solve(s, t, p, ta * 1e-6, dic * 1e-6)
    ph, co3, oc, oa = expected
    assert st.ph_total == pytest.approx(ph, abs=0.002)
    assert st.co3 * 1e6 == pytest.approx(co3, rel=1e-3)
    assert st.omega_calcite == pytest.approx(oc, rel=1e-3) and st.omega_aragonite == pytest.approx(oa, rel=1e-3)


def test_aragonite_is_half_again_as_soluble():
    c = marine.constants(35, 2.0)
    assert c.ksp_aragonite / c.ksp_calcite == pytest.approx(1.5, abs=0.1)


def test_ceara_stations_bracket_the_calcite_saturation_horizon():
    """The three GLODAP Ceara Rise profiles of scripts/ceara_bracket.py, in order of their calcite saturation
    horizon: station 198 (1994) about 4275 m, 484 (2015) about 4380 m, 39 (1972) about 4470 m."""
    from reaction.app.scenario import (Scenario, WATER_PROFILE_CEARA, WATER_PROFILE_CEARA_1972,
                                       WATER_PROFILE_CEARA_2015)
    from reaction.marine import solve

    def horizon(prof):
        card = Scenario(water_profile=prof)
        for z in range(3000, 4700, 5):
            t, s, ta, dic = card.profile_at(float(z))
            p = 1.01325 + z * 1033.0 * 9.81 / 1e5
            if solve(s, t, p, ta, dic).omega_calcite < 1.0:
                return z
        return None

    zs = [horizon(p) for p in (WATER_PROFILE_CEARA, WATER_PROFILE_CEARA_2015, WATER_PROFILE_CEARA_1972)]
    assert zs == sorted(zs) and 4200 < zs[0] < 4330 < zs[1] < 4430 < zs[2] < 4520

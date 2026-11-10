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

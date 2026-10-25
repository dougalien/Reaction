"""Minerals from the Berman (1988) dataset, as OBIGT's default data use them."""
import pytest
from scipy.optimize import fsolve

from reaction.thermo import berman
from reaction.thermo.eosdb import PRESETS
from reaction.thermo.obigt import load_obigt


@pytest.fixture(scope="module")
def current(obigt_dir):
    if not (obigt_dir / "Berman").is_dir():
        pytest.skip("the Berman parameter files (CHNOSZ extdata/Berman) are not in the OBIGT folder")
    return load_obigt(obigt_dir, PRESETS["current"])


def _g(data, name, t_k, p_bar):
    return data.polymorphs(name)[0].eos.gibbs(t_k, p_bar)


def test_berman_minerals_replace_the_older_silicates(current):
    for name in ("albite", "K-feldspar", "quartz", "clinochlore", "epidote", "daphnite"):
        poly = current.polymorphs(name)
        assert len(poly) == 1 and poly[0].model == "Berman", name
    for old in ("microcline,maximum", "clinochlore,14A", "dolomite,ordered", "titanite"):
        assert not current.polymorphs(old), old
    assert current.polymorphs("laumontite")        # no Berman equivalent: the older entry stays


def test_youngest_file_wins_and_g_comes_from_h_and_s(current):
    ab = current.polymorphs("albite")[0].eos
    assert ab.params.source == "SHD91_1991.csv" and ab.params.h0 == -3941903
    for name in ("albite", "muscovite", "calcite", "quartz", "talc"):
        e = current.polymorphs(name)[0].eos
        assert e.gibbs(298.15, 1.0) == pytest.approx(e.g_ref)
        assert e.g_ref == pytest.approx(e.params.g0, abs=100)          # tabulated G agrees


def test_quartz_agrees_with_supcrt92(obigt_dir, current):
    old = load_obigt(obigt_dir, PRESETS["supcrt92"])
    for t_k, p in ((298.15, 1.0), (623.15, 500.0)):
        g_old = min(sp.eos.gibbs(t_k, p) for sp in old.polymorphs("quartz"))
        assert _g(current, "quartz", t_k, p) == pytest.approx(g_old, abs=150)


def test_aluminosilicate_triple_point(current):
    """Berman (1988) places the andalusite-kyanite-sillimanite triple point near 505 C, 3.7 kbar."""
    def f(x):
        t, p = float(x[0]), float(x[1])
        return [_g(current, "andalusite", t, p) - _g(current, "kyanite", t, p),
                _g(current, "andalusite", t, p) - _g(current, "sillimanite", t, p)]
    t, p = fsolve(f, [800.0, 4000.0])
    assert t - 273.15 == pytest.approx(505.6, abs=1.0) and p == pytest.approx(3737, abs=30)


def test_entropy_and_enthalpy_are_consistent(current):
    """S = -dG/dT through the quartz lambda and K-feldspar disorder terms (above 1 bar the
    disorder volume makes this approximate, as in Berman, 1988, and CHNOSZ)."""
    for name in ("quartz", "K-feldspar", "albite"):
        e = current.polymorphs(name)[0].eos
        for t_k, p in ((500.0, 1.0), (900.0, 2000.0), (1200.0, 5000.0)):
            g, h, s = e.properties(t_k, p)
            dt = 1e-3
            dgdt = (e.gibbs(t_k + dt, p) - e.gibbs(t_k - dt, p)) / (2 * dt)
            assert -dgdt == pytest.approx(s, rel=1e-6 if p == 1.0 else 1e-3), (name, t_k)


def test_element_entropy():
    assert berman.element_entropy({"Si": 1, "O": 2}) == pytest.approx((4.49569789674952 + 49.0325047801147) * 4.184)


def test_without_the_parameter_files_the_older_minerals_stay(obigt_dir, tmp_path):
    from reaction.thermo.eosdb import load_eos_database
    for f in obigt_dir.glob("*.csv"):
        (tmp_path / f.name).write_bytes(f.read_bytes())
    db = load_eos_database(tmp_path)
    assert db.species["albite"].provenance.file == "SUPCRT92.csv"
    assert any("Berman" in n for n in db.data_notes)

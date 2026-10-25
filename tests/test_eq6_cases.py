"""Reaction paths against EQ6's own output (the oracle for the path engine, not a design driver).

Only the cases the engine reproduces are listed.  Known open cases (redox fronts with mineral
reactants: methane, pyrsw; solid solutions and kinetics) are not tested here.
"""
import pytest

from reaction.validation.eq6_run import exceeds, run_case

MATCHING = ["micro", "rwtitr", "swtitr"]


@pytest.mark.parametrize("name", MATCHING)
def test_eq6_case_is_reproduced(name, data0_cmp, eq36_tests):
    path = eq36_tests / "6tlib_cmp" / f"{name}.6i"
    if not path.exists():
        pytest.skip(f"{path} not found")
    report, points = run_case(data0_cmp, path)
    bad = [d for d in report.deviations if exceeds(d.what, d)]
    assert report.n_points == len(points) > 5
    assert len(report.deviations) > 100
    worst = sorted(bad, key=lambda d: -d.delta)[:3]
    assert not bad, [(d.what, d.xi, d.ours, d.reference) for d in worst]


def test_dedolo_alias_minerals(data0_cmp, eq36_tests):
    """Dolomite and Dolomite-ord have the same formula and logK within 2e-3: the assemblage must take
    one of them, never both, and keep the reactant's name.  What remains is EQ6 reporting the reactant
    as unconverted at xi = 0 and a few nanomoles of calcite EQ6's tolerances do not precipitate."""
    path = eq36_tests / "6tlib_cmp" / "dedolo.6i"
    if not path.exists():
        pytest.skip(f"{path} not found")
    report, points = run_case(data0_cmp, path)
    bad = [d for d in report.deviations if exceeds(d.what, d)]
    assert report.n_points == len(points) == 5
    assert all("Dolomite" not in p.minerals or "Dolomite-ord" not in p.minerals for p in points)
    for d in bad:
        assert (d.what.startswith("mineral") and d.xi == 0.0) or \
               (d.what == "mineral Calcite" and d.ours < 1e-8), (d.what, d.xi, d.ours, d.reference)


def test_microft_flow_through(data0_cmp, eq36_tests):
    """Fluid-centered flow-through: diaspore, then kaolinite, are left behind where they would start to
    redissolve, and the saturated reactant stops reacting without joining the system.  What remains is
    pH within 1.5e-3 just past the muscovite boundary, where EQ6 moves kaolinite out in partial shifts
    over several small steps instead of all at once."""
    path = eq36_tests / "6tlib_cmp" / "microft.6i"
    if not path.exists():
        pytest.skip(f"{path} not found")
    report, points = run_case(data0_cmp, path)
    bad = [d for d in report.deviations if exceeds(d.what, d)]
    assert report.n_points == len(points) == 16
    assert points[-1].removed.get("Kaolinite", 0.0) > 2e-5
    for d in bad:
        assert d.what in ("pH", "pmH") and d.delta < 1.5e-3 and d.xi > 1.02e-4, (d.what, d.xi, d.ours)

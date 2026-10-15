"""Reaction paths against EQ6's own output (the oracle for the path engine, not a design driver).

Only the cases the engine reproduces are listed.  Known open cases (redox fronts with mineral
reactants: methane, pyrsw; flow-through, solid solutions and kinetics) are not tested here.
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

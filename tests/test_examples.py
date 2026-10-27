"""The one-button examples on the real engine: each shows the feature it is there to show.  Needs the
thermodynamic data.  The column examples take a minute each and are left out here (test_clays.py runs the
stream shortened)."""
import pytest

from reaction.app import config
from reaction.app.engine import Engine
from reaction.app.scenario import EXAMPLES, example_card


@pytest.fixture(scope="module")
def runs():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    engine = Engine()
    return {k: next(iter(engine.run(example_card(k))["modes"].values())) for k, e in EXAMPLES.items()
            if e.changes.get("kind") != "transport"}


def _formed(run):
    return {m for m, v in run["minerals"].items() if v and max(x or 0 for x in v) > 0}


def test_basalt_alteration_forms_mixtures(runs):
    run = runs["solid_solutions"]
    assert {"chlorite", "epidote-clinozoisite", "actinolite", "albite"} <= _formed(run)
    assert 0.05 < run["composition"]["chlorite:daphnite"][-1] < 0.95


def test_boiling_seawater_leaves_a_saltier_brine(runs):
    phases = runs["phase_separation"]["phases"]
    assert phases["brine_salinity"][-1] > 6 and phases["vapour_salinity"][-1] < 1
    assert sum(v is None for v in runs["phase_separation"]["properties"]["pH"]) <= 1


def test_calcite_dissolves_ever_more_slowly(runs):
    left = runs["kinetics"]["kinetic"]["left:calcite"]
    assert 97 < left[-1] < 99.5
    assert left[0] - left[5] > 10 * (left[-6] - left[-1])            # fast at first, then near saturation


def test_sinking_foram_loses_mass(runs):
    run = runs["forams"]
    ug = [v * 100.0869e6 for v in run["minerals"]["calcite"]]                 # the foram's calcite, µg
    assert ug[0] == pytest.approx(20.0, rel=1e-3) and 14 < ug[-1] < 18
    sat = run["saturation"]["calcite"]
    assert sat[0] < -0.3 and sat[-1] < sat[0]                                # more corrosive with depth
    assert all(a > b for a, b in zip(ug, ug[1:]))
    ca = run["elements"]["Ca"]
    assert abs(ca[-1] - ca[0]) < 1e-4 * ca[0]                                # the water hardly notices


def test_evaporating_lake_lays_down_salts_in_order(runs):
    run = runs["lake"]
    first = {m: next(i for i, x in enumerate(v) if x) for m, v in run["minerals"].items() if v and any(v)}
    assert first["halite"] < first["anhydrite"] < first["epsomite"]
    assert run["isotopes"]["d18O:water"][-1] > 10


def test_boiling_vents_co2_and_scales_calcite(runs):
    run = runs["boiling"]
    assert run["gas"]["CO2(g)"][0] == 0 and run["gas"]["CO2(g)"][-1] > 0.01
    assert run["elements"]["C"][-1] < 0.5 * run["elements"]["C"][0]
    assert run["properties"]["pH"][-1] > run["properties"]["pH"][0] + 0.5
    calcite = run["minerals"]["calcite"]
    assert calcite[-1] > 3 * calcite[0]


def test_serpentinization_makes_iron_serpentine_and_hydrogen(runs):
    run = runs["serpentinization"]
    assert {"serpentine", "brucite", "magnetite"} <= _formed(run)
    assert 0.01 < run["composition"]["serpentine:greenalite"][-1] < 0.1
    assert 0.1 < run["species"]["H2"][-1] < 0.25


def test_each_example_opens_on_its_own_plot():
    """What an example shows first must be in its run (checked on the runs above, by name only)."""
    from reaction.app.view import resolve_item
    engine = Engine() if (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists() else None
    if engine is None:
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    for key, e in EXAMPLES.items():
        assert e.plot, key
        if e.changes.get("kind") == "transport":
            continue
        result = engine.run(example_card(key))
        assert all(resolve_item(t, result) for t in e.plot), key

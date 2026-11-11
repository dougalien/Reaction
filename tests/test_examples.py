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
    out = {}
    for k, e in EXAMPLES.items():
        if e.changes.get("kind") != "transport":
            res = engine.run(example_card(k))
            out[k] = {"axes": res["axes"], "varies": res["varies"]} | next(iter(res["modes"].values()))
    return out


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


def _shell(run, mineral):
    depth = run["axes"]["depth"]["values"]
    ug = [v * 100.087e6 for v in run["minerals"][mineral]]
    return depth, ug, run["saturation"][mineral]


def test_sinking_foram_crosses_the_lysocline(runs):
    card = example_card("forams")
    assert card.mode == "flush" and card.carbonate_rates == "seawater"
    run = runs["forams"]
    assert run["varies"][0] == "depth"
    depth, ug, sat = _shell(run, "calcite")
    assert ug[0] == pytest.approx(20.0, rel=1e-3)
    at = lambda z: min(range(len(depth)), key=lambda i: abs(depth[i] - z))          # noqa: E731
    assert 2500 < depth[next(i for i, x in enumerate(sat) if x < 0)] < 3500         # saturation horizon
    assert ug[at(2500)] > 19.999                                                    # nothing above it
    assert ug[at(2500)] - ug[at(4000)] < 1.5                                        # slow just below it
    assert ug[at(4000)] - ug[at(5000)] > 8 * (ug[at(2500)] - ug[at(4000)])          # the lysocline
    assert all(a >= b for a, b in zip(ug, ug[1:]))
    assert 8 < ug[-1] < 17                                  # about 6 days at the Takahashi and Be speed
    days = run["axes"]["time"]["values"]
    assert 5 < days[-1] < 7 and card.sinking
    ca = run["elements"]["Ca"]
    assert abs(ca[-1] - ca[0]) < 1e-4 * ca[0]                                       # the water hardly notices


def test_pteropod_aragonite_starts_dissolving_shallower(runs):
    run = runs["shells"]
    depth, cal, cal_sat = _shell(run, "calcite")
    _, arag, arag_sat = _shell(run, "aragonite")
    first = lambda v: next(i for i in range(1, len(v)) if v[i] < v[0] - 1e-4)       # noqa: E731
    assert first(arag) < first(cal)
    assert all(a < c for a, c in zip(arag_sat, cal_sat))
    half = lambda v: next(i for i in range(len(v)) if v[i] < 10)                    # noqa: E731
    assert depth[half(arag)] < 5300 and arag[-1] < 0.1 and cal[-1] > 8
    alone = _shell(runs["forams"], "calcite")[1]                      # flush: fresh water at each step, so the
    assert max(abs(a - b) for a, b in zip(cal, alone)) < 0.05          # pteropod cannot shelter the foram


def test_ceara_rise_profile_keeps_the_foram_to_the_saturation_horizon(runs):
    """The measured Ceara Rise column (GLODAPv2.2023): calcite-saturated to about 4.3 km, so the foram keeps
    its mass there and loses only a little by 4.6 km; the temperature follows the profile."""
    run = runs["ceara"]
    depth, ug, sat = _shell(run, "calcite")
    horizon = depth[next(i for i, x in enumerate(sat) if x < 0)]
    assert 4000 < horizon < 4500
    assert all(u > 19.999 for z, u in zip(depth, ug) if z < horizon - 100)
    assert 19.5 < ug[-1] < 20.0
    card = example_card("ceara")
    assert card.tp(0.0)[0] == pytest.approx(card.profile_at(987.0)[0], abs=0.05)


def test_water_profile_interpolates_and_holds_its_ends():
    from reaction.app.scenario import Scenario
    card = Scenario(water_profile=[[0.0, 20.0, 35.0, 2300.0, 2000.0], [1000.0, 4.0, 34.0, 2320.0, 2200.0]])
    assert card.profile_at(500.0) == pytest.approx((12.0, 34.5, 2310e-6, 2100e-6))
    assert card.profile_at(3000.0) == pytest.approx((4.0, 34.0, 2320e-6, 2200e-6))
    assert Scenario(water_profile=[[1000.0, 4, 34, 2320, 2200], [0.0, 20, 35, 2300, 2000]]).problems()
    assert Scenario(water_profile=[[0.0, 20, 35, 2300]]).problems()


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


def test_vent_mixing_keeps_disequilibrium_and_gives_energy(runs):
    run = runs["vent_energy"]
    e = run["energy"]
    assert not _formed(run)                                                  # conservative mixing
    # the vent fluid's sulfide is diluted, not oxidized: about 1 mmol/kg at 75 C (a fifth vent fluid)
    sulfide = run["species"]["H2S"][55] + run["species"]["HS-"][55]
    assert 0.8e-3 < sulfide < 1.6e-3
    peak = {k[2:]: max(v or 0 for v in vals) for k, vals in e.items() if k.startswith("E:")}
    assert max(peak, key=peak.get) == "h2s_ox" and 30 < peak["h2s_ox"] < 60
    hot = {k[2:]: vals[2] for k, vals in e.items() if k.startswith("E:")}     # 340 C: anaerobic ahead
    assert max(hot, key=hot.get) in ("sulfate_red", "aom")
    cold = {k[2:]: vals[-1] for k, vals in e.items() if k.startswith("E:")}   # 5 C: sulfide oxidation
    assert max(cold, key=cold.get) == "h2s_ox"
    assert all(v is None or v < -300 for v in e["dG:h2s_ox"][1:])            # far from equilibrium
    es = e["Es:h2s_ox"][-1] / e["E:h2s_ox"][-1]
    assert es == pytest.approx(116, rel=1e-3)                                # per kg of vent fluid


def test_pyrite_oxidation_makes_acid_mine_drainage(runs):
    run = runs["amd"]
    assert not _formed(run)                                          # Fe(III) solids blocked: iron stays dissolved
    assert 1.9 < run["properties"]["pH"][-1] < 2.1
    assert run["properties"]["ionic_strength"][-1] < 0.1             # inside the B-dot range
    # 10 mmol pyrite fully oxidized and hydrolysed gives 40 meq of acid (2000 mg/kg as CaCO3), less the
    # groundwater's alkalinity; net acidity counts it although only part is free H+
    acidity = run["properties"]["net_acidity"]
    assert acidity[0] < 0 and 1800 < acidity[-1] < 2000
    assert all(a < b for a, b in zip(acidity, acidity[1:]))
    assert run["species"]["Fe+2"][-1] > run["species"]["Fe+3"][-1]   # O2 enough for the sulfur, not the iron

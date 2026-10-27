"""Organic carbon: only the organic acids an analysis gives, each held as its own component with its acid
and metal complexes, and never made from carbonate (network.organic_held_apart, ORGANIC_LIGANDS).
Needs the thermodynamic data."""
import pytest

from reaction.app import config
from reaction.app.engine import Engine
from reaction.app.scenario import Scenario
from reaction.network import ORGANIC_LIGANDS, organic_held_apart


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


def _organic(name):
    return any(f"({abbr})" in name for abbr, _ in ORGANIC_LIGANDS.values()) or name in ORGANIC_LIGANDS or any(
        name in acids for _, acids in ORGANIC_LIGANDS.values())


def test_reduced_vent_fluid_keeps_its_carbon_as_carbonate(engine):
    # the low fO2 of the vent fluid once put most of its carbonate into Na(Pent) and other complexes
    sp, _ = engine.starting_fluid(Scenario(fluid="vent_fluid"))
    organic = [n for n, st in sp.species.items() if _organic(n) and st.molality > 1e-15]
    assert organic == []
    held = organic_held_apart(engine.db)
    assert "Na(Pent)" in held and "CH4" not in held


def test_acetate_in_an_analysis_makes_its_complexes(engine):
    sc = Scenario(fluid="custom", analysis="Na+ 500\nCa+2 10\nCl- 500\nHCO3- 2\nacetate 0.1\nformate 0.05",
                  units="mmolal", ph=6.0, charge_balance="Cl-")
    sp, _ = engine.starting_fluid(sc)
    acetate = sum(st.molality for n, st in sp.species.items()
                  if n in ("acetate", "acetic acid") or "(Ac)" in n)
    assert acetate == pytest.approx(1e-4, rel=1e-4)        # Ca(Ac)2 counted once
    assert any("(Ac)" in n and st.molality > 0 for n, st in sp.species.items())
    assert not any("(Pent)" in n or "(Prop)" in n for n in sp.species)


def test_vent_fluid_heats_through_every_point(engine):
    sc = Scenario(fluid="vent_fluid", t_start_c=25, t_end_c=350, p_start_bar=250, p_end_bar=250, steps=6)
    r = engine.run(sc)
    ph = r["modes"]["static"]["properties"]["pH"]
    assert all(v is not None for v in ph)
    assert not any("No equilibrium" in n for n in r["notes"])

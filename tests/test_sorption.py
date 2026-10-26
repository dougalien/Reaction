"""Ion exchange and surface complexation (sorption.py) against PHREEQC 3.7.3 (llnl.dat, with phreeqc.dat's
surface species), seawater at 25 C, no minerals."""
import pytest

from reaction import sorption
from reaction.app.engine import Engine
from reaction.app.scenario import Scenario
from reaction.equilibrium import equilibrate


@pytest.fixture(scope="module")
def engine():
    return Engine()


def _seawater(engine, **card):
    sc, problems = Scenario().changed({"temperature_c": 25, "pressure_bar": 1, **card})
    assert not problems
    sp0, _, _, _, system, bulk0 = engine._prepare(sc, lambda f, m: None)
    return equilibrate(system, 25.0, bulk0, pressure_bar=1.0, eligible=[], guess=sp0)


def test_reactions_parse():
    rx = {r.species: r for r in sorption.reactions()}
    assert len(rx) == 56                     # phreeqc.dat, plus NiX2 and CoX2 (Bradbury and Baeyens 2009)
    assert rx["CaX2"].z == 2 and rx["CaX2"].reactants == {"Ca+2": 1.0} and rx["CaX2"].master == "X-"
    assert rx["Hfo_wOCa+"].reactants == {"Ca+2": 1.0, "H+": -1.0} and rx["Hfo_wOCa+"].z == 1
    assert rx["Hfo_wSO4-"].reactants == {"SO4-2": 1.0, "H+": 1.0, "H2O": -1.0} and rx["Hfo_wSO4-"].z == -1
    # van 't Hoff for the exchange constants that have an enthalpy
    assert rx["CaX2"].log_k_at(25.0) == pytest.approx(0.8)
    assert rx["CaX2"].log_k_at(75.0) > 0.8 and rx["NaX"].log_k_at(75.0) == 0.0


def test_exchanger_against_phreeqc(engine):
    st = _seawater(engine, exchange_meq=100)
    s, cec = st.sorbed, 0.1
    frac = {"Na": s["NaX"] / cec, "K": s["KX"] / cec, "Mg": 2 * s["MgX2"] / cec, "Ca": 2 * s["CaX2"] / cec}
    assert sum(frac.values()) == pytest.approx(1.0, abs=1e-6)          # the sites hold their capacity
    phreeqc = {"Na": 0.5780, "K": 0.0651, "Mg": 0.2621, "Ca": 0.0948}
    for ion, f in frac.items():
        assert f == pytest.approx(phreeqc[ion], rel=0.06), ion
    # loaded from the fluid without changing it
    assert st.speciation.ph == pytest.approx(8.22, abs=0.005)


def test_ferric_oxide_surface_against_phreeqc(engine):
    st = _seawater(engine, hfo_g=1.0)
    s = st.sorbed
    assert st.speciation.ph == pytest.approx(7.218, abs=0.03)
    assert st.surface_potential_mv == pytest.approx(36.3, abs=2.5)
    assert s["Hfo_sOHCa+2"] == pytest.approx(5.1e-5, rel=0.05)
    assert s["Hfo_wOMg+"] == pytest.approx(7.27e-4, rel=0.05)
    assert s["Hfo_wCO3-"] + s["Hfo_wHCO3"] == pytest.approx(6.07e-4, rel=0.05)
    assert s["Hfo_wSO4-"] + s["Hfo_wOHSO4-2"] == pytest.approx(1.49e-4, rel=0.15)
    weak = sum(v for k, v in s.items() if k.startswith("Hfo_w"))
    assert weak == pytest.approx(1.0 / 89 * 0.2, rel=1e-6)                 # every weak site is somewhere


def test_path_with_sorption(engine):
    sc, problems = Scenario().changed({"t_start_c": 25, "t_end_c": 100, "p_start_bar": 1, "p_end_bar": 5, "steps": 4,
                                       "exchange_meq": 100, "hfo_g": 1.0})
    assert not problems
    r = engine.run(sc)
    run = r["modes"]["static"]
    held = [run["sorbed"]["NaX"][i] + run["sorbed"]["KX"][i] + 2 * (run["sorbed"]["CaX2"][i] + run["sorbed"]["MgX2"][i])
            for i in range(4)]
    assert held == pytest.approx([0.1] * 4, rel=1e-4)
    assert all(v is not None for v in run["properties"]["surface_potential_mV"])
    assert r["catalog"]["sorbed"]["CaX2"]["label"] == "Ca on exchange sites (CaX2)"
    assert any("Dzombak" in n for n in r["notes"])


def test_card_and_chat():
    from reaction.app.assistant import interpret_rules
    sc, _ = Scenario().changed({"exchange_meq": -1})
    assert sc.problems()
    assert interpret_rules("heat seawater to 100 C with 50 meq of exchange sites", Scenario()).scenario_changes[
        "exchange_meq"] == 50.0
    got = interpret_rules("react with 10 g calcite and 1 g ferrihydrite", Scenario()).scenario_changes
    assert got["hfo_g"] == 1.0 and got["reactants"] == "10 g calcite"
    assert interpret_rules("no exchanger", Scenario()).scenario_changes == {"exchange_meq": 0.0}

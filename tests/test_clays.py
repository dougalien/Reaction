"""Clays as exchangers, partition coefficients (Kd) and the contaminated-stream card.  The exchange check is
against PHREEQC 3.7.3 (phreeqc.dat, EXCHANGE -equilibrate) for the stream's major ions and metals at 25 C, pH 6.5."""
import math

import pytest

from reaction import sorption
from reaction.app import config
from reaction.app.engine import Engine, sorption_for
from reaction.app.scenario import Scenario, starting_card
from reaction.equilibrium import equilibrate

METALS = "Ca 15\nMg 4.1\nNa 6.3\nK 2.3\nHCO3 58.4\nSO4 11.2\nCl 7.8\nZn 2\nPb 0.2\nCu 0.3\nCd 0.02"


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


def _card(**changes):
    sc, problems = Scenario().changed({"fluid": "custom", "units": "mg/l", "charge_balance": "none",
                                       "temperature_c": 25, "pressure_bar": 1, **changes})
    assert not problems
    return sc


def _solve(engine, sc, t=25.0):
    sp0, _, _, _, system, bulk0 = engine._prepare(sc, lambda f, m: None)
    return sp0, equilibrate(system, t, bulk0, pressure_bar=1.0, eligible=[], guess=sp0)


def test_tables():
    clays = sorption.clays()
    assert {c: clays[c].cec_eq_per_kg for c in clays} == {"kaolinite": 0.020, "illite": 0.225, "montmorillonite": 0.764}
    assert sorption.clay_name("Smectite") == "montmorillonite" and sorption.clay_name("kaolin") == "kaolinite"
    rows = {(r.element, r.oxidation, r.medium): r for r in sorption.kd_table()}
    assert len(rows) == 27                                   # 9 metals x soil, sediment, suspended
    assert rows[("Pb", None, "suspended")].log_kd == 5.7 and rows[("Cr", 6, "soil")].log_kd == 1.1
    assert rows[("Ni", None, "sediment")].stat == "mean"
    assert {r.element for r in sorption.kd_rows("sediment", skip={"Pb", "Zn"})} == {"As", "Cd", "Cr", "Cu", "Hg", "Ni"}
    with pytest.raises(ValueError):
        sorption.kd_rows("river")


def test_clay_exchange_against_phreeqc(engine):
    """100 g illite = 22.5 meq per kg of water, loaded from the water without changing it."""
    _, st = _solve(engine, _card(analysis=METALS, ph=6.5, redox="none", illite_g=100))
    s, cec = st.sorbed, 0.0225
    ours = {"Ca": 2 * s["CaX2"], "Mg": 2 * s["MgX2"], "Na": s["NaX"], "K": s["KX"], "Zn": 2 * s["ZnX2"],
            "Pb": 2 * s["PbX2"], "Cu": 2 * s["CuX2"], "Cd": 2 * s["CdX2"]}
    phreeqc = {"Ca": 0.7246, "Mg": 0.2053, "Na": 0.004843, "K": 0.005211, "Zn": 0.05524, "Pb": 0.001294,
               "Cu": 0.003169, "Cd": 0.0003333}
    for ion in ("Ca", "Mg", "Na", "K"):
        assert ours[ion] / cec == pytest.approx(phreeqc[ion], rel=0.01), ion
    # the metals differ by the aqueous complexes (OBIGT against phreeqc.dat: Cu and Pb carbonate, hydroxide)
    for ion in ("Zn", "Pb", "Cu", "Cd"):
        assert ours[ion] / cec == pytest.approx(phreeqc[ion], rel=0.3), ion
    assert sum(ours.values()) == pytest.approx(cec, rel=1e-6)
    assert st.speciation.ph == pytest.approx(6.5, abs=1e-4)


def test_partition_coefficients_hold_every_species(engine):
    """Kd M times the dissolved metal (all its species); Cr is chromate under air, so Cr(VI)'s Kd applies.  On its
    own a metal is held exactly at its Kd."""
    sc = _card(analysis="Ca 15\nHCO3 58.4\nCl 7.8\nNi 0.1\nCr 0.05\nAs 0.05\nHg 0.002", ph=7.0, redox="air",
               charge_balance="auto", kd_solid_g=2000, kd_medium="sediment")
    sp0, st = _solve(engine, sc)
    assert st.speciation.ph == pytest.approx(sp0.ph, abs=1e-4)              # loaded without changing the water
    rows = {r.label: 10 ** r.log_kd for r in sorption_for(sc).kd}
    dissolved = {}
    for name, s in st.speciation.species.items():
        sp = engine.db.species.get(name)
        for el, k in (sp.elements.items() if sp else ()):
            dissolved[el] = dissolved.get(el, 0.0) + k * s.molality
    for el, label in (("Ni", "Ni"), ("As", "As"), ("Hg", "Hg"), ("Cr", "Cr(VI)")):
        # complexes of two such metals (NiH2AsO4+, NiCrO4) take the larger Kd, which lifts As by about 5%
        assert st.sorbed[f"{el} (Kd)"] / 2.0 / dissolved[el] == pytest.approx(rows[label], rel=0.08), el
    assert "H2AsO4-" in st.speciation.species and "HAsO4-2" in st.speciation.species   # arsenate ions (SLOP98)
    sc = _card(analysis="Ca 15\nHCO3 58.4\nCl 7.8\nNi 0.1", ph=7.0, redox="none", charge_balance="auto",
               kd_solid_g=2000)
    _, st = _solve(engine, sc)
    ni = sum(s.molality * engine.db.species[n].elements.get("Ni", 0) for n, s in st.speciation.species.items()
             if n in engine.db.species)
    assert st.sorbed["Ni (Kd)"] / 2.0 / ni == pytest.approx(10 ** 3.9, rel=1e-6)


def test_kd_left_to_metals_the_clays_cannot_exchange():
    sc, _ = Scenario().changed({"kd_solid_g": 1000, "illite_g": 10})
    assert {r.element for r in sorption_for(sc).kd} == {"As", "Cr", "Hg", "Ni"}
    assert sorption_for(sc).exchange_eq == pytest.approx(0.00225)
    sc, _ = Scenario().changed({"kd_solid_g": 1000})
    assert {r.element for r in sorption_for(sc).kd} == {"As", "Cd", "Cr", "Cu", "Hg", "Ni", "Pb", "Zn"}


def test_stream_card_shortened(engine):
    sc, problems = starting_card("stream").changed({"pore_volumes": 3, "cells": 5, "column_m": 5.0})
    assert not problems and not sc.problems()
    r = engine.run(sc)
    out, prof = r["modes"]["outlet"], r["modes"]["profile"]
    part = out["partition"]
    # one exchanger: each clay's Kd is in the ratio of its exchange capacity
    assert part["Zn|montmorillonite"][-1] / part["Zn|kaolinite"][-1] == pytest.approx(764 / 20, rel=1e-3)
    assert part["Ni|sediment"][-1] == pytest.approx(10 ** 3.9, rel=1e-3)
    assert r["catalog"]["partition"]["Pb|illite"]["label"] == "Pb on illite"
    assert r["catalog"]["partition"]["Ni|sediment"]["label"] == "Ni on bed sediment (Kd)"
    # Cu comes through before Zn; the Kd metals stay in the first cell
    inflow = r["fluid"]["totals_molal"]
    assert out["elements"]["Cu"][-1] / inflow["Cu²⁺"] > 10 * out["elements"]["Zn"][-1] / inflow["Zn²⁺"]
    assert prof["elements"]["Ni"][0] > 1e3 * prof["elements"]["Ni"][-1]
    assert any("Allison and Allison" in n for n in r["notes"])
    assert any("26.3 meq" in n for n in r["notes"])


def test_chat():
    from reaction.app.assistant import interpret_rules
    got = interpret_rules("contaminated stream through a 10 m column with 100 g kaolinite, 40 g illite and 20 g "
                          "montmorillonite", Scenario()).scenario_changes
    assert got["kaolinite_g"] == 100 and got["illite_g"] == 40 and got["montmorillonite_g"] == 20
    assert "reactants" not in got
    assert interpret_rules("react seawater with 10 g kaolinite at 25 C", Scenario()).scenario_changes[
        "reactants"] == "10 g kaolinite"                     # kaolinite alone is still a mineral to react
    got = interpret_rules("add 4 kg of bed sediment with Kd", Scenario()).scenario_changes
    assert got == {"kd_solid_g": 4000.0, "kd_medium": "sediment"}
    assert interpret_rules("no clays", Scenario()).scenario_changes["illite_g"] == 0.0
    from reaction.app.session import _starting_card_named
    assert _starting_card_named("load the contaminated stream card") == "stream"
    sc, _ = Scenario().changed({"illite_g": 40, "kd_solid_g": 4000})
    assert "40 g illite" in sc.describe() and "Kd" in sc.describe()
    assert not math.isnan(sc.clay_meq) and sc.clay_meq == pytest.approx(9.0)

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


def test_kd_left_to_metals_the_clays_cannot_hold():
    sc, _ = Scenario().changed({"kd_solid_g": 1000, "illite_g": 10})
    assert {r.element for r in sorption_for(sc).kd} == {"As", "Cr"}
    assert sorption_for(sc).edges == (("illite", 0.01),)
    sc, _ = Scenario().changed({"kd_solid_g": 1000, "illite_g": 10, "clay_edges": False})
    assert {r.element for r in sorption_for(sc).kd} == {"As", "Cr", "Hg"}
    assert not sorption_for(sc).edges
    assert sorption_for(sc).exchange_eq == pytest.approx(0.00225)
    sc, _ = Scenario().changed({"kd_solid_g": 1000})
    assert {r.element for r in sorption_for(sc).kd} == {"As", "Cd", "Cr", "Cu", "Hg", "Ni", "Pb", "Zn"}


def test_stream_card_shortened(engine):
    sc, problems = starting_card("stream").changed({"pore_volumes": 3, "cells": 5, "column_m": 5.0})
    assert not problems and not sc.problems()
    r = engine.run(sc)
    out, prof = r["modes"]["outlet"], r["modes"]["profile"]
    part = out["partition"]
    # kaolinite has exchange sites only; illite and montmorillonite add their edge sites
    exchange_only = part["Zn|kaolinite"][-1] * 764 / 20
    assert part["Pb|montmorillonite"][-1] > 3 * part["Pb|kaolinite"][-1] * 764 / 20
    assert part["Zn|montmorillonite"][-1] > exchange_only
    assert part["Hg|illite"][-1] > 10 and part["As|sediment"][-1] == pytest.approx(10 ** 2.2, rel=1e-3)
    assert "Ni|sediment" not in part                     # Ni is held by the clays now, not by its Kd
    assert r["catalog"]["partition"]["Pb|illite"]["label"] == "Pb on illite"
    assert r["catalog"]["partition"]["Cr|sediment"]["label"] == "Cr on bed sediment (Kd)"
    inflow = r["fluid"]["totals_molal"]
    pb_edges = out["elements"]["Pb"][-1] / inflow["Pb²⁺"]
    assert prof["elements"]["Pb"][0] > prof["elements"]["Pb"][-1]
    assert any("Allison and Allison" in n for n in r["notes"])
    assert any("26.3 meq" in n for n in r["notes"]) and any("2SPNE" in n for n in r["notes"])
    # without the edge sites lead comes through sooner
    r = engine.run(sc.changed({"clay_edges": False})[0])
    assert r["modes"]["outlet"]["elements"]["Pb"][-1] / inflow["Pb²⁺"] > pb_edges
    assert any("edge sites are off" in n for n in r["notes"])


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
    assert interpret_rules("turn off the clay edge sites", Scenario()).scenario_changes["clay_edges"] is False
    assert interpret_rules("with edge sites on the clays", Scenario()).scenario_changes["clay_edges"] is True
    from reaction.app.session import _starting_card_named
    assert _starting_card_named("load the contaminated stream card") == "stream"
    sc, _ = Scenario().changed({"illite_g": 40, "kd_solid_g": 4000})
    assert "40 g illite" in sc.describe() and "Kd" in sc.describe()
    assert not math.isnan(sc.clay_meq) and sc.clay_meq == pytest.approx(9.0)


def test_edge_constants():
    """Bradbury and Baeyens: fitted Ni and Co on illite (2009a Table 11), the LFER for the rest."""
    ill = {(ion, y): (k, how) for ion, y, k, how in sorption.edge_metal_constants("illite")}
    assert ill[("Ni+2", 0)] == (0.7, "measured") and ill[("Co+2", 2)] == (-16.5, "measured")
    assert ill[("Zn+2", 0)] == (pytest.approx(7.9 + 0.83 * -8.68, abs=0.005), "LFER")
    mnt = {(ion, y): (k, how) for ion, y, k, how in sorption.edge_metal_constants("montmorillonite")}
    assert mnt[("Ni+2", 1)] == (pytest.approx(8.1 + 0.90 * -18.0, abs=0.005), "LFER")
    assert "kaolinite" not in sorption.edge_sites()
    assert sorption.edge_elements() == {"Ni", "Co", "Zn", "Cd", "Pb", "Cu", "Mn", "Fe", "Hg"}
    names = [r.species for r in sorption.edge_reactions("illite")]
    assert {"Ill_sOH", "Ill_sOH2+", "Ill_w2O-", "Ill_sONi+", "Ill_sONiOH", "Ill_sONi(OH)2-"} <= set(names)


def _rd(engine, clay, ion, ph, g_per_l, edges=True):
    """log Rd (L/kg) of a trace metal on a Na-clay in 0.1 M NaCl at a fixed pH."""
    from reaction.equilibrium import System
    el = ion.rstrip("+2")
    kg = g_per_l / 1000.0
    sorb = sorption.Sorption(exchange_eq=sorption.clays()[clay].cec_eq_per_kg * kg,
                             edges=((clay, kg),) if edges else ())
    system = System(engine.db, ["Na+", "Cl-", ion], sorption=sorb)
    st = equilibrate(system, 25.0, {"H2O": 55.508, "Na+": 0.1, "Cl-": 0.1, ion: 1e-9}, pressure_bar=1.0,
                     fixed={"H+": -ph}, eligible=[])
    held = sum(v for k, v in st.sorbed.items() if el in k)
    dissolved = sum(s.molality * engine.db.species[n].elements.get(el, 0) for n, s in st.speciation.species.items()
                    if n in engine.db.species)
    return math.log10(held / kg / dissolved)


def test_edges_against_measured_sorption_edges(engine):
    """Measured log Rd read from Bradbury and Baeyens: Ni (Fig. 3b, 2 g/L) and Co (Fig. 4, 1.6 g/L) on Na-illite
    in 0.1 M NaClO4 (2009a), and Ni on Na-montmorillonite SWy-1 in 0.1 M NaClO4 at 1.42 g/L (2005, Table 5).
    Illite uses the fitted constants (within 0.4 log units); montmorillonite the LFER (within 0.9)."""
    for ion, g, edge in (("Ni+2", 2.0, ((3.5, 1.8), (6.0, 3.1), (7.0, 3.6), (8.0, 3.8), (9.0, 4.2))),
                         ("Co+2", 1.6, ((3.0, 2.15), (6.0, 2.75), (7.0, 3.2), (8.3, 4.4), (9.3, 5.2)))):
        for ph, measured in edge:
            assert _rd(engine, "illite", ion, ph, g) == pytest.approx(measured, abs=0.4), (ion, ph)
    for ph, measured in ((6.1, 2.4), (7.1, 3.3), (8.1, 3.9)):
        assert _rd(engine, "montmorillonite", "Ni+2", ph, 1.42) == pytest.approx(measured, abs=0.9), ph
    # without the edge sites the clay holds Ni by exchange alone, flat with pH
    assert _rd(engine, "illite", "Ni+2", 8.0, 2.0, edges=False) < _rd(engine, "illite", "Ni+2", 8.0, 2.0) - 1.5


def test_edges_start_in_equilibrium_with_the_fluid(engine):
    """The edge sites are loaded from the starting fluid, so the fluid keeps its pH and metals at the start."""
    sc = _card(analysis=METALS, ph=6.8, redox="none", charge_balance="auto", illite_g=40, montmorillonite_g=20)
    sp0, st = _solve(engine, sc)
    assert st.speciation.ph == pytest.approx(sp0.ph, abs=1e-4)
    assert st.speciation.species["Zn+2"].molality == pytest.approx(sp0.species["Zn+2"].molality, rel=1e-4)
    assert st.sorbed["Mnt_sOPb+"] > 0 and st.sorbed["Ill_w1OH"] > 0
    strong = sum(v for k, v in st.sorbed.items() if k.startswith("Ill_s"))
    assert strong == pytest.approx(0.002 * 0.040, rel=1e-6)          # 2 mmol/kg of strong sites, 40 g of illite

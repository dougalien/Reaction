"""Catabolic energies and redox couples held apart (catabolic.py, equilibrium.redox_intermediates and
the family tie in System).  Needs the thermodynamic data."""
import pytest

from reaction import catabolic
from reaction.app import config
from reaction.app.engine import Engine
from reaction.app.scenario import Scenario
from reaction.database.model import Kind
from reaction.equilibrium import REDOX_APART, System, equilibrate, redox_intermediates
from reaction.path import bulk_over


@pytest.fixture(scope="module")
def engine():
    if not (config.external_dir() / "OBIGT" / "inorganic_aq.csv").exists():
        pytest.skip("thermodynamic data not found; set REACTION_DATA_DIR")
    return Engine()


def test_reactions_balance_and_h2_oxidation_matches_water(engine):
    db = engine.db
    for r in catabolic.REACTIONS:
        catabolic.log_k(db, r.reaction, 25.0, 1.0)                           # raises if unbalanced
    # 2 H2 + O2 = 2 H2O at 25 C: dG0 = 2 dGf(H2O) - 2 dGf(H2,aq) - dGf(O2,aq) = -2*237.14 - 2*17.7 - 16.4
    lk = catabolic.log_k(db, catabolic.REACTIONS[0].reaction, 25.0, 1.0)
    dg0 = -8.314462618e-3 * 298.15 * 2.302585093 * lk
    assert dg0 == pytest.approx(-525.9, abs=2.0)


def test_held_families_and_intermediates(engine):
    db = engine.db
    inter = set(redox_intermediates(db, REDOX_APART))
    assert {"S2O3-2", "Na(Pent)", "NO2-", "pyrite", "magnetite"} <= inter
    assert not {"HS-", "H2S", "SO4-2", "CH4", "O2", "H2", "Fe+3", "FeOH+2", "hematite", "pyrrhotite"} & inter
    system = System(db, ["SO4-2", "HCO3-", "Fe+2", "Cl-", "Na+", "O2(g)"], hold_apart=("H2S", "H2", "CH4", "Fe+3"),
                    suppress=inter)
    net = system.net
    row = dict(zip(net.basis, net.nu_dep[net.dep_names.index("HS-")]))
    assert row["H2S"] == 1 and row["H+"] == -1 and row["SO4-2"] == 0 and row["O2(g)"] == 0
    row = dict(zip(net.basis, net.nu_dep[net.dep_names.index("FeOH+2")]))
    assert row["Fe+3"] == 1 and row["Fe+2"] == 0 and row["O2(g)"] == 0


def test_mixture_keeps_both_sides_of_every_couple(engine):
    """Reduced vent fluid with oxygenated seawater, held apart: totals mix linearly, nothing reacts."""
    db = engine.db
    inter = redox_intermediates(db, REDOX_APART)
    vent, _ = engine.starting_fluid(Scenario(fluid="vent_fluid"), inter)
    sea, _ = engine.starting_fluid(Scenario(fluid="seawater"), inter)
    basis = [b for b in dict.fromkeys(vent.basis_names + sea.basis_names) if b not in ("H2O", "H+")]
    hold = [h for h, p in REDOX_APART.items() if p is None or p in basis]
    system = System(db, [b for b in basis if b not in hold], hold_apart=hold, suppress=inter)
    bv, bs = bulk_over(system, vent), bulk_over(system, sea)
    assert bv["H2S"] == pytest.approx(6.5e-3, rel=0.05) and bs["O2(g)"] == pytest.approx(2.5e-4, rel=0.1)
    bulk = {k: max(0.5 * bv[k] + 0.5 * bs[k], 1e-12) if k not in ("H2O", "H+") else 0.5 * bv[k] + 0.5 * bs[k]
            for k in system.basis}
    st = equilibrate(system, 100.0, bulk, pressure_bar=250.0, eligible=[])
    tot = st.speciation.component_totals
    for k in ("H2S", "H2", "O2(g)", "SO4-2", "HCO3-", "Fe+2"):
        assert tot[k] == pytest.approx(bulk[k] / st.water_kg, rel=1e-6)
    assert st.speciation.species["O2"].molality > 1e-4 and st.speciation.species["H2"].molality > 1e-4
    en = catabolic.energies(db, st.speciation)
    assert en["dG:h2_ox"] < -400 and en["E:h2_ox"] > 0


def test_at_equilibrium_the_energies_vanish(engine):
    """The vent fluid on its own, its redox in equilibrium: every reaction among dissolved species has a ΔG of
    zero (iron oxidation to hematite is left out: no mineral formed)."""
    db = engine.db
    sp, _ = engine.starting_fluid(Scenario(fluid="vent_fluid"))
    en = catabolic.energies(db, sp)
    for k, v in en.items():
        if k.startswith("dG:") and v is not None and k != "dG:fe_ox":
            assert abs(v) < 1.0, k


def test_chat_words_set_the_switches():
    from reaction.app.assistant import _disequilibrium_options as read
    assert read("keep redox apart, no minerals, and compute the metabolic energy available") == \
        {"redox_apart": True, "minerals_form": False, "energy": True}
    assert read("turn off the catabolic energy") == {"energy": False}
    assert read("block calcite") == {}

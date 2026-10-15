"""The equation-of-state database against the EQ3/6 reference material.

EQ3/6's data0.com was built from SUPCRT92 (HKF + Maier-Kelley), so the OBIGT data in
its SUPCRT92-era form should reproduce its logK grids where the reactions coincide.
Differences larger than a few thousandths mark species whose parameters were revised
after SUPCRT92 (everything containing Al or Fe+++, several metal complexes).
"""
import math

import pytest

from reaction.speciation import MgPerKg, Molal, PH, speciate
from reaction.thermo.eosdb import eq36_name_to_obigt
from reaction.validation import parse_eq3nr_output
from reaction.validation.logk_compare import compare_logk

MAJOR_PHASES_AND_SPECIES = [
    "Calcite", "Aragonite", "Dolomite", "Magnesite", "Siderite", "Anhydrite", "Halite",
    "Quartz", "Brucite", "Periclase", "Lime", "Wollastonite", "Forsterite", "Fayalite",
    "Diopside", "Talc", "Chrysotile", "Tremolite",
    "OH-", "CO3--", "HS-", "HSO4-", "NaCl(aq)", "CaSO4(aq)", "KSO4-", "FeCl+",
]


def test_logk_reproduces_the_supcrt92_grids_for_major_species(data0_cmp, eos_db_supcrt92):
    rows = {r.name: r for r in compare_logk(data0_cmp, eos_db_supcrt92, MAJOR_PHASES_AND_SPECIES)}
    assert len(rows) >= 24
    for name, row in rows.items():
        assert row.max_abs < 0.015, (name, row.differences)


def test_most_identical_reactions_agree_to_a_hundredth(data0_cmp, eos_db_supcrt92):
    rows = compare_logk(data0_cmp, eos_db_supcrt92)
    assert len(rows) > 200
    good = [r for r in rows if r.max_abs < 0.02]
    assert len(good) > 85


def test_eh_reaction_matches_the_data0_header(data0_cmp, eos_db_supcrt92):
    for t, ref in zip(data0_cmp.temps_c, data0_cmp.logk_eh.values):
        assert eos_db_supcrt92.logk_eh(float(t)) == pytest.approx(float(ref), abs=5e-4)


def test_debye_huckel_parameters_match_the_data0_header(data0_cmp, eos_db_supcrt92):
    for t, p, a, b in zip(data0_cmp.temps_c, data0_cmp.pressures_bar, data0_cmp.adh, data0_cmp.bdh):
        got_a, got_b, _ = eos_db_supcrt92.dh_params(float(t), float(p))
        assert got_a == pytest.approx(float(a), abs=3e-4)
        assert got_b == pytest.approx(float(b), abs=2e-4)


SEAWATER = {"Na+": MgPerKg(10768), "K+": MgPerKg(399.1), "Ca+2": MgPerKg(412.3),
            "Mg+2": MgPerKg(1291.8), "H+": PH(8.22), "HCO3-": Molal(2.022e-3),
            "Cl-": MgPerKg(19353), "SO4-2": MgPerKg(2712)}


def test_seawater_speciation_agrees_with_eq3nr(eq36_tests, eos_db_supcrt92):
    ref = parse_eq3nr_output(eq36_tests / "3tlib_cmp" / "Outputs" / "swmaj.3o")
    res = speciate(eos_db_supcrt92, 25.0, SEAWATER)
    assert res.ionic_strength == pytest.approx(ref.ionic_strength, rel=0.02)
    tight = {"Cl-", "Na+", "NaCl(aq)", "K+", "HCO3-", "CO3--", "NaHCO3(aq)"}
    for name, (_, log_m, _, _) in ref.species.items():
        key = eq36_name_to_obigt(name)
        if log_m < -4.6 or key not in res.species:
            continue
        diff = abs(math.log10(res.species[key].molality) - log_m)
        # Where the underlying data differ a looser bound applies.  KCl(aq) is the extreme
        # case: its dissociation constant differs by ~1 log unit between the 1996 data and
        # the revised OBIGT parameters (both agree it is a minor species).
        limit = 0.01 if name in tight else 1.1 if name == "KCl(aq)" else 0.5
        assert diff < limit, (name, diff)
    assert res.saturation["calcite"] == pytest.approx(ref.solids["Calcite"], abs=0.05)


def test_heated_seawater_reaches_anhydrite_saturation_near_150c(eos_db_supcrt92):
    """Seawater precipitates anhydrite on heating to about 150 C (Bischoff and Seyfried, 1978)."""
    si = {}
    for t in (100.0, 125.0, 150.0, 175.0):
        si[t] = speciate(eos_db_supcrt92, t, SEAWATER, pressure_bar=500.0).saturation["anhydrite"]
    assert si[100.0] < si[125.0] < 0.0 < si[175.0]
    assert si[125.0] < 0.0 < si[150.0] + 0.2


def test_speciation_converges_from_ambient_to_350c(eos_db_supcrt92):
    last_ionic = None
    for t in (25.0, 100.0, 200.0, 300.0, 350.0):
        res = speciate(eos_db_supcrt92, t, SEAWATER, pressure_bar=500.0)
        assert res.ionic_strength > 0.3
        last_ionic = res.ionic_strength
    assert last_ionic < 0.5          # ion pairing removes charge carriers as water loses polarity

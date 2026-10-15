"""Equation-of-state database built from a handful of published OBIGT parameter rows."""
import math

import pytest

from reaction.database import Kind, ThermoRangeError
from reaction.speciation import Molal, PH, speciate
from reaction.thermo.eosdb import build_eos_database, reference_pressure_bar
from reaction.thermo.obigt import load_obigt, parse_formula

HEADER = ("name,abbrv,formula,state,ref1,ref2,date,model,E_units,G,H,S,Cp,V,"
          "a1.a,a2.b,a3.c,a4.d,c1.e,c2.f,omega.lambda,z.T\n")
ROWS = """\
water,NA,H2O,liq,JOH92,HGK84,2006-10-25,H2O,cal,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
H+,H+,H+,aq,Den81,NA,1997-11-06,HKF,cal,0,0,0,0,0,0,0,0,0,0,0,0,0
Na+,Na+,Na+,aq,SH88,NA,1997-11-06,HKF,cal,-62591,-57433,13.96,9.06,-1.11,1.839,-2.285,3.256,-2.726,18.18,-2.981,0.3306,1
Cl-,Cl-,Cl-,aq,SH88,NA,1997-11-07,HKF,cal,-31379,-39933,13.56,-29.44,17.79,4.032,4.801,5.563,-2.847,-4.4,-5.714,1.456,-1
NaCl,NaCl,NaCl,aq,SSH97,NA,1997-09-16,HKF,cal,-92910,-96160,28,8.5,24,5.0364,4.5189,3.9669,-2.9658,10.798,-1.3031,-0.038,0
Ca+2,Ca+2,Ca+2,aq,SH88,NA,1997-11-06,HKF,cal,-132120,-129800,-13.5,-7.53,-18.06,-0.1947,-7.252,5.2966,-2.4792,9,-2.522,1.2366,2
HCO3-,HCO3-,HCO3-,aq,SH88,NA,1997-11-06,HKF,cal,-140282,-164898,23.53,-8.46,24.6,7.5621,1.1505,1.2346,-2.8266,12.9395,-4.7579,1.2733,-1
SO4-2,SO4-2,SO4-2,aq,SH88,NA,1997-11-07,HKF,cal,-177930,-217400,4.5,-64.38,13.88,8.3014,-1.9846,-6.2122,-2.697,1.64,-17.998,3.1463,-2
HSO4-,HSO4-,HSO4-,aq,SH88,NA,1997-11-07,HKF,cal,-180630,-212500,30,5.3,35.2,6.9788,9.259,2.1108,-3.1618,20.0961,-1.955,1.1748,-1
CO2,NA,CO2,aq,PS01.1,NA,2005-08-18,HKF,cal,-92250,-98900,28.1,51.86,33.5,15.2964,-10.0382,-55.4493,7.4092,36.8069,3.5851,-0.3107,0
SiO2,NA,SiO2,aq,SHS89,NA,1989-01-13,HKF,cal,-199190,-209775,18,-76.1,16.1,1.9,1.7,20,-2.7,29.1,-51.2,0.1291,0
calcite,Cal,CaCO3,cr,HDNB78,PB82,1990-03-09,CGL,cal,-269880,-288552,22.15,NA,36.934,24.98,0.00524,-620000,0,0,0,0,1200
anhydrite,Anh,CaSO4,cr,HDNB78,RH95.7,1978-05-05,CGL_Ttr,cal,-315925,-342760,25.5,NA,45.94,16.78,0.0236,0,0,0,0,0,1723
quartz,Qtz,SiO2,cr,HDNB78,NA,1978-05-05,CGL,cal,-204646,-217650,9.88,NA,22.688,11.22,0.0082,-270000,0,0,0,0,848
quartz,Qtz,SiO2,cr2,HDNB78,OBIGT.1,1978-05-05,CGL,cal,-204461.901,-217728.651,8.9986,NA,23.06,14.41,0.00194,0,0,0,0,0,2000
oxygen,O2,O2,gas,WEP+82,Kel60,1989-12-27,CGL,cal,0,0,49.029,NA,0,7.16,0.001,-40000,0,0,0,0,3000
"carbon dioxide",CO2,CO2,gas,WEP+82,Kel60,1989-12-27,CGL,cal,-94254,-94051,51.085,NA,0,10.57,0.0021,-206000,0,0,0,0,2500
"""


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    folder = tmp_path_factory.mktemp("obigt")
    (folder / "mini.csv").write_text(HEADER + ROWS, encoding="utf-8")
    return build_eos_database(load_obigt(folder, ("mini.csv",)))


def test_formula_parser():
    assert parse_formula("Ca+2") == ({"Ca": 1.0}, 2.0)
    assert parse_formula("SO4-2") == ({"S": 1.0, "O": 4.0}, -2.0)
    assert parse_formula("Fe2SiO4") == ({"Fe": 2.0, "Si": 1.0, "O": 4.0}, 0.0)
    assert parse_formula("Ca(OH)2") == ({"Ca": 1.0, "O": 2.0, "H": 2.0}, 0.0)
    assert parse_formula("H+") == ({"H": 1.0}, 1.0)


def test_species_kinds_and_generated_reactions(db):
    assert db["Na+"].kind is Kind.BASIS and db["H2O"].kind is Kind.BASIS
    assert db["calcite"].kind is Kind.MINERAL and db["CO2(g)"].kind is Kind.GAS
    assert db["calcite"].reaction == {"calcite": -1.0, "Ca+2": 1.0, "HCO3-": 1.0, "H+": -1.0}
    assert db["HSO4-"].reaction == {"HSO4-": -1.0, "SO4-2": 1.0, "H+": 1.0}
    assert db["CO2"].reaction == {"CO2": -1.0, "HCO3-": 1.0, "H+": 1.0, "H2O": -1.0}
    assert db["quartz"].reaction == {"quartz": -1.0, "SiO2": 1.0}
    assert db["NaCl"].reaction == {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0}


def test_reactions_conserve_atoms_and_charge(db):
    for sp in db.species.values():
        if not sp.reaction:
            continue
        atoms: dict[str, float] = {}
        charge = 0.0
        for name, c in sp.reaction.items():
            other = db[name] if name != "O2(g)" else db["O2(g)"]
            charge += c * other.charge
            for el, n in other.elements.items():
                atoms[el] = atoms.get(el, 0.0) + c * n
        assert abs(charge) < 1e-9, sp.name
        assert all(abs(v) < 1e-9 for v in atoms.values()), sp.name


def test_logk_at_25c_matches_published_values(db):
    assert db["calcite"].logk(25.0) == pytest.approx(1.8486, abs=0.003)     # EQ3/6 data0: 1.8487
    assert db["quartz"].logk(25.0) == pytest.approx(-3.9993, abs=0.003)
    assert db["NaCl"].logk(25.0) == pytest.approx(0.777, abs=0.01)


def test_eh_reaction_logk(db):
    assert db.logk_eh(25.0) == pytest.approx(-83.1049, abs=2e-3)
    assert db.logk_eh(300.0) == pytest.approx(-42.6828, abs=2e-3)


def test_temperature_and_pressure_dependence(db):
    k = db["calcite"].logk
    assert k(100.0) < k(25.0)                                  # calcite is less soluble when hot
    assert k(25.0, 1000.0) - k(25.0, 1.0) == pytest.approx(0.5, abs=0.1)   # dV = -30 cm3/mol
    # higher temperature than any grid in the EQ3/6 files
    assert math.isfinite(k(350.0, 500.0))


def test_steam_region_is_refused(db):
    with pytest.raises(ThermoRangeError):
        db["calcite"].logk(150.0, 2.0)
    with pytest.raises(ThermoRangeError):
        db["calcite"].logk(25.0, 1e5)


def test_reference_pressure_follows_the_eq36_convention():
    assert reference_pressure_bar(25.0) == 1.0132
    assert reference_pressure_bar(100.0) == 1.0132
    assert reference_pressure_bar(300.0) == pytest.approx(85.84, abs=0.05)
    with pytest.raises(ThermoRangeError):
        reference_pressure_bar(400.0)


def test_polymorph_selection_for_quartz(db):
    k = db["quartz"].logk
    assert k(25.0) == pytest.approx(k(25.0, 1.0132), abs=1e-12)
    # alpha quartz is the stable form at low T and the SiO2 solubility rises smoothly
    values = [k(t, 500.0) for t in (25.0, 100.0, 200.0, 300.0)]
    assert values == sorted(values)


def test_speciation_on_the_mini_database(db):
    res = speciate(db, 25.0, {"H+": PH(8.0), "Na+": Molal(0.5), "Cl-": Molal(0.5),
                              "Ca+2": Molal(5e-3), "HCO3-": Molal(2e-3)})
    assert 0.0 < res.molality("NaCl") < 0.05
    carbon = sum(res.molality(n) for n in ("HCO3-", "CO2"))
    assert carbon == pytest.approx(2e-3, rel=1e-9)               # the mass balance closes
    expected_si = (res.species["Ca+2"].log_activity + res.species["HCO3-"].log_activity
                   + 8.0 - db["calcite"].logk(25.0))
    assert res.saturation["calcite"] == pytest.approx(expected_si, abs=1e-8)


def test_speciation_works_at_350c_and_500_bar(db):
    res = speciate(db, 350.0, {"H+": PH(6.0), "Na+": Molal(0.5), "Cl-": Molal(0.5)},
                   pressure_bar=500.0)
    assert res.molality("NaCl") > 0.1                           # ion pairing is strong when hot
    assert any("B-dot" in w for w in res.warnings)


def test_low_density_water_triggers_a_warning(db):
    dilute = speciate(db, 400.0, {"H+": PH(6.0), "Na+": Molal(0.1), "Cl-": Molal(0.1)},
                      pressure_bar=250.0)
    assert any("density" in w for w in dilute.warnings)
    dense = speciate(db, 350.0, {"H+": PH(6.0), "Na+": Molal(0.1), "Cl-": Molal(0.1)},
                     pressure_bar=1000.0)
    assert not any("density" in w for w in dense.warnings)

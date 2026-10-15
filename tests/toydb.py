"""A tiny hand-built database for tests that must not depend on external files.

Constants are standard textbook 25 C values (log Kw = -13.995, carbonic acid
pK1 = 6.352, pK2 = 10.329, calcite via HCO3-).  Illustrative only.
"""
import numpy as np

from reaction.database.model import GridLogK, Kind, Provenance, Species, ThermoDatabase

TEMPS = np.array([0.01, 25.0, 60.0, 100.0, 150.0, 200.0, 250.0, 300.0])
ELEMENTS = {"H": 1.00794, "O": 15.9994, "C": 12.011, "Na": 22.98977, "Cl": 35.4527, "Ca": 40.078}


def _sp(name, kind, charge, elements, reaction=None, logk=None, a0=4.0, atype=0):
    return Species(
        name=name, kind=kind, formula=None, charge=charge, elements=elements, mol_wt=None,
        a0=a0, activity_type=atype, reaction=reaction or {},
        logk=None if logk is None else GridLogK(TEMPS, np.full(8, logk)),
        molar_volume=None, ref_state={}, provenance=Provenance("toy", "toydb.py", "tests"),
    )


def make_toy_db() -> ThermoDatabase:
    db = ThermoDatabase(
        title="toy", source="toydb.py", temps_c=TEMPS,
        pressures_bar=np.full(8, 1.0132),
        adh=np.full(8, 0.5114), bdh=np.full(8, 0.3288), bdot=np.full(8, 0.0410),
        cco2=np.array([-1.0312, 0.0012806, 255.9, 0.4445, -0.001606]),
        logk_eh=GridLogK(TEMPS, np.full(8, -83.1049)), elements=ELEMENTS,
    )
    B, A, Q, M = Kind.BASIS, Kind.AUX_BASIS, Kind.AQUEOUS, Kind.MINERAL
    for s in [
        _sp("H2O", B, 0, {"H": 2, "O": 1}, a0=3.0),
        _sp("H+", B, 1, {"H": 1}, a0=9.0),
        _sp("HCO3-", B, -1, {"C": 1, "H": 1, "O": 3}),
        _sp("Na+", B, 1, {"Na": 1}),
        _sp("Cl-", B, -1, {"Cl": 1}, a0=3.0),
        _sp("Ca++", B, 2, {"Ca": 1}, a0=6.0),
        _sp("OH-", Q, -1, {"O": 1, "H": 1}, {"OH-": -1, "H+": -1, "H2O": 1}, 13.995, a0=3.5),
        _sp("CO3--", Q, -2, {"C": 1, "O": 3}, {"CO3--": -1, "H+": -1, "HCO3-": 1}, 10.329),
        _sp("CO2(aq)", Q, 0, {"C": 1, "O": 2}, {"CO2(aq)": -1, "H2O": -1, "H+": 1, "HCO3-": 1},
            -6.352, a0=3.0, atype=-1),
        _sp("Calcite", M, 0, {"Ca": 1, "C": 1, "O": 3},
            {"Calcite": -1, "H+": -1, "Ca++": 1, "HCO3-": 1}, 1.85),
    ]:
        db.species[s.name] = s
    return db

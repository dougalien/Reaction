# Reaction

Open geochemical reaction-path modelling in Python: aqueous speciation,
equilibrium with minerals and gases, reaction paths and kinetics.  Early
stage: speciation (B-dot or Davies activity models, 0-300 C on the database's
grid) works and reproduces all 41 speciation cases of the EQ3/6 test library;
reaction paths, high-temperature/pressure physics and kinetics are next
(see [docs/DESIGN.md](docs/DESIGN.md)).

## Try it

Needs Python 3.10+, numpy and scipy.  The thermodynamic database is not part
of this repository; point the importer at your own copy of an EQ3/6 `data0`
file (`scripts/fetch_references.ps1` downloads EQ3/6 v8.0a for you).

```python
from reaction.database import load_data0
from reaction.speciation import speciate, Molal, MgPerKg, PH, Eh

db = load_data0(r"C:\path\to\Eq3_6v8.0a\DB\data0.cmp")
seawater = speciate(
    db, 25.0,
    {"Na+": MgPerKg(10768), "K+": MgPerKg(399.1), "Ca++": MgPerKg(412.3),
     "Mg++": MgPerKg(1291.8), "Cl-": MgPerKg(19353), "SO4--": MgPerKg(2712),
     "HCO3-": Molal(2.022e-3), "H+": PH(8.22)},
    redox=Eh(0.5),
)
print(seawater.table())
print(seawater.saturation["Calcite"])       # log Q/K
print(seawater.log_fugacity["CO2(g)"])
```

Constraints can also be equilibria and a charge balance.  Rain water in
equilibrium with atmospheric CO2 and calcite (pH and Ca come out of the model):

```python
from reaction.speciation import Equilibrium

water = speciate(
    db, 25.0,
    {"H+": PH(7.0),                          # only a starting guess here
     "Ca++": Equilibrium("Calcite"),
     "HCO3-": Equilibrium("CO2(g)", -3.5)},  # log fCO2 = -3.5
    balance_on="H+",
)
print(water.ph, water.totals["Ca++"])        # 8.27, 4.97e-4 mol/kg
```

Other options: `MgPerL`, `LogActivity`, `Eh` / `Pe` / `LogFO2` / `Couple`,
`switches={"NH3(aq)": "NO3-"}` for basis switches, `suppress=[...]` and
`suppress_phases=[...]`, `activity_model="davies"`.  Auxiliary species given a
total (`"O2(aq)": MgPerL(6.6)`) become independent components, as in EQ3NR.

## Tests

```
python -m pytest                       # self-contained tests only
set REACTION_EQ36_DIR=C:\path\to\Eq3_6v8.0a
set REACTION_EQ36_TESTS=C:\path\to\folder-containing-3tlib_cmp
python -m pytest                       # adds the EQ3/6 regression tests
```

With the references present the suite runs all 41 EQ3/6 speciation cases (43
problems) and compares every species, solid saturation index and gas fugacity
with the stored EQ3NR output.

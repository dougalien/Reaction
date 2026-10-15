# Reaction

Open geochemical reaction-path modelling in Python: aqueous speciation,
equilibrium with minerals and gases, reaction paths and kinetics.  Early
stage: speciation (B-dot or Davies activity models) works to 400 C and 1 kbar
and reproduces all 41 speciation cases of the EQ3/6 test library; reaction
paths and kinetics are next (see [docs/DESIGN.md](docs/DESIGN.md)).

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

## High temperature and pressure

An EQ3/6 database tabulates logK from 0 to 300 C along the water saturation curve.
For other temperatures and pressures Reaction can compute logK from the species'
Gibbs energies instead: IAPWS-95 water, the Johnson-Norton dielectric constant,
the HKF equations for aqueous species (with the g function) and Maier-Kelley heat
capacities for minerals and gases.  The data are the OBIGT csv files of CHNOSZ.

```python
from reaction.thermo.eosdb import load_eos_database

eos = load_eos_database(r"C:\path\to\obigt", preset="supcrt92", eq36_database=db)
vent = speciate(
    eos, 350.0,                      # degC
    {"Na+": MgPerKg(10768), "K+": MgPerKg(399.1), "Ca+2": MgPerKg(412.3),
     "Mg+2": MgPerKg(1291.8), "Cl-": MgPerKg(19353), "SO4-2": MgPerKg(2712),
     "HCO3-": Molal(2.022e-3), "H+": PH(6.0)},
    pressure_bar=500.0,
)
print(vent.saturation["anhydrite"], vent.warnings)
```

Species names follow OBIGT here (`Ca+2`, `SO4-2`, `SiO2`, minerals in lower case, gases as
`CO2(g)`).  `preset="supcrt92"` lets the 1990s SUPCRT92 parameters win, as in EQ3/6's
`data0.com`; `"current"` prefers the maintained OBIGT files.  Passing an EQ3/6 database as
`eq36_database` only borrows its ion-size parameters.  The `warnings` list says when a
calculation leaves the calibrated range (B-dot above 300 C, water below 0.35 g/cm3).

## Tests

```
python -m pytest                       # self-contained tests only
set REACTION_EQ36_DIR=C:\path\to\Eq3_6v8.0a
set REACTION_EQ36_TESTS=C:\path\to\folder-containing-3tlib_cmp
set REACTION_OBIGT_DIR=C:\path\to\obigt
python -m pytest                       # adds the EQ3/6 and OBIGT comparison tests
```

With the references present the suite runs all 41 EQ3/6 speciation cases (43
problems) and compares every species, solid saturation index and gas fugacity
with the stored EQ3NR output.

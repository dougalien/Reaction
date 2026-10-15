# Reaction

Open geochemical reaction-path modelling in Python: aqueous speciation,
equilibrium with minerals and gases, reaction paths and kinetics.  Early
stage: speciation with the B-dot model works and is validated against EQ3/6;
reaction paths, high-temperature physics and kinetics are next
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

## Tests

```
python -m pytest                       # self-contained tests only
set REACTION_EQ36_DIR=C:\path\to\Eq3_6v8.0a
set REACTION_EQ36_TESTS=C:\path\to\folder-containing-3tlib_cmp
python -m pytest                       # adds the EQ3/6 regression tests
```

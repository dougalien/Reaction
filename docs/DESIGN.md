# Reaction: design notes

Reaction is an open geochemical reaction-path code: define a fluid, add
minerals/gases/fluids at a given water/rock ratio, follow temperature and
pressure, and get the equilibrium species distribution, mineralogy and
(later) kinetics.  It is a personal research tool, written in readable Python
with numpy/scipy.

The hard part of this problem is not the solver.  It is (1) the
thermodynamic data, (2) the high-temperature/high-pressure physics needed for
hydrothermal and CO2-storage conditions, and (3) showing that results agree
with codes people trust.  The design follows from that.

## Decisions (agreed)

* Own solver in Python, not a wrapper around Reaktoro.  PHREEQC, EQ3/6 and
  Reaktoro are used as cross-checks only.
* After the 25 C core: seawater-basalt, 25-350 C, is the first system to
  validate.  CO2-brine-rock (Pitzer, CO2 equation of state) comes later.
* First interface: input files + tables + matplotlib plots.  A local web UI
  (FastAPI) comes once the engine is trusted.

## Layout

```
src/reaction/
  constants.py            physical constants
  database/
    model.py              Species, Provenance, GridLogK, ThermoDatabase
    eq36.py               EQ3/6 data0 importer
  activity.py             B-dot, Davies, Bates-Guggenheim, Drummond, water activity
  speciation.py           Network (reactions reduced to a basis) + Newton solver
  validation/eq3nr.py     parser for EQ3NR .3o reference outputs
tests/                    toy-database + synthetic-data0 tests; EQ3/6 regression tests
```

Conventions: reactions use the EQ3/6 dissociation form (species coefficient
-1), so `log a_i = sum_k c_k log a_k - log K`.  Reactions that mention other
non-basis species are rewritten recursively until only basis species remain,
which is also what will make basis switching possible.

## Data policy

* No third-party database is ever committed.  Importers read from a path you
  supply (`REACTION_EQ36_DIR`, ...); `.gitignore` blocks `data0.*`.
* **EQ3/6 v8.0a** (github.com/LLNL/EQ3_6): the software is BSD-3.  The
  `Readme.txt` inside the Windows package (2010 text) still says commercial use
  needs notification to DOE/LLNL, and nothing found states the license of the
  `data0.*` database files.  Treat them as reference material for personal and
  academic use and do not redistribute them.
* **CHNOSZ / OBIGT** (github.com/jedick/CHNOSZ): GPL-3.0.  Provides HKF
  parameters (SUPCRT92, SLOP98, DEW, Berman, inorganic aq/cr/gas CSVs, about
  0.5 MB) for the high-T/P work in M2.  Anything derived and redistributed
  inherits GPL-3.
* **Geochemist's Workbench** databases are proprietary and are not used.
* Every imported species keeps its source database, file, section, revision
  date and literature source tags (`Species.provenance`).

## Validation

Oracle: the EQ3/6 v8.0a test library (Windows package, with executables).

| library | what | cases |
|---|---|---|
| `3tlib_cmp` | EQ3NR speciation, data0.cmp, B-dot | 41 inputs with reference outputs |
| `3tlib_hmw`, `3tlib_fmt`, `3tlib_ymp`, `3tlib_ypf` | Pitzer and other databases | later |
| `6tlib_*` | EQ6 reaction paths, kinetics, mixing | oracle for M3/M4 |

Status: `swmaj` (seawater major ions, 25 C) matches EQ3NR for every printed
species (log molality, log gamma, log activity to 2e-4), 35+ solid saturation
indices (5e-4) and gas fugacities.  Run with:

```
set REACTION_EQ36_DIR=...\Eq3_6v8.0a
set REACTION_EQ36_TESTS=...                (the folder that contains 3tlib_cmp)
python -m pytest
```

Known, understood differences from EQ3/6:
* Eh to pe: EQ3/6 uses a legacy Faraday constant (about 23062.3 cal/V/eq);
  Reaction uses CODATA, so pe differs by 6e-5 relative.  Regression tests fix
  log fO2 at EQ3NR's value to compare everything else exactly.
* mg/kg.sol conversion follows EQ3NR (default TDS = 0 means mg per kg water).

## Roadmap

* **M1 (done)**: data0 import; B-dot speciation with NBS pH scale; redox via
  Eh/pe/fO2; saturation indices and gas fugacities.
* **M1.5**: run all 41 `3tlib_cmp` cases automatically.  Needs a `.3i` input
  reader, charge balancing on a species, mg/L and alkalinity units, mineral/gas
  equilibrium constraints ("hetero. equil."), basis switching, species
  suppression, Davies model, and a check at T other than 25 C (3 cases at 250 C).
* **M2**: temperature and pressure.  Water properties (IAPWS), dielectric
  constant, Debye-Hueckel A and B at T,P, HKF for aqueous species and
  Maier-Kelley/Berman for minerals, so logK(T,P) reaches about 400 C and 1 kbar
  instead of the data0 saturation-curve grid (to 300 C).
* **M3**: reaction paths: titration, water/rock ratio, precipitation and
  dissolution under the phase rule, fixed-fugacity gas buffers, fluid mixing,
  temperature paths.  Oracle: `6tlib_*`.
* **M4**: kinetics: rate laws (affinity-based TST, pH-dependent), integrated in
  time and coupled to the equilibrium solver.
* **M5**: plots (Eh-pH, activity diagrams, SI versus reaction progress) and UI.
* **M6**: Pitzer activity model and a CO2 equation of state for brines and
  carbon-storage systems.

Out of scope for now: solid solutions, surface complexation, reactive transport.

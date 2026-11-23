# Reaction

Geochemical reaction-path modelling in Python, from [Dougalien Systems](https://dougalien.com).

Reaction computes aqueous speciation, equilibrium with minerals, gases and solid solutions, reaction paths
and mineral kinetics at any temperature from 0 to 600 °C and any pressure from 1 to 5000 bar. Equilibrium
constants are calculated from Gibbs energies at every point of a run, so there are no fixed-pressure
databases. Every number it reports comes from its own solver and from published data; the Sources page in
the app lists each data set and equation with its reference.

The [user manual](docs/MANUAL.md) explains every control of the app, the models behind it and their limits,
for students and instructors.

## What it does

- Speciation of any water, from a preset (seawater, deep seawater, rainwater, fresh water, groundwater, karst
  spring, vent fluid, formation brine, Great Salt Lake brine and others) or a pasted analysis with units and
  automatic charge balance. Redox from a couple, air, log fO2, Eh or pe.
- Reaction paths along a straight line in temperature and pressure, with reactants added steadily (minerals,
  species, gases, rocks and basalt glass; water/rock ratio and titration), a second fluid mixed in, or water
  evaporated. Static and flow-through modes.
- Activity models: B-dot (Helgeson 1969, the default), b-gamma (Helgeson, Kirkham and Flowers 1981), Davies,
  and Pitzer with PHREEQC's `pitzer.dat` for brines.
- A free gas phase with Peng-Robinson fugacity coefficients.
- Solid solutions (plagioclase, carbonates, chlorite, epidote, olivine, pyroxenes, Fe-Mg brucite and others).
- Mineral kinetics from published rate laws (Palandri and Kharaka 2004 and others), including calcite and
  aragonite dissolution in seawater (Naviaux et al. 2019; Dong et al. 2019).
- The marine carbonate system from the oceanographic constants, as in PyCO2SYS.
- Phase separation of saline fluids (NaCl-H2O after Driesner and Heinrich 2007).
- Stable isotopes (δ18O, δD, δ13C, δ34S) and 87Sr/86Sr carried through the fluid and the minerals.
- Ion exchange and surface complexation, and one-dimensional reactive transport through a column of cells.
- Eh-pH, log fO2-pH and activity diagrams at any temperature and pressure.
- Uncertainty bands: a run can be repeated with log K values varied by their estimated uncertainty.
- A reliability rating for every point: the engine's water ionization constant is compared with IAPWS R11-24
  (2024), which shows where the aqueous thermodynamics can be trusted.

Worked examples fill the scenario card with one click: basalt alteration, an evaporating lake, a
contaminated stream, boiling, supercritical seawater, flow-through weathering, serpentinization, calcite
dissolving, forams and pteropods sinking through the lysocline, acid mine drainage and a vent fluid mixing
with seawater.

## Running it

### On the web (Streamlit)

    pip install -r requirements.txt
    streamlit run streamlit_app.py

The app has an Input window (the scenario card) and a Results window (`/results`, opened from the "Open
results window" button so it can sit on a second screen). The two share one session through a private `?w=`
id in the address. Each visitor has their own card and runs.

### On your own computer

    pip install -e .
    reaction-app

This starts a local server and opens the control window in the browser; nothing leaves your machine. On
Windows, `scripts\reaction.cmd` does the same. Saved runs go to `workspace/` (override with
`REACTION_WORKSPACE`). `packaging/build_exe.ps1` builds a standalone Windows executable.

The control window can also take requests in words through a local language model
([Ollama](https://ollama.com), default model `qwen2.5:7b`, override with `REACTION_OLLAMA_MODEL`). The
model only fills in the card and the plot request; it never produces numbers. This chat is switched off by
default (see the comment in `src/reaction/app/static/control.html`).

## Thermodynamic data

The data Reaction computes with are bundled:

- the OBIGT database of CHNOSZ 2.3.0 (HKF aqueous species, Berman and Maier-Kelley minerals, gases;
  GPL-3.0) in `src/reaction/data/reference/OBIGT/`;
- PHREEQC's `pitzer.dat` (public domain, U.S. Geological Survey) in `src/reaction/data/`;
- the EQ3/6 ion-size parameters, as numbers only, in `src/reaction/data/ion_sizes_eq36.tsv`.

`src/reaction/data/reference/NOTICE.md` gives the source, version and license of each. The EQ3/6 `data0`
files are not redistributed; Reaction can read them from a path you supply. A data folder named by
`REACTION_DATA_DIR` takes precedence over the bundled copies.

## Using the library

Needs Python 3.10 or later, numpy and scipy.

```python
from reaction.app import config
from reaction.thermo.eosdb import load_eos_database
from reaction.speciation import speciate, Molal, MgPerKg, PH, Equilibrium

db = load_eos_database(config.external_dir() / "OBIGT", "current", eq36_database="table")

seawater = speciate(
    db, 25.0,
    {"Na+": MgPerKg(10768), "K+": MgPerKg(399.1), "Ca+2": MgPerKg(412.3),
     "Mg+2": MgPerKg(1291.8), "Cl-": MgPerKg(19353), "SO4-2": MgPerKg(2712),
     "HCO3-": Molal(2.022e-3), "H+": PH(8.22)},
)
print(seawater.table())
print(seawater.saturation["calcite"])        # log Q/K
print(seawater.log_fugacity["CO2(g)"])
```

Constraints can also be equilibria and a charge balance. Rain water in equilibrium with atmospheric CO2
and calcite, where pH and Ca come out of the model:

```python
water = speciate(
    db, 25.0,
    {"H+": PH(7.0),                            # only a starting guess here
     "Ca+2": Equilibrium("calcite"),
     "HCO3-": Equilibrium("CO2(g)", -3.5)},    # log fCO2 = -3.5
    balance_on="H+",
)
print(water.ph, water.totals["Ca+2"])          # 8.33, 5.65e-4 mol/kg
```

Pass `pressure_bar=` for any pressure; `speciate(db, 350.0, ..., pressure_bar=500.0)` works the same way.
The `warnings` list of the result says when a calculation leaves the calibrated range (B-dot above 300 °C,
water below 0.35 g/cm3). Other options: `MgPerL`, `LogActivity`, `Eh`, `Pe`, `LogFO2`, `Couple`, basis
`switches`, `suppress` and `suppress_phases`, and `activity_model`. Species names follow OBIGT (`Ca+2`,
`SO4-2`, `SiO2`, minerals in lower case, gases as `CO2(g)`). `preset="supcrt92"` lets the SUPCRT92
parameters take precedence, as in EQ3/6's `data0.com`; `"current"` uses the maintained OBIGT files.

[docs/DESIGN.md](docs/DESIGN.md) describes the solver, the high-temperature equations and the checks
against other codes.

## Tests

    pip install -e ".[dev]"
    python -m pytest

The suite runs on the bundled data. Comparisons with EQ3/6 v8.0a (the 41 EQ3NR speciation cases of its
test library and the EQ6 reaction-path cases) run when the EQ3/6 files are in the data folder
(`scripts/fetch_references.ps1` downloads them) and are skipped otherwise. Other tests compare with
PHREEQC 3.7.3, PyCO2SYS and published tables.

## How to cite

If you use Reaction in a study or a class, please cite it, and cite the data and methods it used
for your result (the app's Sources page lists them). GitHub's "Cite this repository" button gives
the reference from `CITATION.cff`.

## License and credits

Reaction is free software under the GNU General Public License, version 3 or later (`LICENSE`).
Anyone may use, copy, change and share it at no cost; copies and changed versions that are passed on
must stay under the same license. It comes with no warranty.

Reaction bundles the OBIGT data of CHNOSZ (Dick 2019, Frontiers in Earth Science 7, 180), which are
also GPL-3.0, and PHREEQC's pitzer.dat (U.S. Geological Survey, public domain). Every data set and
equation is cited on the app's Sources page (`src/reaction/app/static/sources.html`).

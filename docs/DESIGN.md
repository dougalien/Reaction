# Reaction: design

Reaction is a geochemical reaction-path code written in Python with numpy and scipy: define a fluid, add
minerals, gases, rocks or a second fluid, follow temperature and pressure, and get the species
distribution, the mineral assemblage, gases, kinetics and isotopes along the way.  PHREEQC, EQ3/6 and
PyCO2SYS serve as numerical cross-checks.

## Layout

```
src/reaction/
  constants.py            physical constants
  database/
    model.py              Species, Provenance, GridLogK, ThermoDatabase
    eq36.py               EQ3/6 data0 importer
  thermo/                 thermodynamics from equations of state
    water.py              IAPWS-95 water: density, Gibbs energy, saturation
    electrostatics.py     Johnson-Norton dielectric constant, Born functions, Debye-Hueckel A, B
    hkf.py                revised HKF aqueous species, Shock et al. g function
    crystal.py            Maier-Kelley minerals, gases and liquids, SUPCRT92 quartz volumes
    berman.py             Berman (1988) minerals and later additions
    gas.py                Peng-Robinson fugacity coefficients of a gas mixture
    bgamma.py             b-gamma of Helgeson, Kirkham and Flowers (1981)
    obigt.py              reader for the OBIGT csv files (CHNOSZ)
    eosdb.py              EosDatabase: formation reactions from formulas, logK(T, P)
    validity.py           water-ionization reliability check (IAPWS R11-24)
  activity.py             B-dot, Davies, Bates-Guggenheim, Drummond, water activity
  pitzer.py, phreeqc_pitzer.py   Pitzer model; PHREEQC pitzer.dat as one consistent model
  network.py              every species and phase written over a chosen basis
                          (basis switches, suppression, decoupled components)
  speciation.py           constraints + Newton solver (speciate)
  aqueous.py              the aqueous model at one temperature and pressure
  equilibrium.py          closed-system equilibrium with minerals, a gas phase and solid solutions
  solidsolution.py        solid solutions
  path.py                 reaction paths: reactants, T-P paths, flow-through, flush, mixing
  kinetics.py             rate laws
  transport.py            one-dimensional reactive transport
  sorption.py             ion exchange and surface complexation
  marine.py               seawater carbonate system (oceanographic constants)
  porewater.py            steady pore water of a carbonate sediment
  nacl_h2o.py, phasesep.py   NaCl-H2O phase relations and phase separation
  isotopes.py             stable isotopes and 87Sr/86Sr on a solved path
  catabolic.py            catabolic energy available to chemolithotrophs
  diagram.py              predominance and activity diagrams
  validation/             readers and runners for the EQ3NR and EQ6 test libraries
  app/                    the local two-window app (server, scenario card, results)
  web/                    the Streamlit front end
tests/                    unit tests and comparisons with EQ3/6, PHREEQC and PyCO2SYS
```

Conventions: reactions use the EQ3/6 dissociation form (species coefficient
-1), so `log a_i = sum_k c_k log a_k - log K`.  Reactions that mention other
non-basis species are rewritten recursively until only basis species remain.

## How the solver is organised

Unknowns are the ln activities `u_j` of the basis species.  Every unknown gets
exactly one equation, chosen by the constraint supplied for it:

| constraint | equation |
|---|---|
| `Molal`, `MgPerKg`, `MgPerL` | mass balance, written as `ln(total / T) = 0` |
| `PH`, `LogActivity`, `LogFO2` | fixed activity |
| `Equilibrium(phase, value)` | log Q/K of a mineral (or log f of a gas) = value |
| `balance_on="X"` | `ln(positive charge / negative charge) = 0`, replacing X's own equation |
| `Eh`, `Pe` | linear relation between `u(O2)` and `u(H+)` |
| `Couple("NO3-")` | redox from a couple whose total is also given |

Newton steps are shortened (never distorted) to at most one decade per
iteration for mass balances, with a sum-of-squares line search and a
Levenberg-Marquardt fallback.  A pre-Newton phase of successive substitution
on the mass balances (EQ3NR does the same) is what makes 40-component fluids
with strongly complexed pairs converge; plain Newton fails on them because the
Jacobian rows of two components in one complex are nearly identical.  The
residuals are logarithmic on purpose: linear ones saturate at 1 when a
polynuclear species (Al13O4(OH)24 7+) explodes from a bad starting guess.

## EQ3NR behaviours worth knowing (all verified against reference outputs)

* **Suppression follows the reaction text.**  Suppressing a species removes
  every species whose database reaction mentions it: suppress `NH3(aq)` and
  `NO2-`, `NH4+`, `HNO2(aq)` vanish; suppress the master organic species
  `Acetic_acid(aq)` and all organics vanish.  Suppressed *phases* are different:
  their saturation states are still reported, they just cannot be used as
  equilibrium constraints.
* **Basis switch** (`NH3(aq)` -> `NO3-`): the switched-out species is defined
  through its replacement's reaction, so species written through `NO3-` (such as
  `CaNO3+`) stay as they are while species written through `NH3(aq)` are
  dropped once it is suppressed.
* **Decoupling**: giving an auxiliary species (`O2(aq)`, `Fe+++`, `NO2-`) its own
  total makes it an independent component, no longer tied to the default redox
  state.  Species that mention it stay written through it.
* **Balancing on a species** replaces that species' constraint: with `H+` the
  input pH is only a starting guess and the final pH is a result.
* **NBS pH scale**: `log gamma(Cl-)` is set by the Bates-Guggenheim equation and
  every ion is shifted by `-z * (BG - own)` of that; neutral species are
  untouched.
* **Water activity** is quasi-consistent with the activity model (B-dot form
  derived in the EQ3/6 theory manual; the Davies form follows from Gibbs-Duhem
  with the 0.2 coefficient EQ3/6 uses).
* **logK between grid temperatures**: one exact polynomial through the valid
  points of the 0-100 C part of the grid, another through those of the
  100-300 C part (EQPT).  A lone 25 C value becomes a constant.  Outside the
  span of valid points the polynomial is still evaluated, as in EQ3/6;
  `Speciation.extrapolated` lists the species that rely on that.

## High temperature and pressure

`EosDatabase` is a `ThermoDatabase` whose logK comes from Gibbs energies:

    log K = -(sum_i nu_i G_i(T, P)) / (R T ln 10)

* **Water**: IAPWS-95 (56-term Helmholtz energy, own analytic derivatives), solved
  for density at T, P.  The 56 coefficients were checked against the published
  verification table (pressure to 1e-8 at ten state points, saturated-liquid density
  at 300 K to 5e-10).
* **Dielectric constant**: Johnson and Norton (1991), the form the HKF parameters
  were regressed with.  Debye-Hueckel A and B at the eight EQ3/6 grid points agree
  with the SUPCRT92 values in the data0 header to 2e-4.
* **Aqueous species**: revised HKF equations; the effective radius depends on T and
  P through the Shock et al. (1992) g function.  Only G is needed for logK.
* **Minerals, gases, liquids**: Berman (1988) minerals where OBIGT has them (as in CHNOSZ's
  defaults, each held only to its transition temperature); otherwise Maier-Kelley Cp, constant
  volume, polymorphs by minimum G, SUPCRT92 quartz/coesite volume functions.  Gases are 1-bar
  ideal-gas standard states; a free gas phase takes Peng-Robinson fugacity coefficients.
* **Reactions** are generated from formulas over a set of master species (one per
  element, plus H2O, H+ and O2(g)), so any species with data and masters is usable.
* **Eh reaction** 2 H2O = O2 + 4 H+ + 4 e- (e- has G = 0 at all T) reproduces the
  data0 header to 3e-4.
* Pressure defaults to the EQ3/6 curve (1.0132 bar to 100 C, saturation above).
  Steam is refused for aqueous species; up to 1% below saturation is taken as
  metastable liquid.  `Speciation.warnings` flags B-dot above 300 C (its parameters
  stop there; the value at 300 C is held) and water below 0.35 g/cm3 (HKF
  extrapolation).

### How well does it agree?

Computed against the SUPCRT92-derived logK grids of `data0.com` (263 species whose
reactions coincide in both databases, 0-300 C):

* 88 agree to 0.01; the rock-forming minerals and main complexes are within
  0.007: quartz, calcite, aragonite, dolomite, magnesite, siderite, anhydrite,
  halite, brucite, periclase, lime, wollastonite, forsterite, fayalite, diopside,
  talc, chrysotile, tremolite, OH-, CO3--, HS-, HSO4-, NaCl(aq), CaSO4(aq),
  KSO4-, FeCl+.
* Clusters with larger differences share one species whose parameters were revised
  after SUPCRT92 (all Al minerals and aluminate, Fe+++, Tb/rare-earth complexes,
  KCl(aq), CaCl+, MgCl+ ...).  These are data differences, not equation errors:
  the same equations give 0.002 for the species whose data did not change.
* Seawater at 25 C on the OBIGT database: Cl-, Na+, NaCl(aq), K+, HCO3-, CO3--
  within 0.003 log units of EQ3NR; complexes of Mg and Ca within 0.01-0.1;
  KCl(aq) 1.0 (revised dissociation constant).
* Seawater heated at 500 bar crosses anhydrite saturation near 145 C (SI -0.19 at
  125 C, +0.07 at 150 C), as observed.

Known limitations: the Born pressure derivative Q (from JN91 and IAPWS-95) is 11%
above the value SUPCRT92's tabulated volumes imply (6.64e-7 against 5.97e-7 /bar at
25 C), which shifts computed standard volumes of ions by up to 1 cm3/mol and logK
at 1 kbar by about 0.01; Y and X agree.  The Akinfiev-Diamond and DEW species of OBIGT
are not read.

### How far to trust it at a given temperature and pressure

The ionization constant of water is measured over the whole range (IAPWS R11-24, 2024) and depends on
the same dielectric constant, g function and Born coefficients as every charged species, so the
engine's own value, compared with IAPWS, is an independent check at any (T, P)
(`thermo/validity.py`).  The app rates each point good (within 0.15 log units and water density at
least 0.35 g/cm3), approximate (within 0.35) or unreliable, and shows the rating before a run, in the
narrative and as a strip under the plot.  Engine log Kw minus IAPWS R11-24, log units:

| T (C) | 500 bar | 1000 bar | 2000 bar | 5000 bar |
|---|---|---|---|---|
| 25 | +0.02 | +0.04 | +0.07 | +0.16 |
| 100 | -0.01 | -0.01 | -0.01 | 0.00 |
| 200 | -0.01 | -0.03 | -0.05 | -0.09 |
| 300 | -0.04 | -0.06 | -0.09 | -0.14 |
| 350 | -0.09 | -0.08 | -0.13 | -0.17 |
| 400 | +0.06 | -0.12 | -0.17 | -0.19 |
| 500 | -0.68 (low density) | -0.25 | -0.30 | -0.26 |
| 600 | -2.72 (low density) | -0.59 | -0.52 | -0.37 |

The model holds to within 0.15 up to about 350 C and 2 kbar and degrades near the critical point and
wherever water density falls below about 0.35 g/cm3.  `tests/test_pressure.py` holds these checks.

## Data policy

* The OBIGT csv files are bundled in `src/reaction/data/reference/` (GPL-3.0, with a NOTICE) and
  PHREEQC's `pitzer.dat` in `src/reaction/data/` (public domain), so the app carries the data it
  computes with.  The ion sizes read from EQ3/6 `data0.cmp` are bundled as numbers only in
  `src/reaction/data/ion_sizes_eq36.tsv`.  Other databases are read from a path you supply
  (`REACTION_EQ36_DIR`, ...); `.gitignore` blocks `data0.*`.
* **EQ3/6 v8.0a** (github.com/LLNL/EQ3_6): the software is BSD-3.  The
  `Readme.txt` inside the Windows package (2010 text) still says commercial use
  needs notification to DOE/LLNL, and nothing found states the license of the
  `data0.*` database files.  Treat them as reference material for personal and
  academic use and do not redistribute them.
* **CHNOSZ / OBIGT** (github.com/jedick/CHNOSZ): GPL-3.0.  Its csv files (SUPCRT92,
  SLOP98, inorganic and organic aq/cr/gas, Berman) are the HKF, Maier-Kelley and Berman
  parameters.  CHNOSZ's R sources (IAPWS95.R, water.R, hkf.R,
  cgl.R ...) were used as a reference for constants and cross-checks only; Reaction's
  code is written from the published equations.  Anything derived and redistributed
  inherits GPL-3.
* **Geochemist's Workbench** databases are proprietary and are not used.
* Every imported species keeps its source database, file, section, revision
  date and literature source tags (`Species.provenance`).

## Validation

Oracle: the EQ3/6 v8.0a test library (Windows package, with executables).

| library | what | status |
|---|---|---|
| `3tlib_cmp` | EQ3NR speciation, data0.cmp, B-dot/Davies | **all 41 inputs (43 problems) reproduced** |
| `6tlib_cmp` | EQ6 reaction paths | `micro`, `rwtitr` and `swtitr` match to printed precision; `microft` (flow-through) except pH within 1.5e-3 at one boundary; `dedolo` runs (`tests/test_eq6_cases.py`) |

For `3tlib_cmp`, 1,625 species (log molality, log gamma), 2,600 solid
saturation indices, gas fugacities and the bulk properties (ionic strength,
water activity, osmotic coefficient, pH) agree with the stored EQ3NR output to
5.0e-5, which is the rounding of the printed values.  Features the cases cover:
charge balancing on H+, Na+ and Cl-; mg/L and mg/kg.sol; equilibrium with
minerals and gases; a special basis switch; redox from Eh, log fO2 and a
couple; decoupled components; suppression; Davies; 9.5, 10, 25.5, 70 and
250 C; multi-problem files.  The ion-exchange cases are compared on their
aqueous part only (exchangers are computed after the aqueous speciation).

```
set REACTION_EQ36_DIR=...\Eq3_6v8.0a
set REACTION_EQ36_TESTS=...                (the folder that contains 3tlib_cmp)
set REACTION_OBIGT_DIR=...                 (the OBIGT csv files)
python -m pytest
```

Known, understood differences from EQ3/6:
* Eh to pe: EQ3/6 uses R = 1.98726 cal/mol/K and F = 23062.3 cal/V/eq;
  Reaction uses CODATA.  The validation harness converts Eh with EQ3/6's
  constants so the inputs match; Reaction itself keeps the modern values (pe
  differs by 6e-5 relative if you pass Eh directly).
* mg/kg.sol conversion follows EQ3NR (default TDS = 0 means mg per kg water).

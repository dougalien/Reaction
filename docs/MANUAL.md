# Reaction user manual

Reaction models what happens when water meets minerals, rocks, gases or another water: which species are
dissolved, which minerals form or dissolve, and how both change along a path in temperature, pressure,
time or amount of reactant. This manual explains each part of the app, what the models behind it are,
and where their limits lie. It is written for students and instructors; no programming is needed.

Contents

1. [The three windows](#1-the-three-windows)
2. [A first run](#2-a-first-run)
3. [The scenario card](#3-the-scenario-card)
4. [The Results window](#4-the-results-window)
5. [The models](#5-the-models)
6. [Limits](#6-limits)
7. [The worked examples](#7-the-worked-examples)
8. [Using Reaction in a class](#8-using-reaction-in-a-class)
9. [Running Reaction on your own computer](#9-running-reaction-on-your-own-computer)
10. [How to cite](#10-how-to-cite)

## 1. The three windows

**Input.** The scenario card: the fluid, the path, the reactants and the options. The Run button and the
list of earlier runs are at the bottom.

**Results.** Opened with "Open results window" at the top right of the Input window. It shows the run
selected in the Input window and follows each new run live, so it can sit on a second screen or a
projector while you change the card. The two windows share one private workspace, named by the `?w=`
part of the address. Anyone with that address sees the same card and runs, so do not share it unless you
mean to. A workspace that has not been used for 12 hours is deleted, with its runs; download what you
want to keep (section 4).

**Sources.** Linked in the footer of every page. It lists each data set and equation with its reference.
Every number Reaction reports is computed by its own solver from these data; nothing you enter is sent
anywhere else.

## 2. A first run

1. Under **Start from**, press **Seawater**. The card now heats standard seawater from 2 to 350 °C at
   500 bar in 36 steps.
2. Press **Run**. A progress bar shows the points as they are computed.
3. Open the Results window. **What we found** describes the run in sentences. The plot shows species
   along the path; choose **Mineral amounts** under "What to plot" to see anhydrite form as the seawater
   heats, and brucite near the end.
4. Back on the card, change **Temperature end** to 250 °C and run again. In the Results window choose the
   first run under **Compare with**; it is drawn dotted over the new one.

Any field can be changed before a run. A field the app cannot accept shows a red message under the card
saying why, and the run does not start until it is fixed. Yellow messages are advice: the run goes ahead.

## 3. The scenario card

### Start from, and worked examples

**Start from** gives five ready cards: Seawater (heating), Fresh water (heating), Brine (Great Salt Lake
brine evaporating with Pitzer activities), Hydrothermal fluid (hot fluid mixing into cold seawater) and
Contaminated stream (a column of streambed). **Reset** returns to the Seawater card.

**Or load a worked example** gives thirteen cards, each built to show one feature of the model. Loading
one opens a box that says what it shows and lists the fields that set it up, one sentence each. Run it
as it is, then change one of those fields and run again to see what that field does. Section 7 lists them.

**Title** is optional; it labels the run in the history and on the Results window.

### Fluid

**Fluid** picks a starting water. The presets, with where their numbers come from:

| Preset | What it is | Source of the analysis |
|---|---|---|
| Seawater | Standard seawater, major ions only | EQ3/6 test case swmaj; oxygen from air |
| Deep seawater | Bottom water near a ridge, 2 °C, with dissolved O2 | major ions as Seawater; O2, carbonate and pH are rounded typical deep Pacific values, not a sample |
| Rainwater | Dilute rain with atmospheric CO2 | illustrative values, not a sample |
| Vent fluid | Black-smoker end-member: acid, reduced, Mg- and sulfate-free | rounded after Von Damm et al. (1985), East Pacific Rise 21°N; the 0.5 mmol/kg H2 is illustrative |
| Formation brine | Na-Ca-Cl basin brine, about 160 g/kg dissolved solids | illustrative, not a sample |
| Fresh water | World-average river water with air and 420 ppm CO2 | Livingstone (1963) |
| Contaminated stream water | River water with sulfate and Zn, Cu, Pb, Cd, Ni, As, Cr, Hg | majors after Livingstone (1963); metals illustrative |
| Great Salt Lake brine | Na-Cl-SO4 brine near halite saturation | Drever (1997); silica left out |
| Groundwater | Shallow granitic groundwater, 10 °C | illustrative, not a sample |
| Karst spring water | Limestone spring water, 12 °C | illustrative, not a sample |
| Hydrothermal fluid | Seawater after the default heating run (350 °C, 500 bar), with its anhydrite and brucite left behind | computed by Reaction itself |

The caption under the menu gives each preset's description. Editing the analysis of a preset turns it
into a **Custom fluid** that starts from the preset.

**Analysis.** One solute per line: a name, an amount and, if it differs from the default unit, a unit.

```
Na 10768
Ca2+ 10.3 mmol/kg
SiO2 12
```

Names can be element symbols, ions in any common charge notation (`Ca`, `Ca2+`, `Ca++`, `Ca+2`) or
English names (`calcium`, `sulfate`, `silica`). An amount given by mass is weighed as the formula
written: `Si 140` is 140 mg of silicon, `SiO2 300` is 300 mg of silica. **Give carbonate as HCO3 (total
dissolved carbonate), not as alkalinity**: a line with "alkalinity" is refused with that message.

**Default unit.** mg/kg, mmol/kg and mol/kg are per kilogram of solution; "mmol/kg water" and "mol/kg
water" (molal) are per kilogram of water; mg/L and meq/L are per litre, read at 1 kg/L.

**pH, Measured at (°C), Measured at (bar).** The pH and the conditions the analysis and pH refer to. The
pressure can be left empty below 300 °C: it is then 1 bar, or the boiling curve above 100 °C.

**Charge balance.** *Auto* adjusts the main ion on the side that is short; *None* uses the analysis as
given; *On X* adjusts that ion. With *Auto* or an ion chosen, the pH you enter stays fixed.

**Redox state.** *Auto* takes it from a redox species in the analysis (dissolved O2 or H2, Fe(III),
sulfide), otherwise from air. *Not tracked* keeps only the basis oxidation states (Fe(II), sulfate).
*Air* is O2 at 0.21 bar. *log fO2*, *Eh* and *pe* take a value. A redox species named in the analysis can
also be chosen as the couple that sets the state.

The expander under the analysis, "As read: ...", shows how each line was read, its molality, and notes
(for example, which ion the charge balance changed).

**Activity model.** How activity coefficients are computed; see section 5. B-dot is the default.

**Evaporate to (times concentrated).** Water is removed until the fluid is this many times more
concentrated, in steps evenly spaced in the logarithm of the factor. 1 means no evaporation. Evaporation
cannot be combined with kinetics.

### Path

**Kind of run.** *One fluid along a path* or *Flow through a column* (one-dimensional reactive
transport, below).

**Steps** (2 to 200): the number of points along the path.

**Temperature start and end, Pressure start and end.** The path is a straight line between them. Equal
start and end hold that variable fixed. Equilibrium constants are computed from the thermodynamic data at
every point, so any temperature from 0 to 600 °C and any pressure from 1 to 5000 bar can be used. The card
refuses a path on which liquid water would boil (the pressure must stay above the boiling curve, checked
at 201 points along the path, not only at the steps). When a second fluid is mixed in at its own
temperature, the temperature follows the mixing and the end temperature is greyed out.

**Mode** (paths only):

- *Static*: minerals that form stay with the fluid and can dissolve again further along.
- *Flow-through*: minerals that form at a step are taken away, so they cannot redissolve. A gas that
  forms leaves too.
- *Flush*: the solids stay, and after each step fresh starting fluid replaces the fluid around them, so a
  solid always meets new water and never its own products. This is the flush model of Bethke's React
  program (The Geochemist's Workbench); Reaction uses it for shells sinking through the ocean. Flush mode
  cannot be used with exchange sites, ferric oxide or clays yet.
- *Both*: static and flow-through side by side, for comparison.

**Column fields** (Flow through a column): **Cells** (2 to 100), **Pore volumes** (how many column
volumes of fluid pass through, up to 100), **Column length** (m) and **Dispersivity** (m; D =
dispersivity × velocity, 0 gives sharp fronts). The fluid on the card flows in at the inlet; **Column
first filled with** sets the water already in the pores. Temperature and pressure run from the start
values at the inlet to the end values at the outlet. Minerals and sorbed matter stay in their cell, and
every cell is at equilibrium at each shift, except for minerals reacting at measured rates (PHREEQC's
TRANSPORT scheme). A column run is limited to 40,000 equilibrium calculations; the card says when a
column needs more and what to reduce. Evaporation, mixing in a second fluid, following the brine or
vapor, uncertainty bands, holding redox apart, suppressing minerals and catabolic energy are not
available in a column.

### Reactants and mixing

**Mix in a second fluid** (paths only). A preset fluid added steadily until there are **Up to** kg of it
per kg of the first fluid (at most 1000). With **Its temperature** given, the mixture's temperature is the
mass-weighted mean of the two (equal heat capacities) and replaces the path's end temperature. **Steps
even in** *Amount mixed* or *Temperature*; temperature spacing resolves the hot end of a mixture, which
even amounts would cross in one or two steps.

**Reactants.** Minerals, aqueous species, gases or rocks, per kg of the starting water. Type one in the
box and press Enter: `10 g albite`, `0.05 mol HCl`, `20 ug forams at 2 m2/g`. Units are g, mg, ug (µg),
kg, mol, mmol and umol. Up to 12 reactants. Each shows as a button; press it to remove it.

- Without kinetics, each reactant is added steadily and all of it has reacted by the last step.
- With kinetics (below), the amount is what is there at the start, and minerals with a measured rate
  dissolve or grow at that rate. Minerals without a rate are still added steadily.

**Rock buttons.** Put a rock on the card as the reactant, in place of what is there, at the grams per kg
of water set above the buttons. Basalt glass dissolves as a glass, all of its oxides at once in proportion
(congruently); it never precipitates. The others (crystalline basalt, granite, lherzolite, Mars
meteorite, sandstone and three limestones) go on the card as their minerals, so each mineral dissolves on
its own and can precipitate again. Their modes are rounded textbook values (the Mars meteorite follows
Shergotty, Stolper and McSween 1979; basalt glass is close to the mean ridge basalt of Gale et al. 2013).

**How reactants react.** *All of it, steadily* or *At measured rates*. At measured rates the run lasts
**Run for** (in the **Time unit**) and time becomes the horizontal axis. Rates come from Palandri and
Kharaka (2004) for minerals and Gislason and Oelkers (2003) for basalt glass. **Surface area (cm²/g)** is
the reactive surface per gram of each reactant; 100 cm²/g is roughly 0.1 mm grains. A reactant can carry
its own surface: `20 ug forams at 2 m2/g`.

With kinetics on, more fields appear:

- **Calcite and aragonite rates.** *Laboratory*: Palandri and Kharaka (2004), measured far from
  equilibrium. *Seawater, in situ*: the lower bound Naviaux et al. (2019, Marine Chemistry 215)
  recommend for particles sinking through the water column, measured near saturation on the BET surface.
  *Seawater, laboratory*: their filtered-seawater fits (Naviaux et al. 2019, Geochimica et Cosmochimica
  Acta 246), about 4 times faster. *Seawater, aragonite after Dong et al.*: Dong et al.'s (2019) own
  aragonite law. With the seawater choices, aragonite takes the calcite law on its own saturation state
  (except with Dong's law), and the saturation state comes from the marine carbonate system (section 5),
  not from the thermodynamic database.
- **Plot against**: *Time* or *Depth (sinking shell)*. Depth reads the pressure as depth in seawater,
  about 10 m per bar.
- **Duration from**: *The run time* or *Sinking speed (Takahashi and Bé)*. A sinking particle falls from
  the start to the end pressure at the speed its mass gives, 10^2.06 W^0.64 m/day for a test of W µg
  (Takahashi and Bé 1984), recomputed as it dissolves. Sinking needs kinetics, a path, and a pressure
  that rises from start to end.
- **Setting**: *Water or rock*, or *Buried in sediment*. A buried test rests on the seafloor (the start
  pressure read as water depth) and is buried at the **Sedimentation rate** to the **Burial depth**, so
  the run lasts burial depth divided by rate. The pore pressure rises with burial and the temperature is
  the bottom water's. Burial needs kinetics and a path, and is not a sinking run.
- **Pore water** (burial only): *Profile or bottom water*, or *From the sediment*. *From the sediment*
  computes the steady pore water of a carbonate sediment after Hales and Emerson (1997): its calcite
  dissolving, CO2 from oxic respiration, both diffusing to the bottom water. It needs the seawater
  carbonate rates and seven site values with no defaults: **Porosity**, **CaCO₃ (wt %)**, **Calcite k
  (% per day)**, **Order n**, **O₂ flux**, **CO₂ per O₂** and **Respiration depth**. The help text of each
  field (the ⓘ) gives the values Hales and Emerson found at the Ceará Rise. **Mixed layer** and
  **Biodiffusion**, given together, add bioturbation: tests are mixed through the top layer by burrowers
  while they are buried. Without a measured profile or the sediment model, the test meets bottom water all
  the way down, which overstates dissolution below the seafloor.
- With the seawater rates: **Shell surface exponent** p, in A = A0 (mass left / starting mass)^p. 0 (the
  default) keeps a shell's whole surface as it dissolves, so its surface per gram rises as its walls thin;
  0.667 treats it as a shrinking solid grain. Shells (forams, pteropods, coccoliths) default to a BET
  surface of 2 m²/g with these rates, a placeholder of the right order from the pteropod values in Dong
  et al. (2019, Table 2); give your own with `at ... m2/g`. **Salinity**, **Alkalinity** and **DIC**
  set the marine carbonate system; left empty, they come from the water on the card.

**Blocked minerals.** Minerals that may not precipitate, because they form too slowly on the time scale
of the model even when the water is supersaturated. Dolomite, magnesite, huntite and goergeyite are
blocked by default. Press a name to allow it again; type a name and press Enter to block another.
Several examples block more (silica polymorphs in a low-temperature water, for instance) and say why.

### Surfaces and sorption

- **Exchange sites** (meq per kg of water): cation exchange sites in equilibrium with the fluid at the
  start (Gaines-Thomas convention, PHREEQC's constants).
- **Ferric oxide surface** (g per kg of water): hydrous ferric oxide with the surface complexation model
  of Dzombak and Morel (1990).
- **Kaolinite, Illite, Montmorillonite** (g per kg of water): each clay carries exchange sites at its
  cation exchange capacity (kaolinite 20, illite 225, montmorillonite 764 meq per kg of clay).
- **Clay edge sites**: *Bind metals (2SPNE)* lets illite and montmorillonite also bind metals on their
  edges (Bradbury and Baeyens 2005, 2009); *Exchange sites only* switches that off.
- **Kd solid** (g per kg of water) and **Kd values for** (bed sediment, suspended matter or soil): a solid
  that holds metals by the partition coefficients of Allison and Allison (2005, U.S. EPA). With clays or
  exchange sites on the card it holds only the metals they cannot.

### Options

| Option | Choices | What it does |
|---|---|---|
| Solid solutions | Pure minerals only / Can form | Lets ideal mixtures (plagioclase, carbonates, chlorite, epidote, olivine, pyroxenes, Fe-Mg brucite, serpentine and others) form in place of their endmembers. On by default |
| Follow | Whole fluid / Brine / Vapor | Where a saline fluid splits into brine and vapor (NaCl-H2O), which part the run follows; the other part leaves at each step |
| Isotopes | Off / On | Carries δ18O, δD, δ13C, δ34S and 87Sr/86Sr through the fluid and minerals |
| Gas phase | Can form / Suppressed | A free gas forms when the dissolved gases' fugacities add up to the pressure |
| Methane and sulfate | Held apart / Can react | Whether methane and sulfate may reach equilibrium with each other (thermochemical sulfate reduction). Held apart by default |
| Redox couples | Reach equilibrium / Held apart | Held apart keeps H2, sulfide, methane, Fe(III) and nitrate from reacting with O2, sulfate, carbonate, Fe(II) and ammonia, so a mixture keeps its redox disequilibrium |
| Minerals | Can form / Suppressed | Suppressed keeps every mineral from forming (conservative mixing) |
| Catabolic energy | Off / On | Gibbs energy of chemosynthetic reactions and the energy they make available per kg of fluid (McCollom and Shock 1997). Meaningful with redox couples held apart |
| Uncertainty range | 0 to 200 extra runs | Repeats the run with every log K shifted at random by its estimated uncertainty and shades the middle 90% of the results. At least 5 runs; 20 or more for a stable range. Each takes as long as the main run |

### Run and earlier runs

**Run** starts the calculation; a progress bar shows the points. **Earlier runs** lists every run in the
workspace with its title, summary and time; press one to show it in the Results window.

## 4. The Results window

**Run shown** and **Compare with.** Pick any run of the workspace, and a second one to draw dotted over
it on the same axes.

**What we found.** The run described in sentences: what formed or dissolved, where, and what the fluid
did. Every number in it is from the run.

**The plot.** Choose a group under **What to plot**, then pick one or more items from the list (type to
search). One axis holds one kind of quantity, so picking from another group replaces the lines. The
groups that the run produced are offered:

| Group | Plots |
|---|---|
| Species | molality of each aqueous species |
| Elements | total dissolved amount of each element |
| Mineral amounts | each mineral per kg of the original water, drawn from zero |
| Gas phase | amount of each gas in the gas phase |
| Rates | rate of each kinetic reactant (positive dissolves) and the % of it left |
| Energy | ΔG of each chemosynthetic reaction and the energy it makes available |
| Solid solutions | mole fraction of each endmember |
| Brine and vapor | salinities and the share of the fluid in the vapor |
| Isotopes | δ values and 87Sr/86Sr |
| Sorbed | what the exchange sites, surfaces and clays hold |
| Partition coefficients | Kd of each metal |
| Saturation index | log Q/K of each mineral, with a line at 0 (saturation); blocked minerals are marked |
| Properties | pH, ionic strength, solvent mass, log water activity, log fO2, Eh, and others where they apply (surface potential, net acidity) |

Above the plot: **Horizontal axis** (when more than one thing varies: temperature, pressure, time,
depth, amount mixed, concentration factor; a column offers pore volumes and distance from the inlet),
**Modes** (static and flow-through in a *Both* run; outlet, start and end profile in a column), **Axis**
log or linear, and **Units** (mol, mmol, µmol, mg, µg) for amounts.

**Shaded columns** behind the plot mark points where the thermodynamic data are less reliable: amber
approximate, rust unreliable (section 5, Reliability). **Shaded bands** around a line are the middle 90%
of the uncertainty runs.

**Downloads.** *Download plotted data (CSV)* gives the lines on the plot, with the path's temperature and
pressure and the 5% and 95% band where there is one. *Download the whole run (JSON)* gives everything the
run computed.

**End of the run.** For each mode, the solids left (mmol and g per kg of starting water), the fluid's
properties, and the dissolved elements (mmol/kg and mg/kg).

**Diagrams.** Eh-pH, log fO2-pH and activity diagrams at any temperature and pressure. Eh-pH and fO2-pH
take an **Element**, its **log activity** and **With** (other elements to include, such as S C); the
activity diagram takes a **Cation** and plots log a(cation)/a(H+) against log a(SiO2). Temperature and
pressure start at the end of the shown run, and the run's path is drawn from the open square (start) to
the red dot (end). Hatching marks where a boundary is within the uncertainty of log K at that temperature
and pressure; dashed lines are the stability limits of water.

**Assumptions and limits.** Every note the run produced: which data and models were used, which
species were extrapolated, where a model is outside its range, which minerals were blocked. Read it
before using a result.

## 5. The models

**Thermodynamic data.** Aqueous species, minerals and gases come from OBIGT, the database of CHNOSZ
(Dick 2019): revised HKF equations of state for aqueous species, Berman (1988) data for silicate minerals
where OBIGT has them, Maier-Kelley heat capacities for the rest. Water is IAPWS-95 (Wagner and Pruss
2002) with the dielectric constant of Johnson and Norton (1991). Equilibrium constants are computed
from Gibbs energies at the temperature and pressure of every point, so there is no fixed grid of log K
values and no fixed pressure.

**Activity models.**

| Model | Use it for | Range |
|---|---|---|
| B-dot (default) | Most natural waters up to about seawater salinity | Helgeson (1969), as in EQ3/6. Tabulated from 0 to 300 °C along the boiling curve, with no pressure term; above 300 °C it is held at its 300 °C value, and the card suggests b-gamma above 300 °C or 1 kbar |
| b-gamma | NaCl-dominated fluids at high temperature and pressure | Helgeson, Kirkham and Flowers (1981), Table 27, to 500 °C and 5 kbar; 550 to 600 °C from Manning et al. (2013) via CHNOSZ |
| Davies | Dilute waters | Davies (1962), with the 0.2 coefficient EQ3/6 uses |
| Pitzer | Brines and evaporites, to halite saturation and beyond | PHREEQC's pitzer.dat (Parkhurst and Appelo 2013), fitted for about 0 to 200 °C and 1 to 1000 atm. With Pitzer, the species and minerals pitzer.dat lists use its own equilibrium constants, and its evaporite minerals (glauberite, bloedite, polyhalite, kainite and others) and borate complexes are available only with this model |

The extended Debye-Hückel models (B-dot, b-gamma, Davies) are unreliable above an ionic strength of about
1 mol/kg; the run says so in Assumptions and limits when a path goes beyond it.

**Gas phase.** Peng-Robinson fugacity coefficients for the gas mixture.

**Solid solutions.** Ideal mixing on sites, as in EQ3/6. Real carbonates and feldspars unmix at low
temperature, so there the minor endmember is overstated.

**Kinetics.** Rate = surface area × the sum of each mechanism's rate constant (acid, neutral, base,
CO2 or other catalyst, each with its own activation energy) × (1 − Ω), the form of Lasaga (1984) and
Palandri and Kharaka (2004), so a rate is zero at equilibrium and changes sign past it. Basalt glass
follows Gislason and Oelkers (2003) and has no (1 − Ω) term because glass is always undersaturated.
Olivine in the serpentinization example follows the reaction-progress law of Malvoisin et al. (2012).
The seawater calcite and aragonite laws (above) are k(1 − Ω)^n in two regimes, joined smoothly, with an
activation energy of 25 kJ/mol (Naviaux et al. 2019, GCA). They measure dissolution only, so a
supersaturated water neither dissolves nor grows a shell.

**Marine carbonate system.** With the seawater carbonate rates, calcite and aragonite saturation comes
from the oceanographic constants that PyCO2SYS uses by default: K1 and K2 of Lueker et al. (2000),
calcite and aragonite solubility of Mucci (1983), pressure terms of Millero (1995), borate of Dickson
(1990) with total boron after Uppström (1974). These are the constants the seawater rates were measured
with. They hold for marine waters (salinity about 20 to 43, 0 to 40 °C).

**Pore water of a carbonate sediment.** Steady one-dimensional diffusion and reaction after Hales and
Emerson (1997): calcite dissolves at k(1 − Ω)^n per mass of calcite, oxic respiration makes CO2 in
proportion to the O2 consumed, and both diffuse to the bottom water. Left out: bioirrigation, suboxic
respiration (oxygen is not tracked), aragonite in the sediment, and precipitation. The diffusion
coefficient, the tortuosity correction and the alkalinity lost to respiration are taken from standard
sources that have not yet been checked against the originals; see section 6.

**Brine and vapor.** Where a saline fluid crosses the two-phase boundary of NaCl-H2O (Driesner and
Heinrich 2007), it splits into a brine and a vapor whose salinities come from that phase diagram.

**Isotopes.** δ18O, δD, δ13C, δ34S and 87Sr/86Sr are carried through the fluid and the minerals in
isotopic equilibrium at each step, with the fractionation factors listed on the Sources page (oxygen
between minerals and water from Vho et al. 2019; carbon and sulfur from Ohmoto and Rye 1979; water and
vapor from Horita and Wesolowski 1994). Strontium is carried with the fluid only. Starting values not
given on the card are taken as defaults, which the notes list.

**Sorption.** Ion exchange and surface complexation with PHREEQC's constants; clay edge sites after
Bradbury and Baeyens; Kd values of Allison and Allison (2005). The constants are for 25 °C.

**Catabolic energy.** The method of McCollom and Shock (1997): the Gibbs energy of each chemosynthetic
reaction in the mixture and the energy it can release per kg of fluid, limited by its scarcer reactant.

**Reliability.** The ionization constant of water is measured at all temperatures and pressures (IAPWS
R11-24, 2024) and depends on the same equations as every charged species, so Reaction compares its own
log Kw with it at every point. A point is *good* within 0.15 log units with water denser than
0.35 g/cm³, *approximate* within 0.35, and *unreliable* beyond. The card shows the rating of the path
before a run, the narrative repeats it, and the plot shades the points that are not good. The data hold
to within 0.15 up to about 350 °C and 2 kbar and degrade near the critical point of water and wherever
its density falls below about 0.35 g/cm³ (docs/DESIGN.md gives the table).

## 6. Limits

- Temperature 0 to 600 °C and pressure 1 to 5000 bar. Steam is not modelled: the path must keep liquid
  water (or a supercritical fluid).
- B-dot, b-gamma and Davies are for ionic strengths below about 1 mol/kg. Use Pitzer for brines, and
  keep Pitzer below about 200 °C and 1000 atm.
- The database lacks some minerals that form in experiments and nature (for example the magnesium
  hydroxide sulfate hydrate of heated seawater, and ferrihydrite and schwertmannite of mine waters), so
  the elements they would take up stay dissolved. The notes of the presets and examples say where this
  matters.
- Equilibrium means equilibrium. Minerals that form slowly must be blocked by hand, and the result depends
  on that choice; the default list and each example's list say what was blocked and why.
- Several presets (rainwater, groundwater, karst spring, formation brine, contaminated stream, deep
  seawater) are illustrative compositions, not samples. Use a measured analysis for a real site.
- Shell surface areas default to 2 m²/g, a placeholder; forams in particular have not been measured
  for this default.
- The seawater rate laws were fitted near 5 °C and in seawater; the activation energy carries them to
  other temperatures. Naviaux et al.'s in-situ Ω is from alkalinity and pH, which they found 5 to 10%
  below the Ω from alkalinity and DIC that Reaction computes, so with these rates dissolution sets in
  somewhat deeper than it would with theirs.
- The pore-water model uses a single exponential for the O2 consumption, and three of its parameters
  (the diffusion coefficient of HCO3-, Boudreau's tortuosity correction and the Redfield alkalinity
  change of 17 per 106 carbon) have not yet been checked against their sources.
- Flush mode does not yet keep sorbed species apart from the fluid that leaves. Columns cannot evaporate,
  mix a second fluid, follow brine or vapor, hold redox apart, suppress minerals, give catabolic energy
  or uncertainty bands.
- Sorption constants are for 25 °C.
- Evaporation and kinetics cannot be combined.

## 7. The worked examples

Each example's box in the app gives its full description and the card fields that set it up.

| Example | Shows | What to look at |
|---|---|---|
| Basalt alteration | Solid solutions | Seawater and 1 kg of basalt glass at 350 °C and 500 bar; chlorite, epidote and actinolite form as mixtures. Solid solutions tab |
| Evaporating lake | Evaporation, brines, isotopes | Great Salt Lake brine to a tenth of its water with Pitzer activities: halite first, then bloedite, leonite and epsomite, with calcium in glauberite and polyhalite. The water left grows heavier in 18O and D. Isotopes tab |
| Contaminated stream | Metals on clays | A 10 m streambed column with three clays and a Kd solid; metals break through in turn. Sorbed and Partition coefficients tabs |
| Boiling | Boiling and a gas phase | A CO2-bearing geothermal water at 150 °C rising from 25 to 5 bar: it boils near 10 bar, loses carbon to the gas, its pH rises and calcite scale forms. Gas phase tab |
| Supercritical seawater | Phase separation | Seawater at 300 bar from 390 to 410 °C splits into brine and vapor past its critical point; the run follows the brine. Brine and vapor tab |
| Flow-through weathering | Weathering in a column | Soil water through 1 m of albite and a little calcite: the calcite front moves down; albite near the inlet weathers to kaolinite |
| Serpentinization | Rates and solid solutions | Mantle olivine (Fo90) and a reduced salt water at 300 °C for three years: serpentine, Fe-bearing brucite and magnetite, and hydrogen |
| Calcite dissolving | Kinetics | 1 g of calcite in rainwater for 12 hours: fast at first, slowing to almost nothing as the water nears saturation. Rates tab |
| A foram sinking through the lysocline | Kinetics on a porous shell | A 20 µg foram falls from 1 to 5.5 km in flush mode with the in-situ seawater rate: it keeps its mass to the saturation horizon and loses a few percent below the lysocline |
| Foram and pteropod: calcite vs aragonite | Two polymorphs | The same fall with a pteropod: aragonite, the more soluble, dissolves far more |
| Ceará Rise: foram and pteropod sinking | A measured water column | The same two shells through the GLODAPv2.2023 profile at the Ceará Rise (Lauvset et al. 2024), western equatorial Atlantic |
| Acid mine drainage | Acid from pyrite, acid stored in iron | Pyrite oxidized in groundwater: pH falls to about 2, and net acidity (Kirby and Cravotta 2005) shows the acid held by dissolved iron. Properties tab |
| Vent fluid meets seawater | Chemosynthesis energy | A 350 °C vent fluid mixed into 2 °C seawater with redox held apart; which reactions offer the most energy at each temperature. Energy tab |

## 8. Using Reaction in a class

- **Change one thing.** Every example names the fields that make it work. Ask students to predict, change
  one field, run, and compare the two runs with **Compare with**. Good first changes: the mode (static
  against flow-through), the activity model (B-dot against Pitzer in the evaporating lake), the blocked
  minerals (unblock dolomite in a seawater run, or goethite in acid mine drainage), the carbonate rates
  (in-situ against laboratory for the sinking foram), or the surface area.
- **Read the notes.** Assumptions and limits is part of the result. Ask students which notes apply to
  their run and what each means for their conclusion.
- **Use the reliability rating.** Run seawater past 400 °C at low pressure and look at the shaded points
  and the rating; then raise the pressure.
- **Bring real data.** Paste a measured analysis (field or literature) as a custom fluid, check "As
  read", and model it. Give carbonate as HCO3, not alkalinity.
- **Keep the data.** Download the CSV of each run students hand in; workspaces are deleted after 12 hours
  without use.
- **Cite the data.** A result rests on the data and models listed on the Sources page; students should
  cite those along with Reaction (section 10).

## 9. Running Reaction on your own computer

The same app runs locally with Python 3.10 or later:

    pip install -r requirements.txt
    streamlit run streamlit_app.py

or, for the desktop version with its own two windows:

    pip install -e .
    reaction-app

Nothing leaves your computer. The desktop version saves runs to `workspace/` (set `REACTION_WORKSPACE`
to change it). It can also take requests in words through a local language model (Ollama); that chat is
switched off by default, and when on it only fills in the card, never the numbers. The README describes
the Python library for scripted calculations.

## 10. How to cite

If you use Reaction in a study or a class, please cite it:

> Dougalien Systems (2026) Reaction: open geochemical reaction-path modelling, version 0.1.0.
> https://github.com/dougalien/Reaction

and cite the data and methods your result rests on, which the Sources page lists. GitHub's "Cite this
repository" button gives the same reference from `CITATION.cff`.

Reaction is free software under the GNU General Public License, version 3 or later. It comes with no
warranty.

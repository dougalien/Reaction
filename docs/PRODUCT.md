# Reaction: what it is for

Written 2026-10-15 after a design conversation with the user. This supersedes the "input files first"
framing in DESIGN.md. DESIGN.md still describes the physics and data layers; this file describes the
product those layers serve.

## Purpose

A geochemical calculator for one person (a geochemist who models seafloor hydrothermal systems and
carbon-sequestration sites). Not for sale and not built for a general audience. The goal is to ask a
question in plain language, get a result, and learn from it graphically. Setup effort has kept many
questions unresearched; removing that effort is the point.

Not a re-implementation of EQ3/6. EQ3/6 is the example of the bare-bones workflow this replaces. It
survives only as an oracle for checking numbers, never as a design driver. The reference for usability
is Geochemist's Workbench: simple, visual, quick to change.

## Questions it should answer (the user's own examples)

1. A methane-bearing fluid reacts with calcite at 100 C (gas seep).
2. How do mineral reaction rates behave under static versus flow-through conditions?
3. Can aragonite converting to calcite over time be simulated? (kinetics)
4. What happens to dissolved silica if quartz is undersaturated or not allowed to precipitate in a
   mid-ocean-ridge reaction zone? (phase what-ifs)
5. The first end-to-end target: "React seawater and basalt at 350 C", static and flow-through, over a
   range of water/rock ratios.

## How it is used

1. The user types a question. The local LLM turns it into a visible, editable scenario card and asks only
   about what is missing ("what W/R? static or flow-through?").
2. The solver runs. The LLM never produces numbers.
3. A results view shows "here is what we found" plus interactive plots. Any species, mineral, ratio or
   property can go on an axis (CaCl+ instead of Ca++ is one click). Log/linear, overlays, static versus
   flow-through together, fluid versus mineral views.
4. Every run is saved with its scenario and the user's notes so earlier questions can be found again.

A clickable mock-up of this flow was shown in the chat on 2026-10-15 (placeholder numbers) and the user
approved the look.

## Layout: two windows

The user works on a two-screen laptop setup. The chat and scenario card live in one window and the
results in another, one on each screen. So the app is a local server (nothing leaves the machine) with
two browser windows that share one session: a control window (chat, scenario card, run history) and a
results window (narrative, plots, species picker). Changes in either update the other. A second results
window must be possible for comparing runs.

## Rules

- Numbers come only from the solver. The LLM fills the card, asks questions and narrates, and the
  narration quotes values from the result table.
- Everything works by hand with forms and clicks. The LLM is a convenience, never a dependency.
- Every datum has a source (database, paper, fit). Rate laws included; Palandri and Kharaka (2004) is the
  likely starting compilation for dissolution kinetics.
- The basis, switching and suppression details are chosen by the program, never shown to the user unless
  asked.
- Solid solutions are a standard phase type, not an add-on. Isotopes (O, H, C, S, Sr) come as a layer on
  the solved speciation.

## Local LLM

Ollama on the user's machine: `qwen2.5:7b` (4.7 GB), `llama3.2:3b`, `nomic-embed-text`; GPU is an RTX
2060 with 6 GB. The 7B model has to fit beside the calculations. Use constrained (JSON-schema) output for
the scenario card so a small model cannot wander. The embedding model can retrieve presets (seawater,
basalt) and earlier runs.

## What the engine needs, against what exists

| Capability | Status |
|---|---|
| Water properties, HKF, logK from 0 to about 400 C, 1 kbar | done (M2), validated |
| Speciation, B-dot/Davies, redox, pH scales | done (M1/M1.5), validated against 41 EQ3NR cases |
| Closed-system equilibrium with minerals | works for carbonate/sulfate/hydroxide systems and for seawater-basalt silicate assemblages with redox (Berman minerals); EQ6 `methane` and `pyrsw` still differ at redox fronts |
| Reaction path with reactants, T path, fixed fugacity | done in the app: minerals, species, basalt glass; static and flow-through |
| Flow-through (products left behind as they would redissolve) | done, in the general path class; the app uses it (EQ6 `microft` matches) |
| Fluid mixing, titration | titration cases (`rwtitr`, `swtitr`) match EQ6; mixing done in the app |
| Kinetics (rate laws, time axis, surface area) | not started |
| Solid solutions | not started; affects solver design |
| Gas phases and methane chemistry | partial; needs a decision about redox coupling |
| Isotopes | O, C, S and 87Sr/86Sr along every run (P7c); H not yet |

## Open decisions

- Solver core: keep extending the mass-action Newton plus active-set approach, or move to Gibbs-energy
  minimization, which handles phase appearance and solid solutions more systematically. Decide before more
  solver work.
- Which solid solutions and which isotope systems first.
- Redox coupling in the gas-seep case: whether CH4 and sulfate are allowed to equilibrate at 100 C is a
  modelling choice the program should ask about.

## Built so far (the skeleton, 2026-10-15)

`src/reaction/app/` holds the two-window calculator, standard library only on the server side:

- `scenario.py`: the card, a flat set of fields (fluid, start and end temperature, start and end pressure,
  steps, mode, blocked minerals, activity model) with validation (for example water must stay liquid all
  along the path). It is stored with every result, so a run can be reproduced from its card.
- `engine.py`: runs a card. Takes a fluid along a straight line in temperature and pressure and finds the
  closed-system equilibrium at each step. Static keeps minerals; flow-through leaves each step's minerals
  behind. Uses OBIGT/HKF logK at every (T, P); also rates how far to trust the thermodynamics along the path.
- `assistant.py`: the chat. Ollama (JSON-schema output) fills the card and asks about missing details; plain
  rules are the fallback. Numbers and modes from the model are kept only if the user's words contain them,
  because the 7B model invents values otherwise.
- `view.py`: what the results window shows; chat and clicks both go through it.
- `narrate.py`: the "here is what we found" text, every number read from the result tables.
- `session.py`, `server.py`, `static/`: shared state, the local server (loopback only, Host/Origin checked,
  never shares a busy port) and the two pages. The plot is a small dependency-free SVG chart (log axis,
  hover readout, SVG and CSV export), so nothing is loaded from the internet.

Decided while building it: the scenario format is the flat card above; the plot library is our own SVG
chart; the app needs no packages beyond numpy and scipy.

Since 2026-10-25 a run can start from any fluid (preset or a pasted analysis, charge-balanced), add reactants
(water/rock, titration) and mix in a second fluid, along any temperature-pressure line; results can be
compared run against run, in two results windows. Since P3 (also 2026-10-25) the redox state is on the card,
the silicates are the Berman (1988) dataset that OBIGT uses, and basalt glass reacts with seawater at 350 C to
the greenschist assemblage (albite, chlorite, epidote, actinolite, sulfides), static and flow-through, down to
a water/rock ratio of 1 or below. Since P4 the results window draws Eh-pH, log fO2-pH and activity diagrams
(a(cation)/a(H+) against a(SiO2), aluminium conserved) at any temperature and pressure, from the same logK as
the runs, with the open run's path on top and the rest of the fluid held at the run's values; boundaries closer
than the logK uncertainty at that (T, P) are hatched. Since P5 a free gas phase (CO2, CH4, steam, H2, H2S,
N2 ..., Peng-Robinson fugacities) forms and dissolves along any path, is vented in flow-through, and gases
can be added; methane and sulfate are held apart unless the card says they may react. Since P6 a run can
last a given time with minerals and basalt glass reacting at measured rates (Palandri and Kharaka 2004;
Gislason and Oelkers 2003) from a reactive surface area: they dissolve toward saturation, grow on their seed
past it, and the results show each rate and how much is left. Since P7a brines use Pitzer's equations and
a fluid can be evaporated (seawater to its bittern salts). Any run can carry an uncertainty range: it is
repeated with the log K values shifted at random by their uncertainty (larger where the reliability check is
worse) and the plots shade the middle 90%, something GWB does not offer. Since P7b solid solutions (plagioclase, carbonate, chlorite, epidote, olivine, pyroxenes and others) can form with their compositions followed along the run. Since P5b every run reports where a saline fluid splits into brine and vapour (NaCl-H2O from Driesner and Heinrich's correlations, within 3 bar of Bischoff and Rosenbauer's measured boundary for seawater), with the vapour share and both salinities, and a run can follow the brine as the vapour boils off, or the vapour as it leaves the brine behind. Since P7c a run can carry δ18O, δ13C, δ34S and 87Sr/86Sr through the fluid and every mineral, from published fractionation factors, and δD of the water and the hydrous minerals (mineral-water D/H from Suzuoki and Epstein, Graham and others, and Sheppard and Gilg). Since P7d ion exchange and surface complexation on ferric oxide (PHREEQC's constants) are solved with the rest, and agree with PHREEQC. Since P7e a fluid can flow through a column of cells (1D reactive transport, PHREEQC's scheme, matching PHREEQC's TRANSPORT), on a temperature gradient from inlet to outlet if wanted, with the outlet's breakthrough curves and the minerals along the column.

## Pressure and temperature are free variables

There are no fixed-pressure databases. Every equilibrium constant is computed from Gibbs energies at the
temperature and pressure of each point (IAPWS-95 water, Johnson-Norton dielectric constant, HKF equations for
aqueous species, SUPCRT92 volume terms for minerals), so any T from 0 to 600 C and any P from 1 to 5000 bar
can be chosen, and either can change along a run. A run follows a straight line from a start (T, P) to an end
(T, P); equal values hold a variable fixed. This replaces the old workflow of generating a database per
pressure (saturation, 300, 400, 500 bar) with SUPCRT92.

How far to trust it. The water ionization constant is measured everywhere (Marshall and Franck, 1981) and
depends on the same dielectric constant, g function and Born coefficients as every charged species, so the
engine's own value, compared with the measurement, is an independent check at any (T, P). The app runs it
along every path, flags stretches as good (within 0.15 log units and density at least 0.35 g/cm3), approximate
(within 0.35) or unreliable, and shows this before the run (card advice), in the narrative, in the notes and as
a coloured strip under the plot. Measured on 2026-10-15 (engine minus measurement, log units):

| T (C) | 500 bar | 1000 bar | 2000 bar | 5000 bar |
|---|---|---|---|---|
| 25 | -0.01 | -0.03 | -0.07 | -0.24 |
| 100 | -0.01 | -0.01 | -0.03 | -0.16 |
| 200 | +0.04 | +0.06 | +0.07 | 0.00 |
| 300 | +0.04 | +0.12 | +0.19 | +0.17 |
| 350 | -0.03 | +0.10 | +0.23 | +0.25 |
| 400 | -0.53 | +0.03 | +0.26 | +0.33 |
| 500 | -1.76 (low density) | -0.36 | +0.26 | +0.46 |
| 600 | -0.61 (low density) | -0.84 | +0.22 | +0.58 |

The surface is smooth, so the pattern is the model's accuracy and not a bug: sound to about 350 C and 2 kbar,
poor near the critical point and wherever water density is below about 0.35 g/cm3. The reaction volumes check
out too (water ionization -22 cm3/mol at 25 C; calcite and anhydrite solubility rise with pressure).
`tests/test_pressure.py` holds these checks.

Not covered yet, and the most important for gas-bearing systems: free gas phases. A gas species'
Gibbs energy is its 1-bar ideal-gas standard state, so the fugacity reported at equilibrium is correct at any
pressure, but there is no gas phase in the assemblage and no fugacity coefficients to turn fugacity into
partial pressure (CO2, CH4, H2O vapour), no gas-liquid equilibrium and no CO2 solubility model at high P. Other
limits: minerals have constant molar volume (as in SUPCRT92), the B-dot activity model is tabulated only to
300 C (held above), and ionic-strength effects beyond a few molal need Pitzer.

## Solver status

The closed-system solver (`equilibrium.py`) is used by the app. Against EQ6's own output: `micro`, `rwtitr`
and `swtitr` match to the printed precision (tested in `tests/test_eq6_cases.py`); `methane` (49 of 1203
numbers outside tolerance, mostly EQ6's pre-conversion reporting at xi = 0 and uranium phases near their
boundaries) and `pyrsw` (17 of 935) now run through a 40-unit jump in log fO2 but are not yet exact. Two
changes got there: a sign-split logarithmic mass balance for components held in no mineral (fixed the
titration cases; using it for mineral-bound components broke `micro`, so it is not applied there), and the
solver now raises an error instead of returning a state that is still supersaturated in a mineral it failed
to add, which makes the step subdivision follow the redox front. `dedolo` (alias minerals with identical
formulas, such as Dolomite and Dolomite-ord) runs since the active set swaps a mineral whose saturation
condition depends on the assemblage's instead of adding it.  Flow-through is in the general path class
(EQ6's rule: a mineral is left behind where it would start to redissolve) and the app uses it; `microft`
matches to 3 of 640 values.  Still open: the cases that need special basis switches, solid solutions, ion exchangers or kinetics.

The solver, paths and EQ6 harness were committed in `3af9559`, the app in `a3076e0` and pressure paths in
`3f37964`.

# Start here (written the evening of 2026-10-15)

Read this first when picking the project back up. Product plan and accuracy map: `PRODUCT.md`. Physics and
data layers: `DESIGN.md`. EQ6 survey: `M3_NOTES.md`.

## Where things stand

Reaction is a personal geochemical calculator. Two browser windows share a local server: a control window
(chat with a local language model, a scenario card, run history) and a results window (what we found, an
interactive plot, a picker for any species, element, mineral or property). The language model fills in the
card and the plot request; it never produces numbers.

What works: any fluid (presets for seawater, rainwater, vent fluid and formation brine, or a pasted analysis
with units and automatic charge balance) taken along a straight line in temperature and pressure (0 to 600 C,
1 to 5000 bar, either or both changing), with reactants added steadily (minerals or species in g or mol per kg
of water: water/rock and titration) and a second fluid mixed in (with the temperature following the mixing if
its own temperature is given), static and flow-through, with minerals precipitating and dissolving. No fixed-pressure
databases: logK is computed from Gibbs energies at every point. A reliability rating (water ionization
against Marshall and Franck 1981) is shown before a run, in the narrative and as a strip under the plot.

Commits in the Reaction repo (remote: github.com/dougalien/Reaction, private): M1 `8d77c10`, M1.5 `8d9cd4a`, M2 `e90e9d9`,
M3 solver and EQ6 harness `3af9559`, the app `a3076e0`, pressure paths `3f37964`.

Housekeeping committed on 2026-10-25: this memo, `scripts/eq6_report.py`, tests default to the app's data
folder (`tests/conftest.py`, `tests/test_eq36_cases.py`).

P1 (solver: phase rule, flow-through in the path class, speed) and P2 (custom fluids, reactants, mixing,
run comparison, a second results window, mg/kg units) were done on 2026-10-25 as pull requests #2 to #8;
see the roadmap doc linked in the project memory.

P3 (seawater and basalt at 350 C) was done on 2026-10-25 as pull requests #9 to #11:
- Berman (1988) minerals (`thermo/berman.py`, CHNOSZ `extdata/Berman/*.csv` in `OBIGT/Berman/`) replace the
  Helgeson silicates in the "current" data, as in OBIGT's defaults. Minerals hold only to their transition
  temperature (CHNOSZ `exceed.Ttr = FALSE`).
- Redox on the card (`redox`, `redox_value`): from a redox species in the analysis (O2, H2, Fe(III), sulfide,
  nitrate) as the couple, air, log fO2, Eh or pe. Seawater and rain are air-saturated, the vent fluid has
  H2S and an illustrative 0.5 mmol/kg H2. log fO2 and Eh can be plotted.
- Basalt glass as a reactant (`app/rocks.py`, MORB glass, a tenth of the iron ferric, Ti left out): it
  dissolves congruently; the water/rock axis comes with it.
- Solver: log fO2 starts where the redox budget is matched; a failed warm start retries from the bulk
  composition; a mineral whose addition will not solve replaces the most similar mineral (or the next one
  in the ratio test) instead of failing the point. Seawater + 1 kg basalt glass at 350 C, 11 steps: static
  2 s, flow-through 8 s, no failed points; it gives albite, chlorite, epidote, actinolite, pyrrhotite.

P4 (diagrams) was done on 2026-10-25: `reaction/diagram.py` (predominance on a grid, balanced on the element;
water and quartz-saturation lines; an uncertainty band from the water-ionization check), `app/diagrams.py`
(requests in words: "Eh-pH diagram for Fe at 300 C", "fO2-pH diagram of iron with sulfur", "activity diagram
for K"; the run's path and fluid), a panel in the results window (`static/diagram.js`).

P5 (gas phase) was done on 2026-10-25: `thermo/gas.py` (Peng-Robinson mixture with k_ij), a free gas phase
in `equilibrium.py` (one more unknown, the moles of gas, and one more condition, the mole fractions f/(phi P)
summing to one; it appears and disappears like a mineral), flow-through vents it. Methane (OBIGT
organic_aq/organic_gas, only CH4 is read) is held apart from sulfate by default (`hold_apart`: CH4 is its own
component, CH4(g) is tied to it); the card's `methane_sulfate` = react lets them equilibrate (Douglas chose
this default). Gases are reactants ("2 mol CO2 gas", "methane"). CO2 solubility is within 10% of Duan and
Sun (2003) to 200 bar; at 400 bar the aqueous side runs about 20% low (a Duan-type model is the next step).

P5b (phase separation, NaCl-H2O after Driesner and Heinrich 2007) is on hold: the coefficient tables are in
the paywalled paper and neither Douglas nor I could get them (2026-10-25). Look for an open implementation
with a usable license (xThermal, for example) before writing it.

P6 (kinetics) was done on 2026-10-25: `kinetics.py` (rate laws r = A sum k(T) prod a^n (1 - Omega), the
Palandri and Kharaka 2004 table in `data/rates_palandri_kharaka_2004.tsv`, 60-odd minerals mapped to the
OBIGT names; basalt glass after Gislason and Oelkers 2003), kinetic reactants in `path.py` (reaction progress
is time in seconds; amounts advanced by backward Euler with an error test against the explicit step,
Newton on all kinetic amounts with a bracketed sweep as the fallback; a kinetic mineral grows on its seed
when supersaturated; its other structural states, quartz,beta or albite,low, cannot precipitate in its place).
The card has `kinetics`, `duration`, `duration_unit` and `surface_area_cm2_per_g` (surface shrinks as
(left/start)^(2/3)); reactants without a rate are added steadily over the run. Results have a Rates tab
(rate and % left of each kinetic reactant). Judgment calls: kyanite's acid activation energy is printed
-53.9 in the PHREEQC transcription, +53.9 is used; the PK base mechanism exponent is applied to a(H+).
Next: Hermanska et al. (2022, 2023) rates as an alternative set, and precipitation rates and nucleation.

308 tests pass with the data folder present. The data folder needs `OBIGT/Berman/` (run
`scripts/fetch_references.ps1` again); without it the old minerals are used and every run says so.

## First things to do, in this order

1. **Ask how the app felt.** The user has tried heating runs but not yet the pressure paths. Ask about the
   layout, the wording, the reliability strip, anything missing. Their feedback outranks everything below.
2. ~~Commit the housekeeping~~ (done 2026-10-25, with the GitHub remote).
3. **Check the pressure dependence against real fixed-pressure databases.** They were in the EQ3/6 package
   all along, in `%LOCALAPPDATA%\Reaction\external\EQ36_DB`: `data0.500` (0 to 350 C every 25 C at 500 bar),
   `data0.1kb` (to 600 C), `data0.2kb` (to 900 C), `data0.5kb` (5 to 1000 C). They are the SUPCRT92-derived
   "database per pressure" the user used to build by hand, so they are a species-wide oracle for logK beyond
   the saturation curve, which until now was only checked with water ionization. Their layout differs from
   `data0.cmp`: a "Temperature grid" and "Pressure grid" in the header (up to 27 points, three ranges), and
   each species has a `* Log K` line followed by that many values instead of the `**** logK grid` marker with
   8 values. Plan: teach `database/eq36.py` to read it, generalize `validation/logk_compare.py` to compare
   the engine's logK at those exact (T, P) points, then add a test and update the accuracy section of
   `PRODUCT.md` (including where the engine differs from SUPCRT92 itself, if anywhere).
4. **Choose the next big piece with the user.** Their call. My suggestions, with the reasons:
   - *Free gas phase with fugacity coefficients* (CO2, CH4, H2O vapour; gas-liquid equilibrium; CO2
     solubility at high P). This is the biggest gap for the gas seep and for circulation where P varies. A
     modelling question to ask first: should CH4 and sulfate be allowed to react at 100 C, or held apart as
     they are at low temperature?
   - *Seawater reacting with basalt* (their first stated target): fluid presets with Si, Al, Fe; basalt as a
     reactant; a water/rock axis; redox; many more minerals. It needs the solver hardened first.
   - *Solver core*: keep extending the Newton plus active-set approach, or move to Gibbs-energy minimization
     (cleaner for phase appearance and solid solutions). Decide before the basalt work.
5. Later: ~~kinetics~~ (done in P6; aragonite to calcite needs aragonite's rate, not in the PK table),
   solid solutions, isotopes (O, H, C, S, Sr as a layer on the solved speciation).

## Known technical gaps

- **Solver against EQ6** (`python scripts/eq6_report.py`): `micro`, `rwtitr`, `swtitr` exact. `methane` 49 of
  1203 numbers outside tolerance (mostly EQ6's reporting at xi = 0 before reactants convert, plus uranium
  phases at boundaries), `pyrsw` 17 of 935, `dedolo` 3 of 236 (EQ6 reports the reactant unconverted at
  xi = 0, plus 4 nmol of calcite; alias minerals such as Dolomite and Dolomite-ord now swap by the phase
  rule instead of failing), `microft` 3 of 640 (pH within 1.5e-3 past the muscovite boundary; flow-through now
  lives in the general path class and the app uses it too). Not run: `heatswfl` (it forms a saponite solid
  solution), special basis switches
  (`heatqf`, `j13wtitr`, `j13wsf`), solid solutions (`heatsw`, `rwssdiag`, `pptmins`, `j13wtuff`), ion
  exchangers (`swxrca`, `swxrcaft`), kinetics (`crisqtz`, `pptcal`, `pptqtz`, `pptqtza`).
- **How the solver fails**: where a redox couple moves by tens of log units the active-set search used to pick
  wrong minerals. It now raises an error instead of returning a state that is still supersaturated in a
  mineral it could not add, and the step subdivision follows the front. A sign-split logarithmic balance is
  used only for components held in no mineral (applying it to mineral-bound components broke `micro`).
- **Thermodynamics**: HKF is poor near the critical point and below 0.35 g/cm3 (flagged in the app); B-dot is
  tabulated to 300 C only; no Pitzer for high ionic strength; minerals have constant molar volume (as in
  SUPCRT92); Born Q is 11% above the SUPCRT-implied value; seawater lacks the Mg hydroxide sulfate hydrate seen
  in heating experiments, so Mg removal and the acid it releases are underestimated.
- **App**: the reliability strip is thin (5 px); "Earlier runs" in the user's workspace includes my test
  runs. Only preset fluids can be mixed in (not a second custom analysis). Mixing temperature assumes equal
  heat capacities. Alkalinity is not read from an analysis yet; give total carbonate as HCO3.
- **Minerals**: no smectites (saponite, nontronite, beidellite) or celadonite in OBIGT, so low-temperature
  basalt alteration shows the equilibrium greenschist-like minerals (albite, K-feldspar, pumpellyite,
  andradite) instead of clays and zeolites; block those to mimic kinetics. Endmembers only until solid
  solutions (P7), so iron goes to minnesotaite, annite and daphnite rather than into chlorite and actinolite.
  The zeolites (laumontite, wairakite, analcime) and some others are still the Helgeson (1978) data, so they
  are not strictly consistent with the Berman minerals.

## How to run things

| What | Command |
|---|---|
| Start the app | `scripts\reaction.cmd` (opens http://127.0.0.1:8790/control; the "Open results window" button gives the second window) |
| Run all tests | `python -m pytest -q` from the repo root (no environment variables needed) |
| EQ6 path report | `python scripts/eq6_report.py [case ...]` |
| Language model | Ollama with `qwen2.5:7b` (RTX 2060, 6 GB). Without it the chat falls back to plain rules |

Data (never committed; only copied locally): `%LOCALAPPDATA%\Reaction\external` holds `OBIGT/` (CHNOSZ csv
files), `EQ36_DB/` (all the EQ3/6 `data0.*` files), `EQ36_tests/` (3tlib_* and 6tlib_*), and
`CHNOSZ_R_reference_only/` (GPL-3 R sources, read for reference only). Override with `REACTION_DATA_DIR`.
Saved runs and the chat session are in `workspace/` in the repo (ignored by git).

## Things that will bite

- Port 8765 belongs to another program (Galaxy.exe). The app uses 8790 and refuses to share a busy port;
  Windows lets two Python servers share a port by default, which silently loses requests.
- In this environment, shell heredocs collapse `\\` to `\`, so a `\b` in a regex written through a Python
  heredoc becomes a backspace character. Write code that contains backslashes with the file tools.
- The 7B model invents values (it once set flow-through mode and a pressure nobody mentioned). The app keeps a
  number or mode from the model only if the user's own words contain it; keep that guard when changing prompts.
- The browser tool's `wait` allows at most 10 seconds. A hidden browser tab can fail to focus; calling the
  page's own `api(...)` function from `javascript_tool` does the same thing as the Send button.
- Boiling is checked along the path: at 300 C water needs at least 86 bar, at 350 C at least 165 bar.

## Decisions already made (do not reopen without a reason)

Own Python solver (numpy and scipy), not a wrapper; EQ3/6 is an oracle for numbers only and never shapes the
design; the interface is question first (chat, then card, then run, then plots that can be changed on the
fly); two synced windows for a two-screen setup; everything local; the language model only fills the card and
the plot request; numbers come only from the solver; the plot is our own SVG chart (no CDN); the card is a flat
set of fields stored with every result; no fixed-pressure databases; ask before every commit.

## Where the code is

`src/reaction/` : `database/` (data0 reader), `thermo/` (IAPWS-95 water, dielectric, HKF, crystals, OBIGT loader,
`eosdb.py` logK at any T and P, `validity.py` reliability check), `network.py`, `aqueous.py`, `speciation.py`,
`equilibrium.py` (closed-system solver), `path.py`, `validation/` (EQ3NR and EQ6 harnesses), `app/` (scenario,
engine, assistant, view, narrate, session, server, `static/` pages). Tests are in `tests/`; the pressure checks are
`tests/test_pressure.py`.

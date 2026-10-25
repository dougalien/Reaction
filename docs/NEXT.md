# Start here (written the evening of 2026-10-15)

Read this first when picking the project back up. Product plan and accuracy map: `PRODUCT.md`. Physics and
data layers: `DESIGN.md`. EQ6 survey: `M3_NOTES.md`.

## Where things stand

Reaction is a personal geochemical calculator. Two browser windows share a local server: a control window
(chat with a local language model, a scenario card, run history) and a results window (what we found, an
interactive plot, a picker for any species, element, mineral or property). The language model fills in the
card and the plot request; it never produces numbers.

What works: seawater taken along a straight line in temperature and pressure (0 to 600 C, 1 to 5000 bar,
either or both changing), static and flow-through, with minerals precipitating and dissolving. No fixed-pressure
databases: logK is computed from Gibbs energies at every point. A reliability rating (water ionization
against Marshall and Franck 1981) is shown before a run, in the narrative and as a strip under the plot.

Commits in the Reaction repo (remote: github.com/dougalien/Reaction, private): M1 `8d77c10`, M1.5 `8d9cd4a`, M2 `e90e9d9`,
M3 solver and EQ6 harness `3af9559`, the app `a3076e0`, pressure paths `3f37964`.

Housekeeping committed on 2026-10-25: this memo, `scripts/eq6_report.py`, tests default to the app's data
folder (`tests/conftest.py`, `tests/test_eq36_cases.py`).

218 tests pass with no environment variables set.

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
5. Later: kinetics (aragonite to calcite, rate laws with sources, Palandri and Kharaka 2004 as the start),
   solid solutions, isotopes (O, H, C, S, Sr as a layer on the solved speciation).

## Known technical gaps

- **Solver against EQ6** (`python scripts/eq6_report.py`): `micro`, `rwtitr`, `swtitr` exact. `methane` 49 of
  1203 numbers outside tolerance (mostly EQ6's reporting at xi = 0 before reactants convert, plus uranium
  phases at boundaries), `pyrsw` 17 of 935. Not run: `dedolo` (alias minerals with identical data, e.g.
  Dolomite and Dolomite-ord, end in "could not add the supersaturated mineral"), `microft` and `heatswfl`
  (flow-through in the general path class; the app's heating runner has its own), special basis switches
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
- **App**: the reliability strip is thin (5 px); no run comparison or second results window yet; "Earlier
  runs" in the user's workspace includes my test runs; no units switch beyond mol, mmol and umol per kg.

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

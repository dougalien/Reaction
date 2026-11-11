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
against IAPWS R11-24, 2024; Marshall and Franck 1981 until 2026-11-11) is shown before a run, in the narrative and as a strip under the plot.

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

P5b (phase separation) was unblocked on 2026-10-25 by Oguntola's open lookup table: see below.

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

P7a (Pitzer brines and evaporation) was done on 2026-10-25: `pitzer.py` (Harvie-Moller-Weare equations with
E-theta; parameters from PHREEQC's pitzer.dat in `data/pitzer_phreeqc.tsv`, with their temperature terms),
`activity_model = "pitzer"` on the card. In Pitzer mode the complexes made only of the Pitzer system's
components (NaCl, CaSO4, NaSO4- ...) are dropped, since the parameters carry them; complexes of other
elements stay and get only the Debye-Hueckel term. The gamma iteration is under-relaxed for strong brines.
Checks: NaCl and CaCl2 gamma and osmotic coefficients within 1-2% of Robinson and Stokes to 6 m, seawater
a(H2O) 0.9813, halite 6.15 mol/kg at 25 C and 6.7 at 100 C. New card field `evaporate` (concentration
factor, steps even in its logarithm): seawater gives gypsum at about 4.6x, halite at 11x, epsomite at 42x.
Judgment calls: alpha1 = 1.4 and alpha2 = 12 for 2-2 pairs, otherwise alpha1 = 2 and alpha2 = 50 (the
PHREEQC defaults as I understand them); no pressure term in the Pitzer parameters (HKF carries pressure
in the log K); the NBS pH convention is not applied with Pitzer. A run whose ionic strength passes 1
mol/kg without Pitzer says so.

Uncertainty bands (Douglas's request, done): card field `uncertainty_runs` (0 = off, 5 to 200; "with
uncertainty" in the chat gives 20). The run is repeated with every log K shifted by its own standard-normal
draw times sigma(T, P) = sqrt(0.1^2 + dev^2), dev being the water-ionization deviation of the reliability
check, and at least 0.35 below 0.35 g/cm3 (`thermo.validity.logk_sigma`). Each draw is kept along the whole
path; polymorphs (same formula) share a draw, and H2O(g) is never shifted (steam against water is exact).
Draws are seeded by run number, so a card always gives the same range. The runs go to spawned worker
processes (all processors but one). Plots shade the 5th to 95th percentile, widened to include the main
line; hover shows the range, the CSV gets 5% and 95% columns, and the notes list minerals that appear or
disappear within the uncertainty. Judgment calls: independent draws for every species and mineral (no
covariance from shared data is known to us); 0.1 log units as the base sigma.

P7b solid solutions (done): card switch `solid_solutions` (off by default; "with solid solutions" in the
chat). Eleven (now twelve, with serpentine) ideal site-mixing models in `data/solid_solutions.tsv` (data0's carbonate, plagioclase, alkali
feldspar, chlorite, epidote, garnet, orthopyroxene, biotite, plus olivine, clinopyroxene and actinolite).
In the solver a solid solution is a mixture like the gas phase: its total moles are an unknown and its
condition is sum x_i = 1 with x_i = (Q_i/K_i)^(1/sites); when it is supersaturated it takes over its pure
endmembers, and if the assemblage then has too many phases it replaces the mineral most like its main
endmember. Results get the mixtures under Mineral amounts and a Solid solutions tab of mole fractions.
Seawater-basalt at 350 C now gives albite, chlorite (about 25% daphnite), epidote (about 90% Fe epidote),
actinolite and quartz. The EQ6 heatsw cases still cannot be compared: they form saponite, which OBIGT lacks.
Judgment calls: ideal mixing everywhere (overstates mixing in carbonates and feldspars at low T); blocked
minerals are left out of the mixtures too; magnesite is blocked by default, so it is not in the carbonate.

P5b phase separation (first part, done): `data/nacl_h2o_phases.tsv` holds the coexisting brine and vapour
of NaCl-H2O on a 5 C by 20 bar grid (1 to 996 C, 10 to 490 bar), with halite saturation, extracted from
Oguntola's CC BY 4.0 lookup table (Zenodo record 20394023, from Driesner and Heinrich 2007). `phasesep.py`
gives the phase state of a fluid of given NaCl-equivalent salinity (from its chloride). Every run reports
where the fluid splits, the vapour share (lever rule) and both salinities (Brine and vapour tab, story,
notes), also at points where the aqueous chemistry fails near the critical point. Seawater at 250 bar
splits from about 390 C. Not yet: the chemistry still follows the whole fluid; following the brine or the
vapour (partitioning of every solute and gas) is the next step, then paths below the boiling curve.
Judgment calls: NaCl-equivalent from Cl; bilinear interpolation, nearest corner at field edges.

P5b part two (done): the card's `follow_phase` (whole, brine, vapour; chat "follow the brine", "let the
vapour boil off", "leave the brine behind"; a Follow select on the card). Where the fluid splits, the
path's `separate` hook removes the other phase after each step and the point is solved again. Water goes
by the lever rule, solutes in proportion to the salt in each phase, neutral species with a gas of the same
formula (CO2, H2S, HCl, ...) with the water; minerals stay with the brine and a gas phase goes with the
vapour. Seawater at 300 bar, 405 C: the brine left has 6.6 wt% NaCl, the vapour 0.56 wt%, as the table
says. Limits: brines past about 20 wt% (seawater above about 410 C) are beyond B-dot and points there are
blank; the near-critical points fail as before. A failing step is halved at most 40 times per point
(before, a failure could halve both halves 25 deep). Judgment calls: gases with the water (they partition
roughly evenly per kg water near the critical point); the vapour leaves its minerals behind.

P7c isotopes (part one, done): `isotopes.py` and `data/isotope_fractionation.tsv`, switched on with the card's
`isotopes` (chat "with isotopes"; starting values in `isotope_inputs`, chat "fluid d18O -5", "basalt glass
87Sr/86Sr 0.7030"). At each point the carriers of O, C and S (water, DIC, methane, sulfate, sulfide, each
mineral) are in equilibrium at that temperature, and the inventory Σ n δ is carried along: reactants add
their own δ, matter that leaves takes the δ of the carriers that shrank, evaporation removes vapour at the
mean of the step's start and end (Horita and Wesolowski 1994 liquid-vapour). Oxygen: the internally
consistent database of Vho, Lanari and Rubatto (2020, CC BY), 70 minerals relative to quartz with their
quartz-water; carbon: Ohmoto and Rye (1979); sulfur: Ohmoto and Rye (1979) for sulfides, Ohmoto and Lasaga
(1982) for sulfate. Minerals without a row use a named analogue (listed in a note). 87Sr/86Sr is a tracer
carried with the fluid: starting fluid plus the Sr that dissolving rock adds (basalt glass 115 ppm, 0.7027);
minerals keep none; it leaves in proportion to the chloride when a phase separates. Checked: seawater with
100 g basalt at 350 C gives the analytic Sr mix exactly and keeps the system's δ34S at 21.0; evaporation to a
tenth at 25 C is within 1.5% of Rayleigh. Results window: an Isotopes tab; the story gives the fluid's change.
Judgment calls: H isotopes are left for part two (mineral-water D/H factors are scattered and
composition-dependent and need checking against the papers); defaults for fluids other than seawater are
illustrative and named in a note; aqueous O is all at the water's δ18O; carbonate ion is taken as
bicarbonate; HCO3- (to 290 C) and CO (to 330 C) are extrapolated above their range with a note.

P7d sorption (done): `sorption.py` and `data/sorption.tsv` (PHREEQC's phreeqc.dat EXCHANGE_SPECIES and
SURFACE_SPECIES, USGS public domain). Card `exchange_meq` (meq of sites per kg water; chat "with 50 meq of
exchange sites", "cation exchange capacity of 100 meq") and `hfo_g` (g of hydrous ferric oxide; chat "with 1 g
of ferrihydrite"). Both are solved inside the Newton system (`_Solve.sorbents`): the exchanger adds one unknown
(ln a X-, equivalent fractions add to one; species take their cation's activity coefficient, as PHREEQC's
-gamma), the surface adds ln n(SOH) for strong and weak sites and y = F psi / RT with the Gouy-Chapman charge
relation (no diffuse-layer composition, PHREEQC's default). The exchanger starts loaded from the fluid without
changing it ("equilibrate with solution"); the HFO starts bare and reacts with the fluid at the first point.
Checked against PHREEQC 3.7.3 (phreeqpy wheel, llnl.dat plus phreeqc.dat surface species) for seawater at
25 C: exchanger fractions within 5%, pH with 1 g HFO 7.231 against 7.218, potential 34.9 against 36.3 mV,
sorbed Ca, Mg, carbonate within 5%, sulfate within 10%. Results window: Sorbed tab, surface potential
property, story sentence. Judgment calls: 25 C constants at all temperatures (van 't Hoff only where
phreeqc.dat gives an enthalpy); borate and silicate surface species use B(OH)3 and SiO2(aq); following the
vapour leaves the sorbents with the fluid followed.

P7e 1D reactive transport (done): `transport.py` (class `Column`) and the card's `kind = "transport"` with
`cells`, `column_m`, `dispersivity_m`, `pore_volumes` and `pore_fluid` (empty: the column starts filled with the
fluid that flows in). The card's temperature and pressure run from the inlet (start values) to the outlet (end
values), so a column can sit on a thermal gradient; each cell is solved at its own T and P. The scheme is
PHREEQC's TRANSPORT at Courant number one: each shift moves every cell's solution (water and solutes, from
`moles_in_solution`) one cell on; minerals, gas and sorbed matter stay; dispersion mixes neighbours by
alpha/dx per shift, in ceil(3 f) mixes of at most 1/3, half before and half after the shift, with the cells
re-equilibrated after every mix (doing all the mixing in one go before reacting put the Ca front of example 11
0.24 mmol off; the split brings it within 0.003). Constant-concentration inlet, zero-gradient outlet; the column
tracks what leaves at each end (`out`, `back`), and a tracer balances exactly. Checked against PHREEQC 3.7.3
TRANSPORT (example 11 chemistry, 20 and 40 cells): Na, K, Ca, Br, Cl at the outlet within 0.006 mmol/kg at every
shift. Reactants go into every cell (per kg of pore water) and are brought to equilibrium with the pore fluid
before the flow (basalt glass is replaced by what it alters to); no evaporation, mixing, phase following,
isotopes or uncertainty bands in a column yet (the card says so); kinetics came later (P7f, below). Results: three "modes" with axes of
their own, `outlet` (against pore volumes), `start` and `profile` (the column at the start and at the end,
against distance and temperature); the results window overlays only modes on the same axes. Chat: "seawater
through a 2 m column of 100 g basalt glass", "30 cells", "3 pore volumes", "dispersivity 5 cm", "a column
filled with rainwater", "no column"; "along the column", "at the outlet" pick the view. The story gives the
outlet's breakthroughs (halfway points and peaks) and the mineral zones along the column. A 20-cell, 2 pore
volume seawater-basalt column takes about 10 s. Judgment calls: one equilibration per mix (as PHREEQC);
a cell that fails keeps its last state with its matter counted; at most 40,000 equilibrium solves per run.

P7c part two, hydrogen (done for the water): `dD` joins the isotope systems (card `isotope_inputs` key `dD`,
also "d2H", "deuterium"; chat "fluid dD -40"). The water's δD follows evaporation with the liquid-vapour D/H
fractionation of Horita and Wesolowski (1994), checked against the values they report (about 100 permil at
7 C, 27.9 at 100 C, -2 near 347 C; it changes sign near 220 C), plus mixing and phase separation. Evaporation of
rainwater to a tenth at 25 C is within 1% of Rayleigh. Judgment call: mineral-water D/H factors are not in. The
papers (Suzuoki and Epstein 1976, Graham and others 1980 and 1984, Gilg and Sheppard 1996, Chacko and others
1999) are paywalled or blocked from here, and the general Suzuoki-Epstein equation is only rough below 400 C,
so hydrous minerals take the water's δD and only the water's δD is shown. In a seawater-basalt run the minerals
hold about 1% of the hydrogen, so this moves the water's δD by well under 1 permil. H2 and CH4 are counted with
the water (in nature they are far lighter). Default δD: seawater 0, vent fluid +1, formation brine -20,
rainwater -54 (meteoric line from δ18O -8), basalt glass -70, other hydrous reactants -60; all are named in a
note. To finish: mineral factors from the papers, once Douglas can share them.

Gas-phase isotopes (Douglas's request, 2026-11-01): the gas is a set of isotope carriers in equilibrium with
the liquid. Steam (H2O(g)) takes the liquid-vapour 18O and D/H of Horita and Wesolowski (1994); CO2 gas its 18O
against water (Brenninkmeijer and others 1983) and its 13C against dissolved CO2 (Mook 1986; CO2(aq) is no longer
counted as the reference, the gas is); H2 and CH4, gas and dissolved, their D/H from water (Horibe and Craig 1995,
with δ = (1000 + δ_ref) α − 1000 because H2 is some 700 permil lighter than water when cool). Flow-through keeps
the state before venting (`PathPoint.before_vent`), so the vented gas leaves with its own δ (Rayleigh
degassing). Results show δ13C and δ18O of CO2 gas, δ18O and δD of steam, δD of H2 and methane, δ34S of H2S gas;
the story names the gas's values at the end. Judgment calls: H2S gas and dissolved H2S share δ34S (about 1
permil apart near 25 C); Mook and Brenninkmeijer are extrapolated above 60 and 100 C with a note; CH4 D/H
below 200 C is extrapolated with a note; dissolved and gas H2 (and CH4) share a ΔD. Not yet: salt effects on
liquid-vapour (Chacko and others 2001 table 4) and isotopes of the brine and vapour above the critical point
(no fractionation there in Horita and Wesolowski).

Starting cards (Douglas's request, 2026-10-26): "Start from" buttons on the scenario card and chat ("brine card",
"load the fresh water scenario") give a ready-made card (`scenario.STARTING_CARDS`, `/api/reset` with `start`):
seawater (the default, 2 to 350 C at 500 bar); fresh water (Livingstone 1963 world-average river water at the
atmosphere's 420 ppm CO2 and air O2, heated 25 to 300 C: calcite, anhydrite); brine (Great Salt Lake brine from
Drever 1997, Pitzer, evaporated 5x at 25 C: halite, anhydrite, epsomite); hydrothermal fluid (the last point of the
seawater card, 350 C and 500 bar, mixed with 9 kg of 2 C seawater: anhydrite forms and redissolves). New presets
`freshwater`, `gsl_brine`, `hydrothermal`; a fluid's analysis can now refer to a pressure (`fluid_p_bar`, preset
`measured_p_bar`), so it can be measured above 300 C. `tests/test_starting_cards.py` checks that the hydrothermal
preset still equals the end of the seawater run. Judgment calls: the pH of the fresh water and the brine is the one
that puts them at the atmosphere's CO2 (8.14 and 8.2; the brine's is not part of Drever's analysis); the brine's
silica (48 mg/kg) is left out so that talc and sepiolite do not form at 25 C.

Metals on clays and the contaminated-stream card (Douglas's request, 2026-10-26): clay minerals on the card
(`kaolinite_g`, `illite_g`, `montmorillonite_g`, g per kg of water; chat "with 100 g kaolinite and 20 g
montmorillonite") carry exchange sites at their cation exchange capacity (`data/clays.tsv`: kaolinite KGa-1b 20 and
montmorillonite SWy-1 764 meq/kg, van Olphen and Fripiat 1979 via the CMS Source Clays; illite du Puy 225 meq/kg,
Nagra NTB 23-05) with phreeqc.dat's Gaines-Thomas constants; checked against PHREEQC 3.7.3 for a metal-bearing river
water at pH 6.5: Ca, Mg, Na, K fractions within 0.5%, Zn, Cd within 6%, Pb and Cu within 30% (their carbonate and
hydroxide complexes differ between OBIGT and phreeqc.dat). Partition coefficients (`data/kd_allison_2005.tsv`, the
U.S. EPA compilation of Allison and Allison 2005, soil, bed sediment and suspended matter for As, Cd, Cr(III), Cr(VI),
Cu, Hg, Ni, Pb, Zn): `kd_solid_g` grams of solid per kg of water and `kd_medium`; `sorption.LinearSorber` holds
Kd M times every dissolved species of the metal inside the Newton system (no unknowns of its own), and the solid
starts loaded from the fluid like the exchanger. With exchange sites on the card the Kd holds only the metals without
an exchange constant (Ni, Hg, As, Cr). Results: a Partition coefficients tab (each metal's Kd on each clay and on the
Kd solid, L/kg) and a story sentence. Arsenic (H2AsO4-) and chromium (CrO4-2) are now components; the arsenate ions
come from SLOP98 because OBIGT's newer rows have no HKF parameters (`obigt.load_obigt` keeps an earlier complete
aqueous entry). New fluid preset `contaminated_stream` and starting card "Contaminated stream" (10 m streambed of
clays and bed sediment at 15 C, 25 pore volumes, about a minute): Cu and Pb break through first, then Zn and Cd;
cerussite forms at the inlet; the Kd metals barely move. Judgment calls: one set of exchange constants for all clays
(phreeqc.dat has no clay-specific ones), so per-clay Kd differ only by capacity; edge-site surface complexation on
the clays (Bradbury and Baeyens 2SPNE) is not in: the site capacities and protolysis constants are public (NTB
23-05 Table 2-1) but the metal binding constants are in papers I could not open, so the clays understate metal
uptake above about pH 6; Kd values are the report's medians (means where it gives no median) and are held at all
temperatures; a complex of two Kd metals takes the larger Kd; the stream card blocks quartz and chalcedony (silica
does not precipitate in a stream on this time scale). Next: the 2SPNE edge sites once Douglas can share Bradbury and
Baeyens (2005, 2009) or the ClaySor appendix.

Clay edge sites (2SPNE SC/CE, from the Bradbury and Baeyens papers Douglas added): `sorption.EdgeSites`, a
non-electrostatic surface with strong (2 mmol/kg) and two weak (40 mmol/kg each) amphoteric sites per kg of illite
and montmorillonite (`data/clay_edges.tsv`: illite 2009a Table 9, montmorillonite 2005 Table A1), metals on the
strong sites (`data/clay_edge_metals.tsv`). Ni and Co on illite use the fitted constants (2009a Table 11); every other
metal and clay uses the LFER (illite a 7.9, b 0.83, 2009b; montmorillonite 8.1, 0.90, 2005) with hydrolysis
constants from 2009b Table 6 for Ni and Co and from OBIGT for Zn, Cd, Pb, Cu, Mn, Fe(II), Hg. The edges start loaded
from the starting fluid (fluid unchanged), like the exchanger. New card field `clay_edges` (default true, a select
under the clays, and chat "turn off the edge sites"). Ni and Co exchange constants added from 2009a (log K 1.1,
1.3). With edges on, the Kd solid holds only As and Cr. Tests compare with measured edges: Ni and Co on illite
within 0.4 log units; Ni on montmorillonite (2005 Table 5) within 0.9. Judgment calls: weak sites carry protons only
(their metal constants come from isotherm papers not in the folder), so strong sites saturate above about 2 mmol of
metal per kg of clay; OBIGT hydrolysis for the metals other than Ni and Co (PbOH+ is 1.5 log units stronger than
Baes and Mesmer, so Pb binds strongly); one exchanger for all clays, so montmorillonite's low-pH Ni plateau is about
0.7 log units high; constants at 25 C at all temperatures; kaolinite has no 2SPNE set, so exchange only. The stream
card now holds Pb about twice as long (outlet Pb at 25 pore volumes 0.41 of the inflow against 0.80 without edges).

Also fixed: a run is marked idle only after the session is saved (a test that restarted the session right
after a run sometimes read the old file).

P7c part three, mineral δD (done, 2026-10-26, from the papers Douglas supplied): `isotopes.mineral_water_h`.
Kaolinite, dickite (and halloysite as kaolinite): Sheppard and Gilg (1996), -2.2e6/T^2 - 7.7, 0 to 350 C.
Muscovite: Suzuoki and Epstein (1976) eq. 1. Tremolite (flat -21.7 to 650 C, then -31.0e6/T^2 + 14.9), pargasite
and ferropargasite (-23.1), actinolite and ferro-actinolite (-29, one point at 400 C): Graham and others (1984).
Every other hydroxyl mineral with Al, Mg or Fe: Suzuoki and Epstein's general equation
-22.4e6/T^2 + 28.2 + 2 X_Al - 4 X_Mg - 68 X_Fe over the Al, Mg and Fe of the formula. Zeolites, salt hydrates and
minerals with S, C, P, Cl or F keep the water's δD (listed in a note). Judgment calls: below 400 C the general
equation (and muscovite's) is held at its 400 C value, because its 1/T^2 term runs to -225 permil for a clay at
25 C where measured clays sit near -30 to -60 (Sheppard and Gilg, figs. 3, 4); all Al is counted as octahedral
(at most 2 permil); epidote, chlorites, boehmite, prehnite and the like have hydrogen-bonded OH that the
general equation makes too heavy (Chacko and others 2001), and a note says so. No salt effect on D/H yet
(Horita's NaCl terms are in Chacko and others, table 4). Mineral δD is good to about 20 permil, the spread
between published calibrations.

P5b part three (done, 2026-10-26): `nacl_h2o.py` implements Driesner and Heinrich (2007) part I: critical curve
(eqs. 5, 7), halite liquidus (8), V+L+H pressure (10), liquid branch of V+L (11, with g0 and below 374 C Xcrit
fixed by the two ends) and vapour branch (13-17), halite-saturated vapour (9), NaCl vapour pressure (2, 3);
pure water's vapour pressure from IAPWS-IF97. `phasesep.phase_state` uses them for the phase compositions
(0 to 1000 C, 1 to 5000 bar); the old lookup table now only gives densities. Checked: against the table (from
xThermal) at its 1856 two-phase nodes, liquid within 1e-6 relative, vapour within 1% (95th percentile; the
table's vapour differs slightly), halite liquidus within 2e-6. Against Bischoff and Rosenbauer (1988), table 1,
the seawater (3.2 wt%) boundary from 196 to 450 C is within 3.2 bar everywhere (the grid was up to 10 bar off),
and the critical point of seawater comes out at 408.2 C and 301 bar (measured 407 C, 298.5 bar). Consequence:
seawater at 300 bar splits from 408 C, not 405 C; following the vapour near the critical point still often fails
to solve (B-dot at low density), as before. The phase-separation tests moved to 400 C/305 bar to 410 C/300 bar.

Next: salt effect on liquid-vapour and mineral D/H; Lewis and Lowell (2009) two-phase flow benchmarks once the
column can follow phases. (Kinetics in the column: done, P7f below.)

Examples (Douglas's request, 2026-10-27): a second row of buttons on the card, "Examples, fill the card"
(`scenario.EXAMPLES`, `/api/example`, `Session.load_example`). One press fills the card and the user presses Run
(Douglas, 2026-11-01: buttons fill, they never run); a note under
the buttons says what the example shows and lists the card fields that set it up (each line scrolls to its field),
and those fields are outlined on the card until another card is loaded or the note is closed. The first run opens
the results on the example's own plot (`Example.plot`, in words, read by `view.resolve_item`). Douglas's list
(same day), in his order, plus two: basalt alteration with solid solutions (seawater + 1 kg basalt glass, 350 C,
500 bar); evaporating lake (Great Salt Lake brine, Pitzer, 10x, isotopes on: halite, anhydrite, epsomite, sylvite;
d18O 0 to +13); the contaminated stream (metals on clays); boiling (an illustrative CO2-bearing NaCl geothermal water
at 150 C decompressed 25 to 5 bar, flow-through: CO2 and steam vent from about 10 bar, the water loses two thirds of
its carbon, pH 6.0 to 6.8, calcite scales); supercritical seawater (390 to 410 C at 300 bar, following the brine: 7
wt% NaCl at the end, the point at the critical point is blank); flow-through weathering (soil water at pH 5 through a
1 m, 10-cell column of albite and a little calcite, 12 pore volumes, about a minute: the calcite front leaves after
about 9 pore volumes, the top three cells weather to kaolinite); serpentinization (Douglas's addition: seawater + 1
kg olivine Fo90, 300 C, 300 bar, solid solutions on, antigorite blocked); calcite dissolving in rainwater at its
measured rate over 12 hours. For serpentinization a serpentine solid solution was added (chrysotile-greenalite,
three sites, ideal, as McCollom and Bach 2009): it takes about 4% greenalite and lowers the H2 from 0.29 to 0.17
mol/kg. Low-temperature examples block `SILICA_POLYMORPHS` and `HIGH_T_MINERALS` (lawsonite, paragonite, ...;
the database makes them stable at 25 C). Judgment calls and limits found: boiling of a salt water below the
water's own boiling curve is refused (no steam-saturated liquid), so "boiling" is gas-driven, by dissolved CO2;
the aqueous solver fails at 250 C with more than about 3 g/kg of carbonate at pH 6 (a bug to look at), so the
boiling example runs at 150 C; whole-fluid chemistry fails between 410 and 440 C at 300 bar; closed-system
weathering of albite by plain rainwater goes to pH 9 to 10 at once, so the example uses soil water; no Fe(OH)2 in
OBIGT, so brucite stays pure.
`tests/test_examples.py` checks that each example shows what its note says.

Vent mixing and catabolic energy (Douglas's request, 2026-10-27): new card fields `redox_apart` (each redox couple held
apart: H2 from O2, H2S from sulfate, CH4 from carbonate, Fe+3 from Fe+2, NO3- from NH3; `equilibrium.REDOX_APART`),
`minerals_form` (false: every mineral suppressed), `energy` (catabolic energies, `catabolic.py`) and `mix_spacing`
("temperature": steps even in the mixture's temperature, `Scenario.fractions`). A held species now carries its whole
redox family (`equilibrium._tie_family`: HS-, metal bisulfides and pyrrhotite with H2S, Fe(III) complexes and hematite
with Fe+3), and species in between (thiosulfate, polysulfides, organic acids and their complexes, nitrite, N2, pyrite,
magnetite; `redox_intermediates`) are suppressed, in the starting fluids too; `path.bulk_over` re-counts an
equilibrium speciation over the held columns. Energies: dG = RT ln(Q/K) of seven reactions (H2, sulfide, CH4 and
Fe(II) oxidation by O2, Fe to hematite; methanogenesis; sulfate reduction by H2; anaerobic methane oxidation) and
-dG times the limiting reactant's total, J per kg of fluid and per kg of the starting fluid (McCollom and Shock 1997;
Amend and others 2011); an Energy tab in the results and a story sentence. New preset `deep_seawater` (2 C, 0.11
mmol/kg O2, 2.3 mmol/kg carbonate, pH 7.8; illustrative). Example "Vent fluid meets seawater": the vent preset plus
0.08 mmol/kg CH4 at 350 C, 250 bar, mixed to 5 C (115 kg per kg), 70 steps, 2 s; sulfate reduction leads above about
285 C, then Fe, H2 and, below about 90 C, sulfide oxidation (up to 43 J/kg). Judgment calls: dG held at
its starting value (upper bound, as in the papers); a reactant or product below 1e-10 mol/kg leaves dG blank;
mixed-valence minerals such as magnetite are suppressed when redox is held apart.

462 tests pass with the data folder present. The data folder needs `OBIGT/Berman/` (run
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
   ~~solid solutions~~ (P7b), isotopes (O, C, S, Sr done in P7c; H for the water, minerals still to do), ~~sorption~~ (P7d), ~~1D transport~~ (P7e).

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

Activity models (2026-11-11, decided by Douglas): B-dot (Helgeson 1969) stays the default; b-gamma (Helgeson,
Kirkham and Flowers 1981, Table 27 via CHNOSZ) is a card option and is suggested above 300 C or 1 kbar
(b-gamma made seawater carbonate saturation worse than B-dot at 25 C, so it was not made the default). Pitzer
uses pitzer.dat's own log K and molar volumes for the species and minerals it lists, checked against PHREEQC
3.7.3; with them the brine card gives halite then gypsum (not anhydrite and epsomite as with HKF log K).
Holland and Powell (2011) minerals were not taken up.

## Where the code is

`src/reaction/` : `database/` (data0 reader), `thermo/` (IAPWS-95 water, dielectric, HKF, crystals, OBIGT loader,
`eosdb.py` logK at any T and P, `validity.py` reliability check), `network.py`, `aqueous.py`, `speciation.py`,
`equilibrium.py` (closed-system solver), `path.py`, `validation/` (EQ3NR and EQ6 harnesses), `app/` (scenario,
engine, assistant, view, narrate, session, server, `static/` pages). Tests are in `tests/`; the pressure checks are
`tests/test_pressure.py`.

P7f kinetics in the column (Douglas's request, 2026-10-27: forams dissolving in a column of seawater). The card's
kinetics fields now work in a column: `duration` is how long the flow lasts (the pore volumes pass in that time),
so each shift is `duration / shifts`. Reactants with a rate law are not put into the fill: each cell keeps its own
amount of each (`Column.left`) and after every shift reacts for the shift's time through a one-cell
`ReactionPath` (the static run's implicit steps with error control; `run` takes `guess` and `step_hint` for
this), the rest of the cell at equilibrium; reactants without a rate are brought to equilibrium in the fill as
before. A kinetic mineral and its other structural states never precipitate at equilibrium; it grows on what is
left of it when supersaturated (operator splitting as PHREEQC's TRANSPORT with KINETICS). Results: the Rates tab
works on all three column views (rate and % left per cell, and the kinetic mineral counts in the Minerals tab);
the outlet also gets a Time axis; a note gives the residence time, velocity and time per shift. Checked: with no
dispersion the first cell's first shift equals a batch kinetic run of that length (within 2%); calcium balances
exactly (dissolved = in the column + through the ends). Judgment calls: `kinetic_tol` 0.05 in the column (0.03
in a static run): about 1-2% low on what dissolves at the inlet of a foram-ooze column against 0.03, at about half
the time; the column budget counts a kinetic cell-shift as 15 solves. Chat: "forams", "foraminifera", "foram
calcite", "coccoliths", "chalk" are calcite.

The foram example (Douglas, same day) became one foram sinking below the CCD, as a flow-through kinetic path rather
than a column: the run follows the foram (20 µg of calcite in 1 kg of the water around it, so the water hardly
changes), 4.5 to 5.5 km (455 to 556 bar) at 1.5 C in 5 days (about 200 m a day). It loses about a fifth of its
mass; the example opens on calcite in µg/kg (new plot unit `ug`; `Example.unit`). Reactant amounts take µg ("20
ug forams"); amounts below a mg are shown in µg (`scenario.reactant_text`). The model's bottom water is further
from calcite saturation (SI about -0.4) than measured (Omega about 0.7-0.9) and the Palandri-Kharaka law is linear
in (1 - Omega), so it is an upper bound on dissolution; a seawater-specific calcite law (Subhas et al. 2015,
(1 - Omega)^n) is the next refinement.

Organic acids (2026-10-27): OBIGT's SLOP metal-organic complexes (Na(Pent), Ca(Prop)+ ...) were read without
their free ligands and written over carbonate + O2, so any reduced fluid (the vent preset at 25 C) put most of its
carbonate into them: metastable equilibrium that does not happen on real time scales. Now organic carbon enters
only as the organic acids an analysis gives (acetate, formate, propanoate, butanoate, pentanoate, glycolate,
lactate, oxalate, malonate, succinate; `network.ORGANIC_LIGANDS`, names and synonyms in `analysis.SOLUTES`).
Each is its own component (decoupled), tied to its acid(s) and metal complexes by name (`tie_family`, or
`via_ligand` in `Network` when O2 is not in the basis); every other organic species is suppressed
(`organic_held_apart`). CH4 stays as before. Path runs warm-start their first point from the starting
fluid's speciation (cold start as fallback), which the vent fluid at 25 C needs now.

Rock cards (Douglas, 2026-11-01): a "Rock" row above the reactants (`rocks.ROCK_CARDS`, `rock_cards()` in the
snapshot) with a grams box; a button replaces the reactants with that rock. Basalt glass stays one congruent glass
reactant; the crystalline rocks (`rocks.MINERAL_ROCKS`: crystalline basalt, granite, lherzolite, Mars meteorite as
Shergotty, sandstone, calcitic/dolomitic/aragonitic limestone) are modes that `scenario.parse_reactants` expands into
their minerals, so typing or saying "100 g granite" does the same. Fe-Ti oxides, sulfides and phosphate are left out
(no redox dependence for a rock card). MAX_REACTANTS raised to 12. New fluids: groundwater (granitic), karst spring.

Foram kinetics on a porous surface (Douglas, 2026-11-10). A card can now use seawater carbonate rates
(`carbonate_rates: seawater`): calcite after Naviaux et al. (2019, GCA 246, 363-384; the 5 C BET-normalized fits
either side of Omega 0.75, Ea 25 kJ/mol) and aragonite after Dong et al. (2019, EPSL 515, 1-12; 0.013 g/g/day
(1 - Omega)^1.37 at 5 C on 2.46 m2/g, about 30 kJ/mol). These dissolve only (zero above saturation). A reactant can
carry its own BET surface ("20 ug forams at 3 m2/g", `ssa_m2_g`); with the seawater rates forams and pteropods
default to 2 m2/g (pteropods 1.6-2.2 m2/g in Dong et al. Table 2; the foram value is a placeholder until a
measured one is found, e.g. in Subhas et al. 2018, Mar. Chem. 205). `depth_axis` plots a kinetic run against
depth (about 10 m per bar). The foram example now falls 1 to 5.5 km in 22.5 days through water tuned (pH 7.98
at 25 C) so calcite saturates near 2.7 km and the Omega 0.75 kink, the lysocline, sits near 4 km; aragonite is
blocked there (the shallow water is supersaturated with it). New example "shells" sinks a foram and a pteropod
together. With these two laws the pteropod starts dissolving ~700 m shallower but loses far less by 5.5 km,
because Naviaux's far-from-equilibrium calcite fit is much faster per m2 than Dong's aragonite law; the
database puts aragonite only ~0.07 log units above calcite in solubility at 2 C (seawater is nearer 0.18).
Shell surface evolves (Douglas, same day): `Reactant.area_exponent` (A = A0 (m/m0)^p; 2/3 for grains as before);
card `shell_area_exponent` (default 0) for forams and pteropods with seawater rates, so their surface per gram goes as
m0/m. The foram example ends at 4.4 µg instead of 8.2 µg.
Aragonite rates (Douglas, same day: "very little of the aragonite dissolved while all of the calcite did; that
does not seem correct"): with Dong's 0.013 g/g/day law aragonite is ~150x slower per m2 than Naviaux's calcite
fit at Omega 0.5. Default `seawater` now applies Naviaux's fits to aragonite on its own Omega (Dong et al. find the
same mechanisms, Omega_crit 0.9 vs 0.87, and similar dissolution relative to abundance); Dong's law is
`seawater_dong`. Shells example: pteropod starts near 2 km, half gone 4.4 km, gone by 5 km; foram half gone ~5 km.
Still open: database aragonite/calcite solubility gap 0.07 log at 2 C vs ~0.18 in seawater (Mucci 1983).
Douglas's sinking-test review (2026-11-10): benchmark script scripts/seawater_omega_benchmark.py (PyCO2SYS);
Reaction's Omega_calcite is 0.06-0.16 log low and its calcite/aragonite ratio 1.2 vs Mucci's 1.5; the pressure
effect is ~0.09 log too strong at 400 bar. Fix pending Douglas's choice (seawater constants vs HKF refit).
Calcite seawater regimes joined by a smooth maximum (s = 10; they cross at Omega 0.79). Card `sinking`: progress is
depth, each step lasts dz / w(W) with Takahashi and Be (1984) w = 10^2.06 W^0.64 m/day (Reactant.speed); the foram
falls 1-5.5 km in ~6 days and keeps ~60%. Not yet: flush with a depth profile of DIC/TA (needs the Omega decision).
Marine carbonate system (Douglas, 2026-11-10: "if we are in a carbonate, marine environment and PyCO2SYS is well
demonstrated, then obviously that should be used"): src/reaction/marine.py (Lueker K1/K2, Dickson KB, Millero Kw
and pressure terms, Mucci Ksp, Uppstrom boron), within 0.1% of PyCO2SYS 1.8 (tests/test_marine.py). With seawater
rates the calcite/aragonite rate laws and the reported saturation use it; TA/DIC/S from the card (alkalinity_umol_kg,
dic_umol_kg, salinity) or from the water (engine.marine_water). Next: a measured depth profile of TA/DIC/T/S
(GLODAP or BATS; those hosts are blocked from the container) and flush through it.
Defaults (Douglas, 2026-11-10: "the default should be as close to reality as possible"): solid_solutions and
isotopes on. A column run leaves isotopes out with a note instead of refusing. Tests that named pure endmembers
(clinochlore, tremolite, chrysotile) accept the solid solutions too.
Flush mode (Douglas, 2026-11-11: "We want the foram test to experience the water chemistry and T/P changes as it
sinks ... We want a flush model"). Card `mode: flush` (`flush_fraction`, default 1): Bethke's flush model as in GWB
React (Reaction Modeling Guide, Release 19, s3.3: a fluid "displaces existing fluid from the system", an equal mass
of water). In `path.py` the solids stay and, after each point, `flush` of the fluid (bulk less minerals and solid
solutions) is replaced by the starting fluid; the point shows the fluid the step reacted with. Not yet with sorption
(refused) or isotopes (left out with a note); not for columns. The foram and shells examples now use flush. Why: in
flow-through (GWB s3.2, EQ6) only precipitates leave and one 1 kg parcel keeps every solute the shells give off, the
wrong process for a shell falling through an ocean. Results hardly change for one 20 ug shell per kg, because the
seawater rates take Omega from the card's fixed TA/DIC at each T and P (engine.marine_log_omega), but the shells
example can no longer let the pteropod's products reach the foram. Before this, the shells example's apparently
deeper calcite lysocline was the plot's y axis (it started at zero only below 0.35 x max, so with aragonite going to
0 the calcite loss looked small); mineral plots now always start at zero (plot.js `zero`, Streamlit tozero).
Sinking runs use kinetic_tol 3e-4 (engine.SINKING_KINETIC_TOL) instead of 0.03: foram end mass 14.295 (0.03),
14.493 (3e-3), 14.532 (1e-3), 14.555 (3e-4), 14.567 ug (1e-4); the step error shrinks only about as tol^0.5, so a
higher-order kinetic step would be the better fix. Still to do: a measured TA/DIC/T/S depth profile for the
flushing water (GLODAP or BATS), which flush now makes possible.
Water profile (Douglas, 2026-11-11: "Add the depth profile card"). Card `water_profile`: rows [depth m, T C,
salinity, TA umol/kg, DIC umol/kg] by increasing depth, the pressure read as depth (scenario.depth_m, moved there
from engine). With it the path temperature follows the profile (Scenario.tp) and the marine carbonate system the
seawater rates and reported saturation use takes S, TA and DIC at each depth (engine.marine_log_omega takes a
function of pressure). New example "ceara": GLODAPv2.2023 (Lauvset et al. 2024, ESSD 16, 2047) Atlantic file, cruise
323019940104 station 198 (6.543 N, 42.994 W, March 1994), 31 good-flag samples; calcite horizon ~4.3 km, the foram
keeps its mass to there and ends 19.96 ug at 4.6 km. Extracted files on Douglas's PC: workspace\ceara\. GLODAPv3.zip
there not used yet (no documentation read). A 20 ug test resting in bottom water at 4358 m (ODP 929, Omega 0.957)
is gone within a year with the Naviaux near-equilibrium rate and 2 m2/g: burial and pore water must control what the
cores keep; that rate at Omega 0.9-0.99 is not yet checked against the paper. Next: sediment-pile plan (Claude
Doc "Sediment pile model card: plan"): burial path, pore-water profile by sediment depth, respiration, diagenetic
column, calibration against grayscale and test weights.
Seawater rates now the in-situ lower bound (Douglas, 2026-11-11, "let's go with your recommendation"): Naviaux et al.
(2019), Mar. Chem. 215, 103684 (PDF in the Reaction folder, 10299274.pdf), abstract and Section 5: 10^-14.3
(1-Omega)^0.11 mol cm-2 s-1 for 0.8 < Omega < 1 and 10^-10.8 (1-Omega)^4.7 below, recommended for water-column
dissolution of particles from high-DOC surface water, ~4x slower than their lab fits. The old GCA 246 lab fits
(-13.07/0.34, -10.01/4.81, kink 0.75) are `carbonate_rates: seawater_lab` (not re-checked against GCA 246: that paper
is not in hand). Ea stays 25 kJ/mol (GCA). Caveat: their in-situ Omega is from Alk and pH, 5-10% below Omega(Alk, DIC)
which Reaction computes. Mar. Chem. Table 1 lab 5 C reads -13.1 n 0.11 (Omega > 0.8), the graphical abstract's lab
line ~-13.4: not reconciled. New results: forams 19.11 ug at 5.5 km (95.5%), shells pteropod 9.8 ug at 5.5 km, ceara
19.995 ug at 4.6 km; resting at ODP 929 depth (Omega 0.957): 16.1 ug after 1 yr, 8.2 after 3, gone by 10.
Ceara card (2026-11-11, Douglas: "add pteropods to any Ceara Rise model"): now "Ceara Rise: foram and pteropod
sinking", forams + pteropods 20 ug each, aragonite unblocked. Run: calcite horizon 4297 m, foram 19.995 ug at 4.6 km
(0.02%; the old card text said 0.2%, which was wrong); aragonite saturated only from just below 1 km to ~2.7 km,
pteropod 19.70 ug at 4.6 km (1.5%).
Burial (sediment pile, step 1 of the plan): `burial: true` with kinetics. The test rests at the seafloor (start
pressure = water depth) and is buried at `sedimentation_cm_kyr` to `burial_depth_cm`; run time = depth / rate; pore
pressure hydrostatic; T and S of the water profile at the seafloor; TA and DIC from `porewater_profile` rows
[z cm, TA, DIC umol/kg], else bottom water all the way down (a bound, not a prediction). Axes: depth below the
seafloor (cm) and time since burial (kyr). No bioturbation, compaction, respiration or diagenetic column yet; the
pore-water profile has no field in the web form (card JSON or the assistant only). No Ceara burial example: the
sedimentation rate and the pore-water profile (Martin and Sayles) are not in hand. Test result with bottom water at
4.4 km (Omega calcite 0.94): both shells gone in the first kyr step; the in-situ rate dissolves a 20 ug test in
years, so preservation turns on the pore water in the top cm, not on the burial rate.
Sediment column (Douglas, 2026-11-11: "the tests will be surrounded by other tests and other sediment that will
fix the pore fluid chemistry"; chose the column over a box): porewater.py, the steady diffusion-reaction model of
Hales and Emerson (1997, GCA 61, 501-514; PDF 1-s2.0-S0016703796003663-main.pdf in the Reaction folder):
R = k (1-phi)[C](1-Omega)^4.5, k 25-300 %/day (best 125) at Ceara Rise Sta. G; CO2 = beta x O2 consumption, O2
flux 10-25 umol/cm2/yr, ~70% within a 3 mm scale; CaCO3 ~65% (A, B), 35-50% (G); F 2.5-3, phi = F^(-1/2.3).
Card: sediment_column + porosity, caco3_wt_pct, sediment_k_pct_day, sediment_n, o2_flux_umol_cm2_yr, co2_per_o2,
respiration_depth_cm (no defaults). Checks: grid-converged, DIC mass balance closes to 0.2%; dissolution flux at
4000 m (Sta. A values) 6.8 and 4685 m (Sta. G) 20.8 umol/cm2/yr vs their Table 2 (beta 0.64) 6-9 and 8-23; our
bottom-water Omega is GLODAP's (1.12 at 4000 m) not theirs (1.02). Buried shells take the sediment's k per mass
(an assumption; their calcite is mostly 8-64 um). Result at Sta. A values, 2 cm/kyr (test value): foram 20 ->
11.5 ug by 10 cm (first reported 12.8: marine_log_omega rounded the pressure to 0.001 bar = 1 cm of burial; now 1e-6
bar, and the buried shells read a precomputed log Omega(z) of the column), all loss in the top 2 cm; pteropod gone; at Sta. G values both gone in the top cm. Missing:
bioturbation (tests held at their depth; real tests are mixed through the mixed layer and leave the corrosive
top mm fast), single-exponential respiration, sedimentation rate of the cores, Martin and Sayles (1996).
From memory, unchecked: Boudreau D0 fit and tortuosity, Wolf-Gladrow alkalinity 17/106.
Sedimentation rate bracket (2026-11-11, Douglas: "bracket the likely sedimentation rate"): ODP Leg 154 Initial
Reports site chapters: Pleistocene averages 33 m/m.y. (925, 3040 m; 926, 3598 m), 35-43 (927, 3315 m), 28-33 (928,
4012 m), up to 40 (929, 4356 m), i.e. 2.8-4.3 cm/kyr including glacials (upper bound). Calcite mass balance with
Martin and Sayles' rain of 16.5 umol/cm2/yr (as quoted by Hales and Emerson) minus Hales and Emerson's dissolution
(beta 0.64), burial = w (1-phi) rho_s f / M: Sta. A 1.1-1.8, B 0.7-2.0, C 1.5-2.0, G <= 2.2 cm/kyr (phi 0.62-0.67,
rho_s 2.7 assumed). Working bracket 1-4 cm/kyr, most likely 1-2. Holocene 14C rates for the cores not found.
Bioturbation: porewater.mixed_mass_left, biodiffusion D_b over a mixed layer L with burial; card mixed_layer_cm +
biodiffusion_cm2_yr; output kinetic "mixed:<mineral>" (% left). Limits checked: no mixing exp(-I/w), fast mixing
1/(1 + I/w), I = int lam dz. At D_b >= ~0.05 cm2/yr the result is the fast-mixing limit whatever L (5-10 cm); Teal
et al. 2008 (Aquat. Biol. 2, 207): global mean L 5.75 +- 5.67 cm, D_b decreases with water depth. Calcite left below
the mixed layer (k 125%/day, n 4.5): Sta. A values 47/64/78% at 1/2/4 cm/kyr, Sta. G values 15/27/42%; no mixing
33/57/76% and 0.5/7/26%. Aragonite 0% in every case. BURIAL_KINETIC_TOL 3e-3 (1.3% from exact).

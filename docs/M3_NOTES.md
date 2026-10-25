# M3 working notes: reaction paths

Written 2026-10-15 before any M3 code; the solver and harness described here were committed in 3af9559.
Current status is in NEXT.md.

## The oracle: EQ6 test library (`6tlib_cmp`)

Needs the EQ3/6 files (`scripts/fetch_references.ps1` re-downloads them; the copies from the
first session lived in a temporary scratch folder and may be gone).  21 cases, each with
`name.6i` (input), `Outputs/name.6o` (text output), `.6p` (pickup) and `Outputs/name.csv`
(tables at every print point: this is the machine-readable oracle).

| group | cases |
|---|---|
| kinetics, belongs to M4 | crisqtz, pptcal, pptqtz, pptqtza, j13wtuff |
| solid solutions, out of scope | pptmins, j13wtuff |
| equilibrium paths, the M3 target | micro, microft, methane, pyrsw, dedolo, heatsw, heatswfl, heatqf, j13wtitr, rwtitr, swtitr, swxrca, swxrcaft, rwssdiag, j13wsf |

What each exercises: micro (microcline into pH 4 HCl, closed, ends when saturated),
microft (same, flow-through), methane (redox swing, hematite + calcite present), pyrsw
(pyrite into seawater, step change in fO2), dedolo (gypsum into carbonate water with T
rising linearly to 75 C), heatsw (seawater 25 -> 350 C closed), heatswfl (same,
flow-through), heatqf (jump to in-situ T), *titr (special reactant = titrant fluid),
j13wsf (fixed O2 and CO2 fugacities, spent fuel).

### CSV tables (rows are print points, first column Xi)
B1 Xi, t, T, P, pH, pmH, log fO2, Eh, pe, aw | B2 H2O(kg), Alk, I, TDS, density ... |
C1 dissolved element molalities | D1 basis-species total molalities | E1 log molality of
basis species (log f for O2(g)) | E2 log activities | J moles of reactants destroyed |
K affinities of reactants | P moles of product minerals | Q saturation indices of potential
phases | T fugacities | W mass and volume changes.  Missing values are -99999.

### `.6i` contents that matter
* `Mass Balance Totals (moles)`: the moles of every **component** (H2O, H+, each element's
  basis species, O2(g) for redox) in the equilibrium system and in the aqueous solution
  (their difference is the solids in the system).  This is the initial bulk composition B_j.
* Reactants: type (Pure mineral, Special reactant, Aqueous species, Solid solution), moles
  remaining, and for special reactants a `Reaction` block giving coefficients over the basis
  species including H2O and O2(g): directly usable as a composition vector.  Rate laws:
  forward "Relative rate equation" dXi_r/dXi = rkb (constants), backward "Partial equilibrium".
* Temperature: constant, linear in Xi, linear in time, fluid-mixing tracking.
* iopt(1): closed / titration / fluid-centered flow-through; iopt(15) = no redox.
* Alter/Suppress (nxmod), mineral subset suppression (nxopt), fixed fugacity options
  (nffg: gas, moles to add, log fugacity), Xi range and print intervals.

Observation: EQ6 output at a given Xi in a closed system is an equilibrium state, so the
harness can solve each CSV row independently (feeding T from the CSV).  Flow-through
cases are path dependent and need stepping.

## Design decided so far

Joint Newton, PHREEQC-style, for a closed system with bulk component totals B_j:

* unknowns: u_j = ln a_j of the basis species (not H2O), ln W (kg solvent), n_p moles of
  each mineral in the assemblage;
* equations: component balance for every component j (incl. H2O and O2(g)):
  W * tot_j(u) + sum_p nu_pj n_p = B_j ; log Q/K = 0 for each mineral in the assemblage;
* fixed-fugacity gases / buffers replace the balance of one component (O2(g): fix u_O2;
  CO2(g): replaces the HCO3- balance);
* active-set loop for the assemblage: drop minerals with n_p <= 0, add the most
  supersaturated mineral, swap one-for-one if the Jacobian goes singular (phase rule);
* reuse log residuals for components no mineral involves, scaled linear residuals
  otherwise; keep the pre-Newton successive substitution and the LM globalization
  already in `speciation.py`.
* reactants are composition vectors: a mineral's vector is its row in `Network.nu_phase`;
  special reactants/fluids come from a Speciation (component totals per kg).
* temperature path T(xi) just re-evaluates logK each step (works for the EosDatabase to 400 C).
* flow-through: after each step remove the new minerals' components from B and clear the assemblage.

Refactors needed first:
1. split `speciate()` into a reusable class (network + logK built once per T, `solve()` many
   times with warm starts); keep all 148 tests green.
2. expose component totals for H2O and O2(g) in `Speciation` (not only solute basis species).

Then: `equilibrium.py` (closed-system solver), `path.py` (reactants, T path, buffers,
flow-through, adaptive step with phase-boundary refinement), `.6i`/CSV readers in
`validation/`, harness comparing pH, log fO2, aw, element molalities, mineral moles and
saturation indices at the CSV's Xi points; then the seawater-basalt example (user's goal).

"""Catabolic energy: the Gibbs energy that chemolithotrophic microbes can gain from a fluid out of redox
equilibrium, as where a reduced vent fluid mixes with oxygenated seawater (McCollom and Shock 1997;
McCollom 2007; Amend and others 2011).

For each reaction, at the fluid's temperature, pressure and activities,

    dG_r = RT ln(Q / K)      (kJ per mole of the reaction as written)

and, where dG_r is negative, the energy available per kg of fluid is -dG_r times the moles of reaction
that the limiting reactant allows (its total dissolved amount over its coefficient).  That is the energy
the fluid holds if the reaction ran to completion with dG_r held at its starting value, as in those
papers; dG_r in fact shrinks as the reaction proceeds.  log K comes from the database's Gibbs energies at
every temperature and pressure, like every other log K in the program.

The numbers mean something only when the fluid is out of redox equilibrium: in a run with the redox
couples held apart (see ``equilibrium.REDOX_APART``).  At redox equilibrium every dG_r is zero.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import KELVIN, LN10, R
from .database.model import Kind, ThermoDatabase, ThermoRangeError
from .speciation import Speciation

ABSENT = 1e-10                 # mol/kg; a reactant or product below this is taken as not there


@dataclass(frozen=True)
class Catabolic:
    key: str
    label: str
    equation: str
    reaction: dict[str, float]          # species: coefficient, reactants negative


# Aerobic and anaerobic reactions of the vent-mixing papers.  Iron oxidation makes hematite, as solid
# Fe(III), since dissolved Fe+3 is all but absent in seawater and would make dG_r meaningless (goethite has
# no data above 100 C).
REACTIONS: tuple[Catabolic, ...] = (
    Catabolic("h2_ox", "Hydrogen oxidation", "2 H₂ + O₂ → 2 H₂O", {"H2": -2, "O2": -1, "H2O": 2}),
    Catabolic("h2s_ox", "Sulfide oxidation", "H₂S + 2 O₂ → SO₄²⁻ + 2 H⁺",
              {"H2S": -1, "O2": -2, "SO4-2": 1, "H+": 2}),
    Catabolic("ch4_ox", "Methane oxidation", "CH₄ + 2 O₂ → CO₂ + 2 H₂O", {"CH4": -1, "O2": -2, "CO2": 1, "H2O": 2}),
    Catabolic("fe_ox", "Iron oxidation", "4 Fe²⁺ + O₂ + 4 H₂O → 2 Fe₂O₃ + 8 H⁺",
              {"Fe+2": -4, "O2": -1, "H2O": -4, "hematite": 2, "H+": 8}),
    Catabolic("methanogenesis", "Methanogenesis", "CO₂ + 4 H₂ → CH₄ + 2 H₂O", {"CO2": -1, "H2": -4, "CH4": 1, "H2O": 2}),
    Catabolic("sulfate_red", "Sulfate reduction", "SO₄²⁻ + 2 H⁺ + 4 H₂ → H₂S + 4 H₂O",
              {"SO4-2": -1, "H+": -2, "H2": -4, "H2S": 1, "H2O": 4}),
    Catabolic("aom", "Anaerobic methane oxidation", "CH₄ + SO₄²⁻ → HCO₃⁻ + HS⁻ + H₂O",
              {"CH4": -1, "SO4-2": -1, "HCO3-": 1, "HS-": 1, "H2O": 1}),
)

# what limits each dissolved species: the total of its component (all sulfide, all carbonate, ...), or
# for O2 the species itself
SUPPLY = {"H2": "H2", "H2S": "H2S", "HS-": "H2S", "CH4": "CH4", "Fe+2": "Fe+2", "SO4-2": "SO4-2",
          "CO2": "HCO3-", "HCO3-": "HCO3-", "O2": None}


def log_k(db: ThermoDatabase, reaction: dict[str, float], t_c: float, p_bar: float | None) -> float:
    """log K of ``reaction`` from the formation reactions of its species (all written over the strict
    basis in this database), which must balance."""
    total, basis = 0.0, {}
    for name, c in reaction.items():
        sp = db[name]
        if sp.kind == Kind.BASIS:
            basis[name] = basis.get(name, 0.0) + c
            continue
        # sp.reaction: -1 for the species itself; log a_sp = sum nu log a_basis - log K_sp
        total -= c * sp.logk(t_c, p_bar)
        for k, nu in sp.reaction.items():
            if k != name:
                basis[k] = basis.get(k, 0.0) + c * nu
    off = {k: v for k, v in basis.items() if abs(v) > 1e-9}
    if off:
        raise ValueError(f"reaction does not balance: {off}")
    return total


def supply(spec: Speciation, name: str) -> float:
    comp = SUPPLY.get(name)
    if comp is not None and comp in spec.component_totals:
        return spec.component_totals[comp]
    st = spec.species.get(name)
    return st.molality if st is not None else 0.0


def energies(db: ThermoDatabase, spec: Speciation, water_kg: float = 1.0,
             reactions: tuple[Catabolic, ...] = REACTIONS) -> dict[str, float | None]:
    """For each reaction: "dG:<key>" in kJ per mole of reaction, "E:<key>" the energy available in J per
    kg of water of this fluid, and "Es:<key>" the same per kg of the starting fluid's water (``water_kg``
    kg of water per kg of it: a mixture's).  dG is None where a reactant or product is absent; E is 0
    where a reactant is absent or dG is not negative."""
    t_c, p_bar = spec.temperature_c, spec.pressure_bar
    rt = R * 1e-3 * (t_c + KELVIN)            # kJ/mol
    out: dict[str, float | None] = {}
    for r in reactions:
        reactants = [n for n, c in r.reaction.items() if c < 0 and n in SUPPLY]
        products = [n for n, c in r.reaction.items() if c > 0 and n in SUPPLY]
        have = {n: supply(spec, n) for n in reactants + products}
        dg = energy = None
        if any(have[n] < ABSENT for n in reactants):
            energy = 0.0
        elif all(have[n] >= ABSENT for n in products):
            log_q = 0.0
            for n, c in r.reaction.items():
                if n == "H2O":
                    log_q += c * spec.log_a_water
                elif db[n].kind == Kind.MINERAL:
                    continue
                else:
                    log_q += c * spec.species[n].log_activity
            try:
                dg = rt * LN10 * (log_q - log_k(db, r.reaction, t_c, p_bar))
            except ThermoRangeError:
                dg = None
        if dg is not None:
            moles = min(have[n] / -r.reaction[n] for n in reactants)
            energy = -dg * 1000.0 * moles if dg < 0 else 0.0
        out[f"dG:{r.key}"] = dg
        out[f"E:{r.key}"] = energy
        out[f"Es:{r.key}"] = None if energy is None else energy * water_kg
    return out


def available(values: list[float | None]) -> float:
    """The largest finite value of a series (0 if none)."""
    return max((v for v in values if v is not None and math.isfinite(v)), default=0.0)

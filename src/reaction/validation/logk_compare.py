"""Compare logK computed from equations of state with the grids of an EQ3/6 database.

For each EQ3/6 species whose dissociation reaction, translated to the spelling of
the equation-of-state database, has the *same stoichiometry* as the generated one,
the two logK values are compared on the EQ3/6 temperature grid (0-300 C at the
1.0132 bar / saturation-pressure curve).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from ..database.model import Kind, ThermoDatabase, ThermoRangeError
from ..thermo.eosdb import EosDatabase, eq36_name_to_obigt


@dataclass
class LogKComparison:
    name: str
    eos_name: str
    differences: list[float | None]      # ours - EQ3/6 at each grid temperature

    @property
    def max_abs(self) -> float:
        return max((abs(d) for d in self.differences if d is not None), default=0.0)


def _translate(name: str, kind: Kind) -> str:
    if kind in (Kind.MINERAL, Kind.LIQUID):
        return name.lower()
    if kind == Kind.GAS:
        return name
    if name == "O2(g)":
        return "O2(g)"
    return eq36_name_to_obigt(name)


def same_reaction(eq36: ThermoDatabase, eos: EosDatabase, name: str) -> str | None:
    """Name of the equation-of-state species with the same reaction as ``name``, or None."""
    sp = eq36.species[name]
    if sp.logk is None or not sp.reaction:
        return None
    candidates = [_translate(name, sp.kind)]
    if sp.kind == Kind.MINERAL:
        candidates.append(name.replace("-", "").lower())
    for cand in candidates:
        other = eos.species.get(cand)
        if other is None or other.kind != sp.kind and sp.kind != Kind.AUX_BASIS:
            continue
        translated = {}
        ok = True
        for k, c in sp.reaction.items():
            if k == name:
                translated[cand] = c
                continue
            kind = eq36.species[k].kind if k in eq36.species else Kind.AQUEOUS
            translated[_translate(k, kind) if k not in ("H2O", "H+", "O2(g)") else k] = c
        if set(translated) == set(other.reaction) and all(
                abs(translated[k] - other.reaction[k]) < 1e-6 for k in translated):
            return cand
    return None


def compare_logk(eq36: ThermoDatabase, eos: EosDatabase, names: list[str] | None = None
                 ) -> list[LogKComparison]:
    """Compare logK on the EQ3/6 grid for every species with an identical reaction."""
    out = []
    pool = names if names is not None else list(eq36.species)
    for name in pool:
        if name not in eq36.species:
            continue
        sp = eq36.species[name]
        if sp.kind == Kind.BASIS:
            continue
        cand = same_reaction(eq36, eos, name)
        if cand is None:
            continue
        diffs: list[float | None] = []
        for t, p, ref in zip(eq36.temps_c, eq36.pressures_bar, sp.logk.values):
            if math.isnan(ref):
                diffs.append(None)
                continue
            try:
                diffs.append(eos.species[cand].logk(float(t), float(p)) - float(ref))
            except ThermoRangeError:
                diffs.append(None)
        out.append(LogKComparison(name, cand, diffs))
    return out

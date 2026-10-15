"""Run an EQ3NR input file through Reaction and compare with the stored output."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from ..constants import KELVIN, LN10
from ..database.model import Kind, ThermoDatabase
from ..speciation import (Couple, Equilibrium, LogActivity, LogFO2, MgPerKg, MgPerL,
                          Molal, PH, Pe, Speciation, speciate)
from .eq3nr import Eq3nrReference, parse_eq3nr_outputs
from .eq3nr_input import Eq3nrInput, read_eq3nr_inputs

# EQ3/6 converts Eh to pe with R = 1.98726 cal/(mol K) and F = 23062.3 cal/(V eq);
# Reaction itself uses CODATA values.  Validation uses EQ3/6's so inputs match.
_LEGACY_R = 1.98726
_LEGACY_F = 23062.3

_NO_OP_UNITS = {"Make non-basis", "Homo. equil."}


class Unsupported(RuntimeError):
    """The input uses an EQ3NR feature Reaction does not implement yet."""


def to_speciate_args(db: ThermoDatabase, inp: Eq3nrInput) -> dict:
    if inp.activity_model not in ("bdot", "davies"):
        raise Unsupported(f"activity model {inp.activity_model}")
    if inp.pressure[0] != "reference":
        raise Unsupported(f"pressure option {inp.pressure[0]}")

    basis: dict = {}
    suppress: list[str] = []
    suppress_phases: list[str] = []
    for row in inp.rows:
        u = row.unit
        if u == "Molality":
            basis[row.name] = Molal(row.value)
        elif u == "mg/kg.sol":
            basis[row.name] = MgPerKg(row.value)
        elif u == "mg/L":
            basis[row.name] = MgPerL(row.value)
        elif u == "pH":
            basis[row.name] = PH(row.value)
        elif u == "Log activity":
            basis[row.name] = LogActivity(row.value)
        elif u == "Hetero. equil.":
            basis[row.name] = Equilibrium(row.phase, row.value)
        elif u == "Suppressed":
            suppress.append(row.name)
        elif u in _NO_OP_UNITS:
            continue
        else:
            raise Unsupported(f"constraint {u!r} on {row.name}")

    for name, option, _ in inp.alter_suppress:
        if option != "Suppress":
            raise Unsupported(f"alter/suppress option {option}")
        if name in db and db[name].kind in (Kind.MINERAL, Kind.GAS, Kind.LIQUID):
            suppress_phases.append(name)
        else:
            suppress.append(name)

    kind, value, couple = inp.redox
    t_k = inp.temperature_c + KELVIN
    redox = None
    if kind == "eh":
        redox = Pe(value * _LEGACY_F / (_LEGACY_R * t_k * LN10))
    elif kind == "pe":
        redox = Pe(value)
    elif kind == "logfo2":
        redox = LogFO2(value)
    elif kind == "couple":
        redox = Couple(couple)
    elif kind != "o2line":
        raise Unsupported(f"redox option {kind}")

    tds_kg = inp.tds[1] if inp.tds[0] == "mg/kg.sol" else inp.tds[1] / inp.density
    return dict(
        temperature_c=inp.temperature_c, basis=basis, redox=redox,
        switches={old: new for old, new in inp.switches},
        suppress=suppress, suppress_phases=suppress_phases,
        balance_on=inp.balance_on, activity_model=inp.activity_model,
        ph_scale=inp.ph_scale, density=inp.density, tds_mg_per_kg=tds_kg,
    )


@dataclass
class Deviation:
    what: str
    ours: float
    reference: float

    @property
    def delta(self) -> float:
        return abs(self.ours - self.reference)


@dataclass
class CaseReport:
    name: str
    description: str
    deviations: list[Deviation] = field(default_factory=list)
    missing_species: list[str] = field(default_factory=list)
    missing_solids: list[str] = field(default_factory=list)
    n_species: int = 0
    n_solids: int = 0
    result: Speciation | None = None

    def worst(self, n: int = 5) -> list[Deviation]:
        return sorted(self.deviations, key=lambda d: -d.delta)[:n]

    @property
    def max_delta(self) -> float:
        return max((d.delta for d in self.deviations), default=0.0)


def compare(name: str, description: str, res: Speciation, ref: Eq3nrReference,
            log_m_floor: float = -12.0) -> CaseReport:
    rep = CaseReport(name, description, result=res)
    for what, ours, theirs in [
        ("ionic strength", res.ionic_strength, ref.ionic_strength),
        ("sum molalities", res.sum_molality, ref.sum_molality),
        ("log a(H2O)", res.log_a_water, ref.log_a_water),
        ("osmotic coeff", res.osmotic_coefficient, ref.osmotic_coefficient),
        ("pmH", res.pmh, ref.pmh),
        ("pH (NBS)", res.ph if res.ph_scale == "nbs" else float("nan"), ref.ph_nbs),
    ]:
        if theirs is not None and not math.isnan(ours):
            rel = what in ("ionic strength", "sum molalities")
            rep.deviations.append(Deviation(what, ours / theirs - 1 if rel else ours, 0.0 if rel else theirs))
    for sp, (_, log_m, log_g, log_a) in ref.species.items():
        if log_m < log_m_floor:
            continue
        got = res.species.get(sp)
        if got is None:
            rep.missing_species.append(sp)
            continue
        rep.n_species += 1
        rep.deviations.append(Deviation(f"log m {sp}", math.log10(got.molality), log_m))
        rep.deviations.append(Deviation(f"log gamma {sp}", got.log_gamma, log_g))
    for ph, log_qk in ref.solids.items():
        if ph in res.saturation:
            rep.n_solids += 1
            rep.deviations.append(Deviation(f"SI {ph}", res.saturation[ph], log_qk))
        else:
            rep.missing_solids.append(ph)
    for gas, logf in ref.gases.items():
        if gas in res.log_fugacity:
            rep.deviations.append(Deviation(f"log f {gas}", res.log_fugacity[gas], logf))
    return rep


def run_case(db: ThermoDatabase, input_path: str | Path,
             output_path: str | Path | None = None) -> list[CaseReport]:
    """Run every problem in an input file; one report per problem."""
    input_path = Path(input_path)
    inputs = read_eq3nr_inputs(input_path)
    output_path = Path(output_path) if output_path else \
        input_path.parent / "Outputs" / (input_path.stem + ".3o")
    refs = parse_eq3nr_outputs(output_path)
    if len(refs) != len(inputs):
        raise ValueError(f"{input_path.name}: {len(inputs)} problems but {len(refs)} outputs")
    reports = []
    for k, (inp, ref) in enumerate(zip(inputs, refs), start=1):
        res = speciate(db, **to_speciate_args(db, inp))
        name = input_path.stem if len(inputs) == 1 else f"{input_path.stem}#{k}"
        reports.append(compare(name, inp.description, res, ref))
    return reports

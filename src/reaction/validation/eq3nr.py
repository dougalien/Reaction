"""Parse EQ3NR ``.3o`` output files so Reaction can be checked against them.

Only what validation needs is extracted: scalar summary values, the aqueous
species table, the pure-solid saturation states and the gas fugacities.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

_FLOAT = r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?"


@dataclass
class Eq3nrReference:
    ionic_strength: float | None = None
    sum_molality: float | None = None
    log_a_water: float | None = None
    osmotic_coefficient: float | None = None
    ph_nbs: float | None = None
    pmh: float | None = None
    pe: float | None = None
    eh_volts: float | None = None
    log_fo2: float | None = None
    # name -> (molality, log molality, log gamma, log activity)
    species: dict[str, tuple[float, float, float, float]] = field(default_factory=dict)
    solids: dict[str, float] = field(default_factory=dict)      # name -> log Q/K
    gases: dict[str, float] = field(default_factory=dict)       # name -> log fugacity


def _scalar(text: str, label: str) -> float | None:
    m = re.search(rf"^\s*{re.escape(label)}\s*=\s*({_FLOAT})", text, re.M)
    return float(m.group(1)) if m else None


def parse_eq3nr_output(path: str | Path) -> Eq3nrReference:
    text = Path(path).read_text(encoding="latin-1")
    ref = Eq3nrReference(
        ionic_strength=_scalar(text, "Ionic strength (I)"),
        sum_molality=_scalar(text, "Sum of molalities"),
        log_a_water=_scalar(text, "Log activity of water"),
        osmotic_coefficient=_scalar(text, "Osmotic coefficient"),
        log_fo2=_scalar(text, "Log oxygen fugacity"),
    )
    m = re.search(rf"^ NBS pH scale\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})", text, re.M)
    if m:
        ref.ph_nbs, ref.eh_volts, ref.pe = map(float, m.groups())
    m = re.search(rf"^ Mesmer pH \(pmH\) scale\s+({_FLOAT})", text, re.M)
    if m:
        ref.pmh = float(m.group(1))

    titles = {"Distribution of Aqueous Solute Species": "species",
              "Saturation States of Pure Solids": "solids",
              "Fugacities": "gases"}
    section = None
    for line in text.splitlines():
        header = re.match(r"^\s+---\s+(.*?)\s+---\s*$", line)
        if header:                      # any new titled block ends the previous one
            section = titles.get(header.group(1))
            continue
        if section == "species":
            m = re.match(rf"^ (\S+)\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*$", line)
            if m:
                ref.species[m.group(1)] = tuple(map(float, m.groups()[1:]))
        elif section == "solids":
            m = re.match(rf"^\s+(\S+)\s+({_FLOAT})\s+({_FLOAT})(\s+\w+)?\s*$", line)
            if m and m.group(1) != "Phase":
                ref.solids[m.group(1)] = float(m.group(2))
        elif section == "gases":
            m = re.match(rf"^\s+(\S+)\s+({_FLOAT})\s+\S+\s*$", line)
            if m and m.group(1) != "Gas":
                ref.gases[m.group(1)] = float(m.group(2))
    return ref

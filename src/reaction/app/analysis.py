"""Reading a water analysis as people write it, and turning it into molalities.

An analysis is plain text, one solute per line (commas, semicolons and tabs work too):

    Na 10768
    Ca2+ 10.3 mmol/kg
    SiO2 (mg/L) 12
    alkalinity is not read: give HCO3 (total carbonate) instead

Names may be element symbols, ions in any common charge notation (``Ca``, ``Ca2+``, ``Ca++``, ``Ca+2``)
or English names.  A unit after the number applies to that line; otherwise the card's default unit does.
Units follow The Geochemist's Workbench: ``mg/kg``, ``mmol/kg`` and ``mol/kg`` are per kilogram of
*solution*; ``molal``, ``mmolal`` (or ``mmol/kgw``) are per kilogram of *water*; ``mg/L`` and ``mmol/L``
are per litre, read with a density of 1 kg/L.  An element given by mass is weighed as that element
(``Si 140 mg/kg`` is 140 mg of silicon, ``SiO2 300 mg/kg`` is 300 mg of silica).

Everything here is stoichiometry and arithmetic; no thermodynamic data are needed, so the card can
show what it read before anything is loaded.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from ..thermo.eosdb import ATOMIC_WEIGHTS

# basis species -> (charge, formula weighed when a mass is given, other names)
SOLUTES: dict[str, tuple[int, str, tuple[str, ...]]] = {
    "Li+": (1, "Li", ("lithium",)),
    "Na+": (1, "Na", ("sodium",)),
    "K+": (1, "K", ("potassium",)),
    "Rb+": (1, "Rb", ("rubidium",)),
    "Cs+": (1, "Cs", ("cesium", "caesium")),
    "Be+2": (2, "Be", ("beryllium",)),
    "Mg+2": (2, "Mg", ("magnesium",)),
    "Ca+2": (2, "Ca", ("calcium",)),
    "Sr+2": (2, "Sr", ("strontium",)),
    "Ba+2": (2, "Ba", ("barium",)),
    "Mn+2": (2, "Mn", ("manganese",)),
    "Fe+2": (2, "Fe", ("iron", "ferrous iron")),
    "Zn+2": (2, "Zn", ("zinc",)),
    "Cu+2": (2, "Cu", ("copper",)),
    "Pb+2": (2, "Pb", ("lead",)),
    "Cd+2": (2, "Cd", ("cadmium",)),
    "Ni+2": (2, "Ni", ("nickel",)),
    "Co+2": (2, "Co", ("cobalt",)),
    "Hg+2": (2, "Hg", ("mercury",)),
    "Ag+": (1, "Ag", ("silver",)),
    "Al+3": (3, "Al", ("aluminium", "aluminum")),
    "Cl-": (-1, "Cl", ("chloride", "chlorine")),
    "F-": (-1, "F", ("fluoride", "fluorine")),
    "Br-": (-1, "Br", ("bromide", "bromine")),
    "I-": (-1, "I", ("iodide", "iodine")),
    "SO4-2": (-2, "SO4", ("sulfate", "sulphate")),
    "HCO3-": (-1, "HCO3", ("bicarbonate", "hydrogen carbonate")),
    "SiO2": (0, "SiO2", ("silica",)),
    "B(OH)3": (0, "BO3H3", ("boric acid",)),
    "NH3": (0, "NH3", ("ammonia",)),
    "HPO4-2": (-2, "HPO4", ("phosphate",)),
}

# other ways of giving the same component, weighed as written: name -> (species, formula)
EXTRA_NAMES: dict[str, tuple[str, str]] = {
    "Si": ("SiO2", "Si"), "B": ("B(OH)3", "B"), "N": ("NH3", "N"), "NH4": ("NH3", "NH4"),
    "NH4+": ("NH3", "NH4"), "ammonium": ("NH3", "NH4"), "P": ("HPO4-2", "P"),
    "PO4": ("HPO4-2", "PO4"), "PO4-3": ("HPO4-2", "PO4"), "PO43-": ("HPO4-2", "PO4"),
    "CO2": ("HCO3-", "CO2"), "DIC": ("HCO3-", "C"), "TIC": ("HCO3-", "C"),
    "total carbonate": ("HCO3-", "C"), "carbonate": ("HCO3-", "CO3"), "CO3": ("HCO3-", "CO3"),
    "CO3-2": ("HCO3-", "CO3"), "CO32-": ("HCO3-", "CO3"), "CO3--": ("HCO3-", "CO3"),
}

# unit -> (amount: "mass" in mg, "mole" in mol or "eq" in equivalents; scale; per "solution", "water" or "litre")
UNITS: dict[str, tuple[str, float, str]] = {
    "mg/kg": ("mass", 1.0, "solution"), "ppm": ("mass", 1.0, "solution"),
    "g/kg": ("mass", 1e3, "solution"), "ug/kg": ("mass", 1e-3, "solution"),
    "ppb": ("mass", 1e-3, "solution"),
    "mg/l": ("mass", 1.0, "litre"), "g/l": ("mass", 1e3, "litre"), "ug/l": ("mass", 1e-3, "litre"),
    "mol/kg": ("mole", 1.0, "solution"), "mmol/kg": ("mole", 1e-3, "solution"),
    "umol/kg": ("mole", 1e-6, "solution"),
    "molal": ("mole", 1.0, "water"), "mmolal": ("mole", 1e-3, "water"),
    "umolal": ("mole", 1e-6, "water"), "mol/kgw": ("mole", 1.0, "water"),
    "mmol/kgw": ("mole", 1e-3, "water"), "umol/kgw": ("mole", 1e-6, "water"),
    "mol/l": ("mole", 1.0, "litre"), "mmol/l": ("mole", 1e-3, "litre"),
    "umol/l": ("mole", 1e-6, "litre"),
    "meq/kg": ("eq", 1e-3, "solution"), "meq/l": ("eq", 1e-3, "litre"),
}
UNIT_LABELS = {"mg/kg": "mg/kg", "mg/l": "mg/L", "mmol/kg": "mmol/kg", "mol/kg": "mol/kg",
               "mmolal": "mmol/kg water", "molal": "mol/kg water", "meq/l": "meq/L"}
DEFAULT_UNIT = "mg/kg"

_SUPERSCRIPTS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻₀₁₂₃₄₅₆₇₈₉µμ−", "0123456789+-0123456789uu-")


def _formula_weight(formula: str) -> float:
    total = 0.0
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        total += ATOMIC_WEIGHTS[el] * (int(n) if n else 1)
    return total


def _names() -> dict[str, tuple[str, str]]:
    """Every accepted spelling -> (basis species, formula to weigh)."""
    out: dict[str, tuple[str, str]] = {}
    for species, (z, formula, words) in SOLUTES.items():
        base = re.sub(r"[+-]\d*$", "", species)
        spellings = {species, base}
        if z:
            sign, n = ("+" if z > 0 else "-"), abs(z)
            spellings |= {base + sign * n, f"{base}{n}{sign}" if n > 1 else base + sign,
                          f"{base}{sign}{n}"}
            if z == -2:
                spellings.add(base + "=")
        if species == "Fe+2":
            spellings |= {"Fe(II)", "Fe2"}
        for s in spellings:
            out[s] = (species, formula)
        for w in words:
            out[w] = (species, formula)
    out.update(EXTRA_NAMES)
    return out


NAMES = _names()
_LOWER = {k.lower(): v for k, v in NAMES.items()}

_UNIT_RE = "|".join(re.escape(u) for u in sorted(
    list(UNITS) + ["mg/kg.sol", "mmol/kg h2o", "mol/kg h2o", "mg/kgw"], key=len, reverse=True))
_NUMBER = r"\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d*\.?\d+(?:[eE][-+]?\d+)?"
_ENTRY = re.compile(
    r"(?P<name>[A-Za-z][A-Za-z0-9()]*(?:[+\-=]+\d?)?(?:\s+(?:carbonate|iron|acid))?)"
    r"\s*(?:\(\s*(?P<u1>[^)]{1,12})\s*\))?\s*[:=]?\s*"
    r"(?P<lt><\s*)?(?P<num>" + _NUMBER + r")"
    r"(?:\s*(?P<u2>" + _UNIT_RE + r")(?![a-z]))?", re.I)
_PH = re.compile(r"\bpH\s*(?:[:=]|of|is)?\s*(?P<num>\d+(?:\.\d+)?)", re.I)


def unit_key(text: str | None) -> str | None:
    """The canonical key of a unit as written (``mg/L`` -> ``mg/l``), or None."""
    if not text:
        return None
    u = text.translate(_SUPERSCRIPTS).strip().lower().replace(" ", "")
    u = {"mg/kg.sol": "mg/kg", "mmol/kgh2o": "mmol/kgw", "mol/kgh2o": "mol/kgw",
         "mg/kgw": "mg/kg"}.get(u, u)
    return u if u in UNITS else None


def lookup(name: str) -> tuple[str, str] | None:
    """(basis species, formula weighed) for a solute name as written, or None."""
    n = name.translate(_SUPERSCRIPTS).strip()
    return NAMES.get(n) or NAMES.get(n.replace(" ", "")) or _LOWER.get(n.lower())


@dataclass
class Entry:
    species: str                 # basis species
    value: float
    unit: str                    # key of UNITS
    written: str                 # the name as the user wrote it
    formula: str                 # what a mass refers to

    @property
    def charge(self) -> int:
        return SOLUTES[self.species][0]

    def line(self) -> str:
        return f"{self.written} {self.value:g} {UNIT_LABELS.get(self.unit, self.unit)}"


@dataclass
class Analysis:
    entries: list[Entry] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    ph: float | None = None
    unit: str | None = None      # a unit the text names on its own (a table heading)

    @property
    def species(self) -> list[str]:
        return [e.species for e in self.entries]


def parse(text: str, default_unit: str = DEFAULT_UNIT, *, strict: bool = True) -> Analysis:
    """Read an analysis.  ``strict`` (the card's own text) reports words that are not solutes;
    the chat reads leniently and keeps only what it recognises."""
    out = Analysis()
    text = (text or "").translate(_SUPERSCRIPTS)
    m = _PH.search(text)
    if m:
        out.ph = float(m.group("num"))
        text = text[:m.start()] + " " + text[m.end():]
    default = unit_key(default_unit) or DEFAULT_UNIT
    lone = re.search(r"(?<![\w/])(" + _UNIT_RE + r")(?![a-z/])", text, re.I)
    seen: dict[str, str] = {}
    last_unit = None if strict or not lone else unit_key(lone.group(1))   # a chat table: one unit for all
    for m in _ENTRY.finditer(text):
        written = m.group("name").strip()
        hit = lookup(written)
        if hit is None:
            if strict and not re.fullmatch(r"(?i)t|temp(erature)?|ph|eh|density|tds", written):
                out.problems.append(f"'{written}' is not a solute I know")
            continue
        species, formula = hit
        unit_text = m.group("u2") or m.group("u1")
        unit = unit_key(unit_text) if unit_text else (default if strict else last_unit or default)
        if unit_text and unit is not None:
            last_unit = unit
        if unit is None:
            out.problems.append(f"{written}: unknown unit '{unit_text}'")
            continue
        value = float(m.group("num").replace(",", ""))
        if m.group("lt") or value <= 0:
            out.notes.append(f"{written} is below detection or zero, so it is left out")
            continue
        if UNITS[unit][0] == "eq" and SOLUTES[species][0] == 0:
            out.problems.append(f"{written} has no charge, so it cannot be given in equivalents")
            continue
        if species in seen:
            out.problems.append(f"{SOLUTES_LABEL(species)} is given twice ({seen[species]} and {written})")
            continue
        seen[species] = written
        out.entries.append(Entry(species, value, unit, written, formula))
    if lone and not any(m.group("u1") or m.group("u2") for m in _ENTRY.finditer(text)):
        out.unit = unit_key(lone.group(1))
    if re.search(r"(?i)alkalinity", text):
        out.problems.append("alkalinity is not read; give HCO3 (total dissolved carbonate) instead")
    return out


def SOLUTES_LABEL(species: str) -> str:
    return re.sub(r"[+-]\d*$", "", species) if species not in ("SiO2", "B(OH)3", "NH3") else species


def molalities(entries: list[Entry], density: float = 1.0) -> tuple[dict[str, float], float]:
    """Moles per kilogram of water of each basis species, and the dissolved solids (mg/kg of solution)
    used to convert from per kilogram of solution.  Per-litre values are read with ``density``."""
    tds = 0.0
    for _ in range(4):                         # amounts per kg of water need the solids, and back
        per_solution: dict[str, float] = {}        # moles per kg of solution
        mass = 0.0
        water_fraction = 1.0 - tds / 1e6
        for e in entries:
            kind, scale, per = UNITS[e.unit]
            fw = _formula_weight(e.formula)
            amount = e.value * scale
            if kind == "mass":
                moles = amount / 1000.0 / fw
            elif kind == "eq":
                moles = amount / abs(e.charge)
            else:
                moles = amount
            if per == "litre":
                moles /= density
            elif per == "water":
                moles *= water_fraction
            per_solution[e.species] = per_solution.get(e.species, 0.0) + moles
            mass += moles * fw * 1000.0
        if abs(mass - tds) < 1e-6 * max(mass, 1.0):
            tds = mass
            break
        tds = mass
    water_fraction = 1.0 - tds / 1e6
    return {k: v / water_fraction for k, v in per_solution.items()}, tds


def balance_ion(molal: dict[str, float], excess: float | None = None) -> str | None:
    """The ion to adjust for charge balance: the most abundant one (in equivalents) of the sign that
    is short, so balancing only ever adds to it.  ``excess`` is the speciated charge imbalance
    (eq/kg, positive for excess cations); without it the nominal charges are used."""
    if excess is None:
        excess = sum(SOLUTES[s][0] * m for s, m in molal.items())
    want = -1 if excess > 0 else 1
    ions = [(abs(SOLUTES[s][0]) * m, s) for s, m in molal.items() if SOLUTES[s][0] * want > 0]
    return max(ions)[1] if ions else None


def preview(text: str, default_unit: str = DEFAULT_UNIT) -> dict:
    """What the card shows under the analysis: each line as read, with its molality."""
    got = parse(text, default_unit)
    molal, tds = molalities(got.entries) if got.entries else ({}, 0.0)
    return {"lines": [{"written": e.written, "species": e.species, "value": e.value,
                       "unit": UNIT_LABELS.get(e.unit, e.unit), "molal": molal.get(e.species),
                       "charge": e.charge}
                      for e in got.entries],
            "tds_mg_per_kg": round(tds, 1), "notes": got.notes, "problems": got.problems}

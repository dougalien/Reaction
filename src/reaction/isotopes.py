"""Stable isotopes (δ18O, δ13C, δ34S) and 87Sr/86Sr as a layer on a solved reaction path.

At each point every carrier of an element (water, dissolved carbon or sulfur species, each mineral, the gas)
is in isotopic equilibrium with the others at that temperature: δ_carrier = δ_ref + Δ_carrier(T), with
Δ = 1000 ln α against a reference (water for O, CO2 for C, H2S for S) from ``data/isotope_fractionation.tsv``.
The system's inventory Σ n δ is carried from point to point: reactants add their own δ, matter that leaves
(flow-through minerals and gas, a phase that separates) takes the δ it had, and evaporation takes vapour in
equilibrium with the water (Horita and Wesolowski 1994).  δ_ref then follows from the mass balance
Σ n (δ_ref + Δ) = inventory.  This is the δ-notation approximation GWB and EQ6 also use.

Strontium does not fractionate: 87Sr/86Sr is a tracer carried with the fluid, mixed from the starting fluid and
the Sr that dissolving rock (or a reactant given Sr) adds; secondary minerals are taken to keep none.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path


DATA = Path(__file__).parent / "data" / "isotope_fractionation.tsv"
SYSTEMS = {"O": "d18O", "C": "d13C", "S": "d34S"}
SR_GRAMS_PER_MOLE = 87.62
STANDARDS = {"d18O": "VSMOW", "d13C": "VPDB", "d34S": "VCDT"}


@dataclass(frozen=True)
class Factor:
    system: str
    phase: str
    a: float
    b: float
    c: float
    d: float
    t_max_c: float | None
    source: str

    def delta(self, t_c: float) -> float:
        """1000 ln α (phase − reference) at t_c."""
        t = t_c + 273.15
        return self.a * 1e8 / t ** 3 + self.b * 1e6 / t ** 2 + self.c * 1e3 / t + self.d


@lru_cache(maxsize=1)
def table() -> dict[tuple[str, str], Factor]:
    rows = [line for line in DATA.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    out = {}
    for r in csv.DictReader(rows, delimiter="\t"):
        out[(r["system"], r["phase"])] = Factor(r["system"], r["phase"], float(r["a"]), float(r["b"]), float(r["c"]),
                                                float(r["d"]), float(r["t_max_c"]) if r["t_max_c"] else None,
                                                r["source"])
    return out


def liquid_vapour_o(t_c: float) -> float:
    """1000 ln α (liquid water − vapour) for 18O, Horita and Wesolowski (1994), 0 at and above 374 °C."""
    if t_c >= 374.0:
        return 0.0
    t = t_c + 273.15
    return max(-7.685 + 6.7123e3 / t - 1.6664e6 / t ** 2 + 0.35041e9 / t ** 3, 0.0)


# minerals of the database without a row of their own: the closest phase that has one
_SAME = {
    "albite,high": "albite", "albite,low": "albite", "K-feldspar,high": "K-feldspar", "K-feldspar,low": "K-feldspar",
    "quartz,beta": "quartz", "chalcedony": "quartz", "amorphous silica": "quartz", "coesite": "quartz",
    "cristobalite,alpha": "quartz", "cristobalite,beta": "quartz", "tridymite,high": "quartz",
    "tridymite,low": "quartz", "antigorite": "lizardite", "clinozoisite": "zoisite", "orthoepidote": "epidote",
    "daphnite": "chamosite", "ferro-actinolite": "actinolite", "enstatite,clino": "enstatite",
    "enstatite,proto": "enstatite", "dickite": "kaolinite", "halloysite": "kaolinite", "celestite": "barite",
    "gypsum": "anhydrite", "laumontite": "wairakite", "pumpellyite": "epidote", "margarite": "muscovite",
    "fluormuscovite": "muscovite", "pargasite": "hornblende", "ferropargasite": "hornblende",
    "anthophyllite": "actinolite", "grunerite": "actinolite", "fluortremolite": "tremolite",
    "fluorphlogopite": "phlogopite", "chlorphlogopite": "phlogopite", "fluorannite": "annite",
    "chlorannite": "annite", "fluortalc": "talc", "minnesotaite": "talc", "greenalite": "lizardite",
    "sepiolite": "talc", "hydroxyapatite": "apatite", "fluorapatite": "apatite", "chlorapatite": "apatite",
    "maghemite": "magnetite", "lepidocrocite": "goethite", "sudoite": "clinochlore", "cookeite": "clinochlore",
    "Fe-amesite": "chamosite", "Mg-amesite": "clinochlore", "huntite": "dolomite", "wurtzite": "sphalerite",
    "bornite": "chalcopyrite", "pseudowollastonite": "wollastonite", "meionite": "anorthite",
}


@dataclass
class Carriers:
    """Moles of one element in each carrier and its Δ against the reference, at one point."""
    moles: dict[str, float] = field(default_factory=dict)
    delta: dict[str, float] = field(default_factory=dict)

    def add(self, name: str, moles: float, delta: float) -> None:
        if moles > 0:
            self.moles[name] = self.moles.get(name, 0.0) + moles
            self.delta[name] = delta

    @property
    def total(self) -> float:
        return sum(self.moles.values())


class Fractionation:
    """Δ of each carrier at a temperature, with the analogues used where a phase has no data of its own."""

    def __init__(self, db):
        self.db = db
        self.analogues: dict[str, set[str]] = {}       # system -> "mineral (as analogue)"
        self.beyond: set[str] = set()                  # factors used above their calibrated range

    def _factor(self, system: str, phase: str, t_c: float) -> float:
        f = table()[(system, phase)]
        if f.t_max_c is not None and t_c > f.t_max_c + 1:
            self.beyond.add(f"{phase} ({system}, to {f.t_max_c:.0f} °C)")
        return f.delta(t_c)

    def mineral(self, system: str, name: str, t_c: float) -> float:
        rows = table()
        key = {"S": self._sulfur_phase, "C": self._carbon_phase}.get(system, self._oxygen_phase)(name)
        if key != name and not (system == "S" and key == "sulfate" and name in ("anhydrite", "barite")):
            self.analogues.setdefault(system, set()).add(f"{name} as {key}")
        return self._factor(system, key, t_c) if (system, key) in rows else 0.0

    def _elements(self, name: str) -> dict[str, float]:
        return self.db.species[name].elements

    def _oxygen_phase(self, name: str) -> str:
        if ("O", name) in table():
            return name
        if name in _SAME:
            return _SAME[name]
        el = self._elements(name)
        if el.get("C"):
            return "calcite"
        if el.get("S"):
            return "anhydrite"
        if el.get("Si"):
            return "albite"
        return "goethite" if el.get("H") else "hematite"

    def _carbon_phase(self, name: str) -> str:
        if ("C", name) in table():
            return name
        if name in ("huntite", "ankerite"):
            return "dolomite"
        return "graphite" if self._elements(name).get("O", 0) == 0 else "calcite"

    def _sulfur_phase(self, name: str) -> str:
        if ("S", name) in table():
            return name
        el = self._elements(name)
        if el.get("O", 0) >= 4 * el.get("S", 1):
            return "sulfate"
        return _SAME.get(name, "H2S")

    def aqueous(self, system: str, name: str, el: dict[str, float], t_c: float) -> float:
        """Δ of a dissolved species (or gas) carrying C or S."""
        if system == "C":
            o, h, c = el.get("O", 0.0), el.get("H", 0.0), el["C"]
            if o >= 3 * c:
                return self._factor("C", "HCO3-", t_c)      # bicarbonate, carbonate and their complexes
            if o >= 2 * c:
                return 0.0                                  # CO2 (aqueous or gas), the reference
            if o == 0 and h == 4 * c:
                return self._factor("C", "CH4", t_c)
            return self._factor("C", "CO" if o >= c else "CH4", t_c)   # CO; organic species like methane
        o, s = el.get("O", 0.0), el["S"]
        if o >= 4 * s:
            return self._factor("S", "sulfate", t_c)
        if o == 0:
            return 0.0                                      # H2S, HS-, and sulfide complexes
        return self._factor("S", "SO2", t_c)                # sulfite, thiosulfate and other intermediates


def element_totals(db, bulk: dict[str, float]) -> dict[str, float]:
    """Moles of O, C and S in a system given as moles of basis components."""
    out = {e: 0.0 for e in SYSTEMS}
    for comp, moles in bulk.items():
        sp = db.species.get(comp)
        if sp is None:
            continue
        for e in SYSTEMS:
            out[e] += moles * sp.elements.get(e, 0.0)
    return out


def carriers(db, state, t_c: float, frac: Fractionation) -> dict[str, Carriers]:
    """The carriers of O, C and S in an equilibrium state.  Oxygen not in a mineral is in the water (dissolved
    species and gas take the water's δ18O); C and S in species without a class of their own go to the reference."""
    totals = element_totals(db, state.bulk)
    out = {e: Carriers() for e in SYSTEMS}
    solids = dict(state.minerals)
    for members in state.solid_solutions.values():
        for m, n in members.items():
            solids[m] = solids.get(m, 0.0) + n
    for name, n in solids.items():
        el = db.species[name].elements
        for e in SYSTEMS:
            if el.get(e):
                out[e].add(name, n * el[e], frac.mineral(e, name, t_c))
    kg = state.water_kg
    fluid = list((name, st.molality * kg) for name, st in state.speciation.species.items()) + list(state.gas.items())
    for name, n in fluid:
        sp = db.species.get(name)
        if sp is None or n <= 0:
            continue
        for e in ("C", "S"):
            if sp.elements.get(e):
                out[e].add(_fluid_label(e, name, sp.elements), n * sp.elements[e],
                           frac.aqueous(e, name, sp.elements, t_c))
    out["O"].add("water", totals["O"] - out["O"].total, 0.0)
    for e in ("C", "S"):
        rest = totals[e] - out[e].total
        if rest > 1e-12 * max(totals[e], 1e-300):
            out[e].add("CO2" if e == "C" else "H2S", rest, 0.0)
    return out


def _fluid_label(e: str, name: str, el: dict[str, float]) -> str:
    """Dissolved carbon and sulfur grouped into DIC, methane and organics, sulfate and sulfide."""
    if name.endswith("(g)"):
        return name
    if e == "C":
        o, c = el.get("O", 0.0), el["C"]
        return "DIC" if o >= 2 * c else "CH4" if name.startswith("CH4") else "organic C"
    o, s = el.get("O", 0.0), el["S"]
    return "sulfate" if o >= 4 * s else "sulfide" if o == 0 else "other S"


@dataclass
class Source:
    """The isotopic make-up of something added to the system: δ per system and its Sr."""
    deltas: dict[str, float]                 # d18O, d13C, d34S
    sr_moles_per_unit: float = 0.0           # mol Sr per unit of the reactant (or per kg of the fluid)
    sr_ratio: float | None = None            # 87Sr/86Sr


class Tracker:
    """Carries the isotope inventory along a path of equilibrium states (see the module text)."""

    def __init__(self, db, fluid: Source, bulk0: dict[str, float]):
        self.db = db
        self.frac = Fractionation(db)
        self.inventory: dict[str, float] = {}
        self.totals = element_totals(db, bulk0)
        for e, key in SYSTEMS.items():
            self.inventory[e] = self.totals[e] * fluid.deltas.get(key, 0.0)
        self.prev: dict[str, Carriers] | None = None
        self.prev_deltas: dict[str, dict[str, float]] = {}
        self.sr = fluid.sr_moles_per_unit            # the fluid starts as 1 kg of water
        self.sr_inv = self.sr * (fluid.sr_ratio or 0.0)
        self.sr_ratio = fluid.sr_ratio
        self._cl = None
        self._vapour = (0.0, 0.0, 0.0)               # this step's evaporation: moles, Σ m δ_start, Σ m Δ(l-v)

    def add(self, comp_moles: dict[str, float], source: Source, units: float) -> None:
        """``units`` of a reactant whose components per unit are ``comp_moles`` join the system."""
        for e in SYSTEMS:
            n = units * sum(m * self.db.species[c].elements.get(e, 0.0) for c, m in comp_moles.items()
                            if c in self.db.species)
            self.totals[e] += n
            self.inventory[e] += n * source.deltas.get(SYSTEMS[e], 0.0)
        sr = units * source.sr_moles_per_unit
        if sr > 0 and source.sr_ratio is not None:
            self.sr += sr
            self.sr_inv += sr * source.sr_ratio

    def evaporate(self, water_moles: float, t_c: float) -> None:
        """Water vapour in equilibrium with the water leaves.  The vapour's δ18O is the mean of the water's at
        the start and at the end of the step (settled in ``point``), so coarse steps follow Rayleigh closely."""
        d_water = self.prev_deltas.get("O", {}).get("water")
        if d_water is None:
            d_water = self.inventory["O"] / self.totals["O"] if self.totals["O"] else 0.0
        self.totals["O"] -= water_moles
        m, start, lv = self._vapour
        self._vapour = (m + water_moles, start + water_moles * d_water, lv + water_moles * liquid_vapour_o(t_c))

    def point(self, state, t_c: float, chloride: float) -> dict[str, dict[str, float]]:
        """δ of every carrier at this point, after matter that left since the last one is taken out."""
        now = carriers(self.db, state, t_c, self.frac)
        totals = element_totals(self.db, state.bulk)
        for e in SYSTEMS:
            gone = self.totals[e] - totals[e]
            if gone > 1e-9 * max(self.totals[e], 1e-300):
                self.inventory[e] -= gone * self._leaving_delta(e, now[e])
            elif gone < 0:                                   # gained without a reactant (rounding)
                mean = self.inventory[e] / self.totals[e] if self.totals[e] else 0.0
                self.inventory[e] -= gone * mean
            self.totals[e] = totals[e]
        if self._cl is not None and self._cl > 0 and chloride < self._cl:   # salt (and Sr) left with a phase
            self.sr_inv *= chloride / self._cl
            self.sr *= chloride / self._cl
        self._cl = chloride
        out: dict[str, dict[str, float]] = {}
        for e, c in now.items():
            total = c.total
            if total <= 0:
                out[e] = {}
                continue
            inv, weight = self.inventory[e], total
            m, start, lv = self._vapour
            if e == "O" and m > 0:                    # vapour at the mean of the water's δ18O before and after
                inv -= 0.5 * start - lv
                weight += 0.5 * m
            ref = (inv - sum(n * c.delta[k] for k, n in c.moles.items())) / weight
            if e == "O" and m > 0:
                self.inventory[e] = inv - 0.5 * m * (ref + c.delta.get("water", 0.0))
            out[e] = {k: ref + c.delta[k] for k in c.moles}
        self.prev, self.prev_deltas = now, out
        self._vapour = (0.0, 0.0, 0.0)
        return out

    @property
    def sr_fluid_ratio(self) -> float | None:
        return self.sr_inv / self.sr if self.sr > 0 else None

    def _leaving_delta(self, e: str, now: Carriers) -> float:
        """δ of what left between the last point and this one: the carriers that shrank, as they were."""
        if self.prev is None:
            return self.inventory[e] / self.totals[e] if self.totals[e] else 0.0
        shrank = {k: n - now.moles.get(k, 0.0) for k, n in self.prev[e].moles.items() if n > now.moles.get(k, 0.0)}
        weight = sum(shrank.values())
        if weight <= 0:
            return self.inventory[e] / self.totals[e] if self.totals[e] else 0.0
        return sum(w * self.prev_deltas[e][k] for k, w in shrank.items()) / weight

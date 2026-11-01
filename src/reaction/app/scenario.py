"""The scenario card: everything needed to describe one model run, in plain terms.

A card is a flat set of fields.  It is what the chat fills in, what the user can edit by hand, and
what is stored next to every result, so any run can be reproduced from its card alone.

A run follows a straight line in temperature and pressure from a start point to an end point, in
``steps`` equal increments.  Either variable (or both) may change; equal start and end values hold a
variable fixed.  Nothing about pressure is baked into a database: the equilibrium constants are
computed from Gibbs energies at every point on the path.
"""
from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field, fields
from typing import Any

from ..thermo import water
from . import analysis, rocks
from .analysis import DEFAULT_UNIT
from ..sorption import KD_MEDIA
from .fluids import CUSTOM, CUSTOM_LABEL, PRESETS, REDOX_CHOICES, FluidPreset, label as fluid_label

MODES = ("static", "flow_through", "both")
KINDS = ("path", "transport")
MAX_COLUMN_SOLVES = 40000             # cells x shifts x (1 + dispersion mixes) in one column run
KINETIC_STEP_SOLVES = 15              # equilibrium calculations in a cell's kinetic step, typically
MODE_LABELS = {"static": "Static", "flow_through": "Flow-through", "both": "Both"}
ACTIVITY_MODELS = ("bdot", "davies", "pitzer")
FLUID_FIELDS = ("analysis", "units", "ph", "fluid_t_c", "fluid_p_bar", "charge_balance", "redox",
                "redox_value")                                       # what a custom fluid is
REACTANT_UNITS = {"g": ("g", 1.0), "gram": ("g", 1.0), "grams": ("g", 1.0), "mg": ("g", 1e-3),
                  "ug": ("g", 1e-6), "µg": ("g", 1e-6), "microgram": ("g", 1e-6), "micrograms": ("g", 1e-6),
                  "kg": ("g", 1e3), "mol": ("mol", 1.0), "mole": ("mol", 1.0), "moles": ("mol", 1.0),
                  "mmol": ("mol", 1e-3), "umol": ("mol", 1e-6)}
MAX_REACTANTS = 12


def reactant_text(r: dict) -> str:
    """'10 g albite', with masses below a milligram in µg: '20 µg forams'."""
    amount, unit = r["amount"], r["unit"]
    if unit == "g" and 0 < amount < 1e-3:
        amount, unit = amount * 1e6, "µg"
    return f"{float(f'{amount:.6g}'):g} {unit} {r['name']}"
MAX_UNCERTAINTY_RUNS = 200
FOLLOW_PHASES = ("whole", "brine", "vapour")
ISOTOPE_KEYS = {"d18o": "d18O", "18o": "d18O", "o18": "d18O", "d13c": "d13C", "13c": "d13C", "c13": "d13C",
                "d34s": "d34S", "34s": "d34S", "s34": "d34S", "sr": "Sr", "sr_ppm": "Sr", "srppm": "Sr",
                "sr_mg_per_kg": "Sr", "87sr/86sr": "87Sr/86Sr", "87sr": "87Sr/86Sr", "sr87": "87Sr/86Sr",
                "sr_ratio": "87Sr/86Sr", "dd": "dD", "d2h": "dD", "2h": "dD", "deuterium": "dD", "dh": "dD"}
DEFAULT_UNCERTAINTY_RUNS = 20
TIME_UNITS = {"seconds": 1.0, "minutes": 60.0, "hours": 3600.0, "days": 86400.0,
              "years": 365.25 * 86400.0}
TIME_WORDS = {"s": "seconds", "sec": "seconds", "second": "seconds", "min": "minutes",
              "minute": "minutes", "h": "hours", "hr": "hours", "hour": "hours", "d": "days",
              "day": "days", "y": "years", "yr": "years", "year": "years", "a": "years"}

T_MIN, T_MAX = 0.0, 600.0            # degrees C
P_MIN, P_MAX = 1.0, 5000.0           # bar (the range the HKF equations are calibrated for)

# carbonates that supersaturated seawater "wants" to form but that do not form on laboratory and
# field time scales at these temperatures (the usual practice is to block them)
DEFAULT_BLOCKED = ["dolomite", "magnesite", "huntite"]

# oxides, elements and ice that the database lists but that are not sensible products here
ALWAYS_EXCLUDED = ("ice", "sulfur", "S_liq", "graphite", "potassium oxide", "sodium oxide",
                   "lime", "periclase")


@dataclass
class Scenario:
    title: str = ""
    kind: str = "path"
    fluid: str = "seawater"
    t_start_c: float = 2.0
    t_end_c: float = 350.0
    p_start_bar: float = 500.0
    p_end_bar: float = 500.0
    steps: int = 36
    mode: str = "static"
    blocked_minerals: list[str] = field(default_factory=lambda: list(DEFAULT_BLOCKED))
    activity_model: str = "bdot"
    # a custom fluid (fluid = "custom"): the analysis as typed, its default unit, the pH and the
    # temperature they were measured at, and how to close the charge balance
    analysis: str = ""
    units: str = DEFAULT_UNIT
    ph: float | None = None
    fluid_t_c: float = 25.0
    fluid_p_bar: float | None = None      # pressure of the measurement; None = 1 bar or the boiling curve
    charge_balance: str = "auto"
    # the redox state of a custom fluid: one of fluids.REDOX_CHOICES or a redox species of the analysis
    # (O2, H2, Fe+3, H2S, NO3-) used as the couple; ``redox_value`` for log_fo2, eh and pe
    redox: str = "auto"
    redox_value: float | None = None
    # minerals or aqueous species added steadily along the path, per kg of the starting fluid's water:
    # [{"name": "albite", "amount": 10.0, "unit": "g"}]; reacted in full by the last step
    reactants: list[dict] = field(default_factory=list)
    # mixing: a second (preset) fluid added steadily until there are ``mix_ratio`` kg of its water per kg
    # of the first fluid's water.  With ``mix_t_c`` the temperature follows the mixing (equal heat
    # capacities: the mass-weighted mean of the start temperature and ``mix_t_c``) instead of the path.
    mix_with: str = ""
    mix_ratio: float = 1.0
    mix_t_c: float | None = None
    # where the steps fall: "ratio" evenly in the kg of second fluid added; "temperature" evenly in the
    # mixture's temperature (only with ``mix_t_c``), so that the hot end, reached by the first few
    # percent of the second fluid, is not skipped
    mix_spacing: str = "ratio"
    # ``redox_apart`` holds every redox couple apart (H2 from O2, sulfide from sulfate, methane from carbonate,
    # Fe(III) from Fe(II), nitrate from ammonia), so a mixture keeps the disequilibrium its fluids bring;
    # ``minerals_form`` false keeps every mineral from forming; ``energy`` gives the catabolic energy of the
    # chemosynthetic reactions at each step (see catabolic.py)
    redox_apart: bool = False
    minerals_form: bool = True
    energy: bool = False
    # a free gas phase (CO2, CH4, H2O, H2, H2S, N2 ...) forms when the dissolved gases' fugacities add
    # up to the pressure; methane and sulfate are held apart (no reaction between them) unless asked
    gas_phase: bool = True
    methane_sulfate: str = "apart"
    # kinetics: the run lasts ``duration`` (in ``duration_unit``) and minerals with measured rates
    # dissolve or grow at those rates instead of all at once; the reactant amounts are then what is
    # there at the start, with ``surface_area_cm2_per_g`` of reactive surface per gram
    kinetics: bool = False
    # evaporation: water is removed until the fluid is ``evaporate`` times more concentrated (1 = none),
    # in steps evenly spaced in the logarithm of the concentration factor
    evaporate: float = 1.0
    duration: float = 1.0
    duration_unit: str = "years"
    surface_area_cm2_per_g: float = 100.0
    # uncertainty bands: this many extra runs, each with every log K shifted by its own random draw
    # (one standard deviation from the reliability check, see thermo.validity.logk_sigma); 0 = none
    uncertainty_runs: int = 0
    # solid solutions (plagioclase, carbonate, chlorite, olivine ...; see solidsolution.py) may form
    # instead of their pure endmembers; a name in blocked_minerals leaves that one out
    solid_solutions: bool = False
    # phase separation (NaCl-H2O): "whole" follows the whole fluid; "brine" lets the vapour leave where the
    # fluid splits (boiling off), "vapour" lets the brine stay behind (the vapour rises on)
    follow_phase: str = "whole"
    # isotopes (see isotopes.py): d18O, d13C, d34S and 87Sr/86Sr along the run; ``isotope_inputs`` sets the
    # starting values of "fluid", "mix" (the second fluid) or a reactant by its card name:
    # {"fluid": {"d18O": -5.0}, "basalt glass": {"Sr": 120.0, "87Sr/86Sr": 0.7030}}; the rest are defaults
    isotopes: bool = False
    isotope_inputs: dict = field(default_factory=dict)
    # sorption (see sorption.py), per kg of the starting water: cation exchange sites (meq, Gaines-Thomas) that
    # start in equilibrium with the fluid, and hydrous ferric oxide (g) with Dzombak and Morel surface sites
    exchange_meq: float = 0.0
    hfo_g: float = 0.0
    # clay minerals (g per kg of water), each with the exchange sites of its cation exchange capacity
    # (sorption.clays: kaolinite 20, illite 225, montmorillonite 764 meq/kg); and a solid with the partition
    # coefficients (Kd) of Allison and Allison (2005) for ``kd_medium`` (sediment, suspended or soil), g per kg
    # of water (0 = none).  With exchange sites the Kd holds only the metals that have no exchange constant.
    kaolinite_g: float = 0.0
    illite_g: float = 0.0
    montmorillonite_g: float = 0.0
    # illite and montmorillonite also bind metals on their edge sites (2SPNE, Bradbury and Baeyens); false = exchange only
    clay_edges: bool = True
    kd_medium: str = "sediment"
    kd_solid_g: float = 0.0
    # a column (kind = "transport", see transport.py): the card's fluid flows into ``cells`` cells filled with
    # ``pore_fluid`` (empty: the same fluid) and, in each, the reactants per kg of pore water brought to
    # equilibrium first; temperature and pressure run from the start values at the inlet to the end values at
    # the outlet.  ``pore_volumes`` of fluid pass through; ``dispersivity_m`` spreads the fronts.
    cells: int = 20
    column_m: float = 1.0
    dispersivity_m: float = 0.05
    pore_volumes: float = 2.0
    pore_fluid: str = ""

    # ---- editing --------------------------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Scenario":
        data = dict(data)
        if "pressure_bar" in data:                       # cards saved before pressure could vary
            p = data.pop("pressure_bar")
            data.setdefault("p_start_bar", p)
            data.setdefault("p_end_bar", p)
        if data.get("kind") == "heating":
            data["kind"] = "path"
        if data.get("title") in ("Seawater heating",):    # the old default title
            data["title"] = ""
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in known})

    def changed(self, changes: dict[str, Any]) -> tuple["Scenario", list[str]]:
        """Copy with ``changes`` applied.  Returns the new card and a list of problems; a problem
        in one field never blocks the others.  ``pressure_bar`` and ``temperature_c`` set both ends."""
        changes = dict(changes)
        for single, (lo, hi) in {"pressure_bar": ("p_start_bar", "p_end_bar"),
                                 "temperature_c": ("t_start_c", "t_end_c")}.items():
            if single in changes:
                value = changes.pop(single)
                changes.setdefault(lo, value)
                changes.setdefault(hi, value)
        problems: list[str] = []
        new = Scenario.from_dict(self.to_dict())
        known = {f.name for f in fields(Scenario)}
        # editing any part of a preset fluid turns it into a custom fluid that starts from the preset
        asked = _coerce("fluid", changes["fluid"]) if "fluid" in changes else self.fluid
        edits_fluid = any(k in changes for k in FLUID_FIELDS)
        start = asked if asked in PRESETS and edits_fluid else self.fluid
        if (asked == CUSTOM or edits_fluid) and start in PRESETS:
            changes.pop("fluid", None)
            p = PRESETS[start]
            new.analysis, new.units, new.ph = p.analysis, p.units, p.ph
            new.fluid_t_c, new.charge_balance = p.measured_t_c, p.charge_balance
            new.fluid_p_bar = p.measured_p_bar
            new.redox, new.redox_value = p.redox, p.redox_value
            new.fluid = CUSTOM
        for key, value in changes.items():
            if key not in known:
                problems.append(f"'{key}' is not something the card can change")
                continue
            try:
                setattr(new, key, _coerce(key, value))
            except (TypeError, ValueError) as exc:
                problems.append(f"{key}: {exc}")
        return new, problems

    # ---- the path -------------------------------------------------------------------------
    @property
    def varies(self) -> list[str]:
        out = []
        t_end = self.tp(1.0)[0]
        if abs(t_end - self.t_start_c) > 1e-9:
            out.append("temperature_c")
        if abs(self.p_end_bar - self.p_start_bar) > 1e-9:
            out.append("pressure_bar")
        return out

    @property
    def mixing(self) -> bool:
        return bool(self.mix_with)

    @property
    def evaporating(self) -> bool:
        return self.evaporate > 1.0

    @property
    def mix_follows_temperature(self) -> bool:
        return self.mixing and self.mix_t_c is not None

    @property
    def mix_even_in_t(self) -> bool:
        return self.mix_follows_temperature and self.mix_spacing == "temperature"

    def fractions(self) -> list[float]:
        """The fraction (0 to 1) of the way along the path at each step: even, or, for mixing evenly
        spaced in temperature, at the mixing ratios that give those temperatures."""
        n = max(int(self.steps), 2)
        even = [i / (n - 1) for i in range(n)]
        if not self.mix_even_in_t:
            return even
        t0, t_mix = self.t_start_c, self.mix_t_c
        t1 = self.tp(1.0)[0]
        out = []
        for f in even:
            t = t0 + (t1 - t0) * f
            out.append(min(max((t0 - t) / (t - t_mix) / self.mix_ratio, 0.0), 1.0) if f > 0 else 0.0)
        out[-1] = 1.0
        return out

    def tp(self, f: float) -> tuple[float, float]:
        """(temperature in C, pressure in bar) at the fraction ``f`` (0 to 1) of the way along the path."""
        p = self.p_start_bar + (self.p_end_bar - self.p_start_bar) * f
        if self.mix_follows_temperature:
            m = self.mix_ratio * f
            return (self.t_start_c + self.mix_t_c * m) / (1.0 + m), p
        return self.t_start_c + (self.t_end_c - self.t_start_c) * f, p

    def points(self) -> list[tuple[float, float]]:
        """(temperature in C, pressure in bar) at each step, evenly spaced along the path."""
        return [self.tp(f) for f in self.fractions()]

    @property
    def display_title(self) -> str:
        if self.title.strip():
            return self.title.strip()
        return f"{fluid_label(self.fluid)} {'column' if self.transport else 'run'}"

    def fluid_spec(self) -> FluidPreset:
        """The starting fluid, preset or custom, in one shape."""
        if self.fluid in PRESETS:
            return PRESETS[self.fluid]
        return FluidPreset(CUSTOM, CUSTOM_LABEL, "A fluid entered as an analysis.", self.fluid_t_c,
                           self.analysis, self.ph if self.ph is not None else float("nan"),
                           "entered on the card", units=self.units, charge_balance=self.charge_balance,
                           redox=self.redox, redox_value=self.redox_value,
                           measured_p_bar=self.fluid_p_bar)

    def fluid_analysis(self) -> analysis.Analysis:
        spec = self.fluid_spec()
        return analysis.parse(spec.analysis, spec.units)

    def redox_setting(self) -> tuple[str, Any]:
        """How the starting fluid's redox state is set: ("none", None), ("air", None or "assumed"),
        ("log_fo2" | "eh" | "pe", value) or ("couple", species).  Raises ValueError when the
        card cannot say."""
        return resolve_redox(self.fluid_spec(), self.fluid_analysis())

    @property
    def duration_s(self) -> float:
        return self.duration * TIME_UNITS.get(self.duration_unit, 1.0)

    @property
    def hold_methane_apart(self) -> bool:
        return self.methane_sulfate == "apart"

    # ---- checking -------------------------------------------------------------------------
    def problems(self) -> list[str]:
        out: list[str] = []
        if self.kind not in KINDS:
            out.append(f"kind must be path or transport, not '{self.kind}'")
        elif self.transport:
            out.extend(self._column_problems())
        if self.fluid not in PRESETS and self.fluid != CUSTOM:
            out.append(f"unknown fluid '{self.fluid}' (known: {', '.join(PRESETS)}, {CUSTOM})")
        elif self.fluid == CUSTOM:
            out.extend(self._fluid_problems())
        ranges_ok = True
        for name, vals, lo, hi, unit in (
                ("temperature", (self.t_start_c, self.t_end_c), T_MIN, T_MAX, "C"),
                ("pressure", (self.p_start_bar, self.p_end_bar), P_MIN, P_MAX, "bar")):
            if not all(lo <= v <= hi for v in vals):
                out.append(f"{name} must be between {lo:g} and {hi:g} {unit}")
                ranges_ok = False
        if not self.transport and not (2 <= self.steps <= 200):
            out.append("steps must be between 2 and 200")
        for r in self.reactants:
            if not (r.get("amount", 0) > 0):
                out.append(f"the amount of {r.get('name', 'a reactant')} must be positive")
        if len(self.reactants) > MAX_REACTANTS:
            out.append(f"at most {MAX_REACTANTS} reactants")
        if self.mixing:
            if self.mix_with not in PRESETS:
                out.append(f"can only mix with a named fluid ({', '.join(PRESETS)}) so far")
            if not (0.0 < self.mix_ratio <= 1000.0):
                out.append("the mixing ratio must be above 0 and at most 1000 kg per kg")
            if self.mix_t_c is not None and not (T_MIN <= self.mix_t_c <= T_MAX):
                out.append(f"the second fluid's temperature must be between {T_MIN:g} and {T_MAX:g} C")
            if self.mix_spacing not in MIX_SPACINGS:
                out.append(f"mix_spacing must be one of {', '.join(MIX_SPACINGS)}")
            elif self.mix_spacing == "temperature" and (self.mix_t_c is None or self.mix_t_c == self.t_start_c):
                out.append("steps even in temperature need a second fluid at another temperature")
        if ranges_ok:
            if not self.varies and not self.reactants and not self.mixing and not self.kinetics \
                    and not self.evaporating and not self.transport:
                out.append("nothing changes along this path: give different start and end values "
                           "for temperature or pressure")
            boil = self._first_boiling_point()
            if boil:
                t, p, psat = boil
                out.append(f"at {t:.0f} C water boils below {psat:.0f} bar but the path has {p:.0f} bar "
                           f"there; keep the pressure at or above {math.ceil(psat)} bar")
        if self.kinetics:
            if not self.reactants:
                out.append("kinetics needs a reactant: give the minerals or rock the fluid reacts with")
            if not (self.duration > 0 and math.isfinite(self.duration)):
                out.append("the duration must be positive")
            elif self.duration_s > 1e6 * TIME_UNITS["years"]:
                out.append("the duration must be at most a million years")
            if self.duration_unit not in TIME_UNITS:
                out.append(f"the duration unit must be one of {', '.join(TIME_UNITS)}")
            if not (0.0 < self.surface_area_cm2_per_g <= 1e7):
                out.append("the surface area must be above 0 and at most 10,000,000 cm2 per g")
        if not (1.0 <= self.evaporate <= 10000.0):
            out.append("the evaporation factor must be between 1 (none) and 10000")
        elif self.evaporating and self.kinetics:
            out.append("evaporation and kinetics cannot be combined yet; turn one of them off")
        if not (0 <= self.uncertainty_runs <= MAX_UNCERTAINTY_RUNS):
            out.append(f"the number of uncertainty runs must be between 0 (none) and {MAX_UNCERTAINTY_RUNS}")
        elif 0 < self.uncertainty_runs < 5:
            out.append("uncertainty bands need at least 5 runs (20 or more for a stable range)")
        if not 0.0 <= self.exchange_meq <= 100000.0:
            out.append("the exchange capacity must be between 0 and 100000 meq per kg of water")
        if not 0.0 <= self.hfo_g <= 1000.0:
            out.append("the hydrous ferric oxide must be between 0 and 1000 g per kg of water")
        for name in ("kaolinite", "illite", "montmorillonite"):
            if not 0.0 <= getattr(self, f"{name}_g") <= 100000.0:
                out.append(f"the {name} must be between 0 and 100000 g per kg of water")
        if not 0.0 <= self.kd_solid_g <= 100000.0:
            out.append("the solid with partition coefficients must be between 0 and 100000 g per kg of water")
        if self.kd_medium not in KD_MEDIA:
            out.append(f"kd_medium must be one of {', '.join(KD_MEDIA)}")
        if self.follow_phase not in FOLLOW_PHASES:
            out.append(f"follow_phase must be one of {', '.join(FOLLOW_PHASES)}")
        if self.mode not in MODES:
            out.append(f"mode must be one of {', '.join(MODES)}")
        if self.activity_model not in ACTIVITY_MODELS:
            out.append(f"activity_model must be one of {', '.join(ACTIVITY_MODELS)}")
        return out

    def clay_grams(self) -> dict[str, float]:
        """Grams of each clay per kg of water (only those on the card)."""
        out = {"kaolinite": self.kaolinite_g, "illite": self.illite_g, "montmorillonite": self.montmorillonite_g}
        return {k: v for k, v in out.items() if v > 0}

    @property
    def clay_meq(self) -> float:
        """Exchange sites the clays carry, meq per kg of water."""
        from ..sorption import clays
        return sum(g * clays()[c].cec_eq_per_kg for c, g in self.clay_grams().items())

    @property
    def transport(self) -> bool:
        return self.kind == "transport"

    @property
    def shifts(self) -> int:
        """Shifts of a column run: each moves the fluid one cell along."""
        return max(1, int(round(self.pore_volumes * self.cells)))

    @property
    def mixing_fraction(self) -> float:
        """Dispersion per shift: the dispersivity over the cell length (D dt / dx^2 with D = alpha v)."""
        return self.dispersivity_m / (self.column_m / max(self.cells, 1))

    def cell_points(self) -> list[tuple[float, float]]:
        """(temperature C, pressure bar) of each cell, from the inlet's (the start values) to the outlet's."""
        n = max(int(self.cells), 2)
        return [self.tp(i / (n - 1)) for i in range(n)]

    def _column_problems(self) -> list[str]:
        out: list[str] = []
        if not (2 <= self.cells <= 100):
            out.append("a column needs between 2 and 100 cells")
        if not (0.0 < self.column_m <= 1e5):
            out.append("the column length must be above 0 and at most 100 km")
        if not (0.0 <= self.dispersivity_m <= self.column_m):
            out.append("the dispersivity must be between 0 and the column length")
        if not (0.0 < self.pore_volumes <= 100.0):
            out.append("pore volumes must be above 0 and at most 100")
        if self.pore_fluid and self.pore_fluid not in PRESETS:
            out.append(f"the column can start filled with a named fluid ({', '.join(PRESETS)}) or the same fluid")
        if 2 <= self.cells <= 100 and 0 < self.column_m and 0 < self.pore_volumes <= 100:
            mixes = math.ceil(3.0 * self.mixing_fraction - 1e-12) if self.mixing_fraction > 0 else 0
            # a cell's kinetic step takes about fifteen equilibrium calculations
            solves = self.cells * self.shifts * ((KINETIC_STEP_SOLVES if self.kinetics else 1) + mixes)
            if solves > MAX_COLUMN_SOLVES:
                out.append(f"this column needs {solves:,} equilibrium calculations, more than the "
                           f"{MAX_COLUMN_SOLVES:,} allowed: use fewer cells or pore volumes, or a smaller dispersivity")
        for flag, what in ((self.evaporating, "evaporation"),
                           (self.mixing, "mixing in a second fluid"), (self.follow_phase != "whole", "following the brine or vapour"),
                           (self.isotopes, "isotopes"), (self.uncertainty_runs > 0, "uncertainty bands"),
                           (self.redox_apart, "holding redox apart"), (not self.minerals_form, "keeping minerals from forming"),
                           (self.energy, "catabolic energy")):
            if flag:
                out.append(f"{what} is not available in a column yet; turn it off for a transport run")
        return out

    def _fluid_problems(self) -> list[str]:
        out: list[str] = []
        if analysis.unit_key(self.units) is None:
            out.append(f"unknown unit '{self.units}'")
        got = self.fluid_analysis()
        out.extend(got.problems)
        if not got.entries:
            out.append("the analysis is empty: give at least one solute, for example \"Na 100, Cl 150\"")
        if self.ph is None or not math.isfinite(self.ph):
            out.append("give the pH of the fluid")
        elif not (0.0 <= self.ph <= 14.0):
            out.append("pH must be between 0 and 14")
        if self.fluid_p_bar is None:
            if not (0.0 <= self.fluid_t_c <= 300.0):
                out.append("the analysis temperature must be between 0 and 300 C (or give the pressure "
                           "it was measured at)")
        elif not (P_MIN <= self.fluid_p_bar <= P_MAX):
            out.append(f"the analysis pressure must be between {P_MIN:g} and {P_MAX:g} bar")
        elif not (T_MIN <= self.fluid_t_c <= T_MAX):
            out.append(f"the analysis temperature must be between {T_MIN:g} and {T_MAX:g} C")
        elif self.fluid_t_c + 273.15 < water.T_CRIT and self.fluid_p_bar < 0.99 * (
                water.saturation_pressure_mpa(self.fluid_t_c + 273.15) * water.BAR_PER_MPA):
            out.append(f"at {self.fluid_t_c:g} C and {self.fluid_p_bar:g} bar water is steam; give a higher "
                       "pressure for the analysis")
        try:
            self.redox_setting()
        except ValueError as exc:
            out.append(str(exc))
        cb = self.charge_balance
        if cb not in ("auto", "none") and cb not in got.species:
            names = ", ".join(got.species) or "none yet"
            out.append(f"charge balance on '{cb}' needs that ion in the analysis (ions given: {names})")
        elif cb not in ("auto", "none") and analysis.SOLUTES[cb][0] == 0:
            out.append(f"{cb} has no charge, so it cannot close the charge balance")
        return out

    def _first_boiling_point(self) -> tuple[float, float, float] | None:
        """First place along the path (checked finely, not just at the steps) where water is steam."""
        n = 201
        for i in range(n):
            t, p = self.tp(i / (n - 1))
            t_k = t + 273.15
            if t_k < water.T_CRIT:
                psat = water.saturation_pressure_mpa(t_k) * water.BAR_PER_MPA
                if p < 0.99 * psat:
                    return t, p, psat
        return None

    def describe(self) -> str:
        if self.transport:
            return self._describe_column()
        text = self._describe_path()
        if self.reactants:
            added = ", ".join(reactant_text(r) for r in self.reactants)
            verb = (f"reacting with {added} for {self.duration:g} {self.duration_unit} at measured rates"
                    if self.kinetics else f"adding {added}")
            text = text.replace(f" ({self._how()}", f", {verb} ({self._how()}")
        if self.evaporating:
            text = text.replace(f" ({self._how()}", f", evaporated to {self.evaporate:g} times its concentration "
                                                   f"({self._how()}")
        if self.solid_solutions:
            text = text.replace(f" ({self._how()}", f", with solid solutions ({self._how()}")
        if self.follow_phase != "whole":
            what = "the brine as vapour boils off" if self.follow_phase == "brine" else "the vapour, leaving the brine"
            text = text.replace(f" ({self._how()}", f", following {what} ({self._how()}")
        if self.isotopes:
            text = text.replace(f" ({self._how()}", f", with isotopes ({self._how()}")
        if self.redox_apart:
            text = text.replace(f" ({self._how()}", f", redox couples held apart ({self._how()}")
        if not self.minerals_form:
            text = text.replace(f" ({self._how()}", f", no minerals forming ({self._how()}")
        if self.energy:
            text = text.replace(f" ({self._how()}", f", with catabolic energies ({self._how()}")
        if self.exchange_meq > 0:
            text = text.replace(f" ({self._how()}", f", with {self.exchange_meq:g} meq of exchange sites ({self._how()}")
        if self.hfo_g > 0:
            text = text.replace(f" ({self._how()}", f", with {self.hfo_g:g} g of ferric oxyhydroxide surface "
                                                   f"({self._how()}")
        for bit in self._sorbent_bits():
            text = text.replace(f" ({self._how()}", f", with {bit} ({self._how()}")
        if self.uncertainty_runs:
            text = text.replace(f" ({self._how()}", f", with a log K uncertainty range from {self.uncertainty_runs} "
                                                   f"runs ({self._how()}")
        if self.mixing:
            other = fluid_label(self.mix_with).lower()
            hot = f" at {self.mix_t_c:g} °C" if self.mix_t_c is not None else ""
            text = text.replace(f" ({self._how()}",
                                f", mixed with up to {self.mix_ratio:g} kg of {other}{hot} per kg ({self._how()}")
        return text

    def _describe_column(self) -> str:
        label = fluid_label(self.fluid)
        t0, t1, p0, p1 = self.t_start_c, self.t_end_c, self.p_start_bar, self.p_end_bar
        temp = f"{t0:g} °C" if t0 == t1 else f"{t0:g} °C at the inlet to {t1:g} °C at the outlet"
        pres = f"{p0:g} bar" if p0 == p1 else f"{p0:g} to {p1:g} bar"
        text = (f"{label} flowing through a {self.column_m:g} m column of {self.cells} cells, "
                f"{self.pore_volumes:g} pore volumes, at {temp} and {pres}")
        if self.pore_fluid and self.pore_fluid != self.fluid:
            text += f", first filled with {fluid_label(self.pore_fluid).lower()}"
        bits = []
        if self.reactants:
            bits.append(", ".join(reactant_text(r) for r in self.reactants) + " per kg of pore water")
        if self.exchange_meq > 0:
            bits.append(f"{self.exchange_meq:g} meq of exchange sites")
        if self.hfo_g > 0:
            bits.append(f"{self.hfo_g:g} g of ferric oxyhydroxide surface")
        bits.extend(self._sorbent_bits())
        if self.solid_solutions:
            bits.append("solid solutions")
        if bits:
            text += ", with " + "; ".join(bits)
        if self.kinetics:
            text += f", reacting at measured rates over {self.duration:g} {self.duration_unit}"
        return text + f" (dispersivity {self.dispersivity_m:g} m)"

    def _sorbent_bits(self) -> list[str]:
        out = []
        clays = self.clay_grams()
        if clays:
            out.append(" and ".join(", ".join(f"{g:g} g {c}" for c, g in clays.items()).rsplit(", ", 1))
                       + ("" if self.clay_edges or not ({"illite", "montmorillonite"} & set(clays))
                          else " (exchange sites only)"))
        if self.kd_solid_g > 0:
            what = {"sediment": "bed sediment", "suspended": "suspended matter"}.get(self.kd_medium, self.kd_medium)
            out.append(f"{self.kd_solid_g:g} g of {what} holding metals by partition coefficients (Kd)")
        return out

    def _how(self) -> str:
        return {"static": "static", "flow_through": "flow-through",
                "both": "static and flow-through"}.get(self.mode, self.mode)

    def _describe_path(self) -> str:
        label = fluid_label(self.fluid)
        if self.mix_follows_temperature:
            p0, p1 = self.p_start_bar, self.p_end_bar
            where = f"{p0:g} bar" if p0 == p1 else f"{p0:g} to {p1:g} bar"
            return f"{label} at {self.t_start_c:g} °C and {where} ({self._how()}, {self.steps} steps)"
        how = {"static": "static", "flow_through": "flow-through",
               "both": "static and flow-through"}.get(self.mode, self.mode)
        tail = f"({how}, {self.steps} steps)"
        t0, t1, p0, p1 = self.t_start_c, self.t_end_c, self.p_start_bar, self.p_end_bar
        vary = self.varies
        if vary == ["temperature_c"]:
            verb = "heated" if t1 > t0 else "cooled"
            return f"{label} {verb} from {t0:g} to {t1:g} °C at {p0:g} bar {tail}"
        if vary == ["pressure_bar"]:
            verb = "compressed" if p1 > p0 else "decompressed"
            return f"{label} at {t0:g} °C {verb} from {p0:g} to {p1:g} bar {tail}"
        if vary:
            return (f"{label} taken from {t0:g} °C and {p0:g} bar to {t1:g} °C and "
                    f"{p1:g} bar {tail}")
        return f"{label} at {t0:g} °C and {p0:g} bar {tail}"


def resolve_redox(spec: FluidPreset, got: analysis.Analysis) -> tuple[str, Any]:
    """See :meth:`Scenario.redox_setting`."""
    present = [s for s in got.species if s in analysis.REDOX]
    # a redox species can be the couple only when its partner (if any) is in the analysis too
    couples = [s for s in ("O2", "H2", "Fe+3", "H2S", "NO3-", "CH4")
               if s in present and (analysis.REDOX[s] is None or analysis.REDOX[s] in got.species)]
    choice = spec.redox
    if choice == "auto":
        if couples:
            return "couple", couples[0]
        if present:
            names = ", ".join(analysis.SOLUTES_LABEL(s) for s in present)
            raise ValueError(f"{names} needs the redox state of the fluid: give its partner too "
                             f"(Fe(II) with Fe(III), sulfate with sulfide), dissolved O2 or H2, or set "
                             f"the redox state (air, log fO2, Eh or pe)")
        # no redox species at all: the fluid is taken to be in equilibrium with air, so the redox state
        # is always known and iron oxides, sulfides and ferric iron can react; the run says so
        return "air", "assumed"
    if choice == "none":
        if present:
            names = ", ".join(analysis.SOLUTES_LABEL(s) for s in present)
            raise ValueError(f"{names} needs the redox state of the fluid, which is set to not tracked")
        return "none", None
    if choice == "air":
        return "air", None
    if choice in ("log_fo2", "eh", "pe"):
        if spec.redox_value is None or not math.isfinite(spec.redox_value):
            raise ValueError(f"give the {REDOX_LABELS[choice]} of the fluid")
        return choice, spec.redox_value
    if choice in analysis.REDOX:
        if choice not in couples:
            partner = analysis.REDOX[choice]
            need = analysis.SOLUTES_LABEL(choice) + (
                f" and {analysis.SOLUTES_LABEL(partner)}" if partner else "")
            raise ValueError(f"the redox couple {analysis.SOLUTES_LABEL(choice)} needs {need} "
                             f"in the analysis")
        return "couple", choice
    raise ValueError(f"unknown redox setting '{choice}'")


REDOX_LABELS = {"log_fo2": "log fO2", "eh": "Eh (volts)", "pe": "pe"}
MIX_SPACINGS = ("ratio", "temperature")


def _coerce(key: str, value: Any) -> Any:
    if key == "redox":
        text = str(value).strip()
        low = text.lower().replace(" ", "_").replace("-", "_")
        low = {"fo2": "log_fo2", "logfo2": "log_fo2", "log_f_o2": "log_fo2", "oxygen_fugacity":
               "log_fo2", "not_tracked": "none", "off": "none", "atmosphere": "air", "": "auto"
               }.get(low, low)
        if low in REDOX_CHOICES:
            return low
        hit = analysis.lookup(text)
        if hit is not None and hit[0] in analysis.REDOX:
            return hit[0]
        raise ValueError(f"'{text}' is not a redox setting (use auto, none, air, log_fo2, eh, pe, or "
                         f"a redox species such as O2, H2, Fe+3 or H2S)")
    if key == "redox_value":
        if value is None or (isinstance(value, str) and not value.strip()):
            return None
        return float(value)
    if key in ("t_start_c", "t_end_c", "p_start_bar", "p_end_bar", "fluid_t_c", "mix_ratio", "duration",
               "surface_area_cm2_per_g", "evaporate", "exchange_meq", "hfo_g", "kaolinite_g", "illite_g",
               "montmorillonite_g", "kd_solid_g", "column_m", "dispersivity_m",
               "pore_volumes"):
        return float(value)
    if key == "duration_unit":
        low = str(value).strip().lower()
        low = TIME_WORDS.get(low, TIME_WORDS.get(low.rstrip("s"), low))
        if low not in TIME_UNITS:
            raise ValueError(f"'{value}' is not a time unit (use {', '.join(TIME_UNITS)})")
        return low
    if key in ("mix_t_c", "fluid_p_bar"):
        if value is None or (isinstance(value, str) and not value.strip()):
            return None
        return float(value)
    if key == "mix_with":
        text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
        return "" if text in ("", "none", "nothing") else text
    if key == "ph":
        if value is None or (isinstance(value, str) and not value.strip()):
            return None
        return float(value)
    if key == "analysis":
        return str(value)
    if key == "reactants":
        return parse_reactants(value)
    if key == "units":
        unit = analysis.unit_key(str(value))
        if unit is None:
            raise ValueError(f"unknown unit '{value}' (try {', '.join(analysis.UNIT_LABELS.values())})")
        return unit
    if key == "charge_balance":
        text = str(value).strip()
        if text.lower() in ("auto", "none", ""):
            return text.lower() or "none"
        hit = analysis.lookup(text)
        if hit is None:
            raise ValueError(f"'{text}' is not an ion I know")
        return hit[0]
    if key == "kind":
        low = str(value).strip().lower()
        low = {"column": "transport", "reactive_transport": "transport", "1d": "transport", "heating": "path"}.get(low, low)
        if low not in KINDS:
            raise ValueError(f"'{value}' is not path or transport")
        return low
    if key == "pore_fluid":
        text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
        return "" if text in ("", "same", "none", "the_same") else text
    if key in ("steps", "uncertainty_runs", "cells"):
        if key == "uncertainty_runs" and isinstance(value, str) and value.strip().lower() in ("", "off", "none", "no"):
            return 0
        if key == "uncertainty_runs" and value is True:
            return DEFAULT_UNCERTAINTY_RUNS
        return int(round(float(value)))
    if key == "blocked_minerals":
        if isinstance(value, str):
            value = _rejoin_comma_names(value.replace(";", ",").split(","))
        if not isinstance(value, (list, tuple)):
            raise TypeError("expected a list of mineral names")
        return [str(v).strip().lower() for v in value if str(v).strip()]
    if key == "isotope_inputs":
        if not isinstance(value, dict):
            raise TypeError("expected {source: {d18O: value, ...}}")
        out = {}
        for source, vals in value.items():
            if not isinstance(vals, dict):
                raise TypeError(f"expected values for {source} as {{d18O: value, ...}}")
            here = {}
            for k, v in vals.items():
                key_ = ISOTOPE_KEYS.get(str(k).strip().lower().replace("δ", "d").replace(" ", ""))
                if key_ is None:
                    raise ValueError(f"'{k}' is not one of d18O, dD, d13C, d34S, Sr, 87Sr/86Sr")
                here[key_] = float(v)
                if key_ == "87Sr/86Sr" and not 0.69 < here[key_] < 0.80:
                    raise ValueError(f"87Sr/86Sr {v} is outside 0.69 to 0.80")
                if key_ == "Sr" and here[key_] < 0:
                    raise ValueError("Sr cannot be negative")
            out[str(source).strip().lower()] = here
        return out
    if key == "mix_spacing":
        low = str(value).strip().lower()
        low = {"even": "ratio", "mixing ratio": "ratio", "mass": "ratio", "t": "temperature",
               "temperatures": "temperature", "": "ratio"}.get(low, low)
        if low not in MIX_SPACINGS:
            raise ValueError(f"'{value}' is not ratio or temperature")
        return low
    if key in ("gas_phase", "kinetics", "solid_solutions", "isotopes", "clay_edges", "redox_apart", "minerals_form",
               "energy"):
        if isinstance(value, str):
            low = value.strip().lower()
            if low in ("true", "yes", "on", "1", "allowed", "allow"):
                return True
            if low in ("false", "no", "off", "0", "none", "suppressed"):
                return False
            raise ValueError(f"'{value}' is not yes or no")
        return bool(value)
    if key == "kd_medium":
        low = str(value).strip().lower()
        low = {"bed sediment": "sediment", "sediments": "sediment", "streambed": "sediment",
               "suspended matter": "suspended", "suspended sediment": "suspended", "spm": "suspended",
               "soils": "soil"}.get(low, low)
        if low not in KD_MEDIA:
            raise ValueError(f"'{value}' is not one of {', '.join(KD_MEDIA)}")
        return low
    if key == "follow_phase":
        low = str(value).strip().lower()
        low = {"vapor": "vapour", "liquid": "brine", "all": "whole", "fluid": "whole", "both": "whole",
               "none": "whole", "": "whole"}.get(low, low)
        if low not in FOLLOW_PHASES:
            raise ValueError(f"'{value}' is not whole, brine or vapour")
        return low
    if key == "methane_sulfate":
        low = str(value).strip().lower().replace("-", "_").replace(" ", "_")
        low = {"held_apart": "apart", "separate": "apart", "no": "apart", "equilibrate": "react",
               "equilibrium": "react", "yes": "react", "reacting": "react"}.get(low, low)
        if low not in ("apart", "react"):
            raise ValueError(f"'{value}' is not apart or react")
        return low
    if key in ("title", "fluid", "mode", "activity_model"):
        text = str(value).strip()
        if key == "title":
            return text
        return text.lower().replace("-", "_").replace(" ", "_")
    raise ValueError("unknown field")


def parse_reactants(value: Any) -> list[dict]:
    """Reactants from the card's list or from text such as "10 g albite, 0.1 mol quartz" (either order)."""
    if isinstance(value, str):
        items = []
        for part in re.split(r"[;\n]|,(?!\s*(?:ordered|disordered|low|high|alpha|beta|maximum|"
                             r"intermediate|\d+A)\b)", value):
            part = part.strip()
            if not part:
                continue
            m = (re.fullmatch(r"(?P<amount>\d*\.?\d+(?:e-?\d+)?)\s*(?P<unit>[a-zµ]+)\s+(?:of\s+)?(?P<name>.+)", part, re.I)
                 or re.fullmatch(r"(?P<name>.+?)\s+(?P<amount>\d*\.?\d+(?:e-?\d+)?)\s*(?P<unit>[a-zµ]+)", part, re.I))
            if not m:
                raise ValueError(f"read '{part}' as an amount, a unit and a name, such as 10 g albite")
            items.append({"name": m["name"], "amount": m["amount"], "unit": m["unit"]})
        value = items
    if not isinstance(value, (list, tuple)):
        raise TypeError("expected a list of reactants")
    out = []
    for item in value:
        if not isinstance(item, dict) or not str(item.get("name", "")).strip():
            raise TypeError("each reactant needs a name, an amount and a unit")
        unit_text = str(item.get("unit", "g")).strip().lower()
        if unit_text not in REACTANT_UNITS:
            raise ValueError(f"unknown unit '{unit_text}' for {item['name']} (use g, mg, ug, kg, mol or mmol)")
        unit, scale = REACTANT_UNITS[unit_text]
        name, amount = str(item["name"]).strip(), float(f"{float(item.get('amount', 0)) * scale:.12g}")
        rock = rocks.find_mineral_rock(name)
        if rock is not None:                   # a crystalline rock goes on the card as its minerals
            if unit != "g":
                raise ValueError(f"give {rock.label.lower()} in grams")
            out.extend(rocks.mineral_reactants(rock, amount))
            continue
        out.append({"name": name, "amount": amount, "unit": unit})
    return out


def _rejoin_comma_names(parts: list[str]) -> list[str]:
    """``["dolomite", "ordered"]`` came from splitting one name; join qualifiers back on."""
    qualifiers = {"ordered", "disordered", "low", "high", "alpha", "beta"}
    out: list[str] = []
    for part in parts:
        p = part.strip()
        if p.lower() in qualifiers and out:
            out[-1] = f"{out[-1]},{p}"
        elif p:
            out.append(p)
    return out


# Ready-made cards to start from.  Each is a fluid with a path that suits it; anything on the card can
# be changed afterwards.  The fluids and their sources are in fluids.py.
@dataclass(frozen=True)
class StartingCard:
    key: str
    label: str
    description: str
    changes: dict


STARTING_CARDS: dict[str, StartingCard] = {c.key: c for c in (
    StartingCard("seawater", "Seawater",
                 "Seawater heated from 2 to 350 C at 500 bar: anhydrite from about 150 C, brucite near 350 C.",
                 {}),
    StartingCard("freshwater", "Fresh water",
                 "World-average river water in equilibrium with the atmosphere, heated from 25 to 300 C at "
                 "500 bar: calcite and anhydrite come out as their solubilities fall.",
                 {"fluid": "freshwater", "t_start_c": 25.0, "t_end_c": 300.0, "steps": 28}),
    StartingCard("brine", "Brine",
                 "Great Salt Lake brine evaporated to a fifth of its water at 25 C and 1 bar, with Pitzer "
                 "activities: halite, then anhydrite, then epsomite.",
                 {"fluid": "gsl_brine", "activity_model": "pitzer", "t_start_c": 25.0, "t_end_c": 25.0,
                  "p_start_bar": 1.0, "p_end_bar": 1.0, "evaporate": 5.0, "steps": 30}),
    StartingCard("hydrothermal", "Hydrothermal fluid",
                 "The fluid from the seawater heating run (350 C, 500 bar) mixed with nine times its mass of "
                 "2 C seawater, as where a vent fluid meets the ocean: anhydrite forms in the hot mixture "
                 "and dissolves again as it cools.",
                 {"fluid": "hydrothermal", "t_start_c": 350.0, "t_end_c": 350.0, "mix_with": "seawater",
                  "mix_ratio": 9.0, "mix_t_c": 2.0, "steps": 30}),
    StartingCard("stream", "Contaminated stream",
                 "Stream water carrying metals from a contaminated site flows through 10 m of clay-bearing streambed "
                 "(kaolinite, illite and montmorillonite on exchange sites, the rest of the sediment holding Ni, As, "
                 "Cr and Hg by the EPA partition coefficients) at 15 C: Cu and Pb break through first, then Zn and "
                 "Cd; some Pb stays at the inlet as cerussite, and the Kd metals barely move.",
                 {"kind": "transport", "fluid": "contaminated_stream", "pore_fluid": "freshwater",
                  "t_start_c": 15.0, "t_end_c": 15.0, "p_start_bar": 1.0, "p_end_bar": 1.0,
                  "cells": 10, "column_m": 10.0, "dispersivity_m": 0.3, "pore_volumes": 25.0,
                  "kaolinite_g": 100.0, "illite_g": 40.0, "montmorillonite_g": 20.0,
                  "kd_medium": "sediment", "kd_solid_g": 4000.0,
                  "blocked_minerals": DEFAULT_BLOCKED + ["quartz", "quartz,beta", "chalcedony"]}),
)}


def starting_card(key: str) -> Scenario:
    """A fresh card from one of STARTING_CARDS."""
    try:
        return Scenario(**STARTING_CARDS[key].changes)
    except KeyError:
        raise KeyError(f"unknown starting card {key!r}; known cards: {', '.join(STARTING_CARDS)}") from None


# One-button examples for teaching: each loads a card that shows one feature of the model and runs it.
# ``setup`` names the card fields that set the feature up, with a sentence on each; the control window
# outlines those fields and lists the sentences, so the student sees where the feature lives on the card.
@dataclass(frozen=True)
class Example:
    key: str
    label: str
    feature: str
    description: str
    changes: dict
    setup: tuple[tuple[str, str], ...]
    plot: tuple[str, ...] = ()        # what the results window shows first, in words (view.resolve_item)
    log: bool = True
    unit: str = ""                    # and in which unit (view.UNITS), when not the default


# minerals that a low-temperature example blocks: silica other than amorphous silica forms too slowly, and the
# metamorphic minerals are stable in the database at 25 C but never form in soils
SILICA_POLYMORPHS = ("quartz", "quartz,beta", "chalcedony", "tridymite,low", "tridymite,high", "cristobalite,alpha",
                     "cristobalite,beta", "coesite")
HIGH_T_MINERALS = ("lawsonite", "paragonite", "pumpellyite", "sudoite", "prehnite", "margarite", "zoisite",
                   "clinozoisite", "epidote", "pyrophyllite", "muscovite", "clinochlore", "talc", "tremolite", "diaspore",
                   "boehmite")


EXAMPLES: dict[str, Example] = {e.key: e for e in (
    Example("solid_solutions", "Basalt alteration", "Solid solutions",
            "Seawater reacts with 1 kg of basalt glass at 350 C and 500 bar, as in the hot part of a mid-ocean "
            "ridge. With solid solutions on, chlorite, epidote and actinolite form as Fe-Mg (or Fe-Al) mixtures "
            "instead of pure endmembers: the Solid solutions tab gives their mole fractions along the run.",
            {"title": "Basalt alteration with solid solutions", "t_start_c": 350.0, "t_end_c": 350.0,
             "p_start_bar": 500.0, "p_end_bar": 500.0, "steps": 11, "solid_solutions": True,
             "reactants": [{"name": "basalt glass", "amount": 1000.0, "unit": "g"}]},
            (("solid_solutions", "Solid solutions: Can form. This one switch lets the mixtures (plagioclase, "
                                 "chlorite, epidote, olivine, serpentine ...) form in place of their endmembers."),
             ("reactants", "Reactants: 1000 g basalt glass, added steadily, so the run goes from no rock to a "
                           "water/rock ratio of 1."),
             ("t_start_c", "Temperature and pressure held at 350 C and 500 bar (equal start and end).")),
            ("chlorite", "epidote-clinozoisite", "actinolite", "albite")),
    Example("lake", "Evaporating lake", "Evaporation, brines and isotopes",
            "Great Salt Lake brine evaporates at 25 C until a tenth of its water is left, with Pitzer activities for "
            "the strong brine. Halite comes out first, then anhydrite, epsomite and sylvite, the order of a closing "
            "salt lake; the water left grows heavier in 18O and D (d18O from 0 to about +13 permil), because the "
            "vapour that leaves is light. See Mineral amounts and the Isotopes tab.",
            {"title": "Evaporating salt lake", "fluid": "gsl_brine", "activity_model": "pitzer", "t_start_c": 25.0,
             "t_end_c": 25.0, "p_start_bar": 1.0, "p_end_bar": 1.0, "steps": 30, "evaporate": 10.0,
             "isotopes": True},
            (("evaporate", "Evaporate to 10 times concentrated: water leaves as vapour, in steps even in the "
                           "logarithm of the concentration."),
             ("activity_model", "Activity model: Pitzer, needed for brines (B-dot is good to about 1 mol/kg)."),
             ("isotopes", "Isotopes: On. The vapour carries the liquid-vapour fractionation of Horita and "
                          "Wesolowski (1994), so the lake water grows heavier."),
             ("fluid", "Fluid: Great Salt Lake brine (Drever 1997), already near halite saturation.")),
            ("halite", "anhydrite", "epsomite", "sylvite")),
    Example("clays", "Contaminated stream", "Metals on clays",
            STARTING_CARDS["stream"].description + " The Sorbed and Partition coefficients tabs show what the clays "
            "hold.",
            {**STARTING_CARDS["stream"].changes, "title": "Metals on clays in a streambed"},
            (("kind", "Flow through a column: the stream water enters a 10 m streambed of 10 cells, and 25 pore "
                      "volumes pass through."),
             ("kaolinite_g", "Clays: 100 g kaolinite, 40 g illite and 20 g montmorillonite per kg of water, each "
                             "with exchange sites at its cation exchange capacity."),
             ("clay_edges", "Clay edge sites: illite and montmorillonite also bind Pb, Cu, Zn, Cd and Ni on "
                            "their crystal edges (surface complexation)."),
             ("kd_solid_g", "Kd solid: 4 kg of bed sediment holds As, Cr, Hg and Ni by their measured partition "
                            "coefficients."),
             ("pore_fluid", "The column is first filled with clean fresh water, so the fronts are visible.")),
            ("Cu", "Pb", "Zn", "Cd")),
    Example("boiling", "Boiling", "Boiling and a gas phase",
            "A CO2-bearing geothermal water at 150 C rises from 25 to 5 bar, as up a well or a hot-spring conduit. "
            "Near 10 bar its dissolved CO2 and steam can no longer stay dissolved and it boils: a gas of CO2 and steam "
            "forms and leaves (flow-through), the water loses about two thirds of its carbon, the pH rises from 6.0 "
            "to 6.8, and calcite precipitates, the scale that clogs geothermal wells. See the Gas phase tab.",
            {"title": "Boiling of a rising geothermal water", "fluid": "custom",
             "analysis": "Na+ 1600\nK+ 100\nCa+2 30\nCl- 2000\nSO4-2 30\nSiO2 150\nHCO3- 4000",
             "units": "mg/kg", "ph": 6.0, "fluid_t_c": 150.0, "fluid_p_bar": 25.0, "charge_balance": "auto",
             "redox": "none", "t_start_c": 150.0, "t_end_c": 150.0, "p_start_bar": 25.0, "p_end_bar": 5.0,
             "steps": 41, "mode": "flow_through", "gas_phase": True,
             "blocked_minerals": DEFAULT_BLOCKED + list(SILICA_POLYMORPHS)},
            (("p_end_bar", "Pressure falls from 25 to 5 bar at 150 C (water itself boils at 4.8 bar); the dissolved "
                           "CO2 makes the water boil well before that."),
             ("gas_phase", "Gas phase: Can form. CO2 and steam make the gas once their fugacities add up to the "
                           "pressure."),
             ("mode", "Flow-through: the gas (and the calcite) is left behind at each step, as the steam leaves a "
                      "rising water."),
             ("analysis", "Fluid: an illustrative sodium chloride geothermal water with 66 mmol/kg of carbonate, "
                          "pH 6.0 at 150 C (a custom fluid: edit the analysis to try your own)."),
             ("blocked_minerals", "Blocked: the silica minerals, which are too slow to form in a rising water.")),
            ("CO2 gas", "water gas")),
    Example("phase_separation", "Supercritical seawater", "Phase separation",
            "Seawater heated at 300 bar from 390 to 410 C, just past its critical point (408 C, 301 bar). There it "
            "splits into a brine and a vapour; the run follows the brine and lets the vapour leave, so the brine "
            "grows saltier (about 7 wt% NaCl at 410 C, with about two thirds of the fluid gone as vapour). The "
            "point right at the critical point does not solve and is left blank. See the Brine and vapour tab.",
            {"title": "Seawater boiling past its critical point", "t_start_c": 390.0, "t_end_c": 410.0,
             "p_start_bar": 300.0, "p_end_bar": 300.0, "steps": 21, "follow_phase": "brine"},
            (("follow_phase", "Follow: Brine. Where the fluid splits, the vapour is taken away after each step "
                              "(water and salts by the NaCl-H2O phase diagram of Driesner and Heinrich 2007)."),
             ("t_end_c", "Temperature from 390 to 410 C at a fixed 300 bar crosses the two-phase boundary at "
                         "408 C.")),
            ("brine salinity", "vapour salinity"), False),
    Example("weathering", "Flow-through weathering", "Weathering in a flow-through column",
            "Soil water (pH 5, with soil CO2) flows through 1 m of soil made of albite and a little calcite at 25 C. "
            "The calcite dissolves first and its front moves down the column (it leaves the outlet after about 9 pore "
            "volumes); albite near the inlet weathers to kaolinite. The outlet tab gives the breakthrough; the profile "
            "shows the weathered top and the fresh bottom of the column. About a minute.",
            {"title": "Flow-through weathering of a soil", "kind": "transport", "fluid": "custom",
             "analysis": "Na+ 0.017\nK+ 0.0026\nCa+2 0.0115\nMg+2 0.0041\nCl- 0.02\nSO4-2 0.0156\nHCO3- 0.4",
             "units": "mmolal", "ph": 5.0, "fluid_t_c": 25.0, "charge_balance": "auto", "redox": "air",
             "pore_fluid": "rainwater", "t_start_c": 25.0, "t_end_c": 25.0, "p_start_bar": 1.0, "p_end_bar": 1.0,
             "cells": 10, "column_m": 1.0, "dispersivity_m": 0.05, "pore_volumes": 12.0,
             "reactants": [{"name": "albite", "amount": 5.0, "unit": "g"}, {"name": "calcite", "amount": 0.2, "unit": "g"}],
             "blocked_minerals": DEFAULT_BLOCKED + list(SILICA_POLYMORPHS) + list(HIGH_T_MINERALS)},
            (("kind", "Flow through a column: 10 cells over 1 m, 12 pore volumes of soil water pass through."),
             ("reactants", "Reactants: 5 g albite and 0.2 g calcite per kg of pore water in every cell."),
             ("analysis", "Fluid flowing in: rain that has picked up soil CO2 (0.4 mmol/kg of carbonate, pH 5)."),
             ("pore_fluid", "Column first filled with rainwater."),
             ("blocked_minerals", "Blocked: silica polymorphs and minerals of metamorphic rocks (lawsonite, "
                                  "paragonite, ...) that cannot form in a soil.")),
            ("Ca", "Na", "Si")),
    Example("serpentinization", "Serpentinization", "Mantle and solid solutions",
            "Seawater reacts with 1 kg of mantle olivine (Fo90) at 300 C and 300 bar. Olivine turns to serpentine, "
            "brucite and magnetite, and the iron that goes to magnetite makes hydrogen. With solid solutions on, "
            "the serpentine takes about 4% of its iron sites as greenalite, which leaves less iron for magnetite and "
            "lowers the H2 from about 0.29 to 0.17 mol/kg. Brucite stays pure: the database has no Fe(OH)2.",
            {"title": "Serpentinization of mantle olivine", "t_start_c": 300.0, "t_end_c": 300.0,
             "p_start_bar": 300.0, "p_end_bar": 300.0, "steps": 20, "solid_solutions": True,
             "reactants": [{"name": "olivine fo90", "amount": 1000.0, "unit": "g"}],
             "blocked_minerals": DEFAULT_BLOCKED + ["antigorite"]},
            (("solid_solutions", "Solid solutions: Can form. Serpentine forms as a chrysotile-greenalite mixture "
                                 "(Mg and Fe on its three octahedral sites)."),
             ("reactants", "Reactant: 1000 g olivine Fo90 (90% forsterite), dissolving congruently."),
             ("blocked_minerals", "Blocked: antigorite, the high-temperature serpentine, which is too slow to "
                                  "form here (seafloor serpentinites are lizardite and chrysotile)."),
             ("t_start_c", "Held at 300 C and 300 bar.")),
            ("serpentine", "brucite", "magnetite")),
    Example("kinetics", "Calcite dissolving", "Kinetics",
            "One gram of calcite grains dissolves in rainwater at 25 C over twelve hours at its measured rate "
            "(Palandri and Kharaka 2004). The rate is fast at first and falls to zero as the water approaches "
            "saturation, so only about 1.4% of the calcite dissolves; pH rises from 5.6 to almost 10. The Rates "
            "tab shows the rate and what is left.",
            {"title": "Calcite dissolving at its measured rate", "fluid": "rainwater", "t_start_c": 25.0,
             "t_end_c": 25.0, "p_start_bar": 1.0, "p_end_bar": 1.0, "steps": 24, "kinetics": True,
             "duration": 12.0, "duration_unit": "hours", "surface_area_cm2_per_g": 100.0,
             "reactants": [{"name": "calcite", "amount": 1.0, "unit": "g"}]},
            (("kinetics", "How they react: At measured rates. The reactant is then what is there at the start, "
                          "not an amount added steadily."),
             ("duration", "Run for 12 hours: the x axis is time."),
             ("surface_area_cm2_per_g", "Surface area 100 cm2 per g (grains of about 0.1 mm); the rate is "
                                        "proportional to it."),
             ("reactants", "Reactant: 1 g calcite per kg of water.")),
            ("calcite left",), False),
    Example("forams", "A foram sinking below the CCD", "Kinetics along a fall",
            "One foraminifer test (20 µg of calcite) sinks from 4.5 to 5.5 km, below the calcite compensation depth, "
            "at about 200 m a day, so the fall takes 5 days. The run follows the foram: it sees the water it falls "
            "through, held at 1.5 C and compressed with depth (455 to 556 bar), so ever more corrosive to calcite. Its calcite dissolves at the measured rate (Palandri and Kharaka 2004). "
            "So little dissolves that the water around it hardly changes, as for a grain in the open ocean. The "
            "plot shows the foram's mass in µg; it loses about a fifth on the way down. The model puts this water "
            "further from saturation than measurements do and the rate law is linear in (1 - Omega), so this is "
            "an upper bound. A few seconds.",
            {"title": "A foram sinking below the CCD", "fluid": "custom",
             "analysis": "Na+ 10768\nK+ 399.1\nCa+2 412.3\nMg+2 1291.8\nCl- 19353\nSO4-2 2712\nHCO3- 2.3 mmolal",
             "units": "mg/kg", "ph": 7.75, "fluid_t_c": 25.0, "charge_balance": "none", "redox": "air",
             "t_start_c": 1.5, "t_end_c": 1.5, "p_start_bar": 455.0, "p_end_bar": 556.0, "steps": 26,
             "kinetics": True, "duration": 5.0, "duration_unit": "days", "surface_area_cm2_per_g": 100.0,
             "mode": "flow_through", "reactants": [{"name": "forams", "amount": 2e-5, "unit": "g"}]},
            (("analysis", "The water: deep seawater with the extra CO2 of old bottom water (2.3 mmol/kg carbonate, "
                          "pH 7.75 measured at 25 C)."),
             ("p_start_bar", "Pressure 455 to 556 bar: the foram's depth, 4.5 to 5.5 km (about 10 m per bar), at "
                             "1.5 C."),
             ("mode", "Flow-through: the foram falls through the water, so it keeps meeting new water and "
                      "leaves behind whatever it gives off. (One foram barely changes 1 kg of water, so static "
                      "gives the same numbers.)"),
             ("reactants", "Reactant: one foram, 20 µg of calcite, in 1 kg of the water around it."),
             ("kinetics", "How they react: At measured rates."),
             ("duration", "Run for 5 days: 1 km at about 200 m a day, a typical sinking speed for a test of a few "
                          "hundred µm."),
             ("surface_area_cm2_per_g", "Surface area 100 cm2 per g, about that of grains of 0.1 mm.")),
            ("calcite",), False, "ug"),
    Example("vent_energy", "Vent fluid meets seawater", "Chemosynthesis energy",
            "A 350 C black-smoker fluid mixes into 2 C deep seawater at a mid-ocean ridge (250 bar), from pure vent "
            "fluid down to 5 C (115 kg of seawater per kg). Mixing is faster than the redox reactions, so the fluid's "
            "H2, sulfide, methane and Fe(II) meet the seawater's O2 and sulfate without reacting, and no minerals "
            "form. The Energy tab gives what that disequilibrium is worth to microbes that let it react: the Gibbs "
            "energy of each chemosynthetic reaction and the energy it makes available per kg of mixed fluid (the "
            "calculation of McCollom and Shock 1997). Sulfate reduction offers the most where the mixture is above "
            "about 285 C, where little O2 has arrived; then iron oxidation (to about 190 C), hydrogen oxidation (to "
            "about 95 C) and, in the cool mixture, sulfide oxidation, up to about 43 J per kg near 5 C, limited by "
            "the seawater's O2. A few seconds.",
            {"title": "Vent fluid mixing into seawater: chemosynthesis energy", "fluid": "custom",
             "analysis": "Na+ 432\nK+ 23.2\nLi+ 0.89\nCa+2 15.6\nSr+2 0.08\nMn+2 0.96\nFe+2 1.66\nAl+3 0.005\n"
                         "SiO2 17.3\nCl- 489\nHCO3- 5.7\nH2S 6.5\nH2 0.5\nCH4 0.08",
             "units": "mmol/kg", "ph": 3.4, "fluid_t_c": 25.0, "charge_balance": "auto", "redox": "H2",
             "t_start_c": 350.0, "t_end_c": 350.0, "p_start_bar": 250.0, "p_end_bar": 250.0, "steps": 70,
             "mix_with": "deep_seawater", "mix_ratio": 115.0, "mix_t_c": 2.0, "mix_spacing": "temperature",
             "redox_apart": True, "minerals_form": False, "gas_phase": False, "energy": True},
            (("redox_apart", "Redox couples: Held apart. H2, sulfide, methane and Fe(III) keep their own totals, apart "
                             "from O2, sulfate, carbonate and Fe(II), so the mixture keeps the disequilibrium."),
             ("minerals_form", "Minerals: Suppressed, so the mixing is conservative (anhydrite and the metal "
                               "sulfides that a real chimney makes are left dissolved)."),
             ("energy", "Catabolic energy: On. The Energy tab: ΔG of each reaction (kJ/mol) and the energy available "
                        "per kg of mixed fluid and per kg of vent fluid."),
             ("mix_with", "Mixed with deep seawater at 2 C, up to 115 kg per kg of vent fluid; the temperature "
                          "follows the mixing."),
             ("mix_spacing", "Steps even in temperature, 5 C apart from 350 to 5 C (even amounts of seawater "
                             "would put nearly every step below 50 C)."),
             ("analysis", "Vent fluid: the Vent fluid preset (rounded East Pacific Rise 21 N, quenched pH 3.4 at "
                          "25 C, an illustrative 0.5 mmol/kg H2) with 0.08 mmol/kg of methane."),
             ("p_start_bar", "250 bar, about 2.5 km of water.")),
            ("sulfide oxidation energy", "hydrogen oxidation energy", "iron oxidation energy",
             "methane oxidation energy", "sulfate reduction energy", "anaerobic methane oxidation energy",
             "methanogenesis energy")),
)}


def example_card(key: str) -> Scenario:
    """The card of one of EXAMPLES."""
    try:
        return Scenario(**EXAMPLES[key].changes)
    except KeyError:
        raise KeyError(f"unknown example {key!r}; known examples: {', '.join(EXAMPLES)}") from None


# What the language model is told about the card.
FIELD_HELP = {
    "title": "short name for the run (text); leave out unless the user gives one",
    "fluid": f"starting fluid, one of: {', '.join(PRESETS)}, or {CUSTOM} for an analysis the user gives",
    "ph": "pH of the fluid (number); changing it makes the fluid custom",
    "fluid_t_c": "temperature the fluid's analysis and pH were measured at, degrees C",
    "fluid_p_bar": "pressure the fluid's analysis and pH refer to, bar; leave out for 1 bar (or the boiling "
                   "curve above 100 C); needed above 300 C",
    "charge_balance": "auto (adjust the main ion that is short), none, or an ion name such as Cl-",
    "redox": "redox state of the fluid: auto, none, air, log_fo2, eh, pe, or a redox species of the "
             "analysis used as the couple (O2, H2, Fe+3, H2S)",
    "redox_value": "the log fO2, Eh (volts) or pe when redox is log_fo2, eh or pe",
    "t_start_c": "temperature the run starts at, degrees C (0 to 600)",
    "t_end_c": "temperature the run ends at, degrees C; equal to the start means temperature is held fixed",
    "p_start_bar": "pressure the run starts at, bar (1 to 5000)",
    "p_end_bar": "pressure the run ends at, bar; equal to the start means pressure is held fixed",
    "steps": "number of steps along the path (2 to 200)",
    "mode": "static = minerals stay in contact with the fluid; flow_through = each step's minerals "
            "are removed and cannot redissolve; both = run the two for comparison",
    "blocked_minerals": "list of mineral names that may not precipitate (kinetically blocked)",
    "activity_model": "bdot (default, to about seawater salinity), davies (dilute), or pitzer (brines and "
                      "evaporites, to saturation with halite; parameters to about 200 C)",
    "evaporate": "concentration factor to evaporate the fluid to (1 = no evaporation; 10 = a tenth of the "
                 "water left)",
    "mix_with": f"a second fluid mixed in steadily ({', '.join(PRESETS)}), or empty for none",
    "mix_ratio": "kg of the second fluid per kg of the first at the end of the run",
    "mix_t_c": "temperature of the second fluid, degrees C; the mixture's temperature then follows the mixing",
    "mix_spacing": "ratio (default): steps even in the amount of second fluid; temperature: steps even in the "
                   "mixture's temperature (needs mix_t_c), which resolves the hot end of a vent-seawater mixture",
    "redox_apart": "false (default): redox couples react to equilibrium; true: H2, sulfide, methane, Fe(III) and nitrate "
                   "are each held apart from O2, sulfate, carbonate, Fe(II) and ammonia, so a mixture keeps its redox "
                   "disequilibrium (no electron transfer)",
    "minerals_form": "true (default): supersaturated minerals precipitate; false: no mineral forms (conservative mixing)",
    "energy": "false (default); true: the Gibbs energy of chemosynthetic reactions (H2, sulfide, methane and Fe(II) "
              "oxidation, methanogenesis, sulfate reduction, anaerobic methane oxidation) and the energy they make "
              "available per kg of fluid at each step; meaningful with redox_apart",
    "gas_phase": "true (default): a gas phase forms when dissolved gases exceed the pressure; false: none",
    "methane_sulfate": "apart (default): methane and sulfate do not react with each other; react: "
                       "they reach equilibrium (thermochemical sulfate reduction)",
    "reactants": "minerals or species added steadily along the path, per kg of water, e.g. "
                 "[{\"name\": \"albite\", \"amount\": 10, \"unit\": \"g\"}]; with kinetics, "
                 "what is there at the start",
    "kinetics": "false (default): reactants are added steadily and react completely; true: minerals and "
                "basalt glass react at measured rates over the duration (time becomes the axis)",
    "follow_phase": "where a saline fluid splits into brine and vapour: whole (default) follows the whole "
                    "fluid; brine lets the vapour boil off at each step; vapour follows the vapour and leaves "
                    "the brine behind",
    "exchange_meq": "cation exchange capacity in meq per kg of water (0 = none); the exchanger starts in "
                    "equilibrium with the fluid (Gaines-Thomas, PHREEQC constants)",
    "hfo_g": "grams of hydrous ferric oxide per kg of water (0 = none): surface complexation on its strong and "
             "weak sites (Dzombak and Morel 1990)",
    "kaolinite_g": "grams of kaolinite per kg of water (0 = none): exchange sites at 20 meq/kg of clay",
    "illite_g": "grams of illite per kg of water (0 = none): exchange sites at 225 meq/kg of clay, and edge sites",
    "montmorillonite_g": "grams of montmorillonite (smectite) per kg of water (0 = none): exchange sites at "
                         "764 meq/kg of clay, and edge sites",
    "clay_edges": "true (default): illite and montmorillonite also bind Ni, Co, Zn, Cd, Pb, Cu, Mn, Fe(II) and Hg "
                  "on their edge sites (2SPNE model, Bradbury and Baeyens 2005, 2009); false: exchange sites only",
    "kd_medium": "which partition coefficients (Allison and Allison 2005) the Kd solid uses: sediment (bed sediment, "
                 "default), suspended (suspended matter in a stream) or soil",
    "kd_solid_g": "grams of solid per kg of water holding metals by their partition coefficient Kd (0 = none); "
                  "with clays or exchange sites it holds only the metals they cannot hold (As and Cr with clay edge sites)",
    "isotopes": "false (default); true: d18O, dD, d13C, d34S and 87Sr/86Sr of the fluid and minerals, "
                "along the run",
    "isotope_inputs": "starting isotope values, e.g. {\"fluid\": {\"d18O\": -5}, \"basalt glass\": {\"Sr\": 120, "
                      "\"87Sr/86Sr\": 0.703}}; d values in permil (VSMOW, VPDB, VCDT), Sr in mg/kg of fluid or ppm "
                      "of a solid",
    "solid_solutions": "false (default): minerals are pure; true: ideal solid solutions (plagioclase, alkali "
                       "feldspar, carbonate, chlorite, epidote, garnet, olivine, pyroxenes, biotite, actinolite) "
                       "form in place of their endmembers",
    "uncertainty_runs": "0 (default) for a single run; otherwise how many extra runs (5 to 200, 20 is a good "
                        "start) with every log K shifted at random by its uncertainty, giving a shaded range "
                        "on the plots",
    "kind": "path (default): one fluid along a temperature-pressure path; transport: the fluid flows through a "
            "column of cells (1D reactive transport), from the start temperature and pressure at the inlet to the "
            "end values at the outlet",
    "cells": "transport: number of cells in the column (2 to 100, default 20)",
    "column_m": "transport: length of the column in metres (default 1)",
    "dispersivity_m": "transport: dispersivity in metres (default 0.05; 0 = sharp fronts)",
    "pore_volumes": "transport: how many pore volumes of the fluid flow through (default 2)",
    "pore_fluid": f"transport: the fluid the column starts filled with ({', '.join(PRESETS)}), or empty for the "
                  "same fluid as flows in",
    "duration": "how long a kinetic run lasts (number, in duration_unit); in a column, how long the flow lasts "
                "(the pore volumes pass in that time)",
    "duration_unit": "seconds, minutes, hours, days or years",
    "surface_area_cm2_per_g": "reactive surface of the reactants, cm2 per gram (default 100, about "
                              "0.1 mm grains)",
}

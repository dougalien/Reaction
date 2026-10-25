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
from . import analysis
from .analysis import DEFAULT_UNIT
from .fluids import CUSTOM, CUSTOM_LABEL, PRESETS, REDOX_CHOICES, FluidPreset, label as fluid_label

MODES = ("static", "flow_through", "both")
MODE_LABELS = {"static": "Static", "flow_through": "Flow-through", "both": "Both"}
ACTIVITY_MODELS = ("bdot", "davies", "pitzer")
FLUID_FIELDS = ("analysis", "units", "ph", "fluid_t_c", "charge_balance", "redox",
                "redox_value")                                       # what a custom fluid is
REACTANT_UNITS = {"g": ("g", 1.0), "gram": ("g", 1.0), "grams": ("g", 1.0), "mg": ("g", 1e-3),
                  "kg": ("g", 1e3), "mol": ("mol", 1.0), "mole": ("mol", 1.0), "moles": ("mol", 1.0),
                  "mmol": ("mol", 1e-3), "umol": ("mol", 1e-6)}
MAX_REACTANTS = 6
MAX_UNCERTAINTY_RUNS = 200
FOLLOW_PHASES = ("whole", "brine", "vapour")
ISOTOPE_KEYS = {"d18o": "d18O", "18o": "d18O", "o18": "d18O", "d13c": "d13C", "13c": "d13C", "c13": "d13C",
                "d34s": "d34S", "34s": "d34S", "s34": "d34S", "sr": "Sr", "sr_ppm": "Sr", "srppm": "Sr",
                "sr_mg_per_kg": "Sr", "87sr/86sr": "87Sr/86Sr", "87sr": "87Sr/86Sr", "sr87": "87Sr/86Sr",
                "sr_ratio": "87Sr/86Sr"}
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
            new.redox, new.redox_value = p.redox, p.redox_value
            new.fluid = CUSTOM
        for key, value in changes.items():
            if key not in known or key == "kind":
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

    def tp(self, f: float) -> tuple[float, float]:
        """(temperature in C, pressure in bar) at the fraction ``f`` (0 to 1) of the way along the path."""
        p = self.p_start_bar + (self.p_end_bar - self.p_start_bar) * f
        if self.mix_follows_temperature:
            m = self.mix_ratio * f
            return (self.t_start_c + self.mix_t_c * m) / (1.0 + m), p
        return self.t_start_c + (self.t_end_c - self.t_start_c) * f, p

    def points(self) -> list[tuple[float, float]]:
        """(temperature in C, pressure in bar) at each step, evenly spaced along the path."""
        n = max(int(self.steps), 2)
        return [self.tp(i / (n - 1)) for i in range(n)]

    @property
    def display_title(self) -> str:
        if self.title.strip():
            return self.title.strip()
        return f"{fluid_label(self.fluid)} run"

    def fluid_spec(self) -> FluidPreset:
        """The starting fluid, preset or custom, in one shape."""
        if self.fluid in PRESETS:
            return PRESETS[self.fluid]
        return FluidPreset(CUSTOM, CUSTOM_LABEL, "A fluid entered as an analysis.", self.fluid_t_c,
                           self.analysis, self.ph if self.ph is not None else float("nan"),
                           "entered on the card", units=self.units, charge_balance=self.charge_balance,
                           redox=self.redox, redox_value=self.redox_value)

    def fluid_analysis(self) -> analysis.Analysis:
        spec = self.fluid_spec()
        return analysis.parse(spec.analysis, spec.units)

    def redox_setting(self) -> tuple[str, Any]:
        """How the starting fluid's redox state is set: ("none", None), ("air", None),
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
        if self.kind != "path":
            out.append(f"only temperature-pressure path runs exist so far, not '{self.kind}'")
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
        if not (2 <= self.steps <= 200):
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
        if ranges_ok:
            if not self.varies and not self.reactants and not self.mixing and not self.kinetics \
                    and not self.evaporating:
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
        if self.follow_phase not in FOLLOW_PHASES:
            out.append(f"follow_phase must be one of {', '.join(FOLLOW_PHASES)}")
        if self.mode not in MODES:
            out.append(f"mode must be one of {', '.join(MODES)}")
        if self.activity_model not in ACTIVITY_MODELS:
            out.append(f"activity_model must be one of {', '.join(ACTIVITY_MODELS)}")
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
        if not (0.0 <= self.fluid_t_c <= 300.0):
            out.append("the analysis temperature must be between 0 and 300 C")
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
        text = self._describe_path()
        if self.reactants:
            added = ", ".join(f"{r['amount']:g} {r['unit']} {r['name']}" for r in self.reactants)
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
        if self.uncertainty_runs:
            text = text.replace(f" ({self._how()}", f", with a log K uncertainty range from {self.uncertainty_runs} "
                                                   f"runs ({self._how()}")
        if self.mixing:
            other = fluid_label(self.mix_with).lower()
            hot = f" at {self.mix_t_c:g} °C" if self.mix_t_c is not None else ""
            text = text.replace(f" ({self._how()}",
                                f", mixed with up to {self.mix_ratio:g} kg of {other}{hot} per kg ({self._how()}")
        return text

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
        return "none", None
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
               "surface_area_cm2_per_g", "evaporate"):
        return float(value)
    if key == "duration_unit":
        low = str(value).strip().lower()
        low = TIME_WORDS.get(low, TIME_WORDS.get(low.rstrip("s"), low))
        if low not in TIME_UNITS:
            raise ValueError(f"'{value}' is not a time unit (use {', '.join(TIME_UNITS)})")
        return low
    if key == "mix_t_c":
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
    if key in ("steps", "uncertainty_runs"):
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
                    raise ValueError(f"'{k}' is not one of d18O, d13C, d34S, Sr, 87Sr/86Sr")
                here[key_] = float(v)
                if key_ == "87Sr/86Sr" and not 0.69 < here[key_] < 0.80:
                    raise ValueError(f"87Sr/86Sr {v} is outside 0.69 to 0.80")
                if key_ == "Sr" and here[key_] < 0:
                    raise ValueError("Sr cannot be negative")
            out[str(source).strip().lower()] = here
        return out
    if key in ("gas_phase", "kinetics", "solid_solutions", "isotopes"):
        if isinstance(value, str):
            low = value.strip().lower()
            if low in ("true", "yes", "on", "1", "allowed", "allow"):
                return True
            if low in ("false", "no", "off", "0", "none", "suppressed"):
                return False
            raise ValueError(f"'{value}' is not yes or no")
        return bool(value)
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
            m = (re.fullmatch(r"(?P<amount>\d*\.?\d+(?:e-?\d+)?)\s*(?P<unit>[a-z]+)\s+(?:of\s+)?(?P<name>.+)", part, re.I)
                 or re.fullmatch(r"(?P<name>.+?)\s+(?P<amount>\d*\.?\d+(?:e-?\d+)?)\s*(?P<unit>[a-z]+)", part, re.I))
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
            raise ValueError(f"unknown unit '{unit_text}' for {item['name']} (use g, mg, kg, mol or mmol)")
        unit, scale = REACTANT_UNITS[unit_text]
        out.append({"name": str(item["name"]).strip(), "amount": float(item.get("amount", 0)) * scale,
                    "unit": unit})
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


# What the language model is told about the card.
FIELD_HELP = {
    "title": "short name for the run (text); leave out unless the user gives one",
    "fluid": f"starting fluid, one of: {', '.join(PRESETS)}, or {CUSTOM} for an analysis the user gives",
    "ph": "pH of the fluid (number); changing it makes the fluid custom",
    "fluid_t_c": "temperature the fluid's analysis and pH were measured at, degrees C",
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
    "isotopes": "false (default); true: d18O, d13C, d34S and 87Sr/86Sr of the fluid and minerals along the run",
    "isotope_inputs": "starting isotope values, e.g. {\"fluid\": {\"d18O\": -5}, \"basalt glass\": {\"Sr\": 120, "
                      "\"87Sr/86Sr\": 0.703}}; d values in permil (VSMOW, VPDB, VCDT), Sr in mg/kg of fluid or ppm "
                      "of a solid",
    "solid_solutions": "false (default): minerals are pure; true: ideal solid solutions (plagioclase, alkali "
                       "feldspar, carbonate, chlorite, epidote, garnet, olivine, pyroxenes, biotite, actinolite) "
                       "form in place of their endmembers",
    "uncertainty_runs": "0 (default) for a single run; otherwise how many extra runs (5 to 200, 20 is a good "
                        "start) with every log K shifted at random by its uncertainty, giving a shaded range "
                        "on the plots",
    "duration": "how long a kinetic run lasts (number, in duration_unit)",
    "duration_unit": "seconds, minutes, hours, days or years",
    "surface_area_cm2_per_g": "reactive surface of the reactants, cm2 per gram (default 100, about "
                              "0.1 mm grains)",
}

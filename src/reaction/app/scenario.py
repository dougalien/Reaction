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
from dataclasses import asdict, dataclass, field, fields
from typing import Any

from ..thermo import water
from .fluids import PRESETS

MODES = ("static", "flow_through", "both")
MODE_LABELS = {"static": "Static", "flow_through": "Flow-through", "both": "Both"}
ACTIVITY_MODELS = ("bdot", "davies")

T_MIN, T_MAX = 0.0, 600.0            # degrees C
P_MIN, P_MAX = 1.0, 5000.0           # bar (the range the HKF equations are calibrated for)

# carbonates that supersaturated seawater "wants" to form but that do not form on laboratory and
# field time scales at these temperatures (the usual practice is to block them)
DEFAULT_BLOCKED = ["dolomite", "dolomite,ordered", "dolomite,disordered", "magnesite", "huntite"]

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
        if abs(self.t_end_c - self.t_start_c) > 1e-9:
            out.append("temperature_c")
        if abs(self.p_end_bar - self.p_start_bar) > 1e-9:
            out.append("pressure_bar")
        return out

    def points(self) -> list[tuple[float, float]]:
        """(temperature in C, pressure in bar) at each step, evenly spaced along the path."""
        n = max(int(self.steps), 2)
        return [(self.t_start_c + (self.t_end_c - self.t_start_c) * i / (n - 1),
                 self.p_start_bar + (self.p_end_bar - self.p_start_bar) * i / (n - 1))
                for i in range(n)]

    @property
    def display_title(self) -> str:
        if self.title.strip():
            return self.title.strip()
        label = PRESETS[self.fluid].label if self.fluid in PRESETS else self.fluid
        return f"{label} run"

    # ---- checking -------------------------------------------------------------------------
    def problems(self) -> list[str]:
        out: list[str] = []
        if self.kind != "path":
            out.append(f"only temperature-pressure path runs exist so far, not '{self.kind}'")
        if self.fluid not in PRESETS:
            out.append(f"unknown fluid '{self.fluid}' (known: {', '.join(PRESETS)})")
        ranges_ok = True
        for name, vals, lo, hi, unit in (
                ("temperature", (self.t_start_c, self.t_end_c), T_MIN, T_MAX, "C"),
                ("pressure", (self.p_start_bar, self.p_end_bar), P_MIN, P_MAX, "bar")):
            if not all(lo <= v <= hi for v in vals):
                out.append(f"{name} must be between {lo:g} and {hi:g} {unit}")
                ranges_ok = False
        if not (2 <= self.steps <= 200):
            out.append("steps must be between 2 and 200")
        if ranges_ok:
            if not self.varies:
                out.append("nothing changes along this path: give different start and end values "
                           "for temperature or pressure")
            boil = self._first_boiling_point()
            if boil:
                t, p, psat = boil
                out.append(f"at {t:.0f} C water boils below {psat:.0f} bar but the path has {p:.0f} bar "
                           f"there; keep the pressure at or above {math.ceil(psat)} bar")
        if self.mode not in MODES:
            out.append(f"mode must be one of {', '.join(MODES)}")
        if self.activity_model not in ACTIVITY_MODELS:
            out.append(f"activity_model must be one of {', '.join(ACTIVITY_MODELS)}")
        return out

    def _first_boiling_point(self) -> tuple[float, float, float] | None:
        """First place along the path (checked finely, not just at the steps) where water is steam."""
        n = 201
        for i in range(n):
            f = i / (n - 1)
            t = self.t_start_c + (self.t_end_c - self.t_start_c) * f
            p = self.p_start_bar + (self.p_end_bar - self.p_start_bar) * f
            t_k = t + 273.15
            if t_k < water.T_CRIT:
                psat = water.saturation_pressure_mpa(t_k) * water.BAR_PER_MPA
                if p < 0.99 * psat:
                    return t, p, psat
        return None

    def describe(self) -> str:
        label = PRESETS[self.fluid].label if self.fluid in PRESETS else self.fluid
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


def _coerce(key: str, value: Any) -> Any:
    if key in ("t_start_c", "t_end_c", "p_start_bar", "p_end_bar"):
        return float(value)
    if key == "steps":
        return int(round(float(value)))
    if key == "blocked_minerals":
        if isinstance(value, str):
            value = _rejoin_comma_names(value.replace(";", ",").split(","))
        if not isinstance(value, (list, tuple)):
            raise TypeError("expected a list of mineral names")
        return [str(v).strip().lower() for v in value if str(v).strip()]
    if key in ("title", "fluid", "mode", "activity_model"):
        text = str(value).strip()
        if key == "title":
            return text
        return text.lower().replace("-", "_").replace(" ", "_")
    raise ValueError("unknown field")


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
    "fluid": f"starting fluid, one of: {', '.join(PRESETS)}",
    "t_start_c": "temperature the run starts at, degrees C (0 to 600)",
    "t_end_c": "temperature the run ends at, degrees C; equal to the start means temperature is held fixed",
    "p_start_bar": "pressure the run starts at, bar (1 to 5000)",
    "p_end_bar": "pressure the run ends at, bar; equal to the start means pressure is held fixed",
    "steps": "number of steps along the path (2 to 200)",
    "mode": "static = minerals stay in contact with the fluid; flow_through = each step's minerals "
            "are removed and cannot redissolve; both = run the two for comparison",
    "blocked_minerals": "list of mineral names that may not precipitate (kinetically blocked)",
    "activity_model": "bdot (default, for brines) or davies",
}

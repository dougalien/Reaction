"""The scenario card: everything needed to describe one model run, in plain terms.

A card is a flat set of fields.  It is what the chat fills in, what the user can edit by hand, and
what is stored next to every result, so any run can be reproduced from its card alone.
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

# carbonates that supersaturated seawater "wants" to form but that do not form on laboratory and
# field time scales at these temperatures (the usual practice is to block them)
DEFAULT_BLOCKED = ["dolomite", "dolomite,ordered", "dolomite,disordered", "magnesite", "huntite"]

# oxides, elements and ice that the database lists but that are not sensible products here
ALWAYS_EXCLUDED = ("ice", "sulfur", "S_liq", "graphite", "potassium oxide", "sodium oxide",
                   "lime", "periclase")


@dataclass
class Scenario:
    title: str = "Seawater heating"
    kind: str = "heating"
    fluid: str = "seawater"
    t_start_c: float = 2.0
    t_end_c: float = 350.0
    steps: int = 36
    pressure_bar: float = 500.0
    mode: str = "static"
    blocked_minerals: list[str] = field(default_factory=lambda: list(DEFAULT_BLOCKED))
    activity_model: str = "bdot"

    # ---- editing --------------------------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Scenario":
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in known})

    def changed(self, changes: dict[str, Any]) -> tuple["Scenario", list[str]]:
        """Copy with ``changes`` applied.  Returns the new card and a list of problems; a problem
        in one field never blocks the others."""
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

    # ---- checking -------------------------------------------------------------------------
    def problems(self) -> list[str]:
        out: list[str] = []
        if self.kind != "heating":
            out.append(f"only 'heating' runs exist so far, not '{self.kind}'")
        if self.fluid not in PRESETS:
            out.append(f"unknown fluid '{self.fluid}' (known: {', '.join(PRESETS)})")
        temps_ok = 0.0 <= self.t_start_c < self.t_end_c <= 370.0
        if not temps_ok:
            out.append("temperatures must satisfy 0 <= start < end <= 370 C")
        if not (2 <= self.steps <= 200):
            out.append("steps must be between 2 and 200")
        if not (1.0 <= self.pressure_bar <= 1500.0):
            out.append("pressure must be between 1 and 1500 bar")
        elif temps_ok:
            t_k = self.t_end_c + 273.15
            if t_k < water.T_CRIT:
                psat = water.saturation_pressure_mpa(t_k) * water.BAR_PER_MPA
                if self.pressure_bar < psat:
                    out.append(f"at {self.t_end_c:g} C water boils below {psat:.0f} bar; use at least "
                               f"{math.ceil(psat)} bar or lower the end temperature")
        if self.mode not in MODES:
            out.append(f"mode must be one of {', '.join(MODES)}")
        if self.activity_model not in ACTIVITY_MODELS:
            out.append(f"activity_model must be one of {', '.join(ACTIVITY_MODELS)}")
        return out

    def describe(self) -> str:
        label = PRESETS[self.fluid].label if self.fluid in PRESETS else self.fluid
        how = {"static": "static", "flow_through": "flow-through",
               "both": "static and flow-through"}.get(self.mode, self.mode)
        return (f"{label} heated from {self.t_start_c:g} to {self.t_end_c:g} °C at "
                f"{self.pressure_bar:g} bar ({how}, {self.steps} steps)")


def _coerce(key: str, value: Any) -> Any:
    if key in ("t_start_c", "t_end_c", "pressure_bar"):
        return float(value)
    if key == "steps":
        return int(round(float(value)))
    if key == "blocked_minerals":
        if isinstance(value, str):
            value = [v for v in value.replace(";", ",").split(",")]
            # "dolomite,ordered" is one mineral name in the database
            value = _rejoin_comma_names(value)
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
    "title": "short name for the run (text)",
    "fluid": f"starting fluid, one of: {', '.join(PRESETS)}",
    "t_start_c": "temperature the fluid starts at, degrees C (0 to 370)",
    "t_end_c": "temperature the fluid is heated to, degrees C (up to 370)",
    "steps": "number of temperature steps (2 to 200)",
    "pressure_bar": "pressure in bar, held constant (must keep water liquid at the end temperature)",
    "mode": "static = minerals stay in contact with the fluid; flow_through = each step's minerals "
            "are removed and cannot redissolve; both = run the two for comparison",
    "blocked_minerals": "list of mineral names that may not precipitate (kinetically blocked)",
    "activity_model": "bdot (default, for brines) or davies",
}

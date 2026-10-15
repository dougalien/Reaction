"""Named starting fluids.  Each is a set of speciation constraints at a measurement temperature."""
from __future__ import annotations

from dataclasses import dataclass

from ..speciation import MgPerKg, Molal, PH


@dataclass(frozen=True)
class FluidPreset:
    key: str
    label: str
    description: str
    measured_t_c: float
    constraints: dict
    source: str
    caveats: tuple[str, ...] = ()


PRESETS: dict[str, FluidPreset] = {
    "seawater": FluidPreset(
        "seawater", "Seawater",
        "Standard seawater, major ions only (no silica, iron or trace metals).",
        25.0,
        {"Na+": MgPerKg(10768), "K+": MgPerKg(399.1), "Ca+2": MgPerKg(412.3),
         "Mg+2": MgPerKg(1291.8), "H+": PH(8.22), "HCO3-": Molal(2.022e-3),
         "Cl-": MgPerKg(19353), "SO4-2": MgPerKg(2712)},
        "composition of the EQ3/6 test case swmaj (pH and alkalinity measured at 25 C)",
        caveats=(
            "The magnesium hydroxide sulfate hydrate that forms when seawater is heated in "
            "experiments is not in the database, so magnesium removal and the acid it releases "
            "are underestimated here.",
            "The fluid contains no silica, iron or aluminium and the redox state is not tracked, "
            "so only carbonate, sulfate and hydroxide minerals can form.",
        )),
}


def preset(key: str) -> FluidPreset:
    try:
        return PRESETS[key]
    except KeyError:
        raise KeyError(f"unknown fluid {key!r}; known fluids: {', '.join(PRESETS)}") from None

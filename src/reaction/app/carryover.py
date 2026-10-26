"""Carrying the fluid at the end of one run into the next: "take that fluid and cool it back to 25 C".

The fluid is written onto the card as a custom analysis: each element's total dissolved molality on its
basis species, the pH, and the oxygen fugacity, all at the run's last point, with the temperature and
pressure of that point as where it was "measured".  Speciating that analysis there gives back the same
fluid.  Minerals and gas are left behind: only the fluid is carried.
"""
from __future__ import annotations

import math

from ..thermo.eosdb import DEFAULT_MASTERS
from .fluids import CUSTOM

# what the previous run used up or did to its fluid, which the next run must not do again
CONSUMED = {"reactants": [], "mix_with": "", "evaporate": 1.0, "kinetics": False, "kind": "path",
            "pore_fluid": ""}
TRACE = 1e-16                  # mol/kg; below this an element is left out of the analysis


def end_point(result: dict, mode: str | None = None) -> tuple[str, int, float, float] | None:
    """(mode, index, temperature C, pressure bar) of the last solved point of ``result``, or None."""
    modes = result.get("modes", {})
    if mode not in modes:
        mode = next((m for m in ("static", "outlet", "flow_through") if m in modes), next(iter(modes), None))
    if mode is None:
        return None
    run = modes[mode]
    ph = run["properties"]["pH"]
    solved = [i for i, v in enumerate(ph) if v is not None]
    if not solved:
        return None
    i = solved[-1]
    sc = result.get("scenario", {})
    axes = result.get("axes", {})
    if mode in ("outlet", "start", "profile") or "temperature_c" not in axes:
        t, p = sc.get("t_end_c"), sc.get("p_end_bar")    # a column: the fluid leaves at the outlet
    else:
        t, p = axes["temperature_c"]["values"][i], axes["pressure_bar"]["values"][i]
    return mode, i, float(t), float(p)


def fluid_from_result(result: dict, mode: str | None = None) -> tuple[dict, str]:
    """Card changes that make the fluid at the end of ``result`` the starting fluid, held at that point's
    temperature and pressure, and a sentence saying what was taken.  Raises ValueError when the result
    has no solved point."""
    end = end_point(result, mode)
    if end is None:
        raise ValueError("that run has no solved point to take the fluid from")
    mode, i, t, p = end
    run = result["modes"][mode]
    elements = {el: values[i] for el, values in run["elements"].items() if values[i] is not None}
    held_apart = result.get("scenario", {}).get("methane_sulfate", "apart") == "apart"
    methane = ((run["species"].get("CH4") or [0.0] * (i + 1))[i] or 0.0) if held_apart else 0.0
    lines, left_out = [], []
    for el, molal in sorted(elements.items()):
        master = DEFAULT_MASTERS.get(el)
        if el == "C" and methane > 0.0:
            molal -= methane                             # methane stays its own component
        if master is None:
            if molal > TRACE:
                left_out.append(el)
            continue
        if molal > TRACE:
            lines.append(f"{master} {molal:.6g}")
    if methane > TRACE:
        lines.append(f"CH4 {methane:.6g}")
    props = run["properties"]
    log_fo2 = (props.get("log_fO2") or [None] * (i + 1))[i]
    redox = {"redox": "log_fo2", "redox_value": float(log_fo2)} \
        if log_fo2 is not None and math.isfinite(log_fo2) else {"redox": "none", "redox_value": None}
    changes = {"fluid": CUSTOM, "analysis": "\n".join(lines), "units": "molal", "ph": float(props["pH"][i]),
               "fluid_t_c": t, "fluid_p_bar": p, "charge_balance": "none", **redox,
               "t_start_c": t, "t_end_c": t, "p_start_bar": p, "p_end_bar": p, **CONSUMED}
    where = {"static": "", "flow_through": " (flow-through)", "outlet": " at the column outlet"}.get(mode, "")
    note = (f"Took the fluid at the end of the last run{where} ({t:g} °C, {p:g} bar, pH {changes['ph']:.2f}) "
            f"as the starting fluid; the minerals that formed stay behind.")
    if left_out:
        note += " Left out (no basis species): " + ", ".join(left_out) + "."
    return changes, note

"""Plain-language summary of a result.  Every number in it is read from the result's own tables."""
from __future__ import annotations

from .names import pretty_mineral

ELEMENT_WORDS = {"Mg": "magnesium", "Ca": "calcium", "S": "sulfur", "C": "carbon", "K": "potassium",
                 "Na": "sodium", "Cl": "chlorine", "Si": "silicon", "Fe": "iron"}
APPEARS = 1e-6          # mol per kg of original water: below this a mineral is not "present"


def _first(values, predicate):
    for i, v in enumerate(values):
        if v is not None and predicate(v):
            return i
    return None


def _last_valid(values):
    for i in range(len(values) - 1, -1, -1):
        if values[i] is not None:
            return i
    return None


def _first_valid(values):
    for i, v in enumerate(values):
        if v is not None:
            return i
    return None


def _amount(mol: float) -> str:
    mmol = mol * 1000.0
    if mmol >= 100:
        return f"{mmol:.0f} mmol"
    if mmol >= 1:
        return f"{mmol:.1f} mmol"
    return f"{mmol:.2g} mmol"


def _retained(result: dict, run: dict, element: str) -> float | None:
    """Fraction of an element, relative to the starting fluid, still dissolved at the end."""
    series = run["elements"].get(element)
    water = run["properties"]["water_kg"]
    start = result.get("start_elements", {}).get(element)
    if not series or not start:
        return None
    i1 = _last_valid(series)
    if i1 is None or water[i1] is None:
        return None
    return (series[i1] * water[i1]) / start


def _axes(result: dict):
    """(temperatures, pressures, varying axes) of the run's path."""
    axes = result["axes"]
    return axes["temperature_c"]["values"], axes["pressure_bar"]["values"], result["varies"]


def _added(result: dict, i: int) -> str:
    axis = result["axes"]["mixed" if "mixed" in result["varies"] else "reacted"]
    if axis["key"] == "mixed":
        return f"{axis['values'][i]:.3g} kg of {axis['label'].lower().removesuffix(' added')} per kg"
    unit = axis["unit"].replace("per kg water", "per kg of water")
    return f"{axis['values'][i]:.3g} {unit} added"


def _time(result: dict, i: int) -> str:
    axis = result["axes"]["time"]
    value, unit = axis["values"][i], axis["unit"]
    return f"{value:.3g} {unit[:-1] if value == 1 else unit}"


def point_text(result: dict, i: int) -> str:
    """Where on the path point ``i`` is, mentioning only what varies."""
    temps, pressures, varies = _axes(result)
    if "reacted" in varies or "mixed" in varies or "time" in varies or "evaporation" in varies:
        if "evaporation" in varies:
            where = [f"{result['axes']['evaporation']['values'][i]:.3g} times concentrated"]
        else:
            where = [_time(result, i) if "time" in varies else _added(result, i)]
        if "temperature_c" in varies:
            where.append(f"{temps[i]:.0f} °C")
        if "pressure_bar" in varies:
            where.append(f"{pressures[i]:.0f} bar")
        return ", ".join(where)
    if "temperature_c" in varies and "pressure_bar" in varies:
        return f"{temps[i]:.0f} °C and {pressures[i]:.0f} bar"
    if "pressure_bar" in varies:
        return f"{pressures[i]:.0f} bar"
    return f"{temps[i]:.0f} °C"


def _when(result: dict, i: int) -> str:
    temps, pressures, varies = _axes(result)
    if i == 0:
        return "from the start"
    if "time" in varies:
        return f"after about {point_text(result, i)}"
    if "evaporation" in varies:
        return f"from about {point_text(result, i)}"
    if "reacted" in varies or "mixed" in varies:
        return f"from about {point_text(result, i)}"
    if "temperature_c" in varies and "pressure_bar" in varies:
        return f"from about {temps[i]:.0f} °C ({pressures[i]:.0f} bar)"
    if "pressure_bar" in varies:
        return f"at about {pressures[i]:.0f} bar"
    return f"from about {temps[i]:.0f} °C"


def _mineral_clauses(result: dict, run: dict, blocked: set[str]) -> list[str]:
    found = []
    kinetic = {k.split(":", 1)[1] for k in run.get("kinetic", {})}
    for name, series in run["minerals"].items():
        if name.lower() in blocked or name in kinetic:          # reactants, told of separately
            continue
        i = _first(series, lambda v: v > APPEARS)
        if i is None:
            continue
        end = series[_last_valid(series)]
        peak = max(v for v in series if v is not None)
        clause = f"{pretty_mineral(name).lower()} {_when(result, i)}"
        if end > APPEARS:
            clause += f" ({_amount(end)} per kg of original water at the end"
            clause += ")" if peak <= end * 1.05 else f", up to {_amount(peak)} on the way)"
        else:
            clause += f" (up to {_amount(peak)} per kg of original water, gone by the end)"
        mix = _composition_text(run, name, end > APPEARS)
        if mix:
            clause += f", {mix}"
        found.append((i, clause))
    return [c for _, c in sorted(found)]


def _composition_text(run: dict, name: str, at_end: bool) -> str:
    """'mostly anorthite (62%) at the end' for a solid solution; empty for a pure mineral."""
    parts = {k.split(":", 1)[1]: v for k, v in run.get("composition", {}).items() if k.split(":", 1)[0] == name}
    if not parts:
        return ""
    series = next(iter(parts.values()))
    idx = [i for i, v in enumerate(series) if v is not None]
    if not at_end:                            # where there was most of it
        amounts = run["minerals"][name]
        idx = [max(idx, key=lambda i: amounts[i] or 0.0)]
    i = idx[-1]
    shares = sorted(((v[i] or 0.0, k) for k, v in parts.items()), reverse=True)
    text = ", ".join(f"{pretty_mineral(k).lower()} {100 * x:.0f}%" for x, k in shares if x >= 0.005)
    return f"{text} {'at the end' if at_end else 'at its largest'}"


def _kinetic_clause(result: dict, run: dict) -> str:
    """How much of each kinetic reactant reacted, and when a mineral reached saturation."""
    parts = []
    names = result.get("catalog", {}).get("kinetic", {})
    for key, series in run.get("kinetic", {}).items():
        if not key.startswith("left:"):
            continue
        name = key.split(":", 1)[1]
        label = (names.get(key) or {}).get("label", name).removesuffix(" left")
        end = series[_last_valid(series)]
        change = 100.0 - end
        if abs(change) < 0.005:
            text = f"the {label.lower()} barely reacted"
        elif change > 0:
            text = f"{change:.2g}% of the {label.lower()} dissolved"
        else:
            text = f"the {label.lower()} grew by {-change:.2g}%"
        si = run["saturation"].get(name)
        if si:
            i = _first(si[1:], lambda v: abs(v) < 0.02)
            if i is not None:
                text += f" (at saturation by {point_text(result, i + 1)})"
        parts.append(text)
    return ("At the measured rates, " + "; ".join(parts) + ".") if parts else ""


def _gas_clause(result: dict, run: dict, mode: str) -> str:
    """When a gas phase forms, and what it is made of at the end."""
    gas = run.get("gas") or {}
    if not gas:
        return ""
    n = len(next(iter(gas.values())))
    totals = [sum((series[i] or 0.0) for series in gas.values()) for i in range(n)]
    i = _first(totals, lambda v: v > APPEARS)
    if i is None:
        return ""
    end = _last_valid(totals)
    gone = totals[end] <= APPEARS
    at = max(range(n), key=lambda k: totals[k]) if gone else end
    total = totals[at]
    parts = sorted(((series[at] or 0.0) / total, name) for name, series in gas.items())[::-1]
    names = result.get("catalog", {}).get("gas", {})
    what = ", ".join(f"{(names.get(k) or {}).get('label', k)} {100 * x:.0f}%" for x, k in parts[:3] if x > 0.005)
    verb = "is vented" if mode == "flow_through" else "forms"
    amount = (f"up to {_amount(total)} per kg of original water, gone by the end" if gone else
              f"{_amount(total)} per kg of original water by the end")
    return f"A gas phase {verb} {_when(result, i)}: {amount}; by moles {what}."


def _split_clause(result: dict, run: dict) -> str:
    """Boiling or condensation of the fluid into brine and vapour (NaCl-H2O)."""
    state = run.get("phase_state") or []
    split = [i for i, s in enumerate(state) if s == "liquid + vapour"]
    if not split:
        return ""
    ph = run["phases"]
    i0, i1 = split[0], split[-1]
    text = (f"The fluid splits into brine and vapour {_when(result, i0)}: there "
            f"{ph['vapour_share'][i0]:.0f}% of its mass is vapour with {ph['vapour_salinity'][i0]:.2g} wt% NaCl, "
            f"the brine has {ph['brine_salinity'][i0]:.3g} wt%")
    if i1 > i0:
        text += (f"; at {point_text(result, i1)} the vapour is {ph['vapour_share'][i1]:.0f}% and the brine "
                 f"{ph['brine_salinity'][i1]:.3g} wt%")
    if i1 < len(state) - 1 and state[i1 + 1] == "vapour + halite":
        text += "; beyond that the brine dries out to vapour and halite"
    elif i1 < len(state) - 1 and state[i1 + 1] == "one phase":
        text += f"; it is one fluid again from {point_text(result, i1 + 1)}"
    return text + "."


def story(result: dict) -> list[str]:
    sc = result["scenario"]
    blocked = {b.lower() for b in sc["blocked_minerals"]}
    modes = result["modes"]
    paras: list[str] = []

    first_key = next(iter(modes))
    ph = modes[first_key]["properties"]["pH"]
    i0, i1 = _first_valid(ph), _last_valid(ph)
    paras.append(
        f"{result['summary'][0].upper()}{result['summary'][1:]}. The in situ pH goes from "
        f"{ph[i0]:.1f} at {point_text(result, i0)} to {ph[i1]:.1f} at {point_text(result, i1)}.")

    labels = {"static": "Static", "flow_through": "Flow-through"}
    for key, run in modes.items():
        clauses = _mineral_clauses(result, run, blocked)
        lost = []
        for el in ("S", "Mg", "Ca", "C"):
            kept = _retained(result, run, el)
            if kept is not None and kept < 0.95:
                lost.append(f"{(1 - kept) * 100:.0f}% of the {ELEMENT_WORDS[el]}")
        sentence = f"{labels[key]}: "
        kinetic = _kinetic_clause(result, run)
        if kinetic:
            sentence += kinetic + " "
        if clauses:
            sentence += ("Minerals" if kinetic else "minerals") + " form in this order: " + "; ".join(clauses) + ". "
        else:
            sentence += ("No" if kinetic else "no") + " minerals form. "
        if lost:
            sentence += "By the end the fluid has lost " + ", ".join(lost) + ". "
        gas = _gas_clause(result, run, key)
        if gas:
            sentence += gas
        split = _split_clause(result, run)
        if split:
            sentence = sentence.rstrip() + " " + split
        paras.append(sentence.strip())

    if len(modes) == 2:
        a, b = modes["static"], modes["flow_through"]
        diffs = []
        for el in ("S", "Mg", "Ca", "C"):
            ka, kb = _retained(result, a, el), _retained(result, b, el)
            if ka is not None and kb is not None and abs(ka - kb) > 0.05:
                diffs.append(f"{ELEMENT_WORDS[el]} ({ka * 100:.0f}% left when static, {kb * 100:.0f}% "
                             f"when flow-through)")
        pa, pb = a["properties"]["pH"], b["properties"]["pH"]
        ea, eb = _last_valid(pa), _last_valid(pb)
        ph_note = ""
        if ea is not None and eb is not None and abs(pa[ea] - pb[eb]) > 0.1:
            ph_note = f" The final pH differs: {pa[ea]:.1f} static, {pb[eb]:.1f} flow-through."
        if diffs:
            paras.append("Where the two modes differ: " + "; ".join(diffs) + "." + ph_note)
        else:
            paras.append("Static and flow-through leave almost the same fluid at the end." + ph_note)

    held = []
    for name in sc["blocked_minerals"]:
        for run in modes.values():
            sat = run["saturation"].get(name)
            if sat and max((v for v in sat if v is not None), default=-9) > 0.1:
                held.append(pretty_mineral(name).lower())
                break
    if held:
        paras.append("Held out as kinetically blocked although the fluid is supersaturated in them: "
                     + ", ".join(held) + ".")

    levels = result.get("reliability", {}).get("levels", [])
    if any(level != "good" for level in levels):
        worst = "poor" if "poor" in levels else "caution"
        where = "some" if worst == "caution" else "part"
        paras.append(f"{where.capitalize()} of this path lies where the thermodynamics of charged species is "
                     f"less certain ({'unreliable' if worst == 'poor' else 'approximate'} there): see "
                     f"Assumptions and limits, and the coloured strip under the plot.")
    return paras

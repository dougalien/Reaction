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


def _mineral_clauses(run: dict, temps: list, blocked: set[str]) -> list[str]:
    found = []
    for name, series in run["minerals"].items():
        if name.lower() in blocked:
            continue
        i = _first(series, lambda v: v > APPEARS)
        if i is None:
            continue
        end = series[_last_valid(series)]
        peak = max(v for v in series if v is not None)
        when = "from the start" if i == 0 else f"from about {temps[i]:.0f} °C"
        clause = f"{pretty_mineral(name).lower()} {when}"
        if end > APPEARS:
            clause += f" ({_amount(end)} per kg of original water at the end"
            clause += ")" if peak <= end * 1.05 else f", up to {_amount(peak)} on the way)"
        else:
            clause += f" (up to {_amount(peak)} per kg of original water, gone by the end)"
        found.append((i, clause))
    return [c for _, c in sorted(found)]


def story(result: dict) -> list[str]:
    sc = result["scenario"]
    temps = result["x"]["values"]
    blocked = {b.lower() for b in sc["blocked_minerals"]}
    modes = result["modes"]
    paras: list[str] = []

    first_key = next(iter(modes))
    props = modes[first_key]["properties"]["pH"]
    i0, i1 = _first_valid(props), _last_valid(props)
    paras.append(
        f"{result['summary'][0].upper()}{result['summary'][1:]}. The in situ pH goes from "
        f"{props[i0]:.1f} at {temps[i0]:.0f} °C to {props[i1]:.1f} at {temps[i1]:.0f} °C.")

    labels = {"static": "Static", "flow_through": "Flow-through"}
    for key, run in modes.items():
        clauses = _mineral_clauses(run, temps, blocked)
        lost = []
        for el in ("S", "Mg", "Ca", "C"):
            kept = _retained(result, run, el)
            if kept is not None and kept < 0.95:
                lost.append(f"{(1 - kept) * 100:.0f}% of the {ELEMENT_WORDS[el]}")
        sentence = f"{labels[key]}: "
        if clauses:
            sentence += "minerals form in this order: " + "; ".join(clauses) + ". "
        else:
            sentence += "no minerals form. "
        if lost:
            sentence += "By the end the fluid has lost " + ", ".join(lost) + "."
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
    return paras

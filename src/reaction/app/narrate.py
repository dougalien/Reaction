"""Plain-language summary of a result.  Every number in it is read from the result's own tables."""
from __future__ import annotations

import re

from ..thermo.eosdb import ATOMIC_WEIGHTS
from .names import pretty_mineral
from .scenario import Scenario

ELEMENT_WORDS = {"Mg": "magnesium", "Ca": "calcium", "S": "sulfur", "C": "carbon", "K": "potassium",
                 "Na": "sodium", "Cl": "chlorine", "Si": "silicon", "Fe": "iron", "Al": "aluminium",
                 "Mn": "manganese", "P": "phosphorus", "Br": "bromine", "Sr": "strontium", "B": "boron",
                 "N": "nitrogen", "F": "fluorine", "Ba": "barium", "Li": "lithium", "Zn": "zinc", "Cu": "copper"}
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


NUMBER_WORDS = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten"]


def _count(n: int, word: str) -> str:
    return f"{NUMBER_WORDS[n] if n < len(NUMBER_WORDS) else n} {word}{'' if n == 1 else 's'}"


def _listing(items: list[str]) -> str:
    """'a', 'a and b', 'a, b and c'."""
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def mineral_word(name: str) -> str:
    """A mineral's name for a sentence: 'albite,low' -> 'low albite', no brackets."""
    if "," in name:
        base, rest = name.split(",", 1)
        return f"{rest.replace(',', ' ')} {base}".lower()
    return name.lower()


def plain(text: str) -> str:
    """A sentence without brackets: 'heat seawater (static)' -> 'heat seawater, static'."""
    text = re.sub(r"\s*\(([^()]*)\)", lambda m: ", " + m.group(1) if m.group(1) else "", text)
    return re.sub(r",\s*,", ",", text)


def _sentence(text: str) -> str:
    return text[:1].upper() + text[1:]


def _mineral_sentences(result: dict, run: dict, blocked: set[str]) -> tuple[list[str], list[str]]:
    """(names in the order they appear, one or two plain sentences about each)."""
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
        word = mineral_word(name)
        text = f"{_sentence(word)} forms {_when(result, i)}."
        if end > APPEARS:
            text += f" At the end there is {_amount(end)} of it per kg of original water"
            text += "." if peak <= end * 1.05 else f", after reaching {_amount(peak)} on the way."
        else:
            text += f" It reaches {_amount(peak)} per kg of original water and has dissolved again by the end."
        mix = _composition_text(run, name, end > APPEARS)
        if mix:
            text += f" Its make-up is {mix}."
        found.append((i, word, text))
    found.sort(key=lambda f: f[0])
    return [w for _, w, _ in found], [t for _, _, t in found]


def _composition_text(run: dict, name: str, at_end: bool) -> str:
    """'anorthite 62%, albite 38% at the end' for a solid solution; empty for a pure mineral."""
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
    text = ", ".join(f"{mineral_word(k)} {100 * x:.0f}%" for x, k in shares if x >= 0.005)
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
                text += f", reaching saturation by {point_text(result, i + 1)}"
        parts.append(text)
    return ("At the measured rates, " + "; ".join(parts) + ".") if parts else ""


def _energy_clause(result: dict, run: dict) -> str:
    """Which catabolic reaction offers the most energy where, and its largest value."""
    energy = run.get("energy") or {}
    names = result.get("catalog", {}).get("energy", {})
    series = {k[2:]: v for k, v in energy.items() if k.startswith("E:")}
    if not series:
        return ""
    label = {k: (names.get(f"E:{k}") or {}).get("label", k).removesuffix(" energy").lower() for k in series}
    n = len(next(iter(series.values())))
    lead = []
    for i in range(n):
        vals = {k: v[i] for k, v in series.items() if v[i] is not None and v[i] > 0}
        lead.append(max(vals, key=vals.get) if vals else None)
    runs = []                                        # stretches with one reaction ahead
    for i, k in enumerate(lead):
        if k is None:
            continue
        if runs and runs[-1][0] == k and runs[-1][2] == i - 1:
            runs[-1][2] = i
        else:
            runs.append([k, i, i])
    if not runs:
        return "None of the catabolic reactions has energy to offer anywhere along this path."
    temps, _, varies = _axes(result)

    def at(i):
        return f"{temps[i]:.0f} °C" if "temperature_c" in varies else point_text(result, i)

    parts = []
    for k, a, b in runs:
        parts.append(f"{label[k]} at {at(a)}" if a == b else f"{label[k]} from {at(a)} to {at(b)}")
    top = max(series, key=lambda k: max((v or 0.0) for v in series[k]))
    i = max(range(n), key=lambda j: series[top][j] or 0.0)
    return ("Catabolic energy: the reaction offering the most is " + "; then ".join(parts) + ". The largest is "
            f"{label[top]}, {series[top][i]:.3g} J per kg of fluid at {at(i)}.")


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
    follow = result["scenario"].get("follow_phase", "whole")
    if follow == "brine":
        text += ". From there the calculation follows the brine as the vapour boils off"
    elif follow == "vapour":
        text += ". From there the calculation follows the vapour and leaves the brine behind"
    return text + "."


def _sorption_clause(result: dict, run: dict) -> str:
    """The exchanger's make-up at the end (equivalent fractions) and the surface potential of the ferric oxide."""
    sc = result["scenario"]
    sorbed = run.get("sorbed") or {}
    bits = []
    cec = (sc.get("exchange_meq", 0.0) + Scenario.from_dict(sc).clay_meq) / 1000.0
    if cec > 0:
        last = max((i for i in range(len(run["properties"]["pH"])) if run["properties"]["pH"][i] is not None),
                   default=None)
        if last is not None:
            shares = []
            for k, line in sorbed.items():
                m = re.fullmatch(r"(.+?)X(\d*)", k)
                if m and not k.startswith("Hfo_") and line[last]:
                    z = int(m.group(2) or 1)
                    shares.append((line[last] * z / cec, m.group(1)))
            shares = [(f, ion) for f, ion in sorted(shares, reverse=True) if f >= 0.005]
            if shares:
                bits.append("the exchange sites end up holding " + ", ".join(
                    f"{100 * f:.0f}% {ion}" for f, ion in shares) + " by charge")
    part = run.get("partition") or {}
    metals = ("Pb", "Cu", "Zn", "Cd", "Ni", "Hg", "As", "Cr")
    kd = {}
    for k, line in part.items():
        el, where = k.split("|", 1)
        vals = [v for v in line if v is not None]
        if el in metals and vals:
            kd.setdefault(el, []).append((vals[-1], where))
    if kd:
        bits.append("partition coefficients at the end (L/kg): " + "; ".join(
            f"{el} " + ", ".join(f"{v:.3g} on {'bed sediment' if w == 'sediment' else w}" for v, w in sorted(vals, reverse=True))
            for el, vals in sorted(kd.items(), key=lambda kv: -max(v for v, _ in kv[1]))))
    pot = [v for v in run["properties"].get("surface_potential_mV", []) if v is not None]
    if pot:
        bits.append(f"the ferric oxide surface potential goes from {pot[0]:+.0f} to {pot[-1]:+.0f} mV")
    if not bits:
        return ""
    text = "; ".join(bits)
    return text[0].upper() + text[1:] + "."


def _isotope_clause(run: dict) -> str:
    """How the fluid's δ18O, δ34S and 87Sr/86Sr change from the first solved point to the last."""
    iso = run.get("isotopes") or {}

    def ends(key):
        vals = [v for v in iso.get(key, []) if v is not None]
        return (vals[0], vals[-1]) if vals else None

    bits = []
    for key, what, fmt in (("d18O:water", "δ18O of the water", "{:+.1f}‰"), ("dD:water", "δD of the water", "{:+.1f}‰"),
                           ("d34S:sulfate", "δ34S of the dissolved sulfate", "{:+.1f}‰"),
                           ("d34S:sulfide", "δ34S of the dissolved sulfide", "{:+.1f}‰"),
                           ("87Sr/86Sr:fluid", "87Sr/86Sr of the fluid", "{:.5f}")):
        e = ends(key)
        if e is None:
            continue
        a, b = e
        if fmt.endswith("‰") and abs(b - a) < 0.05 or not fmt.endswith("‰") and abs(b - a) < 5e-6:
            continue
        bits.append(f"{what} goes from {fmt.format(a + 0.0)} to {fmt.format(b + 0.0)}".replace("-0.0‰", "0.0‰")
                    .replace("+0.0‰", "0.0‰"))
    water = iso.get("dD:water") or []
    gaps = []                                   # how much lighter than the water the hydrous minerals end up
    for key, vals in iso.items():
        if key.startswith("dD:") and key != "dD:water" and not key.endswith(("(g)", ":H2", ":CH4")):
            if vals and vals[-1] is not None and len(water) == len(vals) and water[-1] is not None:
                gaps.append(water[-1] - vals[-1])       # minerals still there at the end
    if gaps:
        lo, hi = round(min(gaps)), round(max(gaps))
        size = f"{lo}‰" if lo == hi else f"{lo} to {hi}‰"
        bits.append(f"the hydrous minerals end up {size} {'lighter' if lo >= 0 else 'different'} in δD than the water")
    gas = []                                    # the gas where it was last present
    for key, what in (("d13C:CO2(g)", "CO₂ δ13C"), ("d18O:CO2(g)", "CO₂ δ18O"), ("d18O:H2O(g)", "steam δ18O"),
                      ("dD:H2O(g)", "steam δD"), ("d13C:CH4(g)", "methane δ13C"), ("dD:CH4(g)", "methane δD"),
                      ("dD:H2(g)", "H₂ δD"), ("d34S:H2S(g)", "H₂S δ34S")):
        e = ends(key)
        if e is not None:
            gas.append(f"{what} {e[1] + 0.0:+.1f}‰")
    if gas:
        bits.append("the gas, where last present, has " + ", ".join(gas))
    if not bits:
        return ""
    return "Isotopes: " + "; ".join(bits) + "."


def _span(axis: list, cells: list[int], unit: str, digits: int) -> str:
    """The stretches of the column that ``cells`` cover, e.g. "0.05 to 0.25 and 0.65 m"."""
    runs: list[list[int]] = []
    for i in cells:
        if runs and i == runs[-1][-1] + 1:
            runs[-1].append(i)
        else:
            runs.append([i])
    parts = [f"{axis[r[0]]:.{digits}f}" if len(r) == 1 else f"{axis[r[0]]:.{digits}f} to {axis[r[-1]]:.{digits}f}"
             for r in runs]
    return (", ".join(parts[:-1]) + " and " + parts[-1] if len(parts) > 1 else parts[0]) + f" {unit}"


def _pv(v: float) -> str:
    v = round(v, 2)
    return f"{v:g} pore volume" + ("" if v == 1 else "s")


def _per_kg(mol: float) -> str:
    return "almost none" if mol < 1e-6 else f"{_amount(mol)}/kg"


def _column_kinetic_clauses(result: dict, end: dict) -> list[str]:
    """How much of each kinetic reactant has reacted at the inlet and the outlet by the end, and how near
    saturation the pore water is there."""
    out = []
    names = result.get("catalog", {}).get("kinetic", {})
    for key, line in end.get("kinetic", {}).items():
        if not key.startswith("left:") or _first_valid(line) is None:
            continue
        name = key.split(":", 1)[1]
        label = (names.get(key) or {}).get("label", name).removesuffix(" left")
        label = label[0].lower() + label[1:] if label[1:2].islower() else label
        a, b = 100.0 - line[_first_valid(line)], 100.0 - line[_last_valid(line)]

        def how(x):
            return ("lost", f"{x:.2g}%") if x > 0.005 else (("grown", f"{-x:.2g}%") if x < -0.005 else ("barely changed", ""))
        (va, xa), (vb, xb) = how(a), how(b)
        second = xb if vb == va and xb else f"{vb} {xb}".strip()
        text = (f"At the measured rates, by the end the {label} has barely changed anywhere in the column"
                if not xa and not xb else
                f"At the measured rates, by the end the {label} has {va} {xa} at the inlet".replace("  ", " ")
                + f" and {second} at the outlet")
        si = end.get("saturation", {}).get(name)
        if si and si[_first_valid(si)] is not None and si[_last_valid(si)] is not None:
            text += (f"; its saturation index in the pore water goes from {si[_first_valid(si)]:.2f} at the inlet "
                     f"to {si[_last_valid(si)]:.2f} at the outlet")
        out.append(text + ".")
    return out


def _redox_sentence(props: dict, start: str, end: str) -> str:
    """How the redox state moves: Eh and log fO2 from the first to the last point (empty without redox)."""
    eh, fo2 = props.get("Eh"), props.get("log_fO2")
    if not eh or _first_valid(eh) is None:
        return ""
    a, b = _first_valid(eh), _last_valid(eh)
    text = f"Eh goes from {eh[a]:+.2f} V {start} to {eh[b]:+.2f} V {end}"
    if fo2 and _first_valid(fo2) is not None:
        fa, fb = fo2[_first_valid(fo2)], fo2[_last_valid(fo2)]
        text += f" (log fO2 {fa:.1f} to {fb:.1f})"
    word = "oxidizing" if eh[b] > 0.2 else "reducing" if eh[b] < 0.0 else "mildly oxidizing"
    return f" {text}, so the fluid ends {word}."


def _column_story(result: dict) -> list[str]:
    """A column run: what comes out at the outlet, and where the minerals sit along the column."""
    modes = result["modes"]
    out, start, end = modes["outlet"], modes["start"], modes["profile"]
    pv = out["axes"]["pore_volumes"]["values"]
    ph = out["properties"]["pH"]
    i0, i1 = _first_valid(ph), _last_valid(ph)
    paras = [f"{result['summary'][0].upper()}{result['summary'][1:]}. At the outlet the pH goes from {ph[i0]:.1f} "
             f"to {ph[i1]:.1f} after {_pv(pv[i1])}." + _redox_sentence(out["properties"], "at first", f"after {_pv(pv[i1])}")]

    moves = []
    for el, line in out["elements"].items():
        vals = [v for v in line if v is not None]
        if not vals:
            continue
        a, b, top = line[_first_valid(line)], line[_last_valid(line)], max(vals)
        if top < 1e-6:
            continue
        name = ELEMENT_WORDS.get(el, el)
        if top > 1.1 * max(a, b) and top - max(a, b) > 1e-6:
            k = line.index(top)
            moves.append((top - min(a, b), f"{name} peaks at {_per_kg(top)} after {_pv(pv[k])}, "
                                           f"above both its starting and its final level"))
        elif abs(b - a) > 0.1 * max(a, b):
            half = (a + b) / 2.0
            k = _first(line, (lambda v: v >= half) if b > a else (lambda v: v <= half))
            verb = "rises" if b > a else "falls"
            moves.append((abs(b - a), f"{name} {verb} from {_per_kg(a)} to {_per_kg(b)}, halfway "
                                      f"after {_pv(pv[k])}"))
    if moves:
        moves.sort(key=lambda m: -m[0])
        text = "; ".join(m[1] for m in moves[:6])
        paras.append("At the outlet " + text + ".")

    dist = end["axes"]["distance"]["values"]
    temps = end["axes"].get("temperature_c", {}).get("values")
    zones, gone = [], []
    names = sorted(set(start["minerals"]) | set(end["minerals"]))
    peak = {m: max((v or 0.0 for run in (start, end) for v in run["minerals"].get(m, [])), default=0.0) for m in names}
    major = 0.01 * max(peak.values(), default=0.0)          # minerals below 1% of the largest are left out
    for m in names:
        if peak[m] < major:
            continue
        s_line = start["minerals"].get(m, [0.0] * len(dist))
        e_line = end["minerals"].get(m, [0.0] * len(dist))
        here = [i for i, v in enumerate(e_line) if v is not None and v > APPEARS]
        before = [i for i, v in enumerate(s_line) if v is not None and v > APPEARS]
        label = pretty_mineral(m)
        label = label[0].lower() + label[1:] if label[1:2].islower() else label
        if here:
            where = _span(dist, here, "m", 2)
            if temps:
                where += f", {_span(temps, here, '°C', 0)}"
            zones.append((here[0], f"{label} ({where})"))
        lost = [i for i in before if i not in here]
        if lost:
            gone.append(f"{label} ({_span(dist, lost, 'm', 2)})")
    if zones:
        zones.sort()
        paras.append("Along the column at the end: " + "; ".join(z[1] for z in zones[:12]) + "."
                     + (" Gone since the start: " + "; ".join(gone) + "." if gone else ""))
    elif gone:
        paras.append("Gone from the column since the start: " + "; ".join(gone) + ".")
    paras.extend(_column_kinetic_clauses(result, end))
    sorb = _sorption_clause(result, out)
    if sorb:
        paras.append("At the outlet, at the end: " + sorb[0].lower() + sorb[1:])
    levels = end.get("reliability", {}).get("levels", [])
    if any(level != "good" for level in levels):
        paras.append("Part of the column lies where the thermodynamics of charged species is less certain: see "
                     "Assumptions and limits, and the coloured strip under the plot.")
    return paras


MODE_WORDS = {"static": "static", "flow_through": "flow-through"}


def story(result: dict) -> list[str]:
    """What happened during the run, in plain sentences, one paragraph per idea."""
    if result["scenario"].get("kind") == "transport":
        return _column_story(result)
    sc = result["scenario"]
    blocked = {b.lower() for b in sc["blocked_minerals"]}
    modes = result["modes"]
    paras: list[str] = []

    first_key = next(iter(modes))
    ph = modes[first_key]["properties"]["pH"]
    i0, i1 = _first_valid(ph), _last_valid(ph)
    paras.append(
        f"{_sentence(plain(result['summary']))}. The in situ pH goes from "
        f"{ph[i0]:.1f} at {point_text(result, i0)} to {ph[i1]:.1f} at {point_text(result, i1)}."
        + _redox_sentence(modes[first_key]["properties"], "at the start", "at the end"))

    for key, run in modes.items():
        run_name = f"In the {MODE_WORDS[key]} run"
        sentences = []
        kinetic = _kinetic_clause(result, run)
        if kinetic:
            sentences.append(kinetic)
        names, about = _mineral_sentences(result, run, blocked)
        if names:
            sentences.append(f"{run_name}, {_count(len(names), 'mineral')} form{'s' if len(names) == 1 else ''}"
                             f"{', in this order' if len(names) > 1 else ''}: {_listing(names)}.")
            sentences += about
        else:
            sentences.append(f"{run_name}, no minerals form.")
        lost = []
        for el in ("S", "Mg", "Ca", "C"):
            kept = _retained(result, run, el)
            if kept is not None and kept < 0.95:
                lost.append(f"{(1 - kept) * 100:.0f}% of its {ELEMENT_WORDS[el]}")
        if lost:
            sentences.append("By the end the fluid has lost " + _listing(lost) + ".")
        for extra in (_gas_clause(result, run, key), _split_clause(result, run), _sorption_clause(result, run),
                      _isotope_clause(run)):
            if extra:
                sentences.append(extra)
        paras.append(" ".join(s.strip() for s in sentences))
        energy = _energy_clause(result, run)
        if energy:
            paras.append(energy)

    if len(modes) == 2:
        a, b = modes["static"], modes["flow_through"]
        diffs = []
        for el in ("S", "Mg", "Ca", "C"):
            ka, kb = _retained(result, a, el), _retained(result, b, el)
            if ka is not None and kb is not None and abs(ka - kb) > 0.05:
                diffs.append(f"{ka * 100:.0f}% of the {ELEMENT_WORDS[el]} is left when static and "
                             f"{kb * 100:.0f}% when flow-through")
        pa, pb = a["properties"]["pH"], b["properties"]["pH"]
        ea, eb = _last_valid(pa), _last_valid(pb)
        ph_note = ""
        if ea is not None and eb is not None and abs(pa[ea] - pb[eb]) > 0.1:
            ph_note = f" The final pH is {pa[ea]:.1f} static and {pb[eb]:.1f} flow-through."
        if diffs:
            paras.append("Where the two modes differ: " + "; ".join(diffs) + "." + ph_note)
        else:
            paras.append("Static and flow-through leave almost the same fluid at the end." + ph_note)

    held = []
    for name in sc["blocked_minerals"]:
        for run in modes.values():
            sat = run["saturation"].get(name)
            if sat and max((v for v in sat if v is not None), default=-9) > 0.1:
                held.append(mineral_word(name))
                break
    if held:
        paras.append("The fluid becomes supersaturated in " + _listing(held) + ", but "
                     + ("it was" if len(held) == 1 else "they were") + " held out as kinetically blocked.")

    levels = result.get("reliability", {}).get("levels", [])
    if any(level != "good" for level in levels):
        worst = "poor" if "poor" in levels else "caution"
        where = "some" if worst == "caution" else "part"
        paras.append(f"{where.capitalize()} of this path lies where the thermodynamics of charged species is "
                     f"less certain, so the numbers there are {'unreliable' if worst == 'poor' else 'approximate'}. "
                     f"See Assumptions and limits, and the coloured strip under the plot.")
    return paras


FLUID_PROPERTIES = (("temperature", "Temperature", "°C"), ("pressure", "Pressure", "bar"), ("pH", "pH", ""),
                    ("Eh", "Eh", "V"), ("log_fO2", "log fO2", ""), ("ionic_strength", "Ionic strength", "mol/kg"),
                    ("water_kg", "Water left", "kg per kg at the start"))


def _num(v: float, digits: int = 3) -> str:
    if v == 0:
        return "0"
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    return f"{v:.{digits}g}"


def end_lists(result: dict) -> list[dict]:
    """The solids and the fluid at the end of each mode, as rows of plain text for two tables.

    Solids are per kg of the original water; dissolved amounts are per kg of the water then present."""
    if result["scenario"].get("kind") == "transport":
        return []
    cat = result.get("catalog", {})
    temps, pressures, _ = _axes(result)
    out = []
    for key, run in result["modes"].items():
        ph = run["properties"]["pH"]
        i = _last_valid(ph)
        if i is None:
            continue
        solids = []
        for name, series in run["minerals"].items():
            v = series[i] if i < len(series) else None
            if v is None or v <= APPEARS:
                continue
            mw = (cat.get("minerals", {}).get(name) or {}).get("mw")
            solids.append((v, [pretty_mineral(name), _num(v * 1000), _num(v * mw) if mw else ""]))
        solids = [row for _, row in sorted(solids, key=lambda r: -r[0])]
        props = []
        for k, label, unit in FLUID_PROPERTIES:
            if k == "temperature":
                v = temps[i]
            elif k == "pressure":
                v = pressures[i]
            else:
                v = (run["properties"].get(k) or [None] * (i + 1))[i]
            if v is None:
                continue
            text = f"{v:.2f}" if k in ("pH", "Eh") else f"{v:.1f}" if k == "log_fO2" else _num(v)
            props.append([label, f"{text} {unit}".strip()])
        dissolved = []
        for el, series in run["elements"].items():
            v = series[i] if i < len(series) else None
            if not v or v < 1e-12 or el in ("H", "O"):
                continue
            mg = v * ATOMIC_WEIGHTS[el] * 1000 if el in ATOMIC_WEIGHTS else None
            label = ELEMENT_WORDS.get(el, el)
            dissolved.append((v, [f"{label[:1].upper()}{label[1:]}", _num(v * 1000), _num(mg) if mg else ""]))
        dissolved = [row for _, row in sorted(dissolved, key=lambda r: -r[0])]
        out.append({"mode": key, "title": f"End of the {MODE_WORDS[key]} run, at {point_text(result, i)}",
                    "solids": solids, "fluid": props, "dissolved": dissolved})
    return out

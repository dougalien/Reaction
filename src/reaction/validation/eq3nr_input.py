"""Reader for EQ3NR menu-style input files (``.3i``).

This is a validation aid: it lets the EQ3/6 test library drive Reaction so
results can be compared case by case with the stored ``.3o`` outputs.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class BasisRow:
    name: str
    value: float
    unit: str                      # EQ3NR constraint string, e.g. "Molality", "pH"
    phase: str | None = None       # partner phase for "Hetero. equil."


@dataclass
class Eq3nrInput:
    path: str
    description: str
    temperature_c: float
    pressure: tuple[str, float]    # ("reference"|"saturation"|"bars", value)
    density: float
    tds: tuple[str, float]         # ("mg/kg.sol"|"mg/L", value)
    balance_on: str | None
    redox: tuple[str, float, str | None]   # (kind, value, couple)
    switches: list[tuple[str, str]]
    rows: list[BasisRow]
    alter_suppress: list[tuple[str, str, float]]
    activity_model: str            # davies | bdot | pitzer | hc
    ph_scale: str                  # internal | nbs | mesmer
    exchangers: int = 0
    solid_solutions: int = 0
    notes: list[str] = field(default_factory=list)


def _checked(lines: list[str], start: int, stop: int) -> list[tuple[int, str]]:
    """(code, label) of [x] option lines in a range."""
    out = []
    for ln in lines[start:stop]:
        m = re.match(r"\|\s*\[x\]\s*\(\s*(-?\d+)\)\s*(.*?)\s*(\||$)", ln)
        if m:
            out.append((int(m.group(1)), m.group(2)))
    return out


def _find(lines: list[str], pattern: str, start: int = 0) -> int:
    for i in range(start, len(lines)):
        if re.search(pattern, lines[i]):
            return i
    raise ValueError(f"pattern {pattern!r} not found")


def _num(text: str) -> float:
    return float(text.replace("D", "E"))


def read_eq3nr_inputs(path: str | Path) -> list[Eq3nrInput]:
    """All problems in an input file (problems end at ``|End of problem|``)."""
    path = Path(path)
    lines = path.read_text(encoding="latin-1").splitlines()
    chunks, current = [], []
    for ln in lines:
        current.append(ln)
        if ln.startswith("|End of problem"):
            chunks.append(current)
            current = []
    if any("Temperature (C)" in ln for ln in current):      # file without a closing line
        chunks.append(current)
    return [_parse_problem(chunk, path) for chunk in chunks]


def read_eq3nr_input(path: str | Path) -> Eq3nrInput:
    """The first problem of an input file."""
    return read_eq3nr_inputs(path)[0]


def _parse_problem(lines: list[str], path: Path) -> Eq3nrInput:
    m = re.search(r'Description=\s*"([^"]*)"', "\n".join(lines))
    description = m.group(1) if m else ""

    # special basis switches
    switches: list[tuple[str, str]] = []
    i0 = _find(lines, r"Special Basis Switches")
    i1 = _find(lines, r"\|Temperature \(C\)")
    olds = [re.match(r"\|Replace\s*\|(.*?)\s*\|", ln) for ln in lines[i0:i1]]
    news = [re.match(r"\|\s*with\s*\|(.*?)\s*\|", ln) for ln in lines[i0:i1]]
    for a, b in zip([x for x in olds if x], [x for x in news if x]):
        if a.group(1) != "None":
            switches.append((a.group(1), b.group(1)))

    temperature = _num(re.match(r"\|Temperature \(C\)\s*\|\s*(\S+)\|", lines[i1]).group(1))

    ip = _find(lines, r"Pressure option")
    pchoice = _checked(lines, ip, ip + 5)[0][0]
    if pchoice == 2:
        pm = re.search(r"Value \(bars\)\s*\|\s*(\S+)\|", "\n".join(lines[ip:ip + 5]))
        pressure = ("bars", _num(pm.group(1)))
    else:
        pressure = (("reference", 0.0) if pchoice == 0 else ("saturation", 0.0))

    density = _num(re.search(r"\|Density \(g/cm3\)\s*\|\s*(\S+)\|", "\n".join(lines)).group(1))

    it = _find(lines, r"Total dissolved solutes option")
    tchoice = _checked(lines, it, it + 4)[0][0]
    tm = re.findall(r"Value \(mg/(?:kg\.sol|L)\)\s*\|\s*(\S+)\|", "\n".join(lines[it:it + 4]))
    tds = (("mg/kg.sol", _num(tm[0])) if tchoice == 0 else ("mg/L", _num(tm[1])))

    ie = _find(lines, r"Electrical balancing option")
    echoice = _checked(lines, ie, ie + 4)[0][0]
    balance_on = None
    if echoice == 1:
        balance_on = re.search(r"Balance on species\s*\|(.*?)\s*\|",
                               "\n".join(lines[ie:ie + 4])).group(1)

    ir = _find(lines, r"Default redox constraint")
    rchoice = _checked(lines, ir, ir + 7)[0][0]
    block = "\n".join(lines[ir:ir + 7])
    kinds = {-3: "o2line", -2: "pe", -1: "eh", 0: "logfo2", 1: "couple"}
    kind = kinds[rchoice]
    value, couple = 0.0, None
    if kind == "pe":
        value = _num(re.search(r"pe units\)\s*\|\s*(\S+)\|", block).group(1))
    elif kind == "eh":
        value = _num(re.search(r"Eh \(volts\)\s*\|\s*(\S+)\|", block).group(1))
    elif kind == "logfo2":
        value = _num(re.search(r"Log fO2 \(log bars\)\s*\|\s*(\S+)\|", block).group(1))
    elif kind == "couple":
        couple = re.search(r"Couple \(aux\. sp\.\)\s*\|(.*?)\s*\|", block).group(1)

    # basis / constraint rows
    ib = _find(lines, r"Aqueous Basis Species/Constraint Species")
    ie_ = _find(lines, r"^\* Valid jflag strings", ib)
    rows: list[BasisRow] = []
    for ln in lines[ib + 3:ie_]:
        m = re.match(r"\|->\|\s*(\S+)", ln)
        if m and rows:
            rows[-1].phase = m.group(1)
            continue
        m = re.match(r"\|(\S+)\s*\|\s*(\S+)\|(.*?)\s*\|", ln)
        if m:
            rows.append(BasisRow(m.group(1), _num(m.group(2)), m.group(3).strip()))

    # exchangers, solid solutions
    text = "\n".join(lines)
    exchangers = len(re.findall(r"\|Exchanger phase \|(?!None)", text))
    solid_solutions = len(re.findall(r"\|Solid Solution\s+\|(?!None)", text))

    # alter/suppress
    alter: list[tuple[str, str, float]] = []
    ia = _find(lines, r"Alter/Suppress Options")
    ia_end = _find(lines, r"^\* Valid alter/suppress", ia)
    for ln in lines[ia + 4:ia_end]:
        m = re.match(r"\|(\S.*?)\s*\|(\S.*?)\s*\|\s*(\S+)\|", ln)
        if m and m.group(1) != "None":
            alter.append((m.group(1), m.group(2), _num(m.group(3))))

    def option(header: str, span: int):
        i = _find(lines, header)
        return _checked(lines, i, i + span)[0][0]

    act = {-1: "davies", 0: "bdot", 1: "pitzer", 2: "hc"}[option(r"iopg\(1\)", 6)]
    scale = {-1: "internal", 0: "nbs", 1: "mesmer"}[option(r"iopg\(2\)", 5)]

    return Eq3nrInput(
        path=str(path), description=description, temperature_c=temperature,
        pressure=pressure, density=density, tds=tds, balance_on=balance_on,
        redox=(kind, value, couple), switches=switches, rows=rows,
        alter_suppress=alter, activity_model=act, ph_scale=scale,
        exchangers=exchangers, solid_solutions=solid_solutions,
    )

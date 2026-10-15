"""Importer for EQ3/6 ``data0`` thermodynamic database files.

The file is read from a path the user supplies; nothing from it is copied
into this repository (see docs/DESIGN.md, "Data policy").

Layout of a data0 file (EQ3/6 v8):
  header      title, temperature/pressure grid, Debye-Hueckel A, B, bdot,
              Drummond CO2 coefficients, logK of the Eh reaction
  bdot table  per-species ion-size parameter and neutral-species type
  elements    atomic weights
  species     blocks for basis, auxiliary basis, aqueous, solid, liquid, gas
Each species block carries a dissociation reaction and a logK grid on the
header's 8-point temperature grid (500 marks "no data").
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..constants import CAL
from .model import GridLogK, Kind, Provenance, Species, ThermoDatabase

_NUM = re.compile(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?")
_NO_DATA = 500.0

_SECTION_KIND = {
    "basis species": Kind.BASIS,
    "auxiliary basis species": Kind.AUX_BASIS,
    "aqueous species": Kind.AQUEOUS,
    "solids": Kind.MINERAL,
    "liquids": Kind.LIQUID,
    "gases": Kind.GAS,
}
_OTHER_SECTIONS = {"miscellaneous parameters", "bdot parameters", "elements",
                   "solid solutions", "references"}


def _numbers(text: str) -> list[float]:
    return [float(x) for x in _NUM.findall(text)]


def _blocks(lines: list[str]) -> list[list[str]]:
    blocks: list[list[str]] = [[]]
    for ln in lines:
        if ln.startswith("+---"):
            blocks.append([])
        else:
            blocks[-1].append(ln)
    return blocks


def _pairs(tokens: list[str]) -> list[tuple[float, str]]:
    return [(float(tokens[i]), tokens[i + 1]) for i in range(0, len(tokens) - 1, 2)]


def _to_si(value: float, unit: str) -> float:
    u = unit.lower()
    if u.startswith("kcal"):
        return value * CAL * 1000.0
    if u.startswith("kj"):
        return value * 1000.0
    if u.startswith("cal"):
        return value * CAL
    return value  # j/mol, j/(mol*K)


def _parse_header(block: list[str]) -> dict[str, np.ndarray]:
    keys = {
        "temperatures": "temps", "pressures": "press",
        "debye huckel a": "adh", "debye huckel b": "bdh", "bdot": "bdot",
        "cco2": "cco2", "log k for eh reaction": "eh",
    }
    out: dict[str, list[float]] = {}
    current = None
    for ln in block:
        low = ln.strip().lower()
        hit = next((v for k, v in keys.items() if low.startswith(k)), None)
        if hit is not None:
            current = hit
            out[current] = []
        elif current is not None and low and _NUM.fullmatch(low.split()[0]):
            out[current].extend(_numbers(ln))
    return {k: np.array(v) for k, v in out.items()}


def _parse_entry(block: list[str], kind: Kind, section: str, db_title: str,
                 path: str, ntypes: dict[str, int]
                 ) -> tuple[Species, np.ndarray | None] | None:
    """Parse one species block; returns the species and its raw logK grid."""
    block = [ln for ln in block if ln.strip()]
    if not block:
        return None
    head = block[0].split()
    name = head[0]
    formula = head[1] if len(head) > 1 else None

    charge = 0.0
    mol_wt = a0 = vol = None
    revised = None
    subtype = ""
    data_sets: list[str] = []
    elements: dict[str, float] = {}
    reaction: dict[str, float] = {}
    grid: list[float] = []
    ref: dict[str, float] = {}
    sources: list[str] = []

    i = 1
    while i < len(block):
        ln = block[i]
        s = ln.strip().lstrip("*").strip()
        low = s.lower()
        if low.startswith("sp.type"):
            subtype = s.split("=", 1)[1].strip()
        elif low.startswith("eq3/6"):
            data_sets = [x.strip() for x in s.split("=", 1)[1].split(",")]
        elif low.startswith("revised"):
            revised = s.split("=", 1)[1].strip()
        elif low.startswith("mol.wt"):
            mol_wt = _numbers(s.split("=", 1)[1])[0]
        elif low.startswith("dhazero"):
            a0 = _numbers(s.split("=", 1)[1])[0]
        elif low.startswith("charge"):
            charge = _numbers(s.split("=", 1)[1])[0]
        elif low.startswith("v0prtr"):
            vol = _numbers(s.split("=", 1)[1])[0]
        elif re.match(r"\d+ element\(s\)", s):
            i += 1
            while i < len(block) and not block[i].startswith("****"):
                for n, el in _pairs(block[i].split()):
                    elements[el] = n
                i += 1
            continue
        elif "species in aqueous dissociation reaction" in s:
            i += 1
            while i < len(block) and block[i].strip() != "*":
                for n, sp in _pairs(block[i].split()):
                    reaction[sp] = n
                i += 1
            continue
        elif ln.startswith("****") and "logk grid" in low:
            i += 1
            while i < len(block) and block[i].strip() != "*":
                grid.extend(_numbers(block[i]))
                i += 1
            continue
        else:
            m = re.match(r"(delg0f|delh0f|s0prtr)\s*=\s*([-+\d.eE]+)\s*(\S+)\s*\[(\w+)\]", low)
            if m:
                key, val, unit, tag = m.groups()
                # prefer the reported value over a calculated one
                tgt = {"delg0f": "dG_f", "delh0f": "dH_f", "s0prtr": "S"}[key]
                if tag == "reported" or tgt not in ref:
                    ref[tgt] = _to_si(float(val), unit)
            m = re.match(r"t\*\*(-?\d+)\s*=\s*([-+\d.eE]+)", low)
            if m:
                ref[f"cp_T{m.group(1)}"] = float(m.group(2)) * CAL
            m = re.match(r"tlimit\s*=\s*([-+\d.eE]+)", low)
            if m:
                ref["cp_Tlimit_C"] = float(m.group(1))
            for src in re.findall(r"\[source:\s*([^\]]*?)\s*\]", s):
                sources.append(src.strip())
        i += 1

    values = None
    if grid:
        values = np.array(grid, dtype=float)
        values[values == _NO_DATA] = np.nan
    sp = Species(
        name=name, kind=kind, formula=formula, charge=charge, elements=elements,
        mol_wt=mol_wt, a0=a0, activity_type=ntypes.get(name, 0),
        reaction=reaction, logk=None, molar_volume=vol, ref_state=ref,
        provenance=Provenance(database=db_title, file=Path(path).name, section=section,
                              revised=revised, sources=sorted(set(sources)),
                              data_sets=data_sets),
        subtype=subtype,
    )
    return sp, values


def load_data0(path: str | Path) -> ThermoDatabase:
    """Parse an EQ3/6 data0 file into a :class:`ThermoDatabase`."""
    path = Path(path)
    lines = path.read_text(encoding="latin-1").splitlines()
    blocks = _blocks(lines)
    title = blocks[0][0].strip() if blocks[0] else path.name

    header: dict[str, np.ndarray] = {}
    ntypes: dict[str, int] = {}
    elements: dict[str, float] = {}
    entries: list[tuple[list[str], Kind, str]] = []

    section = ""
    for blk in blocks[1:]:
        nonblank = [ln for ln in blk if ln.strip()]
        if len(nonblank) == 1 and nonblank[0].strip().lower() in (
                set(_SECTION_KIND) | _OTHER_SECTIONS):
            section = nonblank[0].strip().lower()
            continue
        if section == "miscellaneous parameters":
            header.update(_parse_header(blk))
        elif section == "bdot parameters":
            for ln in blk:
                tok = ln.split()
                if len(tok) == 3 and not ln.startswith("*"):
                    ntypes[tok[0]] = int(float(tok[2]))
        elif section == "elements":
            for ln in blk:
                tok = ln.split()
                if len(tok) == 2:
                    elements[tok[0]] = float(tok[1])
        elif section in _SECTION_KIND:
            entries.append((blk, _SECTION_KIND[section], section))

    temps = header["temps"]
    db = ThermoDatabase(
        title=title, source=str(path), temps_c=temps, pressures_bar=header["press"],
        adh=header["adh"], bdh=header["bdh"], bdot=header["bdot"], cco2=header["cco2"],
        logk_eh=GridLogK(temps, header["eh"]), elements=elements,
    )
    for blk, kind, sec in entries:
        parsed = _parse_entry(blk, kind, sec, title, str(path), ntypes)
        if parsed is None:
            continue
        sp, values = parsed
        if values is not None and len(values) == len(temps):
            sp.logk = GridLogK(temps, values)
        if sp.name in db.species:
            if set(sp.reaction) == {sp.name}:
                # reference-state identity entry (e.g. gas O2(g) = basis O2(g))
                db.skipped.append(f"{sp.name} ({sec}): identity entry, basis species used")
                continue
            raise ValueError(f"duplicate species name in {path.name}: {sp.name}")
        db.species[sp.name] = sp
    return db

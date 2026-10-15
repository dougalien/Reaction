"""Readable names for species and minerals, and tolerant matching of what the user types."""
from __future__ import annotations

import re
import unicodedata

_SUB = str.maketrans("0123456789", "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089")
_SUP = str.maketrans("0123456789+-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207a\u207b")
_FROM_SUP = {c: a for a, c in zip("0123456789+-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207a\u207b")}
_FROM_SUB = {c: a for a, c in zip("0123456789", "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089")}

ELEMENT_NAMES = {
    "hydrogen": "H", "carbon": "C", "oxygen": "O", "sodium": "Na", "magnesium": "Mg",
    "aluminum": "Al", "aluminium": "Al", "silicon": "Si", "silica": "Si", "sulfur": "S",
    "sulphur": "S", "sulfate": "S", "chlorine": "Cl", "chloride": "Cl", "potassium": "K",
    "calcium": "Ca", "iron": "Fe", "manganese": "Mn", "strontium": "Sr", "barium": "Ba",
    "lithium": "Li", "boron": "B", "bromine": "Br", "bromide": "Br", "fluorine": "F",
    "fluoride": "F", "zinc": "Zn", "copper": "Cu", "lead": "Pb", "nitrogen": "N",
}


def _split_charge(name: str):
    m = re.search(r"([+-])(\d?)$", name)
    if not m:
        return name, 0
    # "SO4-2": the digit after the sign is the charge; a bare sign (HCO3-) is one unit
    n = int(m.group(2)) if m.group(2) else 1
    return name[:m.start()], n if m.group(1) == "+" else -n


def pretty_species(name: str, charge: float | None = None, gas: bool = False) -> str:
    """``Ca+2`` -> ``Ca\u00b2\u207a``, ``HCO3-`` -> ``HCO\u2083\u207b``, neutral ``CaSO4`` -> ``CaSO\u2084(aq)``."""
    body, z = _split_charge(name)
    body_p = re.sub(r"(?<=[A-Za-z)])(\d+)", lambda m: m.group(1).translate(_SUB), body)
    if z:
        mag = str(abs(z)) if abs(z) > 1 else ""
        return body_p + (mag + ("+" if z > 0 else "-")).translate(_SUP)
    if gas or body.endswith("(g)"):
        return body_p
    if charge == 0 and not body.endswith("(aq)"):
        return body_p + "(aq)"
    return body_p


def pretty_mineral(name: str) -> str:
    if "," in name:
        base, rest = name.split(",", 1)
        return f"{base[:1].upper()}{base[1:]} ({rest})"
    return name[:1].upper() + name[1:]


def canonical_key(text: str) -> list[str]:
    """Spellings of ``text`` to try when matching against a known name (most literal first)."""
    s = unicodedata.normalize("NFKC", text.strip())
    out = []
    for ch in s:
        out.append(_FROM_SUP.get(ch, _FROM_SUB.get(ch, ch)))
    s = "".join(out)
    # a run of superscript digits/signs ("2+") written after the body is a charge, "^" is decoration
    s = s.replace("^", "").replace(" ", "").lower()
    s = re.sub(r"\(aq\)$", "", s)
    keys = [s]
    m = re.fullmatch(r"(.*?)(\++|-+)", s)
    if m and len(m.group(2)) > 1:                         # Ca++ -> ca+2
        keys.append(f"{m.group(1)}{m.group(2)[0]}{len(m.group(2))}")
    m = re.fullmatch(r"(.*?)(\d)([+-])", s)
    if m:                                                  # Ca2+ -> ca+2
        keys.append(f"{m.group(1)}{m.group(3)}{m.group(2)}")
    return keys


def build_index(names: dict[str, str]) -> dict[str, str]:
    """Map every accepted spelling to an internal name.  ``names``: internal -> pretty."""
    index: dict[str, str] = {}
    for internal, pretty in names.items():
        for spelling in (internal, pretty):
            for key in canonical_key(spelling):
                index.setdefault(key, internal)
        # spelled the EQ3/6 way: Ca+2 -> Ca++
        body, z = _split_charge(internal)
        if z and abs(z) > 1:
            index.setdefault((body + ("+" if z > 0 else "-") * abs(z)).lower(), internal)
    return index


def resolve(text: str, index: dict[str, str]) -> str | None:
    for key in canonical_key(text):
        if key in index:
            return index[key]
    return None

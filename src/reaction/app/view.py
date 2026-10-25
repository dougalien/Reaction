"""What the results window is showing, and how a request such as "plot CaCl+ instead of Ca++" changes it.

The view is shared state: the chat and the clicks in the results window both go through
``apply_patch``, so they behave the same.  A plotted item is ``{"group": ..., "key": ...}`` where
the group is one of species, elements, minerals, gas, saturation, kinetic (rates and what is left of
kinetic reactants) or property.  Items can share a y axis
only if they have the same unit, so adding an item of another unit replaces the current ones.
"""
from __future__ import annotations

import re
from typing import Any

from .names import ELEMENT_NAMES, build_index, resolve

PROPERTIES = {
    "pH": "pH", "ionic_strength": "Ionic strength (mol/kg)", "water_kg": "Solvent mass (kg)",
    "log_aw": "log water activity", "log_fO2": "log fO2", "Eh": "Eh (V)",
}
PROPERTY_WORDS = {"ph": "pH", "ionic strength": "ionic_strength", "ionic_strength": "ionic_strength",
                  "water activity": "log_aw", "activity of water": "log_aw", "water mass": "water_kg",
                  "solvent mass": "water_kg", "fo2": "log_fO2", "log fo2": "log_fO2",
                  "oxygen fugacity": "log_fO2", "eh": "Eh", "redox potential": "Eh"}
MAX_ITEMS = 8
UNITS = ("mol", "mmol", "umol", "mg")


def unit_group(item: dict) -> str:
    g = item["group"]
    if g in ("species", "elements"):
        return "conc"
    if g in ("minerals", "gas"):
        return "mineral"
    if g == "saturation":
        return "si"
    if g == "kinetic":
        return f"kinetic:{item['key'].split(':')[0]}"
    if g == "composition":
        return "composition"
    if g == "phases":
        return "phases:share" if item["key"] == "vapour_share" else "phases:salinity"
    return f"property:{item['key']}"


def default_view(result: dict) -> dict:
    cat = result["catalog"]
    modes = list(result["modes"])
    first = result["modes"][modes[0]]
    wanted = [("species", "Mg+2"), ("species", "Ca+2"), ("species", "SO4-2")]
    items = [{"group": g, "key": k} for g, k in wanted if k in first["species"]]
    if not items:
        biggest = sorted(first["species"], key=lambda k: -max(v or 0.0 for v in first["species"][k]))[:3]
        items = [{"group": "species", "key": k} for k in biggest]
    return {"items": items, "log": True, "modes": modes, "unit": "mmol", "x": result["x"]["key"],
            "compare": None}


def _present(item: dict, result: dict) -> bool:
    run = next(iter(result["modes"].values()))
    g, k = item["group"], item["key"]
    if g == "property":
        return k in PROPERTIES
    return k in run.get(g, {})


def sanitize(view: dict, result: dict | None) -> dict:
    """Drop anything the result cannot plot and enforce one unit group."""
    out = {"items": [], "log": bool(view.get("log", True)), "modes": [], "unit": view.get("unit", "mmol"),
           "x": view.get("x", "temperature_c"), "compare": view.get("compare") or None}
    if out["compare"] is not None and not isinstance(out["compare"], str):
        out["compare"] = None
    if out["unit"] not in UNITS:
        out["unit"] = "mmol"
    if result is None:
        return out
    # the horizontal axis must be a quantity that actually changes along the path
    if out["x"] not in result.get("varies", []):
        out["x"] = result["x"]["key"]
    avail = list(result["modes"])
    out["modes"] = [m for m in view.get("modes", avail) if m in avail] or avail
    group = None
    for item in view.get("items", []):
        if not _present(item, result) or item in out["items"]:
            continue
        ug = unit_group(item)
        if group is None:
            group = ug
        if ug == group and len(out["items"]) < MAX_ITEMS:
            out["items"].append({"group": item["group"], "key": item["key"]})
    return out


def apply_patch(view: dict, patch: dict[str, Any], result: dict | None) -> dict:
    new = {**view, "items": [dict(i) for i in view.get("items", [])]}
    for key in ("log", "modes", "unit", "x"):
        if key in patch and patch[key] is not None:
            new[key] = patch[key]
    if "compare" in patch:
        new["compare"] = patch["compare"] or None
    if "replace" in patch:
        new["items"] = [dict(i) for i in patch["replace"]]
    for item in patch.get("remove", []):
        new["items"] = [i for i in new["items"] if i != item]
    adds = list(patch.get("add", []))
    if "toggle" in patch:
        t = patch["toggle"]
        if t in new["items"]:
            new["items"] = [i for i in new["items"] if i != t]
        else:
            adds.append(t)
    for item in adds:
        if new["items"] and unit_group(new["items"][0]) != unit_group(item):
            new["items"] = []
        if item not in new["items"]:
            new["items"].append(dict(item))
    return sanitize(new, result)


# ---- finding things by what the user typed ---------------------------------------------------------
def resolve_item(text: str, result: dict, prefer: str | None = None) -> dict | None:
    """Match free text to something plottable in ``result``.  ``prefer='saturation'`` makes a
    mineral name mean its saturation index rather than its amount."""
    cat = result["catalog"]
    t = text.strip().replace("_", " ")
    low = t.lower()
    sat = prefer == "saturation" or bool(re.search(r"saturation|\bsi\b|log\s*q/k", low))
    low_clean = re.sub(r"saturation( index)?( of)?|\bsi\b( of)?|log\s*q/k( of)?", "", low).strip()
    if low_clean in PROPERTY_WORDS:
        return {"group": "property", "key": PROPERTY_WORDS[low_clean]}
    run0 = next(iter(result["modes"].values()))
    if low_clean in ("steam", "water vapor", "water vapour", "vapor", "vapour"):
        low_clean = "water gas"
    g = re.fullmatch(r"(.+?)\s*(?:gas|\(g\))", low_clean)
    if g and run0.get("gas"):
        stem = {"methane": "ch4", "carbon dioxide": "co2", "water": "h2o", "steam": "h2o",
                "hydrogen": "h2", "nitrogen": "n2", "hydrogen sulfide": "h2s"}.get(g.group(1), g.group(1))
        for k in run0["gas"]:
            if k.lower() == f"{stem.replace(' ', '')}(g)":
                return {"group": "gas", "key": k}
    species_index = build_index({k: v["label"] for k, v in cat["species"].items()})
    mineral_index = build_index({k: v["label"] for k, v in cat["minerals"].items()})
    run = next(iter(result["modes"].values()))
    rate = re.fullmatch(r"(?:(?:dissolution|reaction)\s+)?rate(?:s)?\s+of\s+(.+)|(.+?)\s+(?:dissolution\s+)?rate", low_clean)
    left = re.fullmatch(r"(?:how much\s+)?(.+?)\s+(?:left|remaining)", low_clean)
    if (rate or left) and run0.get("kinetic"):
        what, stem = ("rate", rate.group(1) or rate.group(2)) if rate else ("left", left.group(1))
        for k in run0["kinetic"]:
            kind, name = k.split(":", 1)
            label = cat.get("kinetic", {}).get(k, {}).get("label", "").lower()
            if kind == what and (stem == name.lower() or label.startswith(stem)):
                return {"group": "kinetic", "key": k}
    if "vapour_share" in run0.get("phases", {}):      # "vapour share", "brine salinity"
        for k, v in cat.get("phases", {}).items():
            if low_clean in (v["label"].lower(), k.replace("_", " "), k.replace("_", " ").replace("vapour", "vapor")):
                return {"group": "phases", "key": k}
    if run0.get("composition"):                       # "anorthite in plagioclase"
        for k, v in cat.get("composition", {}).items():
            if low_clean in (v["label"].lower(), "{1} in {0}".format(*k.lower().split(":", 1))):
                return {"group": "composition", "key": k}
    hit = resolve(low_clean or t, mineral_index)
    if hit and (sat or hit not in species_index.values()):
        if sat and hit in run["saturation"]:
            return {"group": "saturation", "key": hit}
        if hit in run["minerals"]:
            return {"group": "minerals", "key": hit}
        if hit in run["saturation"]:
            return {"group": "saturation", "key": hit}
    hit = resolve(low_clean or t, species_index)
    if hit:
        return {"group": "species", "key": hit}
    if hit is None:
        word = ELEMENT_NAMES.get(low_clean)
        for cand in (word, t, t.capitalize()):
            if cand and cand in cat["elements"]:
                return {"group": "elements", "key": cand}
    return None

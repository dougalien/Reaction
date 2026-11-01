"""Rocks and glasses as reactants: an oxide analysis that dissolves steadily, all at once in proportion
(congruently), as volcanic glass does.  A rock is not a mineral, so it never saturates or precipitates;
what the fluid cannot hold forms secondary minerals.

Amounts are in grams.  Each gram adds the moles of every element in it, written over the basis as the
formation reactions are: one basis cation per metal atom, H+ for the charge and H2O for the oxygen
(FeO = Fe++ + H2O - 2 H+).  Ferric iron also needs O2(g): when the fluid's redox state is not tracked,
all iron enters as Fe(II) and the run says so.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from ..network import OXYGEN
from ..thermo.eosdb import ATOMIC_WEIGHTS, DEFAULT_MASTERS, _formation, _Master


@dataclass(frozen=True)
class Rock:
    key: str
    label: str
    oxides: dict[str, float]          # weight percent
    ferric_fraction: float            # Fe(III) / total Fe, by moles
    source: str
    names: tuple[str, ...] = ()


ROCKS: dict[str, Rock] = {
    "basalt_glass": Rock(
        "basalt_glass", "Basalt glass",
        {"SiO2": 50.5, "TiO2": 1.7, "Al2O3": 14.7, "FeO": 10.4, "MnO": 0.18, "MgO": 7.6,
         "CaO": 11.4, "Na2O": 2.8, "K2O": 0.16, "P2O5": 0.18},
        0.10,
        "rounded mid-ocean ridge basalt glass, close to the global average of Gale and others "
        "(2013); FeO is total iron, a tenth of it ferric",
        ("basalt glass", "basaltic glass", "morb glass", "morb", "basalt", "glass")),
}

# Crystalline rocks: a mode (weight percent of each mineral) rather than a glass.  Putting one on the card
# writes its minerals as separate reactants, so each dissolves (and with kinetics, at its own rate) and
# each can also precipitate again.  Plagioclase, pyroxene and olivine are written as their endmembers.
# The minor phases are in too: Fe-Ti oxides as magnetite (the database has no titanium, so ilmenite and the
# titanium of titanomagnetite are left out), sulfide as pyrrhotite or pyrite and phosphate as hydroxyapatite.
# They need the fluid's redox state, which every fluid now has (air unless the card says otherwise).
@dataclass(frozen=True)
class MineralRock:
    key: str
    label: str
    minerals: dict[str, float]        # weight percent of each mineral (database names)
    description: str
    source: str
    names: tuple[str, ...] = ()


MINERAL_ROCKS: dict[str, MineralRock] = {r.key: r for r in (
    MineralRock("crystalline_basalt", "Crystalline basalt",
                {"anorthite": 33.0, "albite": 22.0, "diopside": 20.0, "hedenbergite": 10.0, "enstatite": 5.0,
                 "forsterite": 7.0, "fayalite": 3.0, "magnetite": 3.0, "hydroxyapatite": 0.5, "pyrrhotite": 0.2},
                "Diabase-like basalt: plagioclase An60, augite, a little orthopyroxene and olivine Fo80.",
                "rounded mode of a tholeiitic dolerite (55% plagioclase, 35% pyroxene, 10% olivine) with a few "
                "percent of Fe-Ti oxide (as magnetite), apatite and pyrrhotite",
                ("crystalline basalt", "diabase", "dolerite", "basalt crystalline")),
    MineralRock("granite", "Granite",
                {"quartz": 30.0, "K-feldspar": 35.0, "albite": 22.0, "anorthite": 3.0, "muscovite": 5.0,
                 "annite": 3.0, "phlogopite": 2.0, "magnetite": 1.0, "hydroxyapatite": 0.3, "pyrite": 0.1},
                "Granite: quartz, K-feldspar, oligoclase (An12), muscovite and biotite.",
                "rounded mode of a two-mica granite; biotite as annite and phlogopite; accessory magnetite, "
                "apatite and pyrite",
                ("granite", "granitic rock", "granodiorite")),
    MineralRock("lherzolite", "Lherzolite",
                {"forsterite": 52.0, "fayalite": 8.0, "enstatite": 22.0, "ferrosilite": 3.0, "diopside": 12.0,
                 "spinel": 3.0},
                "Spinel lherzolite (upper mantle): olivine Fo90, orthopyroxene En90, diopside and spinel.",
                "rounded mode of a fertile spinel lherzolite (60% olivine, 25% orthopyroxene, 12% clinopyroxene, "
                "3% spinel)",
                ("lherzolite", "peridotite", "mantle peridotite", "mantle rock", "spinel lherzolite")),
    MineralRock("mars_meteorite", "Mars meteorite",
                {"enstatite": 22.0, "ferrosilite": 20.0, "diopside": 14.0, "hedenbergite": 14.0, "anorthite": 12.0,
                 "albite": 11.0, "magnetite": 2.0, "hydroxyapatite": 0.5, "pyrrhotite": 0.3},
                "Basaltic shergottite like Shergotty: iron-rich pigeonite and augite with maskelynite "
                "(shocked plagioclase, An50).",
                "rounded mode of the Shergotty meteorite (about 70% pyroxene and 23% maskelynite; Stolper and "
                "McSween 1979) with its few percent of Fe-Ti oxide (as magnetite), phosphate (as apatite) and "
                "pyrrhotite",
                ("mars meteorite", "martian meteorite", "mars meteor", "martian meteor", "shergottite",
                 "shergotty", "martian basalt", "mars basalt", "mars rock", "martian rock")),
    MineralRock("sandstone", "Sandstone",
                {"quartz": 80.0, "K-feldspar": 8.0, "albite": 4.0, "kaolinite": 4.0, "calcite": 4.0},
                "Subarkose sandstone: quartz with some feldspar, kaolinite and calcite cement.",
                "rounded mode of a subarkosic sandstone",
                ("sandstone", "sand", "arkose", "subarkose", "quartz sandstone")),
    MineralRock("limestone_calcitic", "Limestone (calcitic)",
                {"calcite": 97.0, "quartz": 3.0},
                "Calcitic limestone: calcite with a little quartz.",
                "illustrative composition of a pure limestone",
                ("limestone", "calcitic limestone", "limestone calcitic", "limestone (calcitic)", "calcite limestone")),
    MineralRock("limestone_dolomitic", "Limestone (dolomitic)",
                {"dolomite": 90.0, "calcite": 7.0, "quartz": 3.0},
                "Dolomitic limestone (dolostone): dolomite with some calcite and a little quartz.",
                "illustrative composition of a dolostone",
                ("dolomitic limestone", "limestone dolomitic", "limestone (dolomitic)", "dolostone",
                 "dolomite rock")),
    MineralRock("limestone_aragonitic", "Limestone (aragonitic)",
                {"aragonite": 90.0, "calcite": 10.0},
                "Aragonitic carbonate sediment, as young ooid or coral sand: aragonite with some calcite.",
                "illustrative composition of a modern shallow-water carbonate sediment",
                ("aragonitic limestone", "limestone aragonitic", "limestone (aragonitic)", "aragonite limestone",
                 "ooid sand", "coral sand")),
)}

# The rock cards offered on the scenario card, in order: basalt glass dissolves as a glass, the rest as minerals
ROCK_CARDS = ("basalt_glass", "crystalline_basalt", "granite", "lherzolite", "mars_meteorite", "sandstone",
              "limestone_calcitic", "limestone_dolomitic", "limestone_aragonitic")


def find_mineral_rock(name: str) -> MineralRock | None:
    low = re.sub(r"\s+", " ", name.strip().lower().replace("_", " "))
    for rock in MINERAL_ROCKS.values():
        if low == rock.key.replace("_", " ") or low == rock.label.lower() or low in rock.names:
            return rock
    return None


def mineral_reactants(rock: MineralRock, grams: float) -> list[dict]:
    """The minerals of ``grams`` of a crystalline rock, as reactants."""
    total = sum(rock.minerals.values())
    return [{"name": m, "amount": float(f"{grams * w / total:.6g}"), "unit": "g"} for m, w in rock.minerals.items()]


def rock_cards() -> list[dict]:
    """The rock cards for the control window: key, label and what each puts on the card."""
    out = []
    for key in ROCK_CARDS:
        if key in ROCKS:
            r = ROCKS[key]
            out.append({"key": key, "label": r.label, "description": "Mid-ocean ridge basalt glass, dissolving "
                        "congruently as one reactant (" + r.source + ")."})
        else:
            r = MINERAL_ROCKS[key]
            mode = ", ".join(f"{w:g}% {m}" for m, w in r.minerals.items())
            out.append({"key": key, "label": r.label, "description": f"{r.description} {mode}. Source: {r.source}."})
    return out


_OXIDE = re.compile(r"([A-Z][a-z]?)(\d*)O(\d*)$")


DEFAULT_FO = 90.0
_OLIVINE = re.compile(r"(?:(?:olivine|ol|forsterite)\s*\(?\s*fo\s*(?P<a>\d+(?:\.\d+)?)\s*\)?|fo\s*(?P<b>\d+(?:\.\d+)?)"
                      r"(?:\s+olivine)?|(?P<plain>olivine))")


def olivine(fo: float = DEFAULT_FO) -> Rock:
    """Olivine (Mg,Fe)2SiO4 with ``fo`` mole percent forsterite, as a rock that dissolves congruently."""
    x = fo / 100.0
    w = ATOMIC_WEIGHTS
    oxides = {"MgO": 2 * x * (w["Mg"] + w["O"]), "FeO": 2 * (1 - x) * (w["Fe"] + w["O"]),
              "SiO2": w["Si"] + 2 * w["O"]}
    total = sum(oxides.values())
    return Rock(f"olivine_fo{fo:g}", f"Olivine (Fo{fo:g})", {k: 100.0 * v / total for k, v in oxides.items()},
                0.0, f"stoichiometric olivine, {fo:g} mole percent forsterite; it dissolves congruently and "
                     f"is not itself allowed to precipitate (forsterite, fayalite or the olivine solid "
                     f"solution can form)", (f"olivine fo{fo:g}", f"fo{fo:g}"))


def find(name: str) -> Rock | None:
    low = name.strip().lower().replace("_", " ")
    for rock in ROCKS.values():
        if low == rock.key.replace("_", " ") or low == rock.label.lower() or low in rock.names:
            return rock
    m = _OLIVINE.fullmatch(re.sub(r"\s+", " ", low.replace("(", " ").replace(")", " ")).strip())
    if m:
        fo = float(m.group("a") or m.group("b") or DEFAULT_FO)
        return olivine(fo) if 0.0 <= fo <= 100.0 else None
    return None


def get(key: str) -> Rock:
    """The rock with this key (a named rock, or an olivine such as ``olivine_fo88``)."""
    rock = ROCKS.get(key) or find(key)
    if rock is None:
        raise KeyError(key)
    return rock


def elements_per_gram(rock: Rock) -> tuple[dict[str, float], dict[str, float]]:
    """Moles of each element in one gram, total iron split by ``ferric_fraction`` into "Fe" (ferrous)
    and "Fe3" (ferric); and the oxygen that comes with each of them."""
    total = sum(rock.oxides.values())
    out: dict[str, float] = {}
    oxygen: dict[str, float] = {}
    for oxide, wt in rock.oxides.items():
        m = _OXIDE.fullmatch(oxide)
        el, n_el, n_o = m.group(1), int(m.group(2) or 1), int(m.group(3) or 1)
        moles = wt / total / (ATOMIC_WEIGHTS[el] * n_el + ATOMIC_WEIGHTS["O"] * n_o)
        out[el] = out.get(el, 0.0) + moles * n_el
        oxygen[el] = oxygen.get(el, 0.0) + moles * n_o
    if "Fe" in out and rock.ferric_fraction > 0:
        fe3 = out["Fe"] * rock.ferric_fraction
        out["Fe"] -= fe3
        out["Fe3"] = fe3
        oxygen["Fe"] -= fe3
        oxygen["Fe3"] = 1.5 * fe3                     # Fe2O3
    return out, oxygen


def components_per_gram(rock: Rock, db, basis: list[str]) -> tuple[dict[str, float], list[str]]:
    """Moles of basis components added per gram, and sentences about what was left out or changed."""
    notes: list[str] = []
    el, oxygen = elements_per_gram(rock)
    if "Fe3" in el:
        fe3 = el.pop("Fe3")
        el["Fe"] = el.get("Fe", 0.0) + fe3
        if OXYGEN in basis:
            oxygen["Fe"] += oxygen.pop("Fe3")
        else:
            oxygen["Fe"] += oxygen.pop("Fe3") - 0.5 * fe3       # as FeO
            notes.append(f"{rock.label}: the fluid's redox state is not tracked, so all its iron "
                         f"is added as Fe(II).")
    masters: dict[str, _Master] = {}
    left_out = []
    for e in list(el):
        mname = DEFAULT_MASTERS.get(e)
        if mname is None or mname not in db.species or mname not in basis:
            left_out.append(e)
            el.pop(e)
            oxygen.pop(e, None)
            continue
        sp = db.species[mname]
        masters[e] = _Master(mname, sp.elements, sp.charge)
    if left_out:
        notes.append(f"{rock.label}: {', '.join(left_out)} is left out (no dissolved species for it "
                     f"in the data).")
    el["O"] = sum(oxygen.values())
    nu = _formation(el, 0.0, masters) or {}
    comp = {k: v for k, v in nu.items() if abs(v) > 1e-15}
    if OXYGEN not in basis and abs(comp.pop(OXYGEN, 0.0)) > 1e-9:
        raise ValueError(f"{rock.label} needs the fluid's redox state")
    return comp, notes

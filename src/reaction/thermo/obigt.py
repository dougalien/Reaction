"""Loader for the OBIGT thermodynamic database (CHNOSZ; csv files, GPL-3).

Reads rows of the csv files into HKF (aqueous) or CGL (crystal, gas, liquid)
equation-of-state objects.  Files are read from a user-supplied folder and later
files override earlier ones for the same (name, state), as in CHNOSZ.

Column layout:  name, abbrv, formula, state, ref1, ref2, date, model, E_units,
G, H, S, Cp, V, a1.a, a2.b, a3.c, a4.d, c1.e, c2.f, omega.lambda, z.T
HKF rows carry the SUPCRT scaling a1*1e-1, a2*1e2, a4*1e4, c2*1e4, omega*1e5.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

from . import berman
from .crystal import CrystalSpecies
from .hkf import HKFSpecies

CAL = 4.184
CM3BAR = 0.1            # J per cm3 bar

DEFAULT_FILES = ("inorganic_aq.csv", "inorganic_cr.csv", "inorganic_gas.csv",
                 "inorganic_liq.csv", "H2O_aq.csv")
SUPCRT92_FILES = ("SUPCRT92.csv",)
SLOP98_FILES = ("SLOP98-a.csv", "SLOP98-b.csv")
BERMAN_FILES = ("Berman_cr.csv",)
BERMAN_DIR = "Berman"          # parameter files, next to the OBIGT csv files

# Minerals of the older (Helgeson et al., 1978) dataset that the Berman dataset replaces under
# another name; when Berman minerals are loaded these are dropped so the two cannot compete.
SUPERSEDED_BY_BERMAN = {
    "microcline,maximum": "K-feldspar,low", "sanidine,high": "K-feldspar,high",
    "clinochlore,14A": "clinochlore", "clinochlore,7A": "clinochlore",
    "daphnite,14A": "daphnite", "daphnite,7A": "daphnite",
    "amesite,14A": "Mg-amesite", "amesite,7A": "Mg-amesite", "cordierite,dry": "cordierite",
    "cordierite,hydrous": "cordierite", "titanite": "sphene", "epidote,ordered": "epidote",
    "dolomite,ordered": "dolomite", "dolomite,disordered": "dolomite",
    "cristobalite": "cristobalite,alpha", "ferrotremolite": "ferro-actinolite",
}


class ObigtError(ValueError):
    pass


# --- chemical formulas ----------------------------------------------------------------------
_TOKEN = re.compile(r"([A-Z][a-z]?|\(|\))(\d*\.?\d*)")


def parse_formula(formula: str) -> tuple[dict[str, float], float]:
    """Elements and charge of a formula such as 'Ca+2', 'HCO3-', 'Fe2SiO4', 'Ca(OH)2', or a
    hydrate such as 'NaAlSi2O6*H2O'."""
    f = formula.strip()
    if "*" in f:
        total: dict[str, float] = {}
        for part in f.split("*"):
            m = re.match(r"(\d*\.?\d+)?(.*)$", part)
            n = float(m.group(1)) if m.group(1) else 1.0
            elements, _ = parse_formula(m.group(2))
            for el, c in elements.items():
                total[el] = total.get(el, 0.0) + n * c
        return total, 0.0
    charge = 0.0
    m = re.search(r"([+-])(\d*)$", f)
    if m and f != "H" and len(f) > 1:
        sign = 1.0 if m.group(1) == "+" else -1.0
        charge = sign * (float(m.group(2)) if m.group(2) else 1.0)
        f = f[:m.start()]
    stack: list[dict[str, float]] = [{}]
    pos = 0
    while pos < len(f):
        m = _TOKEN.match(f, pos)
        if not m or m.end() == pos:
            raise ObigtError(f"cannot parse formula {formula!r}")
        tok, num = m.groups()
        pos = m.end()
        n = float(num) if num else 1.0
        if tok == "(":
            stack.append({})
        elif tok == ")":
            inner = stack.pop()
            for el, c in inner.items():
                stack[-1][el] = stack[-1].get(el, 0.0) + c * n
        else:
            stack[-1][tok] = stack[-1].get(tok, 0.0) + n
    return stack[0], charge


# --- species -------------------------------------------------------------------------------------
@dataclass
class ObigtSpecies:
    name: str
    state: str                     # aq, cr, cr2, ..., gas, liq
    formula: str
    elements: dict[str, float]
    charge: float
    model: str
    eos: object                    # HKFSpecies | CrystalSpecies | "water"
    source: str                    # file the row came from
    refs: tuple[str, str] = ("", "")

    @property
    def is_aqueous(self) -> bool:
        return self.state == "aq"


def _num(text: str) -> float | None:
    t = text.strip()
    if t in ("", "NA"):
        return None
    return float(t)


def _row_to_species(row: dict, source: str, berman_params=None) -> ObigtSpecies | None:
    name, state, model = row["name"], row["state"], row["model"]
    if model == "Berman":
        params = (berman_params or {}).get(name)
        elements: dict[str, float] = {}
        try:
            elements, charge = parse_formula(row["formula"])
            eos = berman.BermanSpecies(name, state, params, berman.element_entropy(elements)) \
                if params is not None else None
        except (ObigtError, ValueError):
            eos = None
        return ObigtSpecies(name, state, row["formula"], elements if eos else {}, 0.0, model, eos,
                            source, (row["ref1"], row["ref2"]))
    if model in ("AD", "DEW"):                    # other water models
        return None
    unit = CAL if row["E_units"] == "cal" else 1.0
    g, h, s, cp, v = (_num(row[k]) for k in ("G", "H", "S", "Cp", "V"))
    p = [_num(row[k]) for k in ("a1.a", "a2.b", "a3.c", "a4.d", "c1.e", "c2.f",
                                "omega.lambda", "z.T")]
    refs = (row["ref1"], row["ref2"])
    if model == "H2O":
        return ObigtSpecies("water", "liq", "H2O", {"H": 2.0, "O": 1.0}, 0.0, "H2O", "water",
                            source, refs)
    formula = row["formula"]
    try:
        elements, charge = parse_formula(formula)
    except ObigtError:
        return None                               # e.g. the electron "(Z-1)"
    if model == "HKF":
        z = p[7] if p[7] is not None else charge
        if all(x is None for x in p[:6]) and cp is None and v is None:
            return ObigtSpecies(name, state, formula, elements, z, model, None, source, refs)
        # a1..a4, c1, c2, omega are tabulated in SUPCRT's scaled units
        if all(x is None for x in p[:4]):                 # only V known: a1 carries it
            a = [(v or 0.0) * CM3BAR, 0.0, 0.0, 0.0]
        else:
            a = [(p[0] or 0.0) * 1e-1 * unit, (p[1] or 0.0) * 1e2 * unit,
                 (p[2] or 0.0) * unit, (p[3] or 0.0) * 1e4 * unit]
        if all(x is None for x in p[4:6]):                # only Cp known: c1 carries it
            c1, c2 = (cp or 0.0) * unit, 0.0
        else:
            c1, c2 = (p[4] or 0.0) * unit, (p[5] or 0.0) * 1e4 * unit
        eos = HKFSpecies(name, z, (g or 0.0) * unit, (s or 0.0) * unit, a[0], a[1], a[2], a[3],
                         c1, c2, (p[6] or 0.0) * 1e5 * unit,
                         cp * unit if cp is not None else None, v)
        return ObigtSpecies(name, state, formula, elements, z, model, eos, source, refs)
    if model in ("CGL", "CGL_Ttr"):
        if p[6] not in (None, 0.0):               # lambda term in the heat capacity: not supported
            return ObigtSpecies(name, state, formula, elements, charge, model, None, source, refs)
        if g is None:                             # G from H and S, else no usable data
            if h is None or s is None:
                return ObigtSpecies(name, state, formula, elements, charge, model, None, source, refs)
            try:
                g = h - 298.15 * (s - berman.element_entropy(elements) / unit)
            except ValueError:
                return ObigtSpecies(name, state, formula, elements, charge, model, None, source, refs)
        eos = CrystalSpecies(
            name=name, state=state,
            g0=(g if g is not None else 0.0) * unit,
            s0=(s if s is not None else 0.0) * unit,
            v0=v if v is not None else 0.0,
            a=(p[0] or 0.0) * unit, b=(p[1] or 0.0) * unit, c=(p[2] or 0.0) * unit,
            d=(p[3] or 0.0) * unit, e=(p[4] or 0.0) * unit,
            t_max=p[7] or 0.0,
            cp0=(cp * unit) if cp is not None else None)
        return ObigtSpecies(name, state, formula, elements, charge, model, eos, source, refs)
    return None


@dataclass
class ObigtData:
    """Species by (name, state)."""
    species: dict[tuple[str, str], ObigtSpecies] = field(default_factory=dict)
    skipped: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)      # optional parts of the data not found

    def get(self, name: str, state: str = "aq") -> ObigtSpecies | None:
        return self.species.get((name, state))

    def polymorphs(self, name: str) -> list[ObigtSpecies]:
        """The phases of a mineral in order of temperature (cr, cr2, cr3, ...)."""
        return sorted((sp for (n, st), sp in self.species.items()
                       if n == name and re.fullmatch(r"cr\d?", st)), key=lambda sp: sp.state)


def load_obigt(folder: str | Path, files=DEFAULT_FILES) -> ObigtData:
    """Read OBIGT csv files from a folder; later files override earlier rows."""
    folder = Path(folder)
    data = ObigtData()
    berman_params = None
    for fname in files:
        path = folder / fname
        if fname in BERMAN_FILES:
            pdir = folder / BERMAN_DIR
            if not pdir.is_dir():                 # older data folders: keep the Helgeson minerals
                data.skipped.append(f"{fname}: the Berman parameter files are not in {pdir} "
                                    f"(CHNOSZ inst/extdata/Berman); using the older minerals")
                data.missing.append(BERMAN_DIR)
                continue
            berman_params = berman.load_parameters(pdir)
        if not path.exists():
            raise ObigtError(f"{path} not found")
        with path.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                sp = _row_to_species(row, fname, berman_params)
                if sp is None:
                    data.skipped.append(f"{fname}: {row['name']} ({row['state']}) {row['model']}")
                    continue
                if sp.eos is None:
                    data.skipped.append(f"{fname}: {row['name']} ({row['state']}) has no usable EOS")
                    data.species.pop((sp.name, sp.state), None)
                    continue
                if re.fullmatch(r"cr\d", sp.state) and (sp.name, "cr") not in data.species:
                    data.skipped.append(f"{fname}: {row['name']} ({row['state']}) has no low-"
                                        f"temperature phase")
                    continue                      # its G is relative to a phase without data
                if sp.model == "Berman" and "," in sp.name:
                    base = data.species.get((sp.name.split(",")[0], "cr"))
                    if base is not None and base.model == "Berman" and \
                            _same_params(base.eos.params, sp.eos.params):
                        data.skipped.append(f"{fname}: {sp.name} is the same as {base.name}")
                        continue                  # albite,low is albite: one name, not a tie
                if sp.model == "Berman":          # replaces every polymorph of the older data
                    for key in [k for k in data.species
                                if re.fullmatch(r"cr\d?", k[1])
                                and (k[0] == sp.name or SUPERSEDED_BY_BERMAN.get(k[0]) == sp.name)]:
                        data.species.pop(key)
                data.species[(sp.name, sp.state)] = sp
    return data


def _same_params(a, b) -> bool:
    return (a.h0, a.s0, a.v0, a.k, a.v, a.lam, a.dis) == (b.h0, b.s0, b.v0, b.k, b.v, b.lam, b.dis)

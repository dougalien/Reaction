"""PHREEQC's pitzer.dat as one consistent model: log K, molar volumes and the pressure terms.

Pitzer interaction parameters are fitted together with a particular set of equilibrium constants, so
using them with another database's log K (HKF here) mixes two models.  This module reads PHREEQC's
pitzer.dat (USGS, public domain; ``data/pitzer.dat``, unchanged copy of
github.com/usgs-coupled/phreeqc3 ``database/pitzer.dat`` retrieved 2026-11-11, sha256 dada3271...83ba)
and evaluates, the way PHREEQC 3 does:

* log K(T) of its aqueous species and phases: the analytical expression
  A1 + A2 T + A3/T + A4 log10 T + A5/T^2 + A6 T^2 when given, otherwise log K(25 C) with the
  van 't Hoff term of delta_h (PHREEQC ``k_calc`` and ``select_log_k_expression``);
* the pressure term -delta_V (P - 1 atm) / (RT ln 10), with molar volumes of aqueous species from the
  Redlich-type equation of Appelo, Parkhurst and Post (2014, GCA 125, 49-67) and Appelo (2015, Appl.
  Geochem. 55, 62-71): the SUPCRT92 form at infinite dilution (Johnson et al., 1992) with a Born term,
  plus z^2/2 A_V I^0.5 / (1 + b_Av B I^0.5) and (i1 + i2/(T-228) + i3 (T-228)) I^i4; minerals have a
  constant molar volume and gases none (PHREEQC ``calc_vm``, ``calc_lk_phase``);
* the properties of pure water these need: density along the saturation curve from Wagner and Pruss
  (2002, eq. 2.6) with PHREEQC's pressure polynomial (fitted 0-300 C, to 1000 atm), the dielectric
  constant of Bradley and Pitzer (1979, J. Phys. Chem. 83, 1599), and from them the Debye-Hueckel
  A, B, the Pitzer A-phi, the volume slope A_V and the Born Q (PHREEQC ``calc_rho_0``,
  ``calc_dielectrics``).  PHREEQC holds these at 350 C above 350 C, and so does this module.

The ionic-strength part of each species' volume is a pressure effect on its activity coefficient
(PHREEQC puts it in log K; it is the same thing): it is returned separately by
:meth:`PitzerData.volume_excess` and added to ln gamma by the Pitzer model.

Range, from pitzer.dat and Appelo (2015): 0-200 C and 1-1000 atm; some temperature functions of
log K are fitted to narrower ranges (noted on their lines of pitzer.dat).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).with_name("data") / "pitzer.dat"

# PHREEQC constants (global_structures.h)
R_KJ = 0.00831470            # kJ/K/mol
R_LITER_ATM = 0.0820597      # L atm/K/mol
PA_PER_ATM = 1.01325e5
AVOGADRO = 6.02252e23
LN10 = math.log(10.0)
T_CAP_C = 350.0              # PHREEQC's water properties are held at 350 C above this
T_MAX_C = 200.0              # pitzer.dat is fitted to about this temperature
P_MAX_ATM = 1000.0           # and to this pressure

# Reaction (OBIGT) names for pitzer.dat species written with different water content or names:
# pitzer.dat name -> (Reaction name, moles of H2O to subtract from the pitzer.dat formula)
RENAME = {"H4SiO4": ("SiO2", 2), "H3SiO4-": ("HSiO3-", 1), "B(OH)4-": ("BO2-", 2)}
# pitzer.dat species left out: redox-uncoupled gases and sulfide (separate pseudo-elements there)
SKIP_ELEMENTS = {"Hdg", "Oxg", "Mtg", "Ntg", "Sg"}


@dataclass
class Entry:
    """One aqueous species (``kind`` "aq") or phase ("mineral"/"gas") of pitzer.dat."""
    name: str
    kind: str
    reaction: dict[str, float]          # PHREEQC species -> coefficient (+ products, - reactants)
    defined: str                        # the species/phase the reaction defines
    log_k: float = 0.0
    delta_h: float = 0.0                # kJ/mol
    analytic: tuple[float, ...] = ()
    vm_aq: tuple[float, ...] = ()       # a1..a4 (cm3 units), wref, b_Av, i1..i4, as PHREEQC stores them
    vm_phase: float = 0.0               # cm3/mol
    formula: str = ""

    def log_k_t(self, t_k: float) -> float:
        if any(self.analytic):
            a = tuple(self.analytic) + (0.0,) * (6 - len(self.analytic))
            return a[0] + a[1] * t_k + a[2] / t_k + a[3] * math.log10(t_k) + a[4] / t_k ** 2 + a[5] * t_k ** 2
        return self.log_k - self.delta_h * (298.15 - t_k) / (LN10 * R_KJ * t_k * 298.15)


_SPECIES_RE = re.compile(r"^\s*(?:(\d*\.?\d+)\s*)?(\S+)\s*$")


def _side(text: str) -> dict[str, float]:
    """'2 H+ + CO3-2' -> {'H+': 2, 'CO3-2': 1}; a leading '- X' term counts as negative."""
    out: dict[str, float] = {}
    tokens = re.split(r"\s+(?=[+-]\s)", " " + text.strip())
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        sign = 1.0
        if tok[0] in "+-" and len(tok) > 1 and tok[1] == " ":
            sign = -1.0 if tok[0] == "-" else 1.0
            tok = tok[2:].strip()
        m = _SPECIES_RE.match(tok)
        if not m:
            raise ValueError(f"cannot read species term {tok!r}")
        coef = float(m.group(1)) if m.group(1) else 1.0
        out[m.group(2)] = out.get(m.group(2), 0.0) + sign * coef
    return out


def _numbers(text: str) -> list[float]:
    vals = []
    for tok in text.split():
        try:
            vals.append(float(tok))
        except ValueError:
            break
    return vals


def _options(line: str):
    """Split an option line on ';' into (keyword, rest) pairs, keywords without the leading '-'."""
    for part in line.split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, rest = part.partition(" ")
        rest = rest.strip()
        if "\t" in key:
            key, _, more = key.partition("\t")
            rest = (more + " " + rest).strip()
        yield key.lstrip("-").lower(), rest


def _apply(entry: Entry, key: str, rest: str) -> None:
    nums = _numbers(rest)
    if key in ("log_k", "logk"):
        entry.log_k = nums[0]
    elif key in ("delta_h", "deltah"):
        factor = 4.184 if "kcal" in rest.lower() else 1.0
        entry.delta_h = nums[0] * factor
    elif key.startswith("analytic") or key == "a_e":
        entry.analytic = tuple(nums)
    elif key == "vm":
        if entry.kind == "aq":
            v = list(nums) + [0.0] * (10 - len(nums))
            if len(nums) < 10:
                v[9] = 1.0
            # PHREEQC read_aq_species_vm_parms: cal units to cm3, wref * 1e5
            v[0] *= 41.84004e-1
            v[1] *= 41.84004e2
            v[2] *= 41.84004
            v[3] *= 41.84004e4
            v[4] *= 1e5
            entry.vm_aq = tuple(v)
        else:
            entry.vm_phase = nums[0]


@lru_cache(maxsize=1)
def read(path: str | Path = DATA) -> dict[str, Entry]:
    """The SOLUTION_SPECIES and PHASES entries of pitzer.dat, by name (species by formula)."""
    entries: dict[str, Entry] = {}
    block = None
    current: Entry | None = None
    pending_phase: str | None = None
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not raw[0].isspace() and re.fullmatch(r"[A-Z][A-Z_]+[A-Z]", line.strip()):
            block, current, pending_phase = line.strip(), None, None
            continue
        if block not in ("SOLUTION_SPECIES", "PHASES"):
            continue
        text = line.strip()
        if "=" in text and not text.startswith("-"):
            lhs, rhs = text.split("=", 1)
            left, right = _side(lhs), _side(rhs)
            if block == "SOLUTION_SPECIES":
                product = next(iter(right))
                if len(right) != 1 or right[product] != 1.0:
                    # PHREEQC writes e.g. "H2O = OH- + H+": the species defined is the first product
                    # that is not a master (OH- here); keep it general
                    product = [s for s in right if s not in ("H+", "H2O", "e-")][0]
                rxn = {s: c for s, c in right.items()}
                for s, c in left.items():
                    rxn[s] = rxn.get(s, 0.0) - c
                if lhs.strip() == rhs.strip():           # a master species, "X = X": no reaction
                    rxn = {}
                current = Entry(product, "aq", rxn, product)
                entries[product] = current
            else:
                name = pending_phase or "?"
                formula = next(iter(left))                 # the phase itself, coefficient 1
                left = dict(left)
                left[formula] -= 1.0
                rxn = {s: c for s, c in right.items()}
                for s, c in left.items():
                    if c:
                        rxn[s] = rxn.get(s, 0.0) - c
                kind = "gas" if name.endswith("(g)") else "mineral"
                current = Entry(name, kind, rxn, name, formula=formula)
                entries[name] = current
                pending_phase = None
            continue
        if block == "PHASES" and not raw[0].isspace() and len(text.split()) == 1:
            pending_phase, current = text, None
            continue
        if current is not None:
            for key, rest in _options(text):
                _apply(current, key, rest)
    return entries


# ---- pure water, as PHREEQC computes it --------------------------------------------------------------
@dataclass(frozen=True)
class WaterState:
    t_c: float               # the temperature used (held at 350 C above that)
    p_atm: float             # the pressure used (at least the saturation pressure)
    rho: float               # g/cm3
    kappa: float             # compressibility d ln rho / dP, 1/atm
    eps: float
    dh_a: float              # Debye-Hueckel A, log10 form, (kg/mol)^0.5
    dh_b: float              # Debye-Hueckel B, 1/angstrom (kg/mol)^0.5
    a_phi: float             # Pitzer A-phi (natural log)
    dh_av: float             # Debye-Hueckel volume slope, cm3/mol (kg/mol)^0.5
    q_born: float            # Born Q * 41.84004


def p_sat_atm(t_c: float) -> float:
    """PHREEQC's saturation pressure, atm (``calc_rho_0``)."""
    t_k = min(t_c, T_CAP_C) + 273.15
    return math.exp(11.6702 - 3816.44 / (t_k - 46.13))


@lru_cache(maxsize=4096)
def water(t_c: float, p_atm: float) -> WaterState:
    tc = min(t_c, T_CAP_C)
    t_k = tc + 273.15
    th = 1.0 - t_k / 647.096
    b = (1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.7469445e5)
    rho_sat = 322.0 * (1.0 + b[0] * th ** (1 / 3) + b[1] * th ** (2 / 3) + b[2] * th ** (5 / 3)
                       + b[3] * th ** (16 / 3) + b[4] * th ** (43 / 3) + b[5] * th ** (110 / 3))
    p0 = 5.1880000E-02 + tc * (-4.1885519E-04 + tc * (6.6780748E-06 + tc * (-3.6648699E-08 + tc * 8.3501912E-11)))
    p1 = -6.0251348E-06 + tc * (3.6696407E-07 + tc * (-9.2056269E-09 + tc * (6.7024182E-11 + tc * -1.5947241E-13)))
    p2 = -2.2983596E-09 + tc * (-4.0133819E-10 + tc * (1.2619821E-11 + tc * (-9.8952363E-14 + tc * 2.3363281E-16)))
    p3 = 7.0517647E-11 + tc * (6.8566831E-12 + tc * (-2.2829750E-13 + tc * (1.8113313E-15 + tc * -4.2475324E-18)))
    psat = math.exp(11.6702 - 3816.44 / (t_k - 46.13))
    pa = max(p_atm, psat)
    dp = pa - (psat - 1e-6)
    rho = max(rho_sat + dp * (p0 + dp * (p1 + dp * (p2 + math.sqrt(dp) * p3))), 0.01)   # kg/m3
    kappa = (p0 + dp * (2 * p1 + dp * (3 * p2 + math.sqrt(dp) * 3.5 * p3))) / rho
    rho_g = rho / 1e3
    # Bradley and Pitzer (1979)
    u = (3.4279e2, -5.0866e-3, 9.469e-7, -2.0525, 3.1159e3, -1.8289e2, -8.0325e3, 4.2142e6, 2.1417)
    d1000 = u[0] * math.exp(t_k * (u[1] + t_k * u[2]))
    c = u[3] + u[4] / (u[5] + t_k)
    bb = u[6] + u[7] / t_k + u[8] * t_k
    pb = pa * 1.01325
    eps = d1000 + c * math.log((bb + pb) / (bb + 1e3))
    if eps <= 0:
        eps = 10.0
    e2_dkt = 1.671008e-3 / (eps * t_k)
    dh_b = math.sqrt(8 * math.pi * AVOGADRO * e2_dkt * rho_g / 1e3)
    dh_a = dh_b * e2_dkt / (2.0 * LN10)
    a_phi = dh_b * e2_dkt / 6.0
    dh_av = dh_b * e2_dkt * R_LITER_ATM * 1e3 * t_k * (c / (bb + pb) * 1.01325 / eps - kappa / 3.0)
    q_born = c / (bb + pb) / eps / eps * 41.84004
    return WaterState(tc, pa, rho_g, kappa, eps, dh_a, dh_b / 1e8, a_phi, dh_av, q_born)


# ---- volumes, potentials and log K -------------------------------------------------------------------
def charge(name: str) -> float:
    m = re.search(r"([+-])(\d*)$", name)
    if not m or name in ("e-",) and False:
        return 0.0
    return (1.0 if m.group(1) == "+" else -1.0) * (int(m.group(2)) if m.group(2) else 1)


def vm0(e: Entry, w: WaterState) -> float:
    """Volume at infinite dilution, cm3/mol (PHREEQC ``calc_vm``, supcrt part)."""
    if not e.vm_aq:
        return 0.0
    a1, a2, a3, a4, wref = e.vm_aq[:5]
    pb_s = 2600.0 + w.p_atm * 1.01325
    tk_s = w.t_c + 45.15
    return a1 + a2 / pb_s + (a3 + a4 / pb_s) / tk_s - wref * w.q_born


def vm_excess(e: Entry, w: WaterState, ionic: float) -> float:
    """The ionic-strength part of the volume, cm3/mol."""
    if not e.vm_aq:
        return 0.0
    z = charge(e.defined)
    if z == 0:
        return 0.0
    b_av, i1, i2, i3, i4 = e.vm_aq[5:10]
    s = math.sqrt(max(ionic, 0.0))
    if b_av < 1e-5:
        v = z * z * 0.5 * w.dh_av * s
    else:
        v = z * z * 0.5 * w.dh_av * s / (1.0 + b_av * w.dh_b * s)
    if i1 or i2 or i3:
        tk_s = w.t_c + 45.15
        bi = i1 + i2 / tk_s + i3 * tk_s
        v += bi * (ionic if i4 == 1.0 else max(ionic, 0.0) ** i4)
    return v


def pressure_factor(t_c: float, p_atm: float) -> float:
    """(P - 1 atm) / (RT ln 10) in mol/cm3, as PHREEQC ``k_calc`` (delta_V in cm3/mol), with P at least
    the saturation pressure.  Uses the actual temperature, not the 350 C cap."""
    t_k = t_c + 273.15
    pa = max(p_atm, p_sat_atm(t_c))
    dp = pa * PA_PER_ATM - 101325.0
    return max(dp, 0.0) * 1e-9 / (LN10 * R_KJ * t_k)


@dataclass
class PitzerData:
    """log K of pitzer.dat's species and phases at one T and P, expressed as log-activity potentials."""
    t_c: float
    p_bar: float
    potentials: dict[str, float] = field(default_factory=dict)     # PHREEQC names (and phases)
    entries: dict[str, Entry] = field(default_factory=dict)
    water_state: WaterState | None = None
    pfac: float = 0.0

    @classmethod
    def at(cls, t_c: float, p_bar: float) -> "PitzerData":
        entries = read()
        p_atm = p_bar / 1.01325
        w = water(t_c, p_atm)
        pfac = pressure_factor(t_c, p_atm)
        t_k = t_c + 273.15
        vol = {"H2O": 18.016 / w.rho}
        for e in entries.values():
            if e.kind == "aq" and e.defined != "H2O":
                vol[e.defined] = vm0(e, w)
        # log-activity potentials p at 1 atm: masters (X = X) are 0, p(product) = log K(T) +
        # sum(coef * p(reactant)), p(phase) = sum(coef * p(product)) - log K(T).  The pressure term is then
        # -V (P - 1 atm) / (RT ln 10) on each species and mineral, so that every reaction gets
        # -delta_V (P - 1 atm) / (RT ln 10) as in PHREEQC; gases have no volume term.
        pot0: dict[str, float] = {}

        def p_of(name: str, stack=()) -> float:
            if name in pot0:
                return pot0[name]
            e = entries.get(name)
            if e is None or e.kind != "aq":
                pot0[name] = 0.0
                return 0.0
            if name in stack:
                raise ValueError(f"circular definition of {name}")
            total = e.log_k_t(t_k) if e.reaction else 0.0
            for s, c in e.reaction.items():
                if s != name:
                    total += -c * p_of(s, stack + (name,))
            pot0[name] = total
            return total

        pot: dict[str, float] = {}
        for name, e in entries.items():
            if e.kind == "aq":
                pot[name] = p_of(name) - vol.get(name, 0.0) * pfac
        pot.setdefault("H2O", -vol["H2O"] * pfac)
        pot.setdefault("H+", 0.0)
        for name, e in entries.items():
            if e.kind != "aq":
                p0 = sum(c * p_of(s) for s, c in e.reaction.items()) - e.log_k_t(t_k)
                pot[name] = p0 - (e.vm_phase if e.kind == "mineral" else 0.0) * pfac
        out = cls(t_c, p_bar, pot, entries, w, pfac)
        return out

    def volume_excess(self, name: str, ionic: float) -> float:
        """Pressure effect on log10 gamma of a pitzer.dat species: V_ex(I) (P - 1 atm) / (RT ln 10)."""
        e = self.entries.get(name)
        if e is None or e.kind != "aq" or not self.pfac:
            return 0.0
        return vm_excess(e, self.water_state, ionic) * self.pfac


# ---- matching pitzer.dat to the database's names ------------------------------------------------------
_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")
PHASE_ALIASES = {"SiO2(a)": "amorphous silica"}
SKIP_PHASES = {"H2O(g)"}           # water vapour stays with the water model


def elements(formula: str) -> dict[str, float]:
    """Element counts of a PHREEQC formula: 'CaSO4:2H2O', 'Na2Mg(SO4)2:4H2O', 'K2Ca5(SO4)6H2O'."""
    out: dict[str, float] = {}
    for part in formula.split(":"):
        m = re.match(r"^(\d*\.?\d+)(?=[A-Z(])", part)
        mult = float(m.group(1)) if m else 1.0
        if m:
            part = part[m.end():]
        for el, n in _parse_group(part).items():
            out[el] = out.get(el, 0.0) + mult * n
    return out


def _parse_group(text: str) -> dict[str, float]:
    stack: list[dict[str, float]] = [{}]
    i = 0
    while i < len(text):
        ch = text[i]
        if ch == "(":
            stack.append({})
            i += 1
        elif ch == ")":
            i += 1
            m = re.match(r"\d*\.?\d*", text[i:])
            n = float(m.group(0)) if m.group(0) else 1.0
            i += m.end()
            group = stack.pop()
            for el, c in group.items():
                stack[-1][el] = stack[-1].get(el, 0.0) + c * n
        else:
            m = _ELEMENT.match(text, i)
            if not m or not m.group(0):
                break                                  # a charge such as +2 ends the formula
            n = float(m.group(2)) if m.group(2) else 1.0
            stack[-1][m.group(1)] = stack[-1].get(m.group(1), 0.0) + n
            i = m.end()
    return stack[0]


def reaction_names(entries: dict[str, Entry], db) -> dict[str, tuple[str, float]]:
    """Database name -> (pitzer.dat name, moles of H2O to subtract) for every species and phase that
    pitzer.dat defines and the database has with the same composition."""
    out: dict[str, tuple[str, float]] = {}
    for name, e in entries.items():
        if e.kind == "aq":
            if any(el in name for el in SKIP_ELEMENTS):
                continue
            target, water_off = RENAME.get(name, (name, 0))
            if target in db:
                out[target] = (name, water_off)
            continue
        if name in SKIP_PHASES or any(el in name for el in SKIP_ELEMENTS):
            continue
        target = PHASE_ALIASES.get(name, name if e.kind == "gas" else name.lower())
        if target not in db:
            continue
        mine = {k: v for k, v in db[target].elements.items() if v}
        theirs = {k: v for k, v in elements(e.formula).items() if v}
        if mine.keys() == theirs.keys() and all(abs(mine[k] - theirs[k]) < 1e-9 for k in mine):
            out[target] = (name, 0)
    out.setdefault("H2O", ("H2O", 0))
    out.setdefault("H+", ("H+", 0))
    return out


def potential(data: "PitzerData", mapping: dict[str, tuple[str, float]], name: str) -> float | None:
    """The log-activity potential of a database species or phase, or None when pitzer.dat lacks it."""
    hit = mapping.get(name)
    if hit is None:
        return None
    pname, water_off = hit
    p = data.potentials.get(pname)
    if p is None:
        return None
    return p - water_off * data.potentials["H2O"]

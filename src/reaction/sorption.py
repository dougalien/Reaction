"""Ion exchange and surface complexation, solved together with the aqueous system.

Exchange (Gaines-Thomas, as PHREEQC): sites X- of fixed capacity; each species M X_z has the activity of its
equivalent fraction, beta = K a_M a_X^z, and the betas add up to one.  The exchanger holds exactly its capacity
in cation charge, so a fluid it is equilibrated with stays neutral.

Surface complexation (Dzombak and Morel 1990, hydrous ferric oxide): strong and weak sites with mass balances,
species n_i = K_i a(reactants) n(SOH) exp(-z_i F psi / RT), and the Gouy-Chapman relation between the surface
charge and the potential psi (the diffuse layer is not given its own composition, as in PHREEQC's default).

Clay minerals (``data/clays.tsv``): kaolinite, illite and montmorillonite carry exchange sites in proportion to
their cation exchange capacity; they share the exchange constants above.

Partition coefficients (``data/kd_allison_2005.tsv``, Allison and Allison 2005): a linear sorber holds Kd M times
the dissolved amount of each metal, M being the mass of solid.  Every dissolved species of the metal is copied
onto the solid in proportion (so the sorbed matter keeps its charge and ligands, and the fluid stays neutral);
Kd is set by the element and, for chromium, its oxidation state.

Data: ``data/sorption.tsv`` (PHREEQC's phreeqc.dat).  The log K are for 25 C; the exchange constants that have an
enthalpy follow van 't Hoff, the rest are held at 25 C.  Partition coefficients are held at their values.
"""
from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .thermo.electrostatics import epsilon

DATA = Path(__file__).parent / "data" / "sorption.tsv"
CLAYS = Path(__file__).parent / "data" / "clays.tsv"
KD_DATA = Path(__file__).parent / "data" / "kd_allison_2005.tsv"
KD_MEDIA = {"sediment": "bed sediment and its pore water", "suspended": "suspended matter in stream water",
            "soil": "soil and soil water"}
LN10 = math.log(10.0)
FARADAY = 96485.33212            # C/mol
GAS_R = 8.314462618              # J/mol/K
EPS0 = 8.8541878128e-12          # F/m
# hydrous ferric oxide, Dzombak and Morel (1990): 89 g per mol of Fe, 600 m2/g, 0.2 weak and 0.005 strong sites per Fe
HFO_G_PER_MOL_FE = 89.0
HFO_M2_PER_G = 600.0
HFO_WEAK_PER_FE = 0.2
HFO_STRONG_PER_FE = 0.005
ALIASES = {"H4SiO4": "SiO2", "H3BO3": "B(OH)3"}        # same species, written without the water


@dataclass(frozen=True)
class Reaction:
    kind: str                     # "exchange" or "surface"
    master: str                   # X-, Hfo_sOH, Hfo_wOH
    species: str                  # the sorbed species formed
    reactants: dict[str, float]   # aqueous species (and H2O) consumed; produced ones negative
    z: float                      # exchange: X- per species; surface: the species' charge
    log_k: float
    delta_h_kj: float | None

    def log_k_at(self, t_c: float) -> float:
        if not self.delta_h_kj:
            return self.log_k
        return self.log_k - self.delta_h_kj * 1000.0 / (GAS_R * LN10) * (1.0 / (t_c + 273.15) - 1.0 / 298.15)


def _terms(side: str) -> list[tuple[float, str]]:
    out = []
    for term in side.split(" + "):
        m = re.match(r"\s*(\d*\.?\d*)\s*(\S+)\s*$", term)
        out.append((float(m.group(1)) if m.group(1) else 1.0, m.group(2)))
    return out


def _charge(name: str) -> float:
    m = re.search(r"([+-])(\d*)$", name)
    if not m:
        return 0.0
    return (1.0 if m.group(1) == "+" else -1.0) * (float(m.group(2)) if m.group(2) else 1.0)


@lru_cache(maxsize=1)
def reactions() -> tuple[Reaction, ...]:
    rows = [line for line in DATA.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    out = []
    for r in csv.DictReader(rows, delimiter="\t"):
        left, right = r["reaction"].split(" = ")
        lt, rt = _terms(left), _terms(right)
        masters = [(c, s) for c, s in lt if s == "X-" or s.startswith("Hfo_")]
        sorbed = [(c, s) for c, s in rt if s.endswith("X") or "X2" in s or "X3" in s or s.startswith("Hfo_")]
        master, species = masters[0][1], sorbed[0][1]
        reactants: dict[str, float] = {}
        for c, s in lt:
            if s != master:
                reactants[s] = reactants.get(s, 0.0) + c
        for c, s in rt:
            if s != species:
                reactants[s] = reactants.get(s, 0.0) - c
        z = masters[0][0] if r["kind"] == "exchange" else _charge(species)
        out.append(Reaction(r["kind"], master, species, reactants, z, float(r["log_k"]),
                            float(r["delta_h_kj"]) if r["delta_h_kj"] else None))
    return tuple(out)


def _table(path: Path) -> list[dict]:
    rows = [line for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


@dataclass(frozen=True)
class Clay:
    name: str
    cec_eq_per_kg: float
    area_m2_per_g: float | None
    sample: str
    aliases: tuple[str, ...]


@lru_cache(maxsize=1)
def clays() -> dict[str, Clay]:
    out = {}
    for r in _table(CLAYS):
        out[r["clay"]] = Clay(r["clay"], float(r["cec_meq_per_kg"]) / 1000.0,
                              float(r["area_m2_per_g"]) if r["area_m2_per_g"] else None, r["sample"],
                              tuple(a.strip() for a in r["aliases"].split(",") if a.strip()))
    return out


def clay_name(word: str) -> str | None:
    """The clay a word names (kaolinite, kaolin, smectite ...), or None."""
    w = word.strip().lower()
    for c in clays().values():
        if w == c.name or w in c.aliases:
            return c.name
    return None


@dataclass(frozen=True)
class KdRow:
    element: str
    oxidation: int | None          # None: every oxidation state
    medium: str
    log_kd: float
    stat: str
    log_min: float
    log_max: float
    n: int | None

    @property
    def label(self) -> str:
        return f"{self.element}({_roman(self.oxidation)})" if self.oxidation else self.element


def _roman(n: int) -> str:
    return {2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI"}.get(n, str(n))


@lru_cache(maxsize=1)
def kd_table() -> tuple[KdRow, ...]:
    return tuple(KdRow(r["element"], int(r["oxidation"]) if r["oxidation"] else None, r["medium"],
                       float(r["log_kd"]), r["stat"], float(r["log_min"]), float(r["log_max"]),
                       int(r["n"]) if r["n"] else None) for r in _table(KD_DATA))


def kd_rows(medium: str, elements=None, skip=()) -> tuple[KdRow, ...]:
    """The partition coefficients of one medium (sediment, suspended, soil), for ``elements`` (default all)
    and leaving out the elements in ``skip``."""
    if medium not in KD_MEDIA:
        raise ValueError(f"unknown Kd medium {medium!r} (known: {', '.join(KD_MEDIA)})")
    return tuple(r for r in kd_table() if r.medium == medium and r.element not in set(skip)
                 and (elements is None or r.element in set(elements)))


@dataclass(frozen=True)
class Sorption:
    """What sorbs in a system: exchange capacity (equivalents), hydrous ferric oxide (grams) and a solid with
    partition coefficients (kg of solid, and its Kd rows), all for the whole system (per kg of the starting
    water)."""
    exchange_eq: float = 0.0
    hfo_g: float = 0.0
    kd_solid_kg: float = 0.0
    kd: tuple[KdRow, ...] = ()

    @property
    def active(self) -> bool:
        return self.exchange_eq > 0 or self.hfo_g > 0 or (self.kd_solid_kg > 0 and bool(self.kd))


def _row(model, name: str):
    """(components vector, log K offset) of one mole of an aqueous species (or H2O): ln a = e.u - LN10 k."""
    name = ALIASES.get(name, name)
    if name in model.ib:
        e = np.zeros(model.nbf)
        e[model.ib[name]] = 1.0
        return e, 0.0
    if name in model.dep_names:
        k = model.dep_names.index(name)
        return np.asarray(model.nu[k], dtype=float), float(model.logk_dep[k])
    return None


class _Species:
    """The sorbed species of one sorbent at one temperature: ln (amount or fraction) = lnk + e.u + ... ."""

    def __init__(self, model, rxns: list[Reaction], t_c: float):
        names, rows, lnk, z, masters, lead = [], [], [], [], [], []
        for r in rxns:
            e, off = np.zeros(model.nbf), 0.0
            ok = True
            for s, c in r.reactants.items():
                got = _row(model, s)
                if got is None:
                    ok = False
                    break
                e += c * got[0]
                off += c * got[1]
            if not ok:
                continue
            names.append(r.species)
            ion = next((ALIASES.get(sp, sp) for sp, c in r.reactants.items() if c > 0 and sp not in ("H2O",)), None)
            lead.append(("b", model.ib[ion]) if ion in model.ib else
                        ("d", model.dep_names.index(ion)) if ion in model.dep_names else None)
            rows.append(e)
            lnk.append(LN10 * (r.log_k_at(t_c) - off))
            z.append(r.z)
            masters.append(r.master)
        self.names, self.masters = names, masters
        self.e = np.array(rows).reshape(-1, model.nbf)        # components bound per mole, d ln n / d u
        self.lnk = np.array(lnk)
        self.z = np.array(z, dtype=float)
        self.lead = lead                  # the aqueous ion of each reaction, for its activity coefficient

    def ln_gamma(self, lg_b_ln, lg_d_ln) -> np.ndarray:
        """ln of the activity coefficient of each species' aqueous ion (0 where there is none)."""
        out = np.zeros(len(self.names))
        for k, ref in enumerate(self.lead):
            if ref is not None:
                out[k] = lg_b_ln[ref[1]] if ref[0] == "b" else lg_d_ln[ref[1]]
        return out


class Exchanger:
    """Gaines-Thomas exchange; one unknown, ln a(X-).  As in PHREEQC (-gamma), an exchange species' activity is
    its equivalent fraction times the activity coefficient of its cation in the solution."""
    nvar = 1

    def __init__(self, model, cec_eq: float, t_c: float):
        self.cec = cec_eq
        self.sp = _Species(model, [r for r in reactions() if r.kind == "exchange"], t_c)
        self.names = list(self.sp.names)
        self.ln_g = np.zeros(len(self.names))

    def set_gammas(self, lg_b_ln, lg_d_ln) -> None:
        self.ln_g = self.sp.ln_gamma(lg_b_ln, lg_d_ln)

    def start(self, u: np.ndarray) -> np.ndarray:
        """ln a(X-) for which the equivalent fractions add up to one at ``u``."""
        a = self.sp.lnk + self.sp.e @ u - self.ln_g
        lo, hi = -200.0, 200.0
        for _ in range(200):
            v = 0.5 * (lo + hi)
            if np.exp(np.clip(a + self.sp.z * v, -700, 700)).sum() > 1.0:
                hi = v
            else:
                lo = v
        return np.array([0.5 * (lo + hi)])

    def evaluate(self, u, v, ionic):
        sp = self.sp
        beta = np.exp(np.clip(sp.lnk + sp.e @ u + sp.z * v[0] - self.ln_g, -700.0, 300.0))
        n = beta * self.cec / sp.z
        bound = n @ sp.e
        d_bound_du = sp.e.T @ (n[:, None] * sp.e)
        d_bound_dv = (sp.e.T @ (n * sp.z))[:, None]
        eq = np.array([beta.sum() - 1.0])
        d_eq_du = (beta @ sp.e)[None, :]
        d_eq_dv = np.array([[float(beta @ sp.z)]])
        return n, bound, d_bound_du, d_bound_dv, eq, d_eq_du, d_eq_dv

    def scale(self) -> float:
        return self.cec


class Surface:
    """Diffuse double layer on hydrous ferric oxide; unknowns ln n(SOH) per site type and y = F psi / RT."""

    def __init__(self, model, grams: float, t_c: float, p_bar: float | None):
        fe = grams / HFO_G_PER_MOL_FE
        self.sites = {"Hfo_sOH": fe * HFO_STRONG_PER_FE, "Hfo_wOH": fe * HFO_WEAK_PER_FE}
        self.area = grams * HFO_M2_PER_G
        rx = [r for r in reactions() if r.kind == "surface"]
        rx += [Reaction("surface", m, m, {}, 0.0, 0.0, None) for m in self.sites]     # the bare sites
        self.sp = _Species(model, rx, t_c)
        self.names = list(self.sp.names)
        self.types = list(self.sites)
        self.member = np.array([[1.0 if m == t else 0.0 for m in self.sp.masters] for t in self.types])
        self.nvar = len(self.types) + 1
        t_k = t_c + 273.15
        eps = epsilon(t_k, p_bar if p_bar else 1.0)
        # sigma = sqrt(8 RT eps eps0 c) sinh(y/2), c in mol/m3 (taken as 1000 I)
        self._gc = math.sqrt(8.0 * GAS_R * t_k * eps * EPS0 * 1000.0)

    def start(self, u: np.ndarray) -> np.ndarray:
        a = self.sp.lnk + self.sp.e @ u
        out = []
        for k, t in enumerate(self.types):
            mask = self.member[k] > 0
            top = float(np.max(a[mask]))
            out.append(math.log(self.sites[t]) - top - math.log(np.exp(a[mask] - top).sum()))
        return np.array(out + [0.0])

    def evaluate(self, u, v, ionic):
        sp = self.sp
        nt = len(self.types)
        s_of = self.member.T @ v[:nt]                        # ln n(SOH) of each species' site type
        y = v[nt]
        n = np.exp(np.clip(sp.lnk + sp.e @ u + s_of - sp.z * y, -700.0, 300.0))
        bound = n @ sp.e
        d_bound_du = sp.e.T @ (n[:, None] * sp.e)
        dn_dv = np.concatenate([(self.member * n).T, (-sp.z * n)[:, None]], axis=1)   # species x vars
        d_bound_dv = sp.e.T @ dn_dv
        eq = np.zeros(self.nvar)
        d_eq_du = np.zeros((self.nvar, len(u)))
        d_eq_dv = np.zeros((self.nvar, self.nvar))
        for k, t in enumerate(self.types):
            w = self.member[k] / self.sites[t]
            eq[k] = w @ n - 1.0
            d_eq_du[k] = (w * n) @ sp.e
            d_eq_dv[k] = w @ dn_dv
        q_scale = sum(self.sites.values())
        g = self.area / FARADAY * self._gc * math.sqrt(max(ionic, 1e-12))
        eq[nt] = (sp.z @ n - g * math.sinh(y / 2.0)) / q_scale
        d_eq_du[nt] = ((sp.z * n) @ sp.e) / q_scale
        d_eq_dv[nt] = (sp.z @ dn_dv) / q_scale
        d_eq_dv[nt, nt] -= g * 0.5 * math.cosh(y / 2.0) / q_scale
        return n, bound, d_bound_du, d_bound_dv, eq, d_eq_du, d_eq_dv

    def scale(self) -> float:
        return sum(self.sites.values())

    def set_gammas(self, lg_b_ln, lg_d_ln) -> None:
        pass

    def potential_mv(self, v: np.ndarray, t_c: float) -> float:
        return 1000.0 * float(v[-1]) * GAS_R * (t_c + 273.15) / FARADAY


def _oxidation(sp) -> float | None:
    """Oxidation state of the one element other than H and O in a species made of it, H and O (CrO4-2: 6)."""
    other = [(el, n) for el, n in sp.elements.items() if el not in ("H", "O")]
    if len(other) != 1:
        return None
    el, n = other[0]
    return (sp.charge - sp.elements.get("H", 0.0) + 2.0 * sp.elements.get("O", 0.0)) / n


class LinearSorber:
    """Partition coefficients: n(sorbed copy of species s) = Kd(s) M m(s), for every dissolved species s of a
    metal with a Kd; no unknowns of its own."""
    nvar = 0

    def __init__(self, model, rows: tuple[KdRow, ...], solid_kg: float):
        db = model.db
        self.solid_kg = solid_kg
        i_o2 = model.ib.get("O2(g)")
        by_el: dict[str, list[KdRow]] = {}
        for r in rows:
            by_el.setdefault(r.element, []).append(r)

        def pick(el, ox):
            cands = by_el[el]
            exact = [r for r in cands if r.oxidation is None or (ox is not None and abs(r.oxidation - ox) < 0.01)]
            if exact:
                return exact[0]
            return cands[0] if ox is None else min(cands, key=lambda r: abs(r.oxidation - ox))

        # each species once: (where, index) -> {element: (atoms of it, its Kd row)}
        found: dict[tuple[str, int], dict[str, tuple[float, KdRow]]] = {}
        for j, bname in enumerate(model.names):
            sp = db.species.get(bname)
            els = [el for el in (sp.elements if sp else ()) if el in by_el]
            if len(els) != 1:
                continue
            el = els[0]
            n_el = sp.elements[el]
            ox_b = _oxidation(sp)
            found.setdefault(("b", j), {})[el] = (n_el, pick(el, ox_b))
            for k in np.nonzero(np.abs(model.nu[:, j]) > 1e-12)[0]:
                c = float(model.nu[k, j])
                ox = None
                if ox_b is not None:          # O2 given off per atom: four electrons each
                    ox = ox_b + (4.0 * float(model.nu[k, i_o2]) / (c * n_el) if i_o2 is not None else 0.0)
                found.setdefault(("d", int(k)), {})[el] = (c * n_el, pick(el, ox))
        names, rows_e, lnk, lead, held = [], [], [], [], []
        self.used: dict[str, KdRow] = {}
        for (where, i), els in found.items():
            r = max((row for _, row in els.values()), key=lambda row: row.log_kd)   # two metals: the larger Kd
            if where == "b":
                name, e, off = model.names[i], np.eye(model.nbf)[i], 0.0
            else:
                name, e, off = model.dep_names[i], np.asarray(model.nu[i], dtype=float), float(model.logk_dep[i])
            names.append(name)
            rows_e.append(e)
            lnk.append(LN10 * (r.log_kd - off) + math.log(solid_kg))
            lead.append((where, i))
            held.append({el: atoms for el, (atoms, _) in els.items()})
            self.used[name] = r
        self.names = names
        self.sp = SimpleNamespace(e=np.array(rows_e).reshape(-1, model.nbf))
        self.lnk = np.array(lnk)
        self.lead, self.held = lead, held
        self._model_names = list(model.names)
        self.ln_g = np.zeros(len(names))

    def set_gammas(self, lg_b_ln, lg_d_ln) -> None:
        self.ln_g = np.array([lg_b_ln[i] if w == "b" else lg_d_ln[i] for w, i in self.lead])

    def start(self, u):
        return np.zeros(0)

    def evaluate(self, u, v, ionic):
        e = self.sp.e
        n = np.exp(np.clip(self.lnk + e @ u - self.ln_g, -700.0, 300.0))
        nbf = e.shape[1]
        return (n, n @ e, e.T @ (n[:, None] * e), np.zeros((nbf, 0)), np.zeros(0), np.zeros((0, nbf)),
                np.zeros((0, 0)))

    def scale(self) -> float:
        return 1.0

    def report(self, n) -> dict[str, float]:
        """Moles of each element held, as "Pb (Kd)"."""
        out: dict[str, float] = {}
        for held, a in zip(self.held, n):
            for el, atoms in held.items():
                out[f"{el} (Kd)"] = out.get(f"{el} (Kd)", 0.0) + float(a * atoms)
        return out

    def loading(self, molality: dict[str, float]) -> dict[str, float]:
        """Moles of each component held in equilibrium with a fluid of the given species molalities (the fluid
        unchanged, as the exchanger's starting load)."""
        tot = np.zeros(self.sp.e.shape[1])
        for name, e in zip(self.names, self.sp.e):
            m = molality.get(name, 0.0)
            if m > 0:
                tot += self.factor(name) * m * e
        return {nm: float(v) for nm, v in zip(self._model_names, tot) if abs(v) > 0}

    def factor(self, name: str) -> float:
        """Kd M of a species (litres of fluid whose dissolved load the solid holds)."""
        return 10.0 ** self.used[name].log_kd * self.solid_kg


def sorbents(sorption: Sorption | None, model, t_c: float, p_bar: float | None) -> list:
    if sorption is None or not sorption.active:
        return []
    out = []
    if sorption.exchange_eq > 0:
        out.append(Exchanger(model, sorption.exchange_eq, t_c))
    if sorption.hfo_g > 0:
        out.append(Surface(model, sorption.hfo_g, t_c, p_bar))
    if sorption.kd_solid_kg > 0 and sorption.kd:
        lin = LinearSorber(model, sorption.kd, sorption.kd_solid_kg)
        if lin.names:
            out.append(lin)
    return out


def exchanger_loading(model, u: np.ndarray, cec_eq: float, t_c: float, ln_gamma=None) -> dict[str, float]:
    """Moles of each component the exchanger holds in equilibrium with a fluid of ln activities ``u`` (the
    exchanger is loaded from the fluid without changing it, as PHREEQC's 'equilibrate with solution');
    ``ln_gamma`` maps aqueous species to their ln activity coefficients in that fluid."""
    ex = Exchanger(model, cec_eq, t_c)
    if ln_gamma:
        ex.ln_g = np.array([ln_gamma.get(_lead_name(model, ref), 0.0) if ref else 0.0 for ref in ex.sp.lead])
    n, bound, *_ = ex.evaluate(u, ex.start(u), 0.0)
    return {name: float(b) for name, b in zip(model.names, bound) if abs(b) > 0}


def _lead_name(model, ref) -> str:
    return model.names[ref[1]] if ref[0] == "b" else model.dep_names[ref[1]]

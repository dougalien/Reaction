"""Ion exchange and surface complexation, solved together with the aqueous system.

Exchange (Gaines-Thomas, as PHREEQC): sites X- of fixed capacity; each species M X_z has the activity of its
equivalent fraction, beta = K a_M a_X^z, and the betas add up to one.  The exchanger holds exactly its capacity
in cation charge, so a fluid it is equilibrated with stays neutral.

Surface complexation (Dzombak and Morel 1990, hydrous ferric oxide): strong and weak sites with mass balances,
species n_i = K_i a(reactants) n(SOH) exp(-z_i F psi / RT), and the Gouy-Chapman relation between the surface
charge and the potential psi (the diffuse layer is not given its own composition, as in PHREEQC's default).

Data: ``data/sorption.tsv`` (PHREEQC's phreeqc.dat).  The log K are for 25 C; the exchange constants that have an
enthalpy follow van 't Hoff, the rest are held at 25 C.
"""
from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from .thermo.electrostatics import epsilon

DATA = Path(__file__).parent / "data" / "sorption.tsv"
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


@dataclass(frozen=True)
class Sorption:
    """What sorbs in a system: exchange capacity (equivalents) and hydrous ferric oxide (grams), both for the
    whole system (per kg of the starting water)."""
    exchange_eq: float = 0.0
    hfo_g: float = 0.0

    @property
    def active(self) -> bool:
        return self.exchange_eq > 0 or self.hfo_g > 0


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


def sorbents(sorption: Sorption | None, model, t_c: float, p_bar: float | None) -> list:
    if sorption is None or not sorption.active:
        return []
    out = []
    if sorption.exchange_eq > 0:
        out.append(Exchanger(model, sorption.exchange_eq, t_c))
    if sorption.hfo_g > 0:
        out.append(Surface(model, sorption.hfo_g, t_c, p_bar))
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

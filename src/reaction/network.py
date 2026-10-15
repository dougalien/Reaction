"""Reaction network: every species and phase expressed over a chosen basis.

Each species of the database is written over the columns of the chosen basis
by recursive substitution of reactions

    log a_i = sum_j nu_ij log a_j - sum_t coef_t logK_t .

Three EQ3NR behaviours are reproduced:

* **Basis switch** (``NH3(aq)`` -> ``NO3-``): the switched-out strict species
  is defined through the reaction of its replacement, so every species whose
  reaction mentions it is re-expressed over the new basis.
* **Suppression** acts on the original reaction text: suppressing a species
  removes every species whose reaction mentions it (suppress NH3(aq) and NO2-,
  NH4+ ... disappear; suppress the organic master species Acetic_acid(aq) and
  all organics disappear).
* **Decoupled components**: an auxiliary species (O2(aq), Fe+++, NO2-) given
  its own total concentration becomes an extra basis column, independent of the
  default redox state.  Reactions that mention it stay written through it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .database.model import Kind, ThermoDatabase, ThermoRangeError

WATER, PROTON, OXYGEN = "H2O", "H+", "O2(g)"


class NetworkError(RuntimeError):
    pass


@dataclass
class Evaluated:
    """Network quantities at one temperature (and pressure)."""
    dep_index: np.ndarray            # indices into Network.dep_names that have data
    logk_dep: np.ndarray
    phase_index: np.ndarray
    logk_phase: np.ndarray
    logk_couple: dict[str, float]
    dropped: list[str]
    extrapolated: list[str]


class _Reducer:
    """Writes species over a fixed set of columns by recursive substitution."""

    def __init__(self, db: ThermoDatabase, columns, suppressed, derived=None):
        self.db, self.cols, self.suppressed = db, set(columns), set(suppressed)
        self.derived = dict(derived or {})        # switched-out species -> its replacement
        self._memo: dict[str, tuple[dict[str, float], list] | None] = {}

    def reduce(self, name: str):
        if name in self.suppressed:
            return None
        if name in self.cols:
            return {name: 1.0}, []
        if name in self._memo:
            return self._memo[name]
        self._memo[name] = None                      # guards against cycles
        result = self._derived(name) if name in self.derived else self._plain(name)
        self._memo[name] = result
        return result

    def _combine(self, parts):
        """Sum of alpha * (reduced species) as (nu, terms)."""
        nu: dict[str, float] = {}
        terms: list = []
        for alpha, other in parts:
            red = self.reduce(other)
            if red is None:
                return None
            nu_k, terms_k = red
            for j, v in nu_k.items():
                nu[j] = nu.get(j, 0.0) + alpha * v
            terms += [(alpha * cf, lk) for cf, lk in terms_k]
        return nu, terms

    def _plain(self, name: str):
        sp = self.db.species.get(name)
        if sp is None or sp.kind == Kind.BASIS or sp.logk is None or not sp.reaction:
            return None
        res = self._combine([(c, k) for k, c in sp.reaction.items() if k != name])
        if res is None:
            return None
        nu, terms = res
        return nu, [(1.0, sp.logk)] + terms

    def _derived(self, name: str):
        """Define a switched-out strict species through its replacement's reaction."""
        new = self.derived[name]
        sp = self.db[new]
        c_x = sp.reaction.get(name)
        if not c_x:
            raise NetworkError(f"cannot switch {name} -> {new}: the reaction of {new} "
                               f"does not involve {name}")
        parts = [(1.0 / c_x, new)] + [(-c / c_x, k) for k, c in sp.reaction.items()
                                       if k not in (name, new)]
        res = self._combine(parts)
        if res is None:
            return None
        nu, terms = res
        return nu, [(-1.0 / c_x, sp.logk)] + terms


class Network:
    def __init__(self, db: ThermoDatabase, basis: list[str], *,
                 switches: dict[str, str] | None = None,
                 decoupled: tuple[str, ...] = (),
                 suppress=(), suppress_phases=()):
        self.db = db
        switches = dict(switches or {})
        if WATER not in basis or PROTON not in basis:
            raise NetworkError("the basis must contain H2O and H+")
        suppress = set(suppress)
        targets = set(switches.values())
        for name in basis:
            if name not in db:
                raise NetworkError(f"unknown species {name!r}")
            if db[name].kind != Kind.BASIS and name not in targets:
                raise NetworkError(f"{name!r} is not a strict basis species in {db.title} "
                                   f"and is not the target of a basis switch")
        for old, new in switches.items():
            if old not in db or db[old].kind != Kind.BASIS or new not in basis:
                raise NetworkError(f"invalid basis switch {old} -> {new}")
        clash = (set(basis) | set(decoupled)) & suppress
        if clash:
            raise NetworkError(f"cannot suppress basis species {sorted(clash)}")

        self.basis = list(basis) + list(decoupled)
        nb = len(self.basis)
        index = {n: i for i, n in enumerate(self.basis)}
        red = _Reducer(db, self.basis, suppress, derived=switches)

        def vector(nu):
            v = np.zeros(nb)
            for col, c in nu.items():
                v[index[col]] = c
            return v

        # --- aqueous species ------------------------------------------------
        self.dep_names: list[str] = []
        self.dep_terms: list[list] = []
        dep_rows: list[np.ndarray] = []
        self.phase_names: list[str] = []
        self.phase_kind: list[Kind] = []
        self.phase_terms: list[list] = []
        phase_rows: list[np.ndarray] = []

        def add_aqueous(name):
            r = red.reduce(name)
            if r is None:
                return
            nu, terms = r
            if name.endswith("(g)"):                      # e.g. O2(g) switched out
                self.phase_names.append(name); self.phase_kind.append(Kind.GAS)
                phase_rows.append(vector(nu)); self.phase_terms.append(terms)
            else:
                self.dep_names.append(name); dep_rows.append(vector(nu))
                self.dep_terms.append(terms)

        for sp in db.of_kind(Kind.AUX_BASIS, Kind.AQUEOUS):
            if sp.name not in index:
                add_aqueous(sp.name)
        for old in switches:                               # switched-out strict species
            add_aqueous(old)

        # Suppressed phases stay in the network: EQ3NR still reports their saturation
        # states, it only keeps them out of equilibrium constraints.
        self.inactive_phases = set(suppress_phases)
        for sp in db.of_kind(Kind.MINERAL, Kind.GAS, Kind.LIQUID):
            if sp.name in index:
                continue
            r = red.reduce(sp.name)
            if r is not None:
                self.phase_names.append(sp.name); self.phase_kind.append(sp.kind)
                phase_rows.append(vector(r[0])); self.phase_terms.append(r[1])

        # --- decoupled components: their own definition, for redox couples ----
        self.couple: dict[str, tuple[np.ndarray, list]] = {}
        if decoupled:
            core = _Reducer(db, list(basis), suppress, derived=switches)
            for d in decoupled:
                r = core.reduce(d)
                if r is None:
                    raise NetworkError(f"cannot decouple {d!r}: it has no reaction over "
                                       f"the basis of this problem")
                self.couple[d] = (vector(r[0]), r[1])

        self.nu_dep = np.array(dep_rows).reshape(-1, nb)
        self.nu_phase = np.array(phase_rows).reshape(-1, nb)
        sp_of = [db[n] for n in self.basis]
        self.z_basis = np.array([s.charge for s in sp_of])
        self.a0_basis = np.array([s.a0 or 0.0 for s in sp_of])
        self.type_basis = np.array([s.activity_type for s in sp_of])
        self.is_solute = np.array([n != WATER and not n.endswith("(g)") for n in self.basis])
        dsp = [db[n] for n in self.dep_names]
        self.z_dep = np.array([s.charge for s in dsp])
        self.a0_dep = np.array([s.a0 or 0.0 for s in dsp])
        self.type_dep = np.array([s.activity_type for s in dsp])

    # ------------------------------------------------------------------
    def evaluate(self, t_c: float, p_bar: float | None = None) -> Evaluated:
        dropped: list[str] = []
        extrapolated: list[str] = []

        def value(terms, name):
            total, flag = 0.0, False
            for cf, lk in terms:
                total += cf * lk(t_c, p_bar)
                flag = flag or lk.extrapolated(t_c)
            if flag:
                extrapolated.append(name)
            return total

        def run(names, term_lists):
            idx, vals = [], []
            for i, (name, terms) in enumerate(zip(names, term_lists)):
                try:
                    vals.append(value(terms, name))
                    idx.append(i)
                except ThermoRangeError:
                    dropped.append(name)
            return np.array(idx, dtype=int), np.array(vals)

        di, dv = run(self.dep_names, self.dep_terms)
        pi, pv = run(self.phase_names, self.phase_terms)
        couples = {name: value(terms, name) for name, (_, terms) in self.couple.items()}
        return Evaluated(di, dv, pi, pv, couples, dropped, extrapolated)

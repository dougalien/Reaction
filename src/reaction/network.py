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
CARBONATE = "HCO3-"

# Carboxylic acid anions with metal complexes in SLOP98 (Shock and Koretsky 1993, 1995; Sverjensky, Shock and
# Helgeson 1997): ligand -> (the abbreviation in its complexes' names, its protonated forms).  Organic acids do
# not form from carbonate or break down to it on laboratory or vent time scales, so each is a component of its
# own, present only when an analysis gives it, with its acids and complexes written through it.
ORGANIC_LIGANDS: dict[str, tuple[str, tuple[str, ...]]] = {
    "formate": ("For", ("formic acid",)),
    "acetate": ("Ac", ("acetic acid",)),
    "propanoate": ("Prop", ("propanoic acid",)),
    "n-butanoate": ("But", ("n-butanoic acid",)),
    "n-pentanoate": ("Pent", ("n-pentanoic acid",)),
    "glycolate": ("Glyc", ("glycolic acid",)),
    "lactate": ("Lac", ("lactic acid",)),
    "oxalate-2": ("Oxal", ("H-oxalate", "oxalic acid")),
    "malonate-2": ("Mal", ("H-malonate", "malonic acid")),
    "succinate-2": ("Succ", ("H-succinate", "succinic acid")),
}


def ligand_member(ligand: str, name: str) -> bool:
    """Whether species ``name`` is ``ligand`` itself, one of its acids or one of its metal complexes."""
    abbr, acids = ORGANIC_LIGANDS[ligand]
    return name == ligand or name in acids or f"({abbr})" in name


def organic_held_apart(db, held=()) -> list[str]:
    """Species with reduced (organic) carbon other than methane that no held ligand in ``held`` accounts for:
    the organic acids, their complexes, cyanide ...  Without a ligand of their own they could only form by
    reducing carbonate (or methane by oxidizing it), which does not happen at these time scales; so they are
    left out, as methane is held apart."""
    ligands = [h for h in held if h in ORGANIC_LIGANDS]
    out = []
    for sp in db.species.values():
        r = sp.reaction or {}
        if sp.name in ("CH4", "CH4(g)") or CARBONATE not in r or abs(r.get(OXYGEN, 0.0)) < 1e-12:
            continue
        if any(ligand_member(lig, sp.name) for lig in ligands):
            continue
        out.append(sp.name)
    return out


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

        # an organic ligand in a basis without O2(g): its acids and complexes are written over carbonate
        # and O2 in the data, so go through O2 and take back the ligand's own reaction, which cancels it
        ligands = [d for d in decoupled if d in ORGANIC_LIGANDS]
        aux = (_Reducer(db, list(basis) + [OXYGEN], suppress, derived=switches)
               if ligands and OXYGEN not in index and CARBONATE in index else None)

        def via_ligand(name):
            for d in ligands:
                if not ligand_member(d, name):
                    continue
                own, member = aux.reduce(d), aux.reduce(name)
                if own is None or member is None or abs(own[0].get(CARBONATE, 0.0)) < 1e-12:
                    return None
                a = member[0].get(CARBONATE, 0.0) / own[0][CARBONATE]
                nu = {k: member[0].get(k, 0.0) - a * own[0].get(k, 0.0) for k in {*member[0], *own[0]}}
                if abs(nu.pop(OXYGEN, 0.0)) > 1e-9 or abs(nu.pop(CARBONATE, 0.0)) > 1e-9:
                    return None
                nu[d] = nu.get(d, 0.0) + a
                return ({k: v for k, v in nu.items() if abs(v) > 1e-12},
                        list(member[1]) + [(-a * cf, lk) for cf, lk in own[1]])
            return None

        def add_aqueous(name):
            r = red.reduce(name)
            if r is None and aux is not None:
                r = via_ligand(name)
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
                if r is None and d in ORGANIC_LIGANDS:
                    continue         # no O2 to write it over carbonate: its family is its own already
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


def tie_family(net: Network, db, species: str) -> None:
    """Rewrite every species and phase in the same oxidation state as a held-apart species over that
    species' own column, so that they exchange with it and not with the rest of the basis: CH4(g) with
    CH4; HS-, metal bisulfide complexes and pyrrhotite with H2S; Fe(III) complexes and hematite with
    Fe+3.  A row qualifies when it is a multiple of the held species' row plus a row with neither its
    element column nor O2(g), i.e. no electrons to trade.  A held species without such an element
    column (H2) takes only rows that are an exact multiple of its own (its gas).  An organic ligand
    (ORGANIC_LIGANDS) takes only its own acids and complexes, by name: lactate and acetate have carbon in
    the same oxidation state, but one is not made from the other."""
    if species not in net.couple:
        return
    own, terms = net.couple[species]
    col = net.basis.index(species)
    skip = {net.basis.index(n) for n in (WATER, PROTON, OXYGEN) if n in net.basis}
    keys = [j for j in range(len(net.basis)) if abs(own[j]) > 1e-12 and j not in skip]
    o2 = net.basis.index(OXYGEN) if OXYGEN in net.basis else None
    key = keys[0] if len(keys) == 1 else None

    def tie(name, row, row_terms):
        if species in ORGANIC_LIGANDS and not ligand_member(species, name):
            return None
        if key is not None:
            if abs(row[key]) < 1e-12:
                return None
            a = row[key] / own[key]
            rest = row - a * own
            if o2 is not None and abs(rest[o2]) > 1e-9:
                return None
        else:
            j = max(range(len(own)), key=lambda i: abs(own[i]))
            a = row[j] / own[j]
            rest = row - a * own
            if abs(a) < 1e-12 or np.max(np.abs(rest)) > 1e-9:
                return None
        new = rest.copy()
        new[col] += a
        return new, list(row_terms) + [(-a * cf, lk) for cf, lk in terms]

    for names, rows, term_lists in ((net.dep_names, net.nu_dep, net.dep_terms),
                                    (net.phase_names, net.nu_phase, net.phase_terms)):
        for k in range(len(rows)):
            if abs(rows[k][col]) > 1e-12:            # already written through a held column
                continue
            got = tie(names[k], rows[k], term_lists[k])
            if got is not None:
                rows[k], term_lists[k] = got

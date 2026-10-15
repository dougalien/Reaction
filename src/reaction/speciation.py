"""Aqueous speciation (the EQ3NR problem).

Given total concentrations of the basis species, temperature and a redox
state, solve the coupled mass-action / mass-balance system:

    log a_i = sum_j nu_ij log a_j - log K_i          (every dependent species i)
    T_j     = m_j + sum_i nu_ij m_i                  (every basis species j)

Newton iteration runs in ln(molality) with activity coefficients frozen; an
outer loop updates activity coefficients, water activity and (for Eh/pe
constraints) the oxygen fugacity until they stop changing.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import activity
from .constants import FARADAY, KELVIN, LN10, R, WATER_MOLES_PER_KG
from .database.model import Kind, ThermoDatabase, ThermoRangeError

WATER, PROTON, OXYGEN = "H2O", "H+", "O2(g)"


# --------------------------------------------------------------------------
# Constraints
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Molal:
    """Total molality (mol per kg H2O)."""
    value: float


@dataclass(frozen=True)
class MgPerKg:
    """mg of the species per kg of solution (EQ3NR "mg/kg.sol")."""
    value: float


@dataclass(frozen=True)
class PH:
    """Fix the H+ activity through pH on the chosen pH scale."""
    value: float


@dataclass(frozen=True)
class LogActivity:
    value: float


@dataclass(frozen=True)
class Eh:
    volts: float


@dataclass(frozen=True)
class Pe:
    value: float


@dataclass(frozen=True)
class LogFO2:
    value: float


class SpeciationError(RuntimeError):
    pass


# --------------------------------------------------------------------------
# Reaction network reduced to the chosen basis
# --------------------------------------------------------------------------
class Network:
    """All species/phases of ``db`` expressible in the chosen basis.

    Reactions that mention other non-basis species are rewritten recursively
    until only basis species remain, so that

        log a_i = sum_j nu_ij log a_j - sum_t coef_t * logK_t(T).
    """

    def __init__(self, db: ThermoDatabase, basis: list[str]):
        self.db = db
        for name in basis:
            if name not in db or db[name].kind != Kind.BASIS:
                raise SpeciationError(f"{name!r} is not a strict basis species in {db.title}")
        self.basis = basis
        self._index = {n: i for i, n in enumerate(basis)}
        self._memo: dict[str, tuple[dict[str, float], list] | None] = {}

        self.dep_names: list[str] = []
        self._dep_red: list = []
        for sp in db.of_kind(Kind.AUX_BASIS, Kind.AQUEOUS):
            red = self._reduce(sp.name)
            if red is not None:
                self.dep_names.append(sp.name)
                self._dep_red.append(red)
        self.phase_names: list[str] = []
        self._phase_red: list = []
        for sp in db.of_kind(Kind.MINERAL, Kind.GAS, Kind.LIQUID):
            red = self._reduce(sp.name)
            if red is not None:
                self.phase_names.append(sp.name)
                self._phase_red.append(red)

        nb = len(basis)
        self.nu = self._matrix(self._dep_red, nb)
        self.nu_phase = self._matrix(self._phase_red, nb)
        self.z_basis = np.array([db[n].charge for n in basis])
        self.z_dep = np.array([db[n].charge for n in self.dep_names])
        self.a0_basis = np.array([self._ion_size(n) for n in basis])
        self.a0_dep = np.array([self._ion_size(n) for n in self.dep_names])
        self.type_basis = np.array([db[n].activity_type for n in basis])
        self.type_dep = np.array([db[n].activity_type for n in self.dep_names])

    def _ion_size(self, name: str) -> float:
        a0 = self.db[name].a0
        return a0 if a0 else activity.DEFAULT_ION_SIZE

    def _reduce(self, name: str):
        if name in self._index:
            return {name: 1.0}, []
        if name in self._memo:
            return self._memo[name]
        self._memo[name] = None                      # guards against cycles
        sp = self.db.species.get(name)
        if sp is None or sp.kind == Kind.BASIS or sp.logk is None or not sp.reaction:
            return None
        nu: dict[str, float] = {}
        terms = [(1.0, sp.logk)]
        for other, c in sp.reaction.items():
            if other == name:
                continue
            red = self._reduce(other)
            if red is None:
                return None
            nu_k, terms_k = red
            for j, v in nu_k.items():
                nu[j] = nu.get(j, 0.0) + c * v
            terms += [(c * cf, lk) for cf, lk in terms_k]
        self._memo[name] = (nu, terms)
        return self._memo[name]

    def _matrix(self, reduced: list, nb: int) -> np.ndarray:
        m = np.zeros((len(reduced), nb))
        for r, (nu, _) in enumerate(reduced):
            for j, v in nu.items():
                m[r, self._index[j]] = v
        return m

    def logk_at(self, t_c: float, p_bar: float | None):
        """Evaluate effective logK at (T, P); species lacking data are dropped.

        Returns (keep_dep, logk_dep, keep_phase, logk_phase, dropped_names).
        """
        def evaluate(reduced, names):
            vals, keep, dropped = [], [], []
            for r, ((_, terms), name) in enumerate(zip(reduced, names)):
                try:
                    vals.append(sum(cf * lk(t_c, p_bar) for cf, lk in terms))
                    keep.append(r)
                except ThermoRangeError:
                    dropped.append(name)
            return np.array(keep, dtype=int), np.array(vals), dropped

        kd, vd, dd = evaluate(self._dep_red, self.dep_names)
        kp, vp, dp = evaluate(self._phase_red, self.phase_names)
        return kd, vd, kp, vp, dd + dp


# --------------------------------------------------------------------------
# Result
# --------------------------------------------------------------------------
@dataclass
class SpeciesState:
    charge: float
    molality: float
    log_gamma: float
    log_activity: float


@dataclass
class Speciation:
    temperature_c: float
    pressure_bar: float | None
    ionic_strength: float
    sum_molality: float
    log_a_water: float
    osmotic_coefficient: float
    ph: float
    ph_scale: str
    pmh: float
    pe: float
    eh_volts: float
    log_fo2: float
    charge_imbalance: float                  # eq / kg H2O (cations + anions)
    species: dict[str, SpeciesState]
    saturation: dict[str, float]             # minerals/liquids: log Q/K
    log_fugacity: dict[str, float]           # gases: log f at equilibrium with the fluid
    dropped: list[str]
    iterations: int

    def molality(self, name: str) -> float:
        return self.species[name].molality

    def table(self, min_log_molality: float = -12.0) -> str:
        rows = sorted(self.species.items(), key=lambda kv: -kv[1].molality)
        out = [f"T = {self.temperature_c:g} C   I = {self.ionic_strength:.5f} m   "
               f"pH({self.ph_scale}) = {self.ph:.4f}   Eh = {self.eh_volts:.4f} V   "
               f"log fO2 = {self.log_fo2:.3f}",
               f"{'species':<22}{'molality':>12}{'log m':>9}{'log gamma':>11}{'log a':>9}"]
        for name, s in rows:
            lm = math.log10(s.molality) if s.molality > 0 else -999
            if lm < min_log_molality:
                break
            out.append(f"{name:<22}{s.molality:12.4e}{lm:9.4f}{s.log_gamma:11.4f}{s.log_activity:9.4f}")
        return "\n".join(out)


# --------------------------------------------------------------------------
# Solver
# --------------------------------------------------------------------------
def speciate(
    db: ThermoDatabase,
    temperature_c: float,
    basis: dict,
    redox=None,
    *,
    pressure_bar: float | None = None,
    activity_model: str = "bdot",
    ph_scale: str = "nbs",
    tds_mg_per_kg: float = 0.0,
    tol: float = 1e-12,
    max_outer: int = 200,
    max_inner: int = 100,
) -> Speciation:
    """Solve the speciation of a fluid.

    ``basis`` maps strict basis-species names to a constraint (:class:`Molal`,
    :class:`MgPerKg`, :class:`PH`, :class:`LogActivity`).  ``redox`` is an
    :class:`Eh`, :class:`Pe` or :class:`LogFO2`; without it, species that
    depend on oxygen fugacity are left out of the model.

    ``tds_mg_per_kg`` converts :class:`MgPerKg` inputs to molality the way
    EQ3NR does: m = (c / 1000 / MW) / (1 - TDS/1e6).  EQ3NR's default of 0
    means "treat mg/kg.sol as mg per kg of water".
    """
    if activity_model not in ("bdot", "ideal"):
        raise ValueError("activity_model must be 'bdot' or 'ideal'")
    if ph_scale not in ("nbs", "internal", "mesmer"):
        raise ValueError("ph_scale must be 'nbs', 'internal' or 'mesmer'")

    t_k = temperature_c + KELVIN
    names = [WATER, PROTON] + [n for n in basis if n not in (WATER, PROTON, OXYGEN)]
    if redox is not None:
        names.append(OXYGEN)
    net = Network(db, names)
    kd, logk_dep, kp, logk_phase, dropped = net.logk_at(temperature_c, pressure_bar)
    nu = net.nu[kd]
    z_dep, a0_dep, type_dep = net.z_dep[kd], net.a0_dep[kd], net.type_dep[kd]
    dep_names = [net.dep_names[i] for i in kd]
    nu_phase = net.nu_phase[kp]
    phase_names = [net.phase_names[i] for i in kp]
    ib = {n: i for i, n in enumerate(names)}
    nb = len(names)

    # ---- classify constraints --------------------------------------------
    A, B, bdot = db.dh_params(temperature_c)
    unknown: list[str] = []
    totals: list[float] = []
    fixed_log_a: dict[str, float] = {}
    ph_value = None
    for name, c in basis.items():
        if name in (WATER, OXYGEN):
            continue
        if isinstance(c, Molal):
            unknown.append(name); totals.append(c.value)
        elif isinstance(c, MgPerKg):
            mw = db.molecular_weight(name)
            unknown.append(name)
            totals.append(c.value / 1000.0 / mw / (1.0 - tds_mg_per_kg / 1e6))
        elif isinstance(c, PH):
            if name != PROTON:
                raise SpeciationError("PH applies to H+ only")
            ph_value = c.value
        elif isinstance(c, LogActivity):
            fixed_log_a[name] = c.value * LN10
        else:
            raise SpeciationError(f"unsupported constraint {c!r} for {name}")
    if PROTON not in basis:
        raise SpeciationError("H+ needs a constraint (PH, LogActivity or Molal)")
    if any(t <= 0 for t in totals):
        raise SpeciationError("total concentrations must be positive")
    iu = np.array([ib[n] for n in unknown], dtype=int)
    nu_u = nu[:, iu]
    totals_arr = np.array(totals)
    x = np.log(totals_arr)

    # ---- redox ------------------------------------------------------------
    logk_eh = db.logk_eh(temperature_c) if redox is not None else None
    if isinstance(redox, Eh):
        pe_fixed = redox.volts * FARADAY / (R * t_k * LN10)
    elif isinstance(redox, Pe):
        pe_fixed = redox.value
    else:
        pe_fixed = None

    # ---- activity-coefficient bookkeeping --------------------------------
    solute_basis = [i for i, n in enumerate(names) if n not in (WATER, OXYGEN)]
    z_b, a0_b, type_b = net.z_basis, net.a0_basis, net.type_basis
    cl_a0 = net._ion_size("Cl-") if "Cl-" in db else 3.0

    def log_gammas(ionic: float):
        lg_b = np.zeros(nb)
        lg_d = np.zeros(len(dep_names))
        if activity_model == "ideal":
            return lg_b, lg_d
        def one(z, a0, typ):
            if z != 0:
                return activity.log_gamma_bdot(z, a0, ionic, A, B, bdot)
            if typ == -1:
                return activity.log_gamma_nonpolar(ionic, temperature_c, db.cco2)
            return 0.0
        for i in solute_basis:
            lg_b[i] = one(z_b[i], a0_b[i], type_b[i])
        for k in range(len(dep_names)):
            lg_d[k] = one(z_dep[k], a0_dep[k], type_dep[k])
        if ph_scale == "nbs":
            f = (activity.log_gamma_chloride_bates_guggenheim(ionic, A)
                 - activity.log_gamma_bdot(-1.0, cl_a0, ionic, A, B, bdot))
            lg_b -= z_b * f
            lg_d -= z_dep * f
        return lg_b, lg_d

    # ---- iterate -----------------------------------------------------------
    ionic = 0.0
    sum_m = 0.0
    log_aw = 0.0
    lg_b, lg_d = log_gammas(ionic)
    log_fo2 = 0.0
    total_inner = 0
    converged = False

    def ln_activities(x_vec, lg_b_, log_aw_, log_fo2_):
        ln_a = np.zeros(nb)
        ln_a[ib[WATER]] = log_aw_ * LN10
        if OXYGEN in ib:
            ln_a[ib[OXYGEN]] = log_fo2_ * LN10
        ln_a[iu] = x_vec + lg_b_[iu] * LN10
        for name, v in fixed_log_a.items():
            ln_a[ib[name]] = v
        if ph_value is not None:
            if ph_scale == "mesmer":
                ln_a[ib[PROTON]] = -ph_value * LN10 + lg_b_[ib[PROTON]] * LN10
            else:
                ln_a[ib[PROTON]] = -ph_value * LN10
        return ln_a

    for outer in range(max_outer):
        lg_b_ln = lg_b * LN10
        lg_d_ln = lg_d * LN10

        for inner in range(max_inner):
            total_inner += 1
            ln_a = ln_activities(x, lg_b, log_aw, log_fo2)
            ln_md = np.clip(nu @ ln_a - LN10 * logk_dep - lg_d_ln, -700.0, 700.0)
            md = np.exp(ln_md)
            mu = np.exp(x)
            tot = mu + nu_u.T @ md
            scale = mu + np.abs(nu_u).T @ md + totals_arr
            resid = (tot - totals_arr) / scale
            if np.max(np.abs(resid)) < tol:
                break
            jac = (np.diag(mu) + (nu_u * md[:, None]).T @ nu_u) / scale[:, None]
            try:
                dx = np.linalg.solve(jac, -resid)
            except np.linalg.LinAlgError as exc:
                raise SpeciationError("singular Jacobian") from exc
            big = np.max(np.abs(dx))
            if big > 2.3:
                dx *= 2.3 / big
            # backtrack if the step makes the residual worse
            norm0 = np.max(np.abs(resid))
            lam = 1.0
            for _ in range(12):
                x_try = x + lam * dx
                ln_a_t = ln_activities(x_try, lg_b, log_aw, log_fo2)
                md_t = np.exp(np.clip(nu @ ln_a_t - LN10 * logk_dep - lg_d_ln, -700, 700))
                mu_t = np.exp(x_try)
                r_t = (mu_t + nu_u.T @ md_t - totals_arr) / (
                    mu_t + np.abs(nu_u).T @ md_t + totals_arr)
                if np.max(np.abs(r_t)) < norm0 or lam < 1e-3:
                    break
                lam *= 0.5
            x = x + lam * dx
        else:
            raise SpeciationError("inner Newton iteration did not converge")

        # ---- update activity coefficients, water activity, fO2 ---------
        ln_a = ln_activities(x, lg_b, log_aw, log_fo2)
        ln_md = np.clip(nu @ ln_a - LN10 * logk_dep - lg_d_ln, -700.0, 700.0)
        md = np.exp(ln_md)
        m_basis = np.zeros(nb)
        m_basis[iu] = np.exp(x)
        ph_species_m = np.exp(ln_a[ib[PROTON]] - lg_b_ln[ib[PROTON]])
        if PROTON not in unknown:
            m_basis[ib[PROTON]] = ph_species_m
        for name, v in fixed_log_a.items():
            m_basis[ib[name]] = math.exp(v - lg_b_ln[ib[name]])
        ionic_new = 0.5 * (np.sum(m_basis * z_b ** 2) + np.sum(md * z_dep ** 2))
        sum_m_new = np.sum(m_basis[solute_basis]) + np.sum(md)
        lg_b_new, lg_d_new = log_gammas(ionic_new)
        if activity_model == "bdot":
            log_aw_new = activity.log_activity_water_bdot(sum_m_new, ionic_new, A, B, bdot)
        else:
            log_aw_new = 0.0
        fo2_new = log_fo2
        if redox is not None:
            if isinstance(redox, LogFO2):
                fo2_new = redox.value
            else:
                log_ah = ln_a[ib[PROTON]] / LN10
                fo2_new = logk_eh + 4.0 * pe_fixed - 4.0 * log_ah + 2.0 * log_aw_new
        delta = max(abs(ionic_new - ionic), np.max(np.abs(lg_b_new - lg_b)),
                    np.max(np.abs(lg_d_new - lg_d)) if len(lg_d) else 0.0,
                    abs(log_aw_new - log_aw), abs(fo2_new - log_fo2))
        ionic, sum_m, log_aw = ionic_new, sum_m_new, log_aw_new
        lg_b, lg_d, log_fo2 = lg_b_new, lg_d_new, fo2_new
        if delta < 1e-11 and outer > 0:
            converged = True
            break
    if not converged:
        raise SpeciationError("activity-coefficient iteration did not converge")

    # ---- final state ----------------------------------------------------
    lg_b_ln, lg_d_ln = lg_b * LN10, lg_d * LN10
    ln_a = ln_activities(x, lg_b, log_aw, log_fo2)
    md = np.exp(np.clip(nu @ ln_a - LN10 * logk_dep - lg_d_ln, -700.0, 700.0))
    m_basis = np.zeros(nb)
    m_basis[iu] = np.exp(x)
    if PROTON not in unknown:
        m_basis[ib[PROTON]] = math.exp(ln_a[ib[PROTON]] - lg_b_ln[ib[PROTON]])
    for name, v in fixed_log_a.items():
        m_basis[ib[name]] = math.exp(v - lg_b_ln[ib[name]])

    states: dict[str, SpeciesState] = {}
    for i in solute_basis:
        states[names[i]] = SpeciesState(z_b[i], m_basis[i], lg_b[i], ln_a[i] / LN10)
    ln_a_dep = np.log(np.maximum(md, 1e-300)) + lg_d_ln
    for k, n in enumerate(dep_names):
        states[n] = SpeciesState(z_dep[k], md[k], lg_d[k], ln_a_dep[k] / LN10)

    log_q = (nu_phase @ ln_a) / LN10 - logk_phase
    saturation, fugacity = {}, {}
    for n, v in zip(phase_names, log_q):
        (fugacity if db[n].kind == Kind.GAS else saturation)[n] = float(v)

    log_ah = ln_a[ib[PROTON]] / LN10
    m_h = m_basis[ib[PROTON]]
    if redox is not None:
        # inverse of  log fO2 = logK_eh + 4 pe - 4 log aH + 2 log aw
        pe = (log_fo2 - logk_eh + 4.0 * log_ah - 2.0 * log_aw) / 4.0
        eh = pe * R * t_k * LN10 / FARADAY
    else:
        pe = eh = float("nan")
    cations = float(np.sum(m_basis * z_b * (z_b > 0)) + np.sum(md * z_dep * (z_dep > 0)))
    anions = float(np.sum(m_basis * z_b * (z_b < 0)) + np.sum(md * z_dep * (z_dep < 0)))
    phi = -log_aw * LN10 * WATER_MOLES_PER_KG / sum_m if sum_m > 0 else 1.0
    return Speciation(
        temperature_c=temperature_c, pressure_bar=pressure_bar, ionic_strength=ionic,
        sum_molality=sum_m, log_a_water=log_aw, osmotic_coefficient=phi,
        ph=-log_ah if ph_scale != "mesmer" else -math.log10(m_h), ph_scale=ph_scale,
        pmh=-math.log10(m_h), pe=pe, eh_volts=eh, log_fo2=log_fo2,
        charge_imbalance=cations + anions, species=states, saturation=saturation,
        log_fugacity=fugacity, dropped=dropped, iterations=total_inner,
    )

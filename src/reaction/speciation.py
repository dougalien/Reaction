"""Aqueous speciation (the EQ3NR problem).

Unknowns are the natural-log activities ``u_j`` of the basis species.  Every
unknown gets exactly one equation, chosen by the constraint the user supplies:

    mass balance        sum_i nu_ij m_i = T_j               (Molal, MgPerKg, MgPerL)
    fixed activity      u_j = const                         (PH, LogActivity, LogFO2)
    phase equilibrium   log Q/K of a mineral/gas = target   (Equilibrium)
    charge balance      sum_i z_i m_i = 0                   (balance_on=...)
    redox               Eh / pe / redox couple

Dependent species follow from ``ln m_i = sum_j nu_ij u_j - ln10 logK_i - ln gamma_i``.
Newton iteration runs with activity coefficients frozen; an outer loop updates
activity coefficients and water activity until they stop changing.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import activity
from .constants import FARADAY, KELVIN, LN10, R, WATER_MOLES_PER_KG
from .database.model import Kind, ThermoDatabase
from .aqueous import AqueousModel
from .network import OXYGEN, PROTON, WATER, Network, NetworkError


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
class MgPerL:
    """mg of the species per litre of solution (needs ``density``)."""
    value: float


@dataclass(frozen=True)
class PH:
    """Fix the H+ activity through pH on the chosen pH scale."""
    value: float


@dataclass(frozen=True)
class LogActivity:
    value: float


@dataclass(frozen=True)
class Equilibrium:
    """Equilibrium with a mineral/liquid (log Q/K = ``value``, usually 0) or a
    gas (log fugacity = ``value``).  Applied to the basis species it determines."""
    phase: str
    value: float = 0.0


@dataclass(frozen=True)
class Eh:
    volts: float


@dataclass(frozen=True)
class Pe:
    value: float


@dataclass(frozen=True)
class LogFO2:
    value: float


@dataclass(frozen=True)
class Couple:
    """Redox state from a redox-couple species whose total is also given."""
    species: str


MASS_CONSTRAINTS = (Molal, MgPerKg, MgPerL)


class SpeciationError(RuntimeError):
    pass


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
    totals: dict[str, float]                 # total molality of each basis species
    saturation: dict[str, float]             # minerals/liquids: log Q/K
    log_fugacity: dict[str, float]           # gases: log f at equilibrium with the fluid
    dropped: list[str]                       # no logK data at this temperature
    extrapolated: list[str]                  # logK outside the span of valid grid points
    iterations: int
    warnings: list[str] = field(default_factory=list)
    # moles of every basis component per kg of solvent (H2O and O2(g) included), the
    # bulk composition a closed-system calculation starts from
    component_totals: dict[str, float] = field(default_factory=dict)
    basis_names: list[str] = field(default_factory=list)
    ln_activities: object = None            # ln a of the basis species (index 0: water)

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
def _globalized_step(u, unk, res, jac, step_cap, residual_only):
    """One damped Newton step on ``u[unk]`` with a sum-of-squares merit function.

    The Newton direction is shortened (never distorted) so no component moves
    farther than ``step_cap``; if backtracking along it cannot reduce the merit
    function, a Levenberg-Marquardt step takes over.
    """
    phi0 = 0.5 * float(res @ res)

    def limited(d):
        ratio = float(np.max(np.abs(d) / step_cap))
        return d / ratio if ratio > 1.0 else d

    def trial_of(step):
        trial = u.copy()
        trial[unk] += step
        r = residual_only(trial)
        return 0.5 * float(r @ r), trial

    try:
        newton = limited(np.linalg.solve(jac, -res))
    except np.linalg.LinAlgError:
        newton = None
    if newton is not None:
        if phi0 < 1e-16:                      # close to the root: take the full step
            return trial_of(newton)[1]
        lam = 1.0
        for _ in range(12):
            phi, trial = trial_of(lam * newton)
            if phi < phi0 * (1.0 - 1e-4 * lam):
                return trial
            lam *= 0.5
    jtj, grad = jac.T @ jac, jac.T @ res
    diag = np.diag(np.maximum(np.diag(jtj), 1e-12))
    mu = 1e-3
    for _ in range(25):
        try:
            step = limited(np.linalg.solve(jtj + mu * diag, -grad))
        except np.linalg.LinAlgError:
            mu *= 10.0
            continue
        phi, trial = trial_of(step)
        if phi < phi0:
            return trial
        mu *= 10.0
    raise SpeciationError("Newton line search failed to reduce the residual")


def make_speciation(model: AqueousModel, u, log_aw: float, lg_b, lg_d, ionic: float,
                    sum_m: float, iterations: int) -> Speciation:
    """Assemble a :class:`Speciation` from a solved state (ln activities ``u``)."""
    names, nbf, ib = model.names, model.nbf, model.ib
    nu, dep_names, logk_dep = model.nu, model.dep_names, model.logk_dep
    z_dep, z_b, solute = model.z_dep, model.z_b, model.solute
    t_k = model.temperature_c + KELVIN
    u = np.array(u, dtype=float)
    lg_b_ln, lg_d_ln = lg_b * LN10, lg_d * LN10
    u[0] = log_aw * LN10
    m_dep, m_b = model.state(u, lg_b_ln, lg_d_ln)
    states: dict[str, SpeciesState] = {}
    totals: dict[str, float] = {}
    for i in range(1, nbf):
        if solute[i]:
            states[names[i]] = SpeciesState(z_b[i], m_b[i], lg_b[i], u[i] / LN10)
            totals[names[i]] = float(m_b[i] + nu[:, i] @ m_dep)
    u_dep = nu @ u - LN10 * logk_dep                  # ln activity of dependent species
    for k, n in enumerate(dep_names):
        states[n] = SpeciesState(z_dep[k], m_dep[k], lg_d[k], u_dep[k] / LN10)

    log_q = (model.nu_ph @ u) / LN10 - model.logk_ph
    saturation, fugacity = {}, {}
    for n, kind, v in zip(model.ph_names, model.ph_kind, log_q):
        (fugacity if kind == Kind.GAS else saturation)[n] = float(v)

    log_ah = u[ib[PROTON]] / LN10
    m_h = m_b[ib[PROTON]]
    if OXYGEN in ib:
        log_fo2 = u[ib[OXYGEN]] / LN10
        pe = (log_fo2 - model.logk_eh + 4.0 * log_ah - 2.0 * log_aw) / 4.0
        eh = pe * R * t_k * LN10 / FARADAY
    else:
        log_fo2 = pe = eh = float("nan")
    cations = float(np.sum(m_b * z_b * (z_b > 0)) + np.sum(m_dep * z_dep * (z_dep > 0)))
    anions = float(np.sum(m_b * z_b * (z_b < 0)) + np.sum(m_dep * z_dep * (z_dep < 0)))
    phi = -log_aw * LN10 * WATER_MOLES_PER_KG / sum_m if sum_m > 0 else 1.0
    components = {WATER: WATER_MOLES_PER_KG + float(nu[:, 0] @ m_dep)}
    for j in range(1, nbf):
        components[names[j]] = float((m_b[j] if solute[j] else 0.0) + nu[:, j] @ m_dep)
    ph_scale = model.ph_scale
    return Speciation(
        temperature_c=model.temperature_c, pressure_bar=model.pressure_bar,
        ionic_strength=float(ionic), sum_molality=float(sum_m), log_a_water=float(log_aw),
        osmotic_coefficient=float(phi),
        ph=float(-log_ah if ph_scale != "mesmer" else -math.log10(m_h)), ph_scale=ph_scale,
        pmh=float(-math.log10(m_h)), pe=float(pe), eh_volts=float(eh), log_fo2=float(log_fo2),
        charge_imbalance=cations + anions, species=states, totals=totals,
        saturation=saturation, log_fugacity=fugacity, dropped=model.ev.dropped,
        extrapolated=model.ev.extrapolated, iterations=iterations, warnings=list(model.notes),
        component_totals=components, basis_names=list(names), ln_activities=u.copy(),
    )


def speciate(
    db: ThermoDatabase,
    temperature_c: float,
    basis: dict,
    redox=None,
    *,
    pressure_bar: float | None = None,
    switches: dict[str, str] | None = None,
    suppress=(),
    suppress_phases=(),
    balance_on: str | None = None,
    activity_model: str = "bdot",
    ph_scale: str = "nbs",
    density: float = 1.0,
    tds_mg_per_kg: float = 0.0,
    tol: float = 1e-12,
    max_outer: int = 200,
    max_inner: int = 200,
    verbose: bool = False,
) -> Speciation:
    """Solve the speciation of a fluid.

    ``basis`` maps species names to constraints.  Keys are strict basis
    species, targets of ``switches`` ({old strict species: new basis species}),
    or auxiliary species given their own total concentration, which become
    independent ("decoupled") components, as in EQ3NR.  ``redox`` is an
    :class:`Eh`, :class:`Pe`, :class:`LogFO2` or :class:`Couple` constraint on
    O2(g); without it, species that depend on oxygen fugacity are omitted.

    ``tds_mg_per_kg`` and ``density`` convert :class:`MgPerKg` / :class:`MgPerL`
    to molality as EQ3NR does: m = (c / 1000 / MW) / (1 - TDS/1e6), with
    mg/L divided by the density first.  ``balance_on`` replaces that species'
    constraint by a charge-balance equation.
    """
    if activity_model not in ("bdot", "davies", "ideal", "pitzer"):
        raise ValueError("activity_model must be 'bdot', 'davies', 'pitzer' or 'ideal'")
    if ph_scale not in ("nbs", "internal", "mesmer"):
        raise ValueError("ph_scale must be 'nbs', 'internal' or 'mesmer'")
    switches = dict(switches or {})
    t_k = temperature_c + KELVIN

    rows = dict(basis)
    if redox is not None:
        if OXYGEN in rows:
            raise SpeciationError("redox given twice (redox= and an O2(g) entry)")
        rows[OXYGEN] = redox
    if PROTON not in rows:
        raise SpeciationError("H+ needs a constraint (PH, LogActivity, Molal ... or balance_on)")

    # ---- decide the basis ------------------------------------------------
    targets = set(switches.values())
    core, decoupled = [WATER, PROTON], []
    for name, c in rows.items():
        if name in (WATER, PROTON, OXYGEN):
            continue
        if name not in db:
            raise SpeciationError(f"unknown species {name!r}")
        if db[name].kind == Kind.BASIS or name in targets:
            core.append(name)
        elif isinstance(c, MASS_CONSTRAINTS):
            decoupled.append(name)
        else:
            raise SpeciationError(f"{name!r} is not a basis species; only a total "
                                  f"concentration can be given for it")
    for old, new in switches.items():
        if new not in rows:
            raise SpeciationError(f"switch {old} -> {new}: no constraint given for {new}")
    if OXYGEN in rows:
        core.append(OXYGEN)
    for old in switches:
        if old in rows:
            raise SpeciationError(f"{old} was switched out of the basis; constrain "
                                  f"{switches[old]} instead")
    try:
        net = Network(db, core, switches=switches, decoupled=tuple(decoupled),
                      suppress=suppress, suppress_phases=suppress_phases)
    except NetworkError as exc:
        raise SpeciationError(str(exc)) from exc
    if pressure_bar is None:
        pressure_bar = db.reference_pressure(temperature_c)
    model = AqueousModel(db, net, temperature_c, pressure_bar, activity_model, ph_scale)
    ev = model.ev
    names, nbf, ib = model.names, model.nbf, model.ib
    nu, dep_names, logk_dep = model.nu, model.dep_names, model.logk_dep
    z_dep, type_dep = model.z_dep, model.type_dep
    nu_ph, ph_names, ph_kind, logk_ph = model.nu_ph, model.ph_names, model.ph_kind, model.logk_ph
    notes = list(model.notes)
    logk_eh = model.logk_eh
    unk = list(range(1, nbf))                    # everything but H2O
    pos = {j: r for r, j in enumerate(unk)}
    nu_unk = nu[:, unk]

    # ---- one equation per unknown ----------------------------------------
    solvent_fraction = 1.0 - tds_mg_per_kg / 1e6
    eqs: list[tuple] = []
    u = np.zeros(nbf)
    for j in unk:
        name = names[j]
        if name not in rows:
            raise SpeciationError(f"no constraint given for {name!r}")
        c = rows[name]
        if isinstance(c, (Molal, MgPerKg, MgPerL)):
            if isinstance(c, Molal):
                total = c.value
            else:
                grams = c.value / 1000.0 / db.molecular_weight(name)
                total = grams / (density if isinstance(c, MgPerL) else 1.0) / solvent_fraction
            if total <= 0:
                raise SpeciationError(f"total concentration of {name} must be positive")
            eqs.append(("mass", total)); u[j] = math.log(total)
        elif isinstance(c, PH):
            if name != PROTON:
                raise SpeciationError("PH applies to H+ only")
            eqs.append(("ph", c.value)); u[j] = -c.value * LN10
        elif isinstance(c, LogActivity):
            eqs.append(("loga", c.value * LN10)); u[j] = c.value * LN10
        elif isinstance(c, LogFO2):
            if name != OXYGEN:
                raise SpeciationError("LogFO2 applies to O2(g) only")
            eqs.append(("loga", c.value * LN10)); u[j] = c.value * LN10
        elif isinstance(c, (Eh, Pe)):
            if name != OXYGEN:
                raise SpeciationError("Eh/Pe apply to O2(g) only")
            pe = c.value if isinstance(c, Pe) else c.volts * FARADAY / (R * t_k * LN10)
            eqs.append(("eh", pe)); u[j] = -30.0 * LN10
        elif isinstance(c, Couple):
            if name != OXYGEN or c.species not in net.couple:
                raise SpeciationError(f"Couple({c.species!r}) needs the species to carry its "
                                      f"own total concentration")
            eqs.append(("couple", c.species)); u[j] = -30.0 * LN10
        elif isinstance(c, Equilibrium):
            if c.phase not in ph_names:
                raise SpeciationError(f"phase {c.phase!r} is unavailable for this problem")
            if c.phase in net.inactive_phases:
                raise SpeciationError(f"phase {c.phase!r} is suppressed; it cannot be used "
                                      f"as an equilibrium constraint")
            eqs.append(("phase", ph_names.index(c.phase), c.value))
            u[j] = -6.0 * LN10 if name != OXYGEN else -30.0 * LN10
        else:
            raise SpeciationError(f"unsupported constraint {c!r} for {name}")
    start_eqs = list(eqs)                        # totals still guide the pre-Newton phase
    if balance_on is not None:
        if balance_on not in ib or balance_on == WATER:
            raise SpeciationError(f"cannot balance charge on {balance_on!r}")
        r = pos[ib[balance_on]]
        eqs[r] = ("charge",)
    for r, j in enumerate(unk):                  # start redox at the value its equation implies
        if eqs[r][0] == "eh":
            u[j] = LN10 * (logk_eh + 4.0 * eqs[r][1]) - 4.0 * u[ib[PROTON]]
    # per-component Newton step limit (ln units): generous for equations linear in u
    step_cap = np.array([{"mass": 2.3, "charge": 4.6}.get(eqs[r][0], 40.0)
                         for r in range(len(unk))])
    if any(e[0] == "couple" for e in eqs):
        couple_species = next(e[1] for e in eqs if e[0] == "couple")
        couple_row, couple_logk = net.couple[couple_species][0], ev.logk_couple[couple_species]
        couple_col = ib[couple_species]

    # ---- activity model (shared with the equilibrium solver) -----------------------------
    solute, z_b = model.solute, model.z_b
    log_gammas, log_a_water, state = model.log_gammas, model.log_a_water, model.state

    def residual(uvec, lg_b_ln, lg_d_ln, log_aw, want_jac=True):
        m_dep, m_b = state(uvec, lg_b_ln, lg_d_ln)
        n = len(unk)
        res = np.zeros(n)
        jac = np.zeros((n, n)) if want_jac else None
        for r, j in enumerate(unk):
            kind = eqs[r][0]
            if kind == "mass":
                # ln(total/T) stays well behaved when polynuclear species explode
                col = nu[:, j]
                total_ = eqs[r][1]
                tot = m_b[j] + col @ m_dep
                if tot > 0:
                    res[r] = math.log(tot / total_)
                    if want_jac:
                        jac[r] = nu_unk.T @ (col * m_dep) / tot
                        jac[r, r] += m_b[j] / tot
                else:                                   # net-negative total (e.g. H+): linear form
                    scale = m_b[j] + np.abs(col) @ m_dep + total_
                    res[r] = (tot - total_) / scale
                    if want_jac:
                        jac[r] = nu_unk.T @ (col * m_dep) / scale
                        jac[r, r] += m_b[j] / scale
            elif kind == "charge":
                # ln(positive charge / negative charge): zero at balance, no saturation
                wp_d, wn_d = np.where(z_dep > 0, z_dep, 0.0) * m_dep, np.where(z_dep < 0, -z_dep, 0.0) * m_dep
                wp_b, wn_b = np.where(z_b > 0, z_b, 0.0) * m_b, np.where(z_b < 0, -z_b, 0.0) * m_b
                c_pos, c_neg = wp_d.sum() + wp_b.sum(), wn_d.sum() + wn_b.sum()
                res[r] = math.log(c_pos / c_neg)
                if want_jac:
                    d_pos = nu_unk.T @ wp_d + wp_b[unk]
                    d_neg = nu_unk.T @ wn_d + wn_b[unk]
                    jac[r] = d_pos / c_pos - d_neg / c_neg
            elif kind == "ph":
                const = -eqs[r][1] * LN10 + (lg_b_ln[j] if ph_scale == "mesmer" else 0.0)
                res[r] = uvec[j] - const
                if want_jac:
                    jac[r, r] = 1.0
            elif kind == "loga":
                res[r] = uvec[j] - eqs[r][1]
                if want_jac:
                    jac[r, r] = 1.0
            elif kind == "phase":
                p = eqs[r][1]
                row = nu_ph[p]
                res[r] = (row[0] * log_aw * LN10 + row[unk] @ uvec[unk]) / LN10 \
                    - logk_ph[p] - eqs[r][2]
                if want_jac:
                    jac[r] = row[unk] / LN10
            elif kind == "eh":
                res[r] = uvec[j] + 4.0 * uvec[ib[PROTON]] - LN10 * (
                    logk_eh + 4.0 * eqs[r][1] + 2.0 * log_aw)
                if want_jac:
                    jac[r, r] += 1.0
                    jac[r, pos[ib[PROTON]]] += 4.0
            elif kind == "couple":
                res[r] = uvec[couple_col] - (couple_row[0] * log_aw * LN10
                                             + couple_row[unk] @ uvec[unk]
                                             - LN10 * couple_logk)
                if want_jac:
                    jac[r] = -couple_row[unk]
                    jac[r, pos[couple_col]] += 1.0
        return res, jac

    def pre_newton(u, lg_b_ln, lg_d_ln, log_aw, rounds=300):
        """Successive substitution on the mass balances (EQ3NR's pre-Newton-Raphson
        optimisation): each component corrects only for its own excess, which
        avoids the near-singular Jacobian of strongly complexed components."""
        mass = [(r, j, start_eqs[r][1]) for r, j in enumerate(unk) if start_eqs[r][0] == "mass"]
        if not mass:
            return u
        for _ in range(rounds):
            m_dep, m_b = state(u, lg_b_ln, lg_d_ln)
            biggest = 0.0
            step = np.zeros(nbf)
            for r, j, total_ in mass:
                tot = m_b[j] + nu[:, j] @ m_dep
                if tot > 0:
                    delta = max(-2.3, min(2.3, math.log(total_ / tot)))
                    step[j] = delta
                    biggest = max(biggest, abs(delta))
            u = u + step
            for r, j in enumerate(unk):
                if eqs[r][0] == "eh":
                    u[j] = LN10 * (logk_eh + 4.0 * eqs[r][1] + 2.0 * log_aw) - 4.0 * u[ib[PROTON]]
            if biggest < 0.02:
                break
        return u

    ionic, sum_m, log_aw = 0.0, 0.0, 0.0
    lg_b, lg_d = log_gammas(ionic)
    u = pre_newton(u, lg_b * LN10, lg_d * LN10, log_aw)
    total_inner, converged = 0, False
    relax, last_delta, outer_tol = 1.0, math.inf, 1e-11
    if model.pitzer is not None:
        max_outer, outer_tol = max_outer * 4, 1e-9
    for outer in range(max_outer):
        lg_b_ln, lg_d_ln = lg_b * LN10, lg_d * LN10
        u[0] = log_aw * LN10                     # ln a(H2O) enters every reaction involving water
        for inner in range(max_inner):
            total_inner += 1
            res, jac = residual(u, lg_b_ln, lg_d_ln, log_aw)
            norm0 = np.max(np.abs(res))
            if verbose:
                worst = int(np.argmax(np.abs(res)))
                print(f"outer {outer} inner {inner}: |res|={norm0:.3e} worst eq "
                      f"{names[unk[worst]]} ({eqs[worst][0]})")
                if verbose and eqs[worst][0] == "mass":
                    md, mb = state(u, lg_b_ln, lg_d_ln)
                    col = nu[:, unk[worst]]
                    top = np.argsort(-np.abs(col * md))[:4]
                    print("    basis m=%.3e  T=%.3e  top:" % (mb[unk[worst]], eqs[worst][1]),
                          [(dep_names[i], f"{col[i]:g}", f"{md[i]:.2e}") for i in top])
            if norm0 < tol:
                break
            u = _globalized_step(u, unk, res, jac, step_cap,
                                 lambda v: residual(v, lg_b_ln, lg_d_ln, log_aw, want_jac=False)[0])
        else:
            raise SpeciationError("inner Newton iteration did not converge")

        m_dep, m_b = state(u, lg_b_ln, lg_d_ln)
        ionic_new = 0.5 * (np.sum(m_b * z_b ** 2) + np.sum(m_dep * z_dep ** 2))
        sum_m_new = np.sum(m_b) + np.sum(m_dep)
        lg_b_new, lg_d_new, log_aw_new = model.activities(m_dep, m_b, ionic_new, sum_m_new)
        delta = max(abs(ionic_new - ionic), np.max(np.abs(lg_b_new - lg_b)),
                    np.max(np.abs(lg_d_new - lg_d)) if len(lg_d) else 0.0,
                    abs(log_aw_new - log_aw))
        if model.pitzer is not None:             # under-relaxed, as in equilibrium.py
            if delta > 0.7 * last_delta:
                relax = max(0.05, 0.5 * relax)
            elif delta < 0.3 * last_delta:
                relax = min(1.0, 1.5 * relax)
            last_delta = delta
            lg_b, lg_d = lg_b + relax * (lg_b_new - lg_b), lg_d + relax * (lg_d_new - lg_d)
            log_aw += relax * (log_aw_new - log_aw)
            ionic, sum_m = ionic_new, sum_m_new
        else:
            ionic, sum_m, log_aw, lg_b, lg_d = ionic_new, sum_m_new, log_aw_new, lg_b_new, lg_d_new
        if delta < outer_tol and outer > 0:
            converged = True
            break
    if not converged:
        raise SpeciationError("activity-coefficient iteration did not converge")

    return make_speciation(model, u, log_aw, lg_b, lg_d, ionic, sum_m, total_inner)

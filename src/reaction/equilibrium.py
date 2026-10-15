"""Closed-system equilibrium: an aqueous solution plus pure minerals.

The system is defined by the *bulk* moles of its components (the basis species,
including H2O for the solvent and O2(g) for the oxidation state).  The solver finds the
mass of solvent W, the ln activities u_j of the basis species and the moles n_p of each
mineral in the stable assemblage, so that

    W * tot_j(u) + sum_p nu_pj n_p = B_j      for every component j (mass balance)
    log Q/K(p) = 0                            for every mineral p present

A joint Newton iteration solves these for a given assemblage; an active-set loop adds
the most supersaturated mineral and drops minerals whose amount goes negative.
Components held at fixed activity (``fixed``) or tied to a gas fugacity (``buffers``)
are open to an external reservoir and not conserved.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np

from .aqueous import AqueousModel
from .constants import LN10, WATER_MOLES_PER_KG
from .database.model import Kind, ThermoDatabase
from .network import PROTON, WATER, Network, NetworkError
from .speciation import Speciation, SpeciationError, _globalized_step, make_speciation


class EquilibriumError(SpeciationError):
    pass


class System:
    """The basis, reaction network and options shared by all conditions of one problem."""

    def __init__(self, db: ThermoDatabase, basis: list[str], *,
                 switches: dict[str, str] | None = None, suppress=(), suppress_phases=(),
                 activity_model: str = "bdot", ph_scale: str = "nbs"):
        basis = [WATER, PROTON] + [b for b in basis if b not in (WATER, PROTON)]
        try:
            self.net = Network(db, basis, switches=switches, suppress=suppress,
                               suppress_phases=suppress_phases)
        except NetworkError as exc:
            raise EquilibriumError(str(exc)) from exc
        self.db = db
        self.activity_model, self.ph_scale = activity_model, ph_scale
        self._models: dict[tuple, AqueousModel] = {}

    @property
    def basis(self) -> list[str]:
        return self.net.basis

    def model(self, temperature_c: float, pressure_bar: float | None = None) -> AqueousModel:
        if pressure_bar is None:
            pressure_bar = self.db.reference_pressure(temperature_c)
        key = (round(temperature_c, 9), None if pressure_bar is None else round(pressure_bar, 9))
        if key not in self._models:
            if len(self._models) > 32:
                self._models.clear()
            self._models[key] = AqueousModel(self.db, self.net, temperature_c, pressure_bar,
                                             self.activity_model, self.ph_scale)
        return self._models[key]


@dataclass
class EquilibriumState:
    speciation: Speciation
    water_kg: float
    minerals: dict[str, float]               # moles of each mineral in the assemblage
    bulk: dict[str, float]                   # moles of each component in the whole system
    iterations: int = 0
    open_components: tuple[str, ...] = ()
    _u: np.ndarray = field(default=None, repr=False)

    @property
    def temperature_c(self) -> float:
        return self.speciation.temperature_c

    @property
    def ph(self) -> float:
        return self.speciation.ph

    def moles_in_solution(self, component: str) -> float:
        """Moles of a basis component held by the aqueous phase (solvent included)."""
        return self.water_kg * self.speciation.component_totals[component]


# --------------------------------------------------------------------------------------------
def _initial_u(model: AqueousModel, B: np.ndarray, fixed: dict[int, float]) -> np.ndarray:
    """Starting ln activities from the bulk composition (fluid only, no minerals)."""
    nb = model.nbf
    w0 = max(B[0], 1e-9) / WATER_MOLES_PER_KG
    u = np.zeros(nb)
    for j in range(1, nb):
        if j in fixed:
            u[j] = fixed[j]
        elif B[j] > 0:
            u[j] = math.log(B[j] / w0)
        else:
            u[j] = -9.0 * LN10 if model.names[j] == PROTON else -30.0 * LN10
    lg0 = np.zeros(nb)
    lgd = np.zeros(len(model.dep_names))
    for _ in range(200):
        m_dep, m_b = model.state(u, lg0, lgd)
        step = np.zeros(nb)
        biggest = 0.0
        for j in range(1, nb):
            if j in fixed or B[j] <= 0 or model.names[j] == PROTON or not model.solute[j]:
                continue
            tot = m_b[j] + model.nu[:, j] @ m_dep
            if tot > 0:
                d = max(-2.3, min(2.3, math.log(B[j] / w0 / tot)))
                step[j] = d
                biggest = max(biggest, abs(d))
        u += step
        if biggest < 0.05:
            break
    return u


class _Solve:
    """Joint Newton solution for a fixed mineral assemblage."""

    def __init__(self, model: AqueousModel, B: np.ndarray, fixed: dict[int, float],
                 buffer_rows: dict[int, tuple[int, float]]):
        self.m, self.B = model, B
        self.fixed, self.buffer_rows = fixed, buffer_rows
        nb = model.nbf
        self.unk = [j for j in range(1, nb) if j not in fixed]
        self.pos = {j: r for r, j in enumerate(self.unk)}
        self.rows = [j for j in range(nb) if j not in fixed]

    # ---- residual and Jacobian of the combined system --------------------------------------
    def residual(self, x, assemblage, u_base, lg_b_ln, lg_d_ln, log_aw, want_jac=True):
        m = self.m
        unk, nu_cnt = self.unk, len(self.unk)
        u = u_base.copy()
        u[unk] = x[:nu_cnt]
        u[0] = log_aw * LN10
        w = math.exp(x[nu_cnt])
        n = x[nu_cnt + 1:]
        nu, nu_ph = m.nu, m.nu_ph
        m_dep, m_b = m.state(u, lg_b_ln, lg_d_ln)
        nu_u = nu[:, unk]
        size = nu_cnt + 1 + len(assemblage)
        res = np.zeros(len(self.rows) + len(assemblage))
        jac = np.zeros((len(res), size)) if want_jac else None
        for r, j in enumerate(self.rows):
            if j in self.buffer_rows:                    # gas fugacity replaces this balance
                p, target = self.buffer_rows[j]
                row = nu_ph[p]
                res[r] = row @ u - LN10 * (m.logk_ph[p] + target)
                if want_jac:
                    jac[r, :nu_cnt] = row[unk]
                continue
            if j == 0:
                tot = WATER_MOLES_PER_KG + nu[:, 0] @ m_dep
            else:
                tot = m_b[j] + nu[:, j] @ m_dep
            dtot = (nu[:, j] * m_dep) @ nu_u
            if j != 0 and j in self.pos:
                dtot[self.pos[j]] += m_b[j]
            mineral = sum(nu_ph[a, j] * n[k] for k, a in enumerate(assemblage))
            aq = w * tot
            b = self.B[j]
            involved = any(abs(nu_ph[a, j]) > 1e-12 for a in assemblage)
            if j != 0 and not involved:
                # Balance written as ln(sources / sinks) = 0: the aqueous terms are split by sign
                # and the bulk amount c = b joins the smaller side.  Far from the root this is
                # linear in u (a redox couple moves by tens of log units per unit of reaction
                # progress); near it, it is the usual mass balance.  Components held in a
                # mineral keep the linear form: the logarithm is singular as the mineral takes
                # up the last of the component.
                col = nu[:, j]
                terms = col * m_dep
                plus, minus = terms * (col > 0), terms * (col < 0)
                p_aq = w * (m_b[j] + float(plus.sum()))
                n_aq = -w * float(minus.sum())
                c = b
                if c >= 0.0 and n_aq + c > 0.0:
                    den = n_aq + c
                    res[r] = math.log(p_aq / den)
                    mode = 1
                elif c < 0.0 and n_aq > 0.0:
                    den = p_aq - c
                    res[r] = math.log(den / n_aq)
                    mode = 2
                else:
                    mode = 0
                if mode:
                    if want_jac:
                        dp = w * (plus @ nu_u)
                        dn = -w * (minus @ nu_u)
                        if j in self.pos:
                            dp[self.pos[j]] += w * m_b[j]
                        if mode == 1:
                            jac[r, :nu_cnt] = dp / p_aq - dn / den
                            jac[r, nu_cnt] = 1.0 - n_aq / den
                        else:
                            jac[r, :nu_cnt] = dp / den - dn / n_aq
                            jac[r, nu_cnt] = p_aq / den - 1.0
                    continue
            scale = abs(aq) + abs(b) + sum(abs(nu_ph[a, j] * n[k]) for k, a in enumerate(assemblage))
            scale = scale if scale > 0 else 1.0
            res[r] = (aq + mineral - b) / scale
            if want_jac:
                jac[r, :nu_cnt] = w * dtot / scale
                jac[r, nu_cnt] = aq / scale
                for k, a in enumerate(assemblage):
                    jac[r, nu_cnt + 1 + k] = nu_ph[a, j] / scale
        base = len(self.rows)
        for k, a in enumerate(assemblage):               # mineral saturation, log Q/K = 0
            row = nu_ph[a]
            res[base + k] = row @ u - LN10 * m.logk_ph[a]
            if want_jac:
                jac[base + k, :nu_cnt] = row[unk]
        return res, jac

    # ---- solve for one assemblage ----------------------------------------------------------------
    def solve(self, assemblage, u, ln_w, n, ionic0, tol=1e-10, max_newton=80, max_outer=60):
        m = self.m
        unk, nu_cnt = self.unk, len(self.unk)
        u = u.copy()
        for j, v in self.fixed.items():
            u[j] = v
        x = np.concatenate([u[unk], [ln_w], n])
        cap = np.concatenate([np.full(nu_cnt, 2.3), [0.5], np.full(len(assemblage), np.inf)])
        everything = np.arange(len(x))
        ionic = ionic0
        lg_b, lg_d = m.log_gammas(ionic)
        log_aw = 0.0                          # refined after the first pass
        iterations = 0
        sum_m = 0.0
        for outer in range(max_outer):
            lg_b_ln, lg_d_ln = lg_b * LN10, lg_d * LN10
            for inner in range(max_newton):
                iterations += 1
                res, jac = self.residual(x, assemblage, u, lg_b_ln, lg_d_ln, log_aw)
                if np.max(np.abs(res)) < tol:
                    break
                x = _globalized_step(
                    x, everything, res, jac, cap,
                    lambda v: self.residual(v, assemblage, u, lg_b_ln, lg_d_ln, log_aw,
                                            want_jac=False)[0])
            else:
                raise EquilibriumError("Newton iteration did not converge")
            u_now = u.copy()
            u_now[unk] = x[:nu_cnt]
            u_now[0] = log_aw * LN10
            m_dep, m_b = m.state(u_now, lg_b_ln, lg_d_ln)
            ionic_new, sum_m = m.ionic_and_sum(m_dep, m_b)
            lg_b_new, lg_d_new = m.log_gammas(ionic_new)
            log_aw_new = m.log_a_water(sum_m, ionic_new)
            delta = max(abs(ionic_new - ionic), float(np.max(np.abs(lg_b_new - lg_b))),
                        float(np.max(np.abs(lg_d_new - lg_d))) if len(lg_d) else 0.0,
                        abs(log_aw_new - log_aw))
            ionic, lg_b, lg_d, log_aw = ionic_new, lg_b_new, lg_d_new, log_aw_new
            if delta < 1e-11 and outer > 0:
                break
        else:
            raise EquilibriumError("activity coefficients did not converge")
        u_out = u.copy()
        u_out[unk] = x[:nu_cnt]
        return dict(u=u_out, ln_w=float(x[nu_cnt]), n=np.array(x[nu_cnt + 1:]), ionic=ionic,
                    sum_m=sum_m, log_aw=log_aw, lg_b=lg_b, lg_d=lg_d, iterations=iterations)


def equilibrate(system: System, temperature_c: float, bulk: dict[str, float], *,
                pressure_bar: float | None = None,
                fixed: dict[str, float] | None = None,
                buffers: Iterable[tuple[str, float]] = (),
                minerals: Iterable[str] = (),
                amounts: dict[str, float] | None = None,
                eligible: Iterable[str] | None = None,
                guess: EquilibriumState | Speciation | None = None,
                tol: float = 1e-10, max_attempts: int = 80) -> EquilibriumState:
    """Equilibrium of a closed system with the given bulk composition.

    ``bulk`` maps basis-species (component) names to total moles in the whole system.
    ``fixed`` holds components at a log activity (log fugacity for O2(g)), open to a
    reservoir.  ``buffers`` are (gas, log fugacity) pairs; each replaces the mass balance
    of the component the gas depends on most.  ``minerals`` are the minerals to start the
    assemblage with and ``eligible`` the minerals allowed to precipitate (default: all
    that are not suppressed).  ``amounts`` suggests starting moles for minerals not in the
    ``guess`` (for example a reactant that has just been converted to an equilibrium mineral).
    """
    model = system.model(temperature_c, pressure_bar)
    names, ib = model.names, model.ib
    nb = model.nbf
    unknown = [k for k in bulk if k not in ib]
    if unknown:
        raise EquilibriumError(f"components {unknown} are not in the basis {names}")
    B = np.array([bulk.get(n, 0.0) for n in names], dtype=float)
    if B[0] <= 0:
        raise EquilibriumError("the system contains no water")

    fixed_idx = {ib[n]: v * LN10 for n, v in (fixed or {}).items()}
    for n in (fixed or {}):
        if n not in ib or n == WATER:
            raise EquilibriumError(f"cannot fix {n!r}")
    buffer_rows: dict[int, tuple[int, float]] = {}
    used = {0, ib[PROTON], *fixed_idx}
    for gas, target in buffers:
        if gas not in model.ph_names:
            raise EquilibriumError(f"buffer gas {gas!r} is unavailable at this temperature")
        p = model.ph_names.index(gas)
        row = np.abs(model.nu_ph[p])
        order = [j for j in np.argsort(-row) if j not in used and row[j] > 1e-12]
        if not order:
            raise EquilibriumError(f"no free component to buffer with {gas}")
        buffer_rows[int(order[0])] = (p, target)
        used.add(int(order[0]))

    phase_index = {n: i for i, n in enumerate(model.ph_names)}
    inactive = system.net.inactive_phases
    allowed = [i for i, (n, k) in enumerate(zip(model.ph_names, model.ph_kind))
               if k == Kind.MINERAL and n not in inactive
               and (eligible is None or n in set(eligible))]
    allowed_set = set(allowed)
    assemblage = [phase_index[n] for n in minerals if n in phase_index]
    for a in assemblage:
        if a not in allowed_set:
            allowed.append(a)
            allowed_set.add(a)

    problem = _Solve(model, B, fixed_idx, buffer_rows)

    # warm start
    if isinstance(guess, EquilibriumState):
        u0 = np.array(guess.speciation.ln_activities, dtype=float)
        ln_w0 = math.log(guess.water_kg)
        n_prev = dict(guess.minerals)
        ionic0 = guess.speciation.ionic_strength
    elif isinstance(guess, Speciation):
        u0 = np.array(guess.ln_activities, dtype=float)
        ln_w0 = math.log(B[0] / WATER_MOLES_PER_KG)
        n_prev, ionic0 = {}, guess.ionic_strength
    else:
        u0 = _initial_u(model, B, fixed_idx)
        ln_w0 = math.log(B[0] / WATER_MOLES_PER_KG)
        n_prev, ionic0 = {}, 0.0
    if len(u0) != nb:
        raise EquilibriumError("the starting guess belongs to a different basis")
    for j, v in fixed_idx.items():
        u0[j] = v

    hints = dict(amounts or {})

    def initial_n(assem):
        return np.array([n_prev.get(model.ph_names[a], hints.get(model.ph_names[a], 0.0))
                         for a in assem], dtype=float)

    scale_b = max(float(np.max(np.abs(B[1:]))) if nb > 1 else 1.0, 1e-30)
    n_floor = 1e-12 * scale_b
    seen: set[frozenset] = set()
    blocked: set[int] = set()                 # minerals whose addition failed
    u_cur, ln_w_cur, ionic_cur = u0, ln_w0, ionic0
    last = None
    total_iter = 0
    for attempt in range(max_attempts):
        key = frozenset(assemblage)
        if key in seen:
            raise EquilibriumError("the phase assemblage cycles: no consistent set of minerals found")
        seen.add(key)
        try:
            last = problem.solve(assemblage, u_cur, ln_w_cur, initial_n(assemblage), ionic_cur,
                                 tol=tol)
        except (SpeciationError, np.linalg.LinAlgError) as exc:
            if not assemblage:
                raise EquilibriumError(f"aqueous solution alone failed: {exc}") from exc
            blocked.add(assemblage[-1])           # the newest mineral could not be added
            assemblage = assemblage[:-1]
            continue
        total_iter += last["iterations"]
        u_cur, ln_w_cur, ionic_cur = last["u"], last["ln_w"], last["ionic"]
        n = last["n"]
        n_prev = {model.ph_names[a]: float(n[k]) for k, a in enumerate(assemblage)}
        negative = [k for k in range(len(assemblage)) if n[k] < -n_floor]
        if negative:
            worst = min(negative, key=lambda k: n[k])
            n_prev.pop(model.ph_names[assemblage[worst]])
            assemblage = [a for k, a in enumerate(assemblage) if k != worst]
            continue
        si = (model.nu_ph @ last["u"]) / LN10 - model.logk_ph
        candidates = [(si[a], a) for a in allowed
                      if a not in assemblage and a not in blocked and si[a] > 1e-7]
        if candidates:
            _, add = max(candidates)
            assemblage = assemblage + [add]
            continue
        stuck = [model.ph_names[a] for a in blocked if a in allowed_set and a not in assemblage
                 and si[a] > 1e-6]
        if stuck:                                     # never return a state that is supersaturated
            raise EquilibriumError(f"could not add the supersaturated mineral(s) {stuck}")
        break
    else:
        raise EquilibriumError("no stable phase assemblage found")

    spec = make_speciation(model, last["u"], last["log_aw"], last["lg_b"], last["lg_d"],
                           last["ionic"], last["sum_m"], total_iter)
    n = last["n"]
    present = {model.ph_names[a]: float(n[k]) for k, a in enumerate(assemblage)}
    open_names = tuple(names[j] for j in list(fixed_idx) + list(buffer_rows))
    return EquilibriumState(spec, math.exp(last["ln_w"]), present, dict(bulk), total_iter,
                            open_names, _u=last["u"])

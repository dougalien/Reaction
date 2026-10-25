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
from .network import OXYGEN, PROTON, WATER, Network, NetworkError
from .speciation import Speciation, SpeciationError, _globalized_step, make_speciation
from .thermo import gas as gas_eos


class EquilibriumError(SpeciationError):
    pass


class System:
    """The basis, reaction network and options shared by all conditions of one problem."""

    def __init__(self, db: ThermoDatabase, basis: list[str], *,
                 switches: dict[str, str] | None = None, suppress=(), suppress_phases=(),
                 activity_model: str = "bdot", ph_scale: str = "nbs", gas_phase: bool = False,
                 hold_apart=()):
        basis = [WATER, PROTON] + [b for b in basis if b not in (WATER, PROTON)]
        # a species in the basis that is not a strict basis species (CH4 from a speciation that held it
        # apart) is a component of its own
        targets = set((switches or {}).values())
        own = [b for b in basis if b in db and db[b].kind == Kind.AQUEOUS and b not in targets]
        basis = [b for b in basis if b not in own]
        hold_apart = tuple(hold_apart) + tuple(b for b in own if b not in hold_apart)
        try:
            self.net = Network(db, basis, switches=switches, suppress=suppress,
                               suppress_phases=suppress_phases, decoupled=tuple(hold_apart))
        except NetworkError as exc:
            raise EquilibriumError(str(exc)) from exc
        for d in hold_apart:                         # its gas follows it, not the rest of the basis
            _tie_gas_to(self.net, db, d)
        self.db = db
        self.gas_phase = gas_phase
        self.hold_apart = tuple(hold_apart)
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
            if len(self._models) > 400:
                self._models.clear()
            self._models[key] = AqueousModel(self.db, self.net, temperature_c, pressure_bar,
                                             self.activity_model, self.ph_scale)
        return self._models[key]


def _tie_gas_to(net: Network, db, species: str) -> None:
    """Rewrite the gas of a held-apart species (CH4(g) for CH4) over that species' own column, so that
    the gas exchanges with it and not with the rest of the basis.  Both are written over the core basis
    with the same formula, so CH4(g) = CH4(aq) with the difference of their logK."""
    gas = f"{db[species].formula}(g)"
    if gas not in net.phase_names:
        return
    k = net.phase_names.index(gas)
    own, terms = net.couple[species]
    if np.max(np.abs(net.nu_phase[k] - own)) > 1e-9:
        return
    col = net.basis.index(species)
    row = np.zeros(len(net.basis))
    row[col] = 1.0
    net.nu_phase[k] = row
    net.phase_terms[k] = list(net.phase_terms[k]) + [(-cf, lk) for cf, lk in terms]


@dataclass
class EquilibriumState:
    speciation: Speciation
    water_kg: float
    minerals: dict[str, float]               # moles of each mineral in the assemblage
    bulk: dict[str, float]                   # moles of each component in the whole system
    iterations: int = 0
    open_components: tuple[str, ...] = ()
    gas: dict[str, float] = field(default_factory=dict)      # moles of each gas in the gas phase
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
        o2 = model.ib.get(OXYGEN)
        if o2 is not None and o2 not in fixed:     # the redox budget can be negative (a reduced
            before = u[o2]                          # fluid): match it by bisection, as it is monotonic
            u[o2] = _match_total(model, u, o2, B[o2] / w0, lg0, lgd)
            biggest = max(biggest, abs(u[o2] - before))
        if biggest < 0.05:
            break
    return u


def _match_total(model: AqueousModel, u: np.ndarray, j: int, target: float, lg0, lgd) -> float:
    """The ln activity of basis species ``j`` at which its component total equals ``target``, the
    others held; the total rises with u_j (its derivative is a sum of nu^2 m)."""
    trial = u.copy()

    def total(v):
        trial[j] = v
        m_dep, m_b = model.state(trial, lg0, lgd)
        return (m_b[j] if model.solute[j] else 0.0) + model.nu[:, j] @ m_dep

    lo, hi = -250.0 * LN10, 10.0 * LN10
    with np.errstate(over="ignore", invalid="ignore"):
        if not total(lo) < target < total(hi):
            return float(u[j])
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if total(mid) < target:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-3:
                break
    return 0.5 * (lo + hi)


class _Gas:
    """The gases of the network as one ideal mixture of real gases: x_i = f_i / (phi_i P), the phase
    exists when the x_i add up to one.  O2(g), a basis column, is one of them."""

    def __init__(self, model: AqueousModel):
        nb = model.nbf
        idx = [i for i, k in enumerate(model.ph_kind) if k == Kind.GAS]
        rows = [model.nu_ph[i] for i in idx]
        logk = [float(model.logk_ph[i]) for i in idx]
        names = [model.ph_names[i] for i in idx]
        if OXYGEN in model.ib:
            row = np.zeros(nb)
            row[model.ib[OXYGEN]] = 1.0
            rows.append(row); logk.append(0.0); names.append(OXYGEN)
        self.names = names
        self.formulas = [gas_eos.formula_of(n) for n in names]
        self.nu = np.array(rows).reshape(-1, nb)
        self.logk = np.array(logk)
        self.t_k = model.temperature_c + 273.15
        self.p = float(model.pressure_bar)
        self.log_p = math.log10(self.p)
        self.log_phi = np.zeros(len(names))

    def fractions(self, u: np.ndarray) -> np.ndarray:
        """x_i at ln activities ``u`` with the current fugacity coefficients."""
        return np.exp(np.clip(self.nu @ u - LN10 * (self.logk + self.log_phi + self.log_p), -600.0, 50.0))

    def update_phi(self, u: np.ndarray) -> float:
        """Fugacity coefficients for the composition at ``u``; returns the largest change."""
        x = self.fractions(u)
        new = gas_eos.log_phi(self.formulas, x, self.t_k, self.p)
        change = float(np.max(np.abs(new - self.log_phi))) if len(new) else 0.0
        self.log_phi = new
        return change


class _Solve:
    """Joint Newton solution for a fixed mineral assemblage (and the gas phase, if present)."""

    def __init__(self, model: AqueousModel, B: np.ndarray, fixed: dict[int, float],
                 buffer_rows: dict[int, tuple[int, float]], gas: _Gas | None = None):
        self.m, self.B = model, B
        self.gas = gas
        self.fixed, self.buffer_rows = fixed, buffer_rows
        nb = model.nbf
        self.unk = [j for j in range(1, nb) if j not in fixed]
        self.pos = {j: r for r, j in enumerate(self.unk)}
        self.rows = [j for j in range(nb) if j not in fixed]
        self.nu_r = model.nu[:, self.rows]                   # species x balance rows, and its signs
        self.nu_pos, self.nu_neg = np.maximum(self.nu_r, 0.0), np.minimum(self.nu_r, 0.0)
        self.nu_u = model.nu[:, self.unk]

    # ---- residual and Jacobian of the combined system --------------------------------------
    def residual(self, x, assemblage, u_base, lg_b_ln, lg_d_ln, log_aw, want_jac=True, gas=False):
        """Every balance at once (matrix products over all components); the rows then pick the
        form that suits them.  This is the inner loop of every line search.  With ``gas`` the last
        unknown is the moles of gas and the last equation says its mole fractions add up to one."""
        m = self.m
        unk, nu_cnt = self.unk, len(self.unk)
        u = u_base.copy()
        u[unk] = x[:nu_cnt]
        u[0] = log_aw * LN10
        w = math.exp(x[nu_cnt])
        n_end = len(x) - (1 if gas else 0)
        n = np.asarray(x[nu_cnt + 1:n_end], dtype=float)
        nu_ph = m.nu_ph
        m_dep, m_b = m.state(u, lg_b_ln, lg_d_ln)
        rows = self.rows
        ph = nu_ph[assemblage][:, rows] if assemblage else np.zeros((0, len(rows)))   # minerals x rows
        mb_r = m_b[rows]
        tot = mb_r + m_dep @ self.nu_r
        if rows and rows[0] == 0:
            tot[0] = WATER_MOLES_PER_KG + m_dep @ self.nu_r[:, 0]
        aq = w * tot
        b = self.B[rows]
        mineral = n @ ph
        involved = (np.abs(ph) > 1e-12).any(axis=0)
        if gas:
            g_amt = float(x[-1])
            y = self.gas.fractions(u)
            g_nu = self.gas.nu[:, rows]
            per_mole = y @ g_nu                              # components per mole of gas
            mineral = mineral + g_amt * per_mole
            involved = involved | (np.abs(g_nu) > 1e-12).any(axis=0)
        p_aq = w * (mb_r + m_dep @ self.nu_pos)
        n_aq = -w * (m_dep @ self.nu_neg)
        rows_arr = np.array(rows)
        log_ok = (rows_arr != 0) & ~involved
        mode1 = log_ok & (b >= 0.0) & (n_aq + b > 0.0)
        mode2 = log_ok & ~mode1 & (b < 0.0) & (n_aq > 0.0)
        linear = ~(mode1 | mode2)
        scale = np.abs(aq) + np.abs(b) + np.abs(n) @ np.abs(ph)
        if gas:
            scale = scale + abs(g_amt) * (y @ np.abs(g_nu))
        scale = np.where(scale > 0, scale, 1.0)
        size = nu_cnt + 1 + len(assemblage) + (1 if gas else 0)
        res = np.zeros(len(rows) + len(assemblage) + (1 if gas else 0))
        den1 = np.where(mode1, n_aq + b, 1.0)
        den2 = np.where(mode2, p_aq - b, 1.0)
        nr = len(rows)
        with np.errstate(divide="ignore", invalid="ignore"):
            res[:nr] = np.where(mode1, np.log(np.where(mode1, p_aq, 1.0) / den1),
                                np.where(mode2, np.log(den2 / np.where(mode2, n_aq, 1.0)),
                                         (aq + mineral - b) / scale))
        for r, j in enumerate(rows):
            if j in self.buffer_rows:                    # gas fugacity replaces this balance
                p, target = self.buffer_rows[j]
                res[r] = nu_ph[p] @ u - LN10 * (m.logk_ph[p] + target)
        base = nr
        na = len(assemblage)
        if assemblage:                                   # mineral saturation, log Q/K = 0
            res[base:base + na] = nu_ph[assemblage] @ u - LN10 * m.logk_ph[assemblage]
        if gas:
            res[-1] = y.sum() - 1.0
        if not want_jac:
            return res, None

        jac = np.zeros((len(res), size))
        weighted = m_dep[:, None] * self.nu_u                # d m_dep / d u, species x unknowns
        own = np.zeros((nr, nu_cnt))                     # d m_b[j] / d u_j on the diagonal
        for r, j in enumerate(rows):
            if j != 0 and j in self.pos:
                own[r, self.pos[j]] = m_b[j]
        dtot = self.nu_r.T @ weighted + own
        dp = w * (self.nu_pos.T @ weighted + own)
        dn = -w * (self.nu_neg.T @ weighted)
        lin = linear[:, None]
        jac[:nr, :nu_cnt] = np.where(
            mode1[:, None], dp / np.where(mode1, p_aq, 1.0)[:, None] - dn / den1[:, None],
            np.where(mode2[:, None], dp / den2[:, None] - dn / np.where(mode2, n_aq, 1.0)[:, None],
                     w * dtot / scale[:, None]))
        jac[:nr, nu_cnt] = np.where(mode1, 1.0 - n_aq / den1,
                                    np.where(mode2, p_aq / den2 - 1.0, aq / scale))
        if assemblage:
            jac[:nr, nu_cnt + 1:nu_cnt + 1 + na] = np.where(lin, ph.T / scale[:, None], 0.0)
        if gas:
            g_unk = self.gas.nu[:, unk]
            d_gas = g_amt * (g_nu.T @ (y[:, None] * g_unk))          # rows x unknowns
            jac[:nr, :nu_cnt] += np.where(lin, d_gas / scale[:, None], 0.0)
            jac[:nr, -1] = np.where(linear, per_mole / scale, 0.0)
            jac[-1, :nu_cnt] = y @ g_unk
        for r, j in enumerate(rows):
            if j in self.buffer_rows:
                p, _ = self.buffer_rows[j]
                jac[r, :] = 0.0
                jac[r, :nu_cnt] = nu_ph[p][unk]
        if assemblage:
            jac[base:base + na, :nu_cnt] = nu_ph[assemblage][:, unk]
        return res, jac

    # ---- solve for one assemblage ----------------------------------------------------------------
    def solve(self, assemblage, u, ln_w, n, ionic0, tol=1e-10, max_newton=80, max_outer=60,
              gas_moles=None):
        m = self.m
        unk, nu_cnt = self.unk, len(self.unk)
        gas = gas_moles is not None
        u = u.copy()
        for j, v in self.fixed.items():
            u[j] = v
        x = np.concatenate([u[unk], [ln_w], n] + ([[gas_moles]] if gas else []))
        cap = np.concatenate([np.full(nu_cnt, 2.3), [0.5], np.full(len(assemblage) + gas, np.inf)])
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
                res, jac = self.residual(x, assemblage, u, lg_b_ln, lg_d_ln, log_aw, gas=gas)
                if np.max(np.abs(res)) < tol:
                    break
                x = _globalized_step(
                    x, everything, res, jac, cap,
                    lambda v: self.residual(v, assemblage, u, lg_b_ln, lg_d_ln, log_aw,
                                            want_jac=False, gas=gas)[0])
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
            if gas:
                delta = max(delta, self.gas.update_phi(u_now))
            ionic, lg_b, lg_d, log_aw = ionic_new, lg_b_new, lg_d_new, log_aw_new
            if delta < 1e-11 and outer > 0:
                break
        else:
            raise EquilibriumError("activity coefficients did not converge")
        u_out = u.copy()
        u_out[unk] = x[:nu_cnt]
        u_out[0] = log_aw * LN10              # the water activity the solution was solved with
        n_end = len(x) - (1 if gas else 0)
        return dict(u=u_out, ln_w=float(x[nu_cnt]), n=np.array(x[nu_cnt + 1:n_end]), ionic=ionic,
                    sum_m=sum_m, log_aw=log_aw, lg_b=lg_b, lg_d=lg_d, iterations=iterations,
                    gas=float(x[-1]) if gas else None, log_phi=self.gas.log_phi.copy() if gas else None)


def _dependence(nu_ph: np.ndarray, assemblage: list[int], phase: int,
                unk: list[int]) -> np.ndarray | None:
    """Coefficients lam with nu(phase) = sum_k lam_k nu(assemblage[k]) over the unknown
    activities, or None when the phase's saturation condition is independent of theirs."""
    if not assemblage:
        return None
    rows = nu_ph[np.ix_(assemblage, unk)]
    target = nu_ph[phase, unk]
    lam, *_ = np.linalg.lstsq(rows.T, target, rcond=None)
    miss = float(np.max(np.abs(rows.T @ lam - target))) if len(unk) else 0.0
    return lam if miss <= 1e-9 * max(1.0, float(np.max(np.abs(target)))) else None


def _phase_rule_swap(nu_ph, assemblage, add, unk, n_prev, ph_names):
    """(mineral that leaves, moles of the newcomer, lam) when ``add`` cannot join the assemblage
    without making the saturation conditions dependent; (None, 0, None) when it can join.  The
    leaving mineral is picked by the ratio test, so every amount stays non-negative."""
    options, lam = _swap_options(nu_ph, assemblage, add, unk, n_prev, ph_names)
    if not options:
        return None, 0.0, None
    leave, theta = options[0]
    return leave, theta, lam


def _swap_options(nu_ph, assemblage, add, unk, n_prev, ph_names):
    """Every mineral that could leave to make room for ``add``, best (the ratio test) first, as
    (mineral, moles of the newcomer); and lam.  Empty when ``add`` can simply join."""
    lam = _dependence(nu_ph, assemblage, add, unk)
    if lam is None or not np.any(np.abs(lam) > 1e-12):      # independent, or no condition at all
        return [], None
    amounts = [max(n_prev.get(ph_names[a], 0.0), 0.0) for a in assemblage]
    ratios = sorted((amounts[k] / lam[k], k) for k in range(len(assemblage)) if lam[k] > 1e-12)
    if not ratios:
        return [(assemblage[int(np.argmax(np.abs(lam)))], 0.0)], lam
    return [(assemblage[k], theta) for theta, k in ratios], lam


def _most_alike(nu_ph, assemblage, add, unk, limit=4):
    """Minerals of the assemblage whose compositions are closest in direction to ``add``'s, as
    (mineral, 0) swap options, most alike first."""
    v = nu_ph[add, unk]
    scored = []
    for a in assemblage:
        w = nu_ph[a, unk]
        denom = float(np.linalg.norm(v) * np.linalg.norm(w))
        if denom > 0:
            scored.append((-abs(float(v @ w)) / denom, a))
    return [(a, 0.0) for _, a in sorted(scored)[:limit]]


def _swapped(assemblage, n_prev, ph_names, lam, leave, add, theta):
    """The assemblage and amounts after ``add`` replaces ``leave`` (a simplex pivot, or without
    ``lam`` a plain exchange)."""
    n_new = dict(n_prev)
    for k, a in enumerate(assemblage if lam is not None else ()):
        n_new[ph_names[a]] = n_new.get(ph_names[a], 0.0) - theta * lam[k]
    n_new.pop(ph_names[leave], None)
    n_new[ph_names[add]] = theta
    return [a for a in assemblage if a != leave] + [add], n_new


def _independent(nu_ph, assemblage: list[int], unk: list[int]) -> list[int]:
    """The assemblage with any mineral whose saturation condition depends on earlier ones left out."""
    kept: list[int] = []
    for a in assemblage:
        if _dependence(nu_ph, kept, a, unk) is None:
            kept.append(a)
    return kept


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

    gas = _Gas(model) if system.gas_phase and model.pressure_bar else None
    if gas is not None and not len(gas.names):
        gas = None
    problem = _Solve(model, B, fixed_idx, buffer_rows, gas)
    # minerals given starting amounts (a reactant that has just been converted) win ties against
    # earlier minerals with the same saturation condition
    hinted = [a for a in assemblage if model.ph_names[a] in (amounts or {})]
    assemblage = _independent(model.nu_ph, hinted + [a for a in assemblage if a not in hinted],
                              problem.unk)

    # warm start
    gas_on = None
    if isinstance(guess, EquilibriumState):
        u0 = np.array(guess.speciation.ln_activities, dtype=float)
        ln_w0 = math.log(guess.water_kg)
        n_prev = dict(guess.minerals)
        ionic0 = guess.speciation.ionic_strength
        if guess.gas and gas is not None:
            gas_on = sum(guess.gas.values())
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
    o2 = ib.get(OXYGEN)
    if guess is not None and o2 is not None and o2 not in fixed_idx and o2 not in buffer_rows:
        # log fO2 moves by many units with temperature: start it where the fluid's share of the
        # redox budget (the minerals keep their amounts) is matched at the new conditions
        held = sum(v * model.nu_ph[model.ph_names.index(k), o2] for k, v in n_prev.items()
                   if k in model.ph_names)
        zeros_b, zeros_d = np.zeros(nb), np.zeros(len(model.dep_names))
        u0[o2] = _match_total(model, u0, o2, (B[o2] - held) / math.exp(ln_w0), zeros_b, zeros_d)

    hints = dict(amounts or {})
    if gas is not None:                       # fugacity coefficients for the starting composition
        for _ in range(3):
            gas.update_phi(u0)

    def initial_n(assem):
        return np.array([n_prev.get(model.ph_names[a], hints.get(model.ph_names[a], 0.0))
                         for a in assem], dtype=float)

    # log Q/K per mole of solute components in the formula: unscaled, a clay or an amphibole looks far
    # more supersaturated than calcite for the same driving force per atom, and is tried first
    solute_cols = [j for j in range(1, nb) if names[j] not in (PROTON, OXYGEN)]
    formula_size = np.maximum(np.abs(model.nu_ph[:, solute_cols]).sum(axis=1), 1.0)
    scale_b = max(float(np.max(np.abs(B[1:]))) if nb > 1 else 1.0, 1e-30)
    n_floor = 1e-12 * scale_b
    # an assemblage reached twice is not solved twice, and from it the next most supersaturated mineral is
    # tried instead of the one that led back to it: the search walks on rather than cycling
    solved: dict[frozenset, dict] = {}
    tried: dict[frozenset, set[int]] = {}
    blocked: set[int] = set()                 # minerals whose addition failed
    swaps: dict[frozenset, tuple] = {}        # assemblage reached by a pivot -> how to try another
    u_cur, ln_w_cur, ionic_cur = u0, ln_w0, ionic0
    last = None
    total_iter = 0
    gas_tried: set = set()
    gas_blocked = False
    for attempt in range(max_attempts):
        key = frozenset(assemblage) | ({"gas"} if gas_on is not None else set())
        if key in solved:
            last = solved[key]
            assemblage = last["assemblage"]
        else:
            try:
                last = problem.solve(assemblage, u_cur, ln_w_cur, initial_n(assemblage), ionic_cur,
                                     tol=tol, gas_moles=gas_on)
            except (SpeciationError, np.linalg.LinAlgError) as exc:
                # across a redox jump the last state can be a poor start: try once from the
                # bulk composition itself before giving up on this assemblage
                try:
                    cold = _initial_u(model, B, fixed_idx)
                    last = problem.solve(assemblage, cold, math.log(B[0] / WATER_MOLES_PER_KG),
                                         initial_n(assemblage), 0.0, tol=tol, gas_moles=gas_on)
                except (SpeciationError, np.linalg.LinAlgError):
                    if gas_on is not None and key not in gas_tried:   # the gas phase would not solve
                        gas_tried.add(key)
                        gas_on = None
                        gas_blocked = True
                        continue
                    if not assemblage:
                        raise EquilibriumError(f"aqueous solution alone failed: {exc}") from exc
                    if key in swaps:                      # a pivot that did not solve: try the
                        prev, n_before, add, lam, options = swaps.pop(key)   # next mineral to leave
                        if options:
                            leave, theta = options.pop(0)
                            assemblage, n_prev = _swapped(prev, n_before, model.ph_names, lam,
                                                          leave, add, theta)
                            swaps[frozenset(assemblage)] = (prev, n_before, add, lam, options)
                            continue
                        blocked.add(add)
                        assemblage, n_prev = prev, n_before
                        continue
                    blocked.add(assemblage[-1])           # the newest mineral could not be added
                    assemblage = assemblage[:-1]
                    continue
            last["assemblage"] = list(assemblage)
            solved[key] = last
        total_iter += last["iterations"]
        u_cur, ln_w_cur, ionic_cur = last["u"], last["ln_w"], last["ionic"]
        n = last["n"]
        n_prev = {model.ph_names[a]: float(n[k]) for k, a in enumerate(assemblage)}
        if gas_on is not None:
            gas_on = last["gas"]
            if gas_on < -n_floor:                     # the gas phase has dissolved
                gas_on = None
                continue
        negative = [k for k in range(len(assemblage)) if n[k] < -n_floor]
        if negative:
            worst = min(negative, key=lambda k: n[k])
            n_prev.pop(model.ph_names[assemblage[worst]])
            assemblage = [a for k, a in enumerate(assemblage) if k != worst]
            continue
        if gas is not None and gas_on is None and not gas_blocked:
            gas.update_phi(last["u"])
            if gas.fractions(last["u"]).sum() > 1.0 + 1e-9:     # bubbles form
                gas_on = 1e-6 * scale_b
                continue
        si = (model.nu_ph @ last["u"]) / LN10 - model.logk_ph
        candidates = [(si[a] / formula_size[a], a) for a in allowed
                      if a not in assemblage and a not in blocked and si[a] > 1e-7]
        untried = [c for c in candidates if c[1] not in tried.setdefault(key, set())]
        if candidates and not untried:
            raise EquilibriumError("the phase assemblage cycles: no consistent set of minerals found")
        if candidates:
            _, add = max(untried)
            tried[key].add(add)
            options, lam = _swap_options(model.nu_ph, assemblage, add, problem.unk, n_prev,
                                         model.ph_names)
            if not options:
                # it can join; if that assemblage will not solve, it may replace the mineral most
                # like it instead (tremolite for diopside when talc-like Mg is short)
                prev, n_before = list(assemblage), dict(n_prev)
                assemblage = assemblage + [add]
                swaps[frozenset(assemblage)] = (prev, n_before, add, None,
                                                _most_alike(model.nu_ph, prev, add, problem.unk))
            else:
                # the new mineral's saturation condition is a combination of those already present
                # (Dolomite against Dolomite-ord, gypsum against anhydrite): adding it would break
                # the phase rule, so it replaces one of them, as in a simplex pivot; if that
                # assemblage cannot be solved, the next candidate leaves instead
                prev, n_before = list(assemblage), dict(n_prev)
                leave, theta = options.pop(0)
                assemblage, n_prev = _swapped(prev, n_before, model.ph_names, lam, leave, add, theta)
                swaps[frozenset(assemblage)] = (prev, n_before, add, lam, options)
            continue
        stuck = [model.ph_names[a] for a in blocked if a in allowed_set and a not in assemblage
                 and si[a] > 1e-6]
        if stuck:                                     # never return a state that is supersaturated
            raise EquilibriumError(f"could not add the supersaturated mineral(s) {stuck}")
        if gas is not None and gas_on is None and gas_blocked and \
                gas.fractions(last["u"]).sum() > 1.0 + 1e-6:
            raise EquilibriumError("could not add the gas phase that the dissolved gases call for")
        break
    else:
        raise EquilibriumError("no stable phase assemblage found")

    spec = make_speciation(model, last["u"], last["log_aw"], last["lg_b"], last["lg_d"],
                           last["ionic"], last["sum_m"], total_iter)
    n = last["n"]
    present = {model.ph_names[a]: float(n[k]) for k, a in enumerate(assemblage)}
    open_names = tuple(names[j] for j in list(fixed_idx) + list(buffer_rows))
    gas_moles = {}
    if gas_on is not None and last.get("gas") is not None:
        gas.log_phi = last["log_phi"]
        x = gas.fractions(last["u"])
        gas_moles = {nm: float(last["gas"] * xi) for nm, xi in zip(gas.names, x) if xi > 1e-12}
    return EquilibriumState(spec, math.exp(last["ln_w"]), present, dict(bulk), total_iter,
                            open_names, gas=gas_moles, _u=last["u"])

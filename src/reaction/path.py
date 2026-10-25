"""Reaction paths: follow a system as reactants are added, temperature changes and
minerals precipitate or dissolve.

The system is *closed*: at every value of reaction progress xi its bulk composition is the
initial composition plus what the reactants have contributed, and the state is the
equilibrium of that bulk (see `equilibrium.py`).  Reactants are irreversible: their
composition is added at a rate dxi_r/dxi; they cannot precipitate while they are reacting.
A mineral reactant that becomes saturated is converted, as in EQ6: its unreacted mass
joins the equilibrium system, where it may dissolve or stay.

In a fluid-centered *flow-through* system minerals are left behind and cannot react again.
As in EQ6, a mineral stays in the equilibrium system while it grows and is moved to the
"physically removed system" at the point where it would start to redissolve: a growing solid
never back-reacts, so the fluid follows exactly the path of continuous removal, and the amount
reported for each mineral is what has formed since it was last removed.

A reactant with a rate law is *kinetic*: reaction progress is then time in seconds, and the moles it
has reacted follow d(reacted)/dt = rate(state), where the rate depends on the fluid it is in (see
`kinetics.py`).  The amounts are advanced implicitly (backward Euler: the rate is taken at the end of
each step), which is stable however fast a reaction is next to the step, so a mineral that reaches
saturation within a step simply sits at saturation.  A kinetic mineral is never part of the
equilibrium assemblage: it dissolves while undersaturated and grows on what is left of it while
supersaturated, at its rate.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from .equilibrium import EquilibriumError, EquilibriumState, System, equilibrate
from .kinetics import RateLaw
from .network import WATER
from .speciation import Speciation, SpeciationError


@dataclass
class Reactant:
    """Something added to the system as reaction progress increases.

    ``composition`` is moles of each basis component per mole of reactant.  For a pure
    mineral build it with :func:`mineral_reactant` so the saturation conversion works.
    ``rates`` = (r1, r2, r3) gives the moles reacted as r1 xi + r2 xi^2/2 + r3 xi^3/6;
    the reactant stops when ``amount`` moles have reacted.

    With a ``rate_law`` the reactant is kinetic instead: ``amount`` is how much of it there is at the
    start, ``area`` its reactive surface in m2 then, and the surface shrinks or grows with the amount
    left as (left / amount)^(2/3), as for a population of grains of one size.
    """
    name: str
    composition: dict[str, float]
    amount: float
    rates: tuple[float, float, float] = (1.0, 0.0, 0.0)
    mineral: str | None = None                # phase name if this reactant is a pure mineral
    rate_law: RateLaw | None = None
    area: float = 0.0

    @property
    def kinetic(self) -> bool:
        return self.rate_law is not None

    def surface(self, reacted: float) -> float:
        left = self.amount - reacted
        return self.area * (left / self.amount) ** (2.0 / 3.0) if left > 0 and self.amount > 0 else 0.0

    def reacted(self, xi: float) -> float:
        r1, r2, r3 = self.rates
        return min(self.amount, max(0.0, r1 * xi + r2 * xi ** 2 / 2.0 + r3 * xi ** 3 / 6.0))


@dataclass
class PathPoint:
    xi: float
    temperature_c: float
    pressure_bar: float | None
    state: EquilibriumState
    reacted: dict[str, float]                  # moles of each reactant destroyed
    transferred: tuple[str, ...] = ()          # mineral reactants now in the equilibrium system
    removed: dict[str, float] = field(default_factory=dict)   # flow-through: moles left behind so far
    rates: dict[str, float] = field(default_factory=dict)     # kinetic reactants: units per second

    @property
    def speciation(self) -> Speciation:
        return self.state.speciation

    @property
    def minerals(self) -> dict[str, float]:
        return self.state.minerals


def bulk_from_speciation(spec: Speciation, kilograms: float = 1.0) -> dict[str, float]:
    """Moles of every basis component in ``kilograms`` of solvent plus the species in it."""
    return {name: kilograms * value for name, value in spec.component_totals.items()}


def _composition_row(system: System, name: str) -> dict[str, float]:
    """Moles of each basis component in one mole of ``name``.  Stoichiometry does not depend on
    temperature or pressure, so it comes from the network itself: building an aqueous model here
    would need a pressure, and there is no default one above the critical temperature."""
    net = system.net
    if name in net.basis:                         # O2(g) left with a gas phase
        return {name: 1.0}
    if name in net.phase_names:
        row = net.nu_phase[net.phase_names.index(name)]
    elif name in net.dep_names:
        row = net.nu_dep[net.dep_names.index(name)]
    else:
        raise EquilibriumError(f"{name!r} is not a phase or species of this basis")
    return {net.basis[j]: float(c) for j, c in enumerate(row) if abs(c) > 1e-12}


def mineral_reactant(system: System, name: str, amount: float, *,
                     rates=(1.0, 0.0, 0.0), rate_law: RateLaw | None = None,
                     area: float = 0.0) -> Reactant:
    """A pure mineral (or gas, or aqueous species) as a reactant, from its formation reaction."""
    return Reactant(name, _composition_row(system, name), amount, rates, name, rate_law, area)


class ReactionPath:
    def __init__(self, system: System, bulk: dict[str, float], *,
                 reactants: Sequence[Reactant] = (),
                 temperature: float | Callable[[float], float] = 25.0,
                 pressure: float | Callable[[float], float] | None = None,
                 fixed: dict[str, float] | None = None,
                 buffers: Sequence[tuple[str, float]] = (),
                 minerals: Sequence[str] = (),
                 suppress_minerals: Sequence[str] = (),
                 flow_through: bool = False,
                 no_precipitation: bool = False,
                 saturation_tol: float = 1e-6,
                 kinetic_tol: float = 0.03,
                 separate: Callable[[float, EquilibriumState], dict[str, float] | None] | None = None):
        self.system, self.bulk0 = system, dict(bulk)
        self.reactants = list(reactants)
        self._temperature = temperature if callable(temperature) else (lambda xi, t=temperature: t)
        self._pressure = pressure if callable(pressure) or pressure is None else (
            lambda xi, p=pressure: p)
        self.fixed, self.buffers = dict(fixed or {}), tuple(buffers)
        self.initial_minerals = tuple(minerals)
        self.suppress_minerals = set(suppress_minerals)
        self.flow_through = flow_through
        self.no_precipitation = no_precipitation
        self.saturation_tol = saturation_tol
        self.kinetic_tol = kinetic_tol
        # phase separation: given (xi, state), the moles of each component of the fluid phase that leaves
        # (the vapour when the brine is followed), or None; taken away like the flow-through gas
        self.separate = separate
        self._split_off: dict[str, float] = {}
        self._transferred: set[str] = set()
        self._removed: dict[str, float] = {}
        self._removed_members: dict[str, float] = {}    # endmembers of the solid solutions left behind
        self._stopped: dict[str, float] = {}      # flow-through: saturated reactants, moles reacted
        self.kinetic = [r for r in self.reactants if r.kinetic]
        # kinetic amounts reacted at the two ends of the current step; in between, a straight line
        self._span: tuple[float, dict, float, dict] = (0.0, {}, 0.0, {})
        self._dt_hint: float | None = None

    # ------------------------------------------------------------------------------------------
    def _solve(self, xi: float, guess: EquilibriumState | None) -> tuple[EquilibriumState, dict]:
        """Equilibrium at one xi with the current set of transferred reactants."""
        t_c = self._temperature(xi)
        p_bar = self._pressure(xi) if self._pressure else None
        model = self.system.model(t_c, p_bar)
        bulk = dict(self.bulk0)
        reacted = {}
        kin = self._kinetic_reacted(xi)
        for r in self.reactants:
            if r.kinetic:
                done = kin.get(r.name, 0.0)
            else:
                done = (r.amount if r.name in self._transferred
                        else self._stopped.get(r.name, r.reacted(xi)))
            reacted[r.name] = done
            for comp, c in r.composition.items():
                bulk[comp] = bulk.get(comp, 0.0) + c * done
        ss_names = {ss.name for ss in self.system.solid_solutions}
        for comp, moles in self._split_off.items():
            bulk[comp] = bulk.get(comp, 0.0) - moles
        for name, moles in [*((k, v) for k, v in self._removed.items() if k not in ss_names),
                            *self._removed_members.items()]:
            for comp, c in _composition_row(self.system, name).items():
                bulk[comp] = bulk.get(comp, 0.0) - c * moles
        reacting = {r.mineral for r in self.reactants
                    if r.mineral and (r.kinetic or (r.name not in self._transferred
                                                    and r.name not in self._stopped))}
        excluded = reacting | self.suppress_minerals
        # a kinetic mineral's other structural states (quartz,beta; albite,low) would grow in its place
        stems = {r.mineral.split(",")[0] for r in self.kinetic if r.mineral}
        eligible = [n for n, k in zip(model.ph_names, model.ph_kind)
                    if k.value == "mineral" and n not in excluded and not self.no_precipitation
                    and n.split(",")[0] not in stems]
        present = list(guess.minerals) if guess is not None else []
        present += [n for n in self.initial_minerals if n not in present]
        converted = {r.mineral: r.amount for r in self.reactants
                     if r.mineral and r.name in self._transferred}
        present += [n for n in converted if n not in present]
        state = equilibrate(self.system, t_c, bulk, pressure_bar=p_bar, fixed=self.fixed,
                            buffers=self.buffers, minerals=[n for n in present if n in eligible],
                            amounts=converted, eligible=eligible, guess=guess)
        return state, reacted

    # ---- kinetic reactants ------------------------------------------------------------------------
    def _kinetic_reacted(self, xi: float) -> dict[str, float]:
        xa, na, xb, nb = self._span
        if not nb or xb <= xa:
            return nb or na
        f = min(max((xi - xa) / (xb - xa), 0.0), 1.0)
        return {k: na.get(k, 0.0) + f * (v - na.get(k, 0.0)) for k, v in nb.items()}

    def _kinetic_rates(self, state: EquilibriumState, reacted: dict[str, float]) -> dict[str, float]:
        sat = state.speciation.saturation
        return {r.name: r.rate_law.rate(state.speciation, r.surface(reacted.get(r.name, 0.0)),
                                        sat.get(r.mineral) if r.mineral else None)
                for r in self.kinetic}

    def _implicit(self, xi0: float, st0: EquilibriumState, n0: dict[str, float],
                  xi1: float) -> tuple[dict[str, float], EquilibriumState, dict[str, float]]:
        """Amounts reacted at xi1 from n = n0 + dt * rate(state at xi1 with n): one root in each amount
        (bracketed: the rate falls as a mineral dissolves toward saturation, and its surface shrinks),
        swept over the kinetic reactants in turn until they agree.  Also returns the explicit step
        (rate taken at the start), for the error estimate."""
        dt = xi1 - xi0
        n = dict(n0)
        last = [st0]

        def state_with(trial):
            self._span = (xi0, n0, xi1, trial)
            st = self._solve(xi1, last[0])[0]
            last[0] = st
            return st

        def residual(r, value):
            trial = {**n, r.name: value}
            return value - n0[r.name] - dt * self._kinetic_rates(state_with(trial), trial)[r.name]

        rates0 = self._kinetic_rates(st0, n0)
        explicit = {name: n0[name] + dt * rate for name, rate in rates0.items()}
        found = self._newton(n0, explicit, state_with, dt)
        if found is not None:
            return found, last[0], explicit
        sweeps = 1 if len(self.kinetic) == 1 else 8
        for sweep in range(sweeps):
            moved = False
            for r in self.kinetic:
                # minus the step at the rate where this one starts and the others are now
                untouched = all(n[k] == n0[k] for k in n if k != r.name)
                f0 = (n0[r.name] - explicit[r.name]) if untouched else residual(r, n0[r.name])
                scale = max(abs(r.amount), abs(n0[r.name]), 1e-300)
                if abs(f0) <= 1e-14 * scale:
                    new = n0[r.name]
                else:
                    far = min(n0[r.name] - f0, r.amount)    # the explicit step; never more than there is
                    f_far = None
                    for _ in range(60):                     # a step too long to solve: shorten it
                        try:
                            f_far = residual(r, far)
                            break
                        except SpeciationError:
                            far = n0[r.name] + 0.5 * (far - n0[r.name])
                    if f_far is None:
                        raise SpeciationError(f"{r.name}: no solution for its kinetic step")
                    near, f_near = n0[r.name], f0
                    for _ in range(12):                     # the rate grew along the step: go further
                        if f_far * f0 <= 0 or far >= r.amount:
                            break
                        try:
                            further = min(n0[r.name] + 2.0 * (far - n0[r.name]), r.amount)
                            f_further = residual(r, further)
                        except SpeciationError:
                            break
                        near, f_near, far, f_far = far, f_far, further, f_further
                    if f_far * f_near > 0:
                        new = far
                    else:
                        lo, hi = sorted((near, far))
                        new = brentq(lambda v, r=r: residual(r, v), lo, hi,
                                     xtol=1e-5 * (hi - lo) + 1e-15 * scale, rtol=1e-12, maxiter=60)
                change = abs(new - n0[r.name])
                moved = moved or abs(new - n[r.name]) > 1e-3 * change + 1e-14 * scale
                n[r.name] = new
            if not moved:
                break
            if len(self.kinetic) > 1:                       # close enough now for Newton to finish
                found = self._newton(n0, explicit, state_with, dt, start=n)
                if found is not None:
                    return found, last[0], explicit
        return n, state_with(n), explicit

    def _newton(self, n0: dict[str, float], explicit: dict[str, float], state_with,
                dt: float, start: dict[str, float] | None = None) -> dict[str, float] | None:
        """The implicit step by Newton's method on all kinetic amounts at once (Jacobian by
        differences), from ``start`` (no change by default); a few solves when the step is smooth, as
        it is once a fast reaction sits at saturation.  None when it does not settle, for the slower,
        safer sweep."""
        names = [r.name for r in self.kinetic]
        amount = np.array([r.amount for r in self.kinetic])
        base = np.array([n0[k] for k in names])
        guess_step = np.minimum(np.abs(np.array([explicit[k] for k in names]) - base), np.abs(amount))

        def resid(x):
            trial = dict(zip(names, x))
            rates = self._kinetic_rates(state_with(trial), trial)
            return x - base - dt * np.array([rates[k] for k in names])

        def tol(x):
            return 1e-3 * np.maximum(np.abs(x - base), 1e-6 * np.abs(x)) + 1e-14 * np.abs(amount)

        x = base.copy() if start is None else np.array([start[k] for k in names])
        try:
            f = resid(x)
            for _ in range(8):
                if np.all(np.abs(f) <= tol(x)):
                    return dict(zip(names, x))
                h = 1e-4 * np.maximum(np.maximum(np.abs(x - base), guess_step), 1e-8 * np.abs(amount))
                jac = np.empty((len(names), len(names)))
                for j in range(len(names)):
                    xp = x.copy()
                    xp[j] += h[j] if xp[j] + h[j] <= amount[j] else -h[j]
                    jac[:, j] = (resid(xp) - f) / (xp[j] - x[j])
                step = np.linalg.solve(jac, -f)
                for _ in range(5):                           # halve a step that makes things worse
                    xn = np.minimum(x + step, amount)
                    try:
                        fn = resid(xn)
                    except SpeciationError:
                        step *= 0.5
                        continue
                    if np.max(np.abs(fn) / tol(xn)) <= max(np.max(np.abs(f) / tol(x)), 1.0):
                        break
                    step *= 0.5
                else:
                    return None
                x, f = xn, fn
        except (SpeciationError, np.linalg.LinAlgError):
            return None
        return None

    def _step_error(self, n0: dict[str, float], n1: dict[str, float], explicit: dict[str, float],
                    st1: EquilibriumState) -> float:
        """The step's error over what is allowed (above 1: too long).  Half the gap between the implicit
        and explicit steps is the error; it may be ``kinetic_tol`` of what reacted in the step plus a
        third of that of all that has reacted so far.  Next to saturation the implicit step is right
        whatever its length, so a reaction there does not count."""
        sat = st1.speciation.saturation
        worst = 0.0
        for r in self.kinetic:
            if r.mineral and abs(sat.get(r.mineral, -99.0)) < 0.02:
                continue
            allowed = self.kinetic_tol * (abs(n1[r.name] - n0[r.name]) + 0.3 * abs(n1[r.name]))
            worst = max(worst, 0.5 * abs(n1[r.name] - explicit[r.name]) / (allowed + 1e-14 * r.amount))
        return worst

    def _kinetic_step(self, xi0: float, st0: EquilibriumState, n0: dict[str, float],
                      xi1: float) -> tuple[EquilibriumState, dict, dict[str, float]]:
        """State and kinetic amounts at xi1, in as many implicit steps as the error allows: a step is
        halved when its error is too large and the next one doubled when it is well inside, starting
        from the last step that worked (or, the first time, from a ten-thousandth of the way, shorter
        if a reactant would change by more than a percent at its starting rate in that time)."""
        span = xi1 - xi0
        h = self._dt_hint
        if h is None:
            rates0 = self._kinetic_rates(st0, n0)
            h = min([0.01 * abs(r.amount) / abs(rates0[r.name]) for r in self.kinetic
                     if rates0[r.name] != 0.0] + [1e-4 * span])
        x, st, n = xi0, st0, n0
        reacted: dict = {}
        while x < xi1 - 1e-12 * span:
            h = min(h, xi1 - x)
            last_try = h <= 1e-9 * span
            try:
                n1, st1, explicit = self._implicit(x, st, n, x + h)
                err = 0.0 if last_try else self._step_error(n, n1, explicit, st1)
            except SpeciationError:
                if last_try:
                    raise
                n1, err = None, math.inf
            if err > 1.0:
                h *= 0.5
                continue
            self._span = (x, n, x + h, n1)
            st, reacted = self._advance(x, st, x + h)
            x, n = x + h, n1
            if err < 0.25:
                h *= 2.0
        self._dt_hint = h
        return st, reacted, n

    def _peak(self, mineral: str, xi0: float, st0: EquilibriumState, xi1: float) -> float:
        """Reaction progress in [xi0, xi1] where the amount of ``mineral`` is largest (Brent's method:
        the amount rises, then falls, smoothly or with a kink where another mineral appears).  A point
        the solver cannot settle counts as no mineral, so the search moves away from it."""
        seen = [(xi0, st0)]

        def loss(xi):                             # warm start from the nearest point solved so far
            guess = min(seen, key=lambda pair: abs(pair[0] - xi))[1]
            try:
                state = self._solve(xi, guess)[0]
            except SpeciationError:
                return 0.0
            seen.append((xi, state))
            return -state.solids.get(mineral, 0.0)

        start = st0.solids.get(mineral, 0.0)
        tol = 1e-6 * max(xi1 - xi0, 1e-300)
        if -loss(xi0 + tol) <= start:                 # already falling where the step begins
            return xi0
        found = minimize_scalar(loss, bounds=(xi0, xi1), method="bounded",
                                options={"xatol": tol, "maxiter": 200})
        return xi0 if start >= -found.fun else float(found.x)

    def _saturated(self, state: EquilibriumState) -> list[str]:
        sat = state.speciation.saturation
        return [r.name for r in self.reactants
                if r.mineral and not r.kinetic and r.name not in self._transferred and r.name not in self._stopped
                and sat.get(r.mineral, -1e9) > -self.saturation_tol]

    def _advance(self, xi0: float | None, st0: EquilibriumState | None, xi1: float,
                 depth: int = 0) -> tuple[EquilibriumState, dict]:
        """State at xi1, starting from the state st0 at xi0; handles reactant conversions."""
        try:
            st1, reacted = self._solve(xi1, st0)
        except SpeciationError:
            # halve the step a limited number of times per point: where the chemistry has no answer,
            # halving both halves again and again would take 2**depth solves
            self._halvings += 1
            if xi0 is None or depth > 25 or self._halvings > 40:
                raise
            mid = 0.5 * (xi0 + xi1)
            st_mid, _ = self._advance(xi0, st0, mid, depth + 1)
            return self._advance(mid, st_mid, xi1, depth + 1)
        if self.flow_through and st0 is not None and depth <= 200:
            # amounts below the solver's resolution (a mineral just removed and still listed) are noise
            floor = 1e-12 * max([abs(v) for k, v in st0.bulk.items() if k != WATER] + [1e-30])
            before, after = st0.solids, st1.solids
            shrinking = [m for m, n in before.items()
                         if n > floor and after.get(m, 0.0) < n * (1.0 - 1e-9) - floor]
            if shrinking:
                # remove the first mineral to start redissolving, where its amount peaks
                peak, name = min((self._peak(m, xi0, st0, xi1), m) for m in shrinking)
                st_p = st0 if peak <= xi0 else self._solve(peak, st0)[0]
                self._removed[name] = self._removed.get(name, 0.0) + st_p.solids.get(name, 0.0)
                for member, moles in st_p.solid_solutions.get(name, {}).items():
                    self._removed_members[member] = self._removed_members.get(member, 0.0) + moles
                st_p, _ = self._solve(peak, st_p)
                if peak >= xi1:
                    return st_p, self._solve(xi1, st_p)[1]
                return self._advance(peak, st_p, xi1, depth + 1)
        newly = self._saturated(st1)
        if not newly:
            return st1, reacted
        # a reacting mineral became saturated somewhere in (xi0, xi1]: locate that point
        if xi0 is None:
            boundary, st_b = xi1, st1
        else:
            st_lo = st0
            if self._saturated(st_lo):                      # already saturated at xi0
                boundary, st_b = xi0, st_lo
            else:
                lo, hi, st_hi = xi0, xi1, st1
                for _ in range(80):
                    if hi - lo <= 1e-10 * max(abs(hi), 1e-300):
                        break
                    mid = 0.5 * (lo + hi)
                    st_mid, _ = self._solve(mid, st_lo)
                    if self._saturated(st_mid):
                        hi, st_hi = mid, st_mid
                    else:
                        lo, st_lo = mid, st_mid
                boundary, st_b = hi, st_hi
        if self.flow_through:
            # the fluid moves on and leaves the unreacted mineral behind: it stops reacting and
            # what is left of it never joins the equilibrium system (EQ6 does the same)
            for name in self._saturated(st_b):
                r = next(r for r in self.reactants if r.name == name)
                self._stopped[name] = r.reacted(boundary)
        else:
            self._transferred.update(self._saturated(st_b))
        st_conv, reacted = self._solve(boundary, st_b)          # same solution, now with the solid
        if boundary >= xi1:
            return st_conv, reacted
        return self._advance(boundary, st_conv, xi1, depth + 1)

    def run(self, xis: Sequence[float], *,
            on_failure: Callable[[float, Exception], None] | None = None,
            progress: Callable[[int, int], None] | None = None) -> list[PathPoint | None]:
        """Solve the system at each reaction-progress value (ascending).

        Without ``on_failure`` a point where no equilibrium is found raises.  With it, the failure
        is reported there, that point is None and the path goes on from the last solved point."""
        self._transferred, self._removed, self._stopped = set(), {}, {}
        self._removed_members, self._split_off = {}, {}
        n_kin = {r.name: 0.0 for r in self.kinetic}
        self._span = (0.0, n_kin, 0.0, n_kin)
        self._dt_hint = None
        self._halvings = 0
        points: list[PathPoint | None] = []
        xi_prev, state = None, None
        for i, xi in enumerate(xis):
            if progress:
                progress(i, len(xis))
            self._halvings = 0
            try:
                if self.kinetic and state is not None:
                    new, reacted, n_kin = self._kinetic_step(xi_prev, state, n_kin, xi)
                else:
                    new, reacted = self._advance(xi_prev, state, xi)
            except SpeciationError as exc:
                if on_failure is None:
                    raise
                on_failure(xi, exc)
                points.append(None)
                continue
            if self.flow_through and new.gas:
                # bubbles leave with the fluid's flow, like the minerals: the gas phase is vented
                for name, moles in new.gas.items():
                    self._removed[name] = self._removed.get(name, 0.0) + moles
                new, reacted = self._solve(xi, new)
            leaving = self.separate(xi, new) if self.separate is not None else None
            if leaving:
                for comp, moles in leaving.items():
                    self._split_off[comp] = self._split_off.get(comp, 0.0) + moles
                try:
                    try:
                        new, reacted = self._solve(xi, new)
                    except SpeciationError:           # much of the fluid has gone: start afresh
                        new, reacted = self._solve(xi, None)
                except SpeciationError as exc:
                    if on_failure is None:
                        raise
                    on_failure(xi, exc)
                    points.append(None)
                    continue
            state, xi_prev = new, xi
            points.append(PathPoint(xi, self._temperature(xi), state.speciation.pressure_bar,
                                    state, reacted, tuple(sorted(self._transferred)),
                                    dict(self._removed), self._kinetic_rates(state, n_kin)))
        if progress:
            progress(len(xis), len(xis))
        return points

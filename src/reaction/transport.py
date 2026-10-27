"""One-dimensional reactive transport: a column of cells that a fluid flows through.

The scheme is PHREEQC's TRANSPORT with a Courant number of one.  Each shift moves the solution of every
cell (its water and everything dissolved in it) one cell downstream; the inflowing fluid enters the first
cell and the last cell's solution leaves at the outlet.  Minerals, the gas phase and sorbed matter stay
where they are.  Dispersion then mixes each cell's solution with its neighbours' by the fraction
alpha / dx per shift (alpha the dispersivity, dx the cell length), which is D dt / dx^2 for D = alpha v;
the shift itself adds no numerical dispersion.  Mixing fractions above 1/3 are taken in several smaller
mixes, as PHREEQC does, so the scheme stays stable.  Then every cell is brought to equilibrium again at
its own temperature and pressure, so a column can sit on a thermal gradient.

The upstream boundary holds the inflow's composition (a constant-concentration inlet) and the
downstream boundary has no gradient (what reaches the last cell leaves).

Minerals with a rate law (``kinetic`` reactants, see kinetics.py) are not at equilibrium: each cell keeps
its own amount of each, and after every shift the cell reacts for the shift's time ``dt`` as a kinetic
reaction path does (path.py: implicit steps with an error test), the rest of the cell staying at
equilibrium.  A kinetic mineral dissolves while the fluid is undersaturated in it and grows on what is left
of it while the fluid is supersaturated; it never precipitates at equilibrium, nor do its other structural
states.  This is the operator splitting of PHREEQC's TRANSPORT with KINETICS.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Sequence

from .equilibrium import EquilibriumError, EquilibriumState, System, equilibrate
from .path import Reactant, ReactionPath
from .speciation import SpeciationError


@dataclass
class TransportResult:
    pore_volumes: list[float]                         # at each recorded outlet state
    outlet: list[EquilibriumState | None]             # the last cell after each shift (None: not solved)
    profile: list[EquilibriumState | None]            # every cell at the end
    snapshots: dict[int, list[EquilibriumState | None]] = field(default_factory=dict)   # shift -> cells
    failures: list[tuple[int, int, str]] = field(default_factory=list)                    # (shift, cell, why)
    out: dict[str, float] = field(default_factory=dict)       # moles of each component that left at the outlet
    back: dict[str, float] = field(default_factory=dict)      # and that left upstream by dispersion (< 0: came in)
    # kinetic reactants: (moles left, rate in units per second) by name, alongside outlet, profile, snapshots
    outlet_kinetic: list[tuple[dict, dict]] = field(default_factory=list)
    profile_kinetic: list[tuple[dict, dict]] = field(default_factory=list)
    snapshot_kinetic: dict[int, list[tuple[dict, dict]]] = field(default_factory=dict)


def mixing_steps(fraction: float) -> tuple[int, float]:
    """How many mixes a shift's dispersion takes, and the fraction of each (at most 1/3)."""
    if fraction <= 0:
        return 0, 0.0
    k = max(1, math.ceil(3.0 * fraction - 1e-12))
    return k, fraction / k


class Column:
    """``cells`` gives (temperature C, pressure bar) of each cell, upstream first.  ``fill`` is each
    cell's starting bulk composition (one dict for all, or one per cell), with any minerals and sorbed
    load already in it; ``inflow`` the moles of each component in the solution that enters per shift
    (its water included).  ``minerals`` starts every cell's assemblage; ``eligible`` limits the minerals
    that may form (default: all that ``suppress`` does not name).  ``kinetic`` reactants (with a rate law)
    sit in every cell in their ``amount`` (not part of ``fill``) and react for ``dt`` seconds per shift, with
    ``kinetic_tol`` the error allowed in each step (see ReactionPath; about 1% of what dissolves at 0.05)."""

    def __init__(self, system: System, cells: Sequence[tuple[float, float]],
                 fill: dict[str, float] | Sequence[dict[str, float]], inflow: dict[str, float], *,
                 dispersion: float = 0.0, minerals: Sequence[str] = (),
                 amounts: dict[str, float] | None = None, suppress: Sequence[str] = (),
                 guess=None, kinetic: Sequence[Reactant] = (), dt: float = 0.0, kinetic_tol: float = 0.05):
        self.system = system
        self.cells = [(float(t), float(p)) for t, p in cells]
        n = len(self.cells)
        if n < 1:
            raise ValueError("a column needs at least one cell")
        fills = [dict(fill) for _ in range(n)] if isinstance(fill, dict) else [dict(f) for f in fill]
        if len(fills) != n:
            raise ValueError("one starting composition per cell")
        self.bulk = fills
        self.inflow = dict(inflow)
        self.dispersion = float(dispersion)
        self.minerals = tuple(minerals)
        self.amounts = dict(amounts or {})
        self.suppress = {s.lower() for s in suppress}
        self.guess = guess
        self.kinetic = [r for r in kinetic if r.kinetic]
        if self.kinetic and not dt > 0:
            raise ValueError("kinetic reactants need the time of a shift")
        self.dt = float(dt)
        self.kinetic_tol = kinetic_tol
        self.left = [{r.name: r.amount for r in self.kinetic} for _ in range(n)]
        self.rates = [{r.name: 0.0 for r in self.kinetic} for _ in range(n)]
        self._hints: list[float | None] = [self.dt] * n      # try the whole shift in one step first
        # a kinetic mineral and its other structural states (quartz,beta; albite,low) never form at equilibrium
        self._stems = {r.mineral.split(",")[0].lower() for r in self.kinetic if r.mineral}
        self.states: list[EquilibriumState | None] = [None] * n
        self.out: dict[str, float] = {}
        self.back: dict[str, float] = {}

    # ------------------------------------------------------------------------------------------
    def _eligible(self, t_c: float, p_bar: float) -> list[str]:
        model = self.system.model(t_c, p_bar)
        return [nm for nm, k in zip(model.ph_names, model.ph_kind)
                if k.value == "mineral" and nm.lower() not in self.suppress
                and nm.split(",")[0].lower() not in self._stems]

    def _solve(self, i: int, guess) -> EquilibriumState:
        t_c, p_bar = self.cells[i]
        eligible = self._eligible(t_c, p_bar)
        present = list(guess.minerals) if isinstance(guess, EquilibriumState) else []
        present += [m for m in self.minerals if m not in present]
        tries = [guess] + [g for g in (self.guess, None) if g is not guess]
        for k, start in enumerate(tries):     # a warm start that leads astray: the starting guess, then none
            try:
                return equilibrate(self.system, t_c, self.bulk[i], pressure_bar=p_bar,
                                   minerals=[m for m in present if m in eligible], amounts=self.amounts,
                                   eligible=eligible, guess=start)
            except (EquilibriumError, SpeciationError):
                if k == len(tries) - 1:
                    raise
            present = list(self.minerals)
        raise AssertionError("unreachable")

    def _solution(self, i: int) -> dict[str, float]:
        """Moles of each component in cell i's solution (the part that moves)."""
        st = self.states[i]
        return {c: st.moles_in_solution(c) for c in self.bulk[i] if c in st.speciation.component_totals}

    def _split(self) -> tuple[list[dict], list[dict]]:
        """Each cell's bulk as (solution, everything that stays: minerals, gas, sorbed)."""
        n = len(self.cells)
        sol = [self._solution(i) if self.states[i] is not None else {} for i in range(n)]
        stay = [{c: v - sol[i].get(c, 0.0) for c, v in self.bulk[i].items()} for i in range(n)]
        return sol, stay

    def _set(self, sol: list[dict], stay: list[dict]) -> None:
        for i in range(len(self.cells)):
            self.bulk[i] = {c: stay[i].get(c, 0.0) + sol[i].get(c, 0.0) for c in set(stay[i]) | set(sol[i])}

    def _shift(self) -> None:
        sol, stay = self._split()
        for c, v in sol[-1].items():
            self.out[c] = self.out.get(c, 0.0) + v
        self._set([dict(self.inflow)] + sol[:-1], stay)          # one cell downstream

    def _mix(self, f: float) -> None:
        sol, stay = self._split()
        up = [dict(self.inflow)] + sol[:-1]
        down = sol[1:] + [sol[-1]]
        for c in set(sol[0]) | set(self.inflow):                 # dispersion across the inlet
            self.back[c] = self.back.get(c, 0.0) - f * (self.inflow.get(c, 0.0) - sol[0].get(c, 0.0))
        mixed = [{c: sol[i].get(c, 0.0) + f * (up[i].get(c, 0.0) - 2.0 * sol[i].get(c, 0.0) + down[i].get(c, 0.0))
                  for c in set(sol[i]) | set(up[i]) | set(down[i])} for i in range(len(self.cells))]
        self._set(mixed, stay)

    def _react(self, s: int, failures: list) -> None:
        # a cell that will not solve keeps its last state; its matter is still counted in its bulk
        for i in range(len(self.cells)):
            self._react_one(s, i, failures)

    def _react_kinetic(self, s: int, failures: list) -> None:
        """Every cell reacts for one shift's time: its kinetic minerals at their rates, the rest at
        equilibrium.  A cell whose kinetic step will not solve is brought to equilibrium without it."""
        for i, (t_c, p_bar) in enumerate(self.cells):
            reactants = []
            for r in self.kinetic:
                left = self.left[i][r.name]
                if left > 1e-12 * r.amount:
                    reactants.append(Reactant(r.name, r.composition, left, mineral=r.mineral,
                                              rate_law=r.rate_law, area=r.surface(r.amount - left)))
            if not reactants:
                self.rates[i] = {r.name: 0.0 for r in self.kinetic}
                self._react_one(s, i, failures)
                continue
            model = self.system.model(t_c, p_bar)
            path = ReactionPath(self.system, self.bulk[i], reactants=reactants, temperature=t_c, pressure=p_bar,
                                suppress_minerals=[nm for nm in model.ph_names if nm.lower() in self.suppress
                                                   or nm.split(",")[0].lower() in self._stems],
                                kinetic_tol=self.kinetic_tol)
            prev = self.states[i]
            try:
                point = path.run([0.0, self.dt], guess=prev if prev is not None else self.guess,
                                 step_hint=self._hints[i])[-1]
            except (EquilibriumError, SpeciationError) as exc:
                failures.append((s, i, f"kinetic step: {exc}"))
                self._hints[i] = self.dt
                self._react_one(s, i, failures)
                continue
            for r in reactants:
                done = point.reacted.get(r.name, 0.0)
                self.left[i][r.name] = max(self.left[i][r.name] - done, 0.0)
                for comp, c in r.composition.items():
                    self.bulk[i][comp] = self.bulk[i].get(comp, 0.0) + c * done
            self.rates[i] = {r.name: point.rates.get(r.name, 0.0) for r in self.kinetic}
            self.states[i] = point.state
            self._hints[i] = path.step_hint

    def _react_one(self, s: int, i: int, failures: list) -> None:
        prev = self.states[i]
        try:
            self.states[i] = self._solve(i, prev if prev is not None else self.guess)
        except (EquilibriumError, SpeciationError) as exc:
            failures.append((s, i, str(exc)))
            if prev is None:
                raise

    def _kinetic_now(self) -> list[tuple[dict, dict]]:
        return [(dict(left), dict(rates)) for left, rates in zip(self.left, self.rates)]

    # ------------------------------------------------------------------------------------------
    def run(self, shifts: int, *, snapshots: Sequence[int] = (),
            progress: Callable[[int, int], None] | None = None) -> TransportResult:
        n = len(self.cells)
        failures: list[tuple[int, int, str]] = []
        for i in range(n):
            self.states[i] = self._solve(i, self.states[i - 1] if i else self.guess)
        outlet = [self.states[-1]]
        out_kin = [self._kinetic_now()[-1]]
        pv = [0.0]
        snaps = {0: list(self.states)} if 0 in snapshots else {}
        snap_kin = {0: self._kinetic_now()} if 0 in snapshots else {}
        k_mix, f_mix = mixing_steps(self.dispersion)
        for s in range(1, shifts + 1):
            for _ in range(k_mix // 2):                         # half the dispersion before the shift
                self._mix(f_mix)
                self._react(s, failures)
            self._shift()
            if self.kinetic:
                self._react_kinetic(s, failures)
            else:
                self._react(s, failures)
            for _ in range(k_mix - k_mix // 2):                 # and half after it
                self._mix(f_mix)
                self._react(s, failures)
            if progress:
                progress(s, shifts)
            outlet.append(self.states[-1])
            out_kin.append(self._kinetic_now()[-1])
            pv.append(s / n)
            if s in snapshots:
                snaps[s] = list(self.states)
                snap_kin[s] = self._kinetic_now()
        return TransportResult(pv, outlet, list(self.states), snaps, failures, dict(self.out), dict(self.back),
                               out_kin, self._kinetic_now(), snap_kin)

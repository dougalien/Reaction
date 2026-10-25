"""Reaction paths: follow a system as reactants are added, temperature changes and
minerals precipitate or dissolve.

The system is *closed*: at every value of reaction progress xi its bulk composition is the
initial composition plus what the reactants have contributed, and the state is the
equilibrium of that bulk (see `equilibrium.py`).  Reactants are irreversible: their
composition is added at a rate dxi_r/dxi; they cannot precipitate while they are reacting.
A mineral reactant that becomes saturated is converted, as in EQ6: its unreacted mass
joins the equilibrium system, where it may dissolve or stay.

In a fluid-centered *flow-through* system the minerals formed at each step are left
behind and cannot react again.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

from .equilibrium import EquilibriumError, EquilibriumState, System, equilibrate
from .speciation import Speciation, SpeciationError


@dataclass
class Reactant:
    """Something added to the system as reaction progress increases.

    ``composition`` is moles of each basis component per mole of reactant.  For a pure
    mineral build it with :func:`mineral_reactant` so the saturation conversion works.
    ``rates`` = (r1, r2, r3) gives the moles reacted as r1 xi + r2 xi^2/2 + r3 xi^3/6;
    the reactant stops when ``amount`` moles have reacted.
    """
    name: str
    composition: dict[str, float]
    amount: float
    rates: tuple[float, float, float] = (1.0, 0.0, 0.0)
    mineral: str | None = None                # phase name if this reactant is a pure mineral

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
    if name in net.phase_names:
        row = net.nu_phase[net.phase_names.index(name)]
    elif name in net.dep_names:
        row = net.nu_dep[net.dep_names.index(name)]
    else:
        raise EquilibriumError(f"{name!r} is not a phase or species of this basis")
    return {net.basis[j]: float(c) for j, c in enumerate(row) if abs(c) > 1e-12}


def mineral_reactant(system: System, name: str, amount: float, *,
                     rates=(1.0, 0.0, 0.0)) -> Reactant:
    """A pure mineral (or gas, or aqueous species) as a reactant, from its formation reaction."""
    return Reactant(name, _composition_row(system, name), amount, rates, name)


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
                 saturation_tol: float = 1e-6):
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
        self._transferred: set[str] = set()

    # ------------------------------------------------------------------------------------------
    def _solve(self, xi: float, guess: EquilibriumState | None) -> tuple[EquilibriumState, dict]:
        """Equilibrium at one xi with the current set of transferred reactants."""
        t_c = self._temperature(xi)
        p_bar = self._pressure(xi) if self._pressure else None
        model = self.system.model(t_c, p_bar)
        bulk = dict(self.bulk0)
        reacted = {}
        for r in self.reactants:
            done = r.amount if r.name in self._transferred else r.reacted(xi)
            reacted[r.name] = done
            for comp, c in r.composition.items():
                bulk[comp] = bulk.get(comp, 0.0) + c * done
        reacting = {r.mineral for r in self.reactants if r.mineral and r.name not in self._transferred}
        excluded = reacting | self.suppress_minerals
        eligible = [n for n, k in zip(model.ph_names, model.ph_kind)
                    if k.value == "mineral" and n not in excluded and not self.no_precipitation]
        present = list(guess.minerals) if guess is not None else []
        present += [n for n in self.initial_minerals if n not in present]
        converted = {r.mineral: r.amount for r in self.reactants
                     if r.mineral and r.name in self._transferred}
        present += [n for n in converted if n not in present]
        state = equilibrate(self.system, t_c, bulk, pressure_bar=p_bar, fixed=self.fixed,
                            buffers=self.buffers, minerals=[n for n in present if n in eligible],
                            amounts=converted, eligible=eligible, guess=guess)
        return state, reacted

    def _saturated(self, state: EquilibriumState) -> list[str]:
        sat = state.speciation.saturation
        return [r.name for r in self.reactants
                if r.mineral and r.name not in self._transferred
                and sat.get(r.mineral, -1e9) > -self.saturation_tol]

    def _advance(self, xi0: float | None, st0: EquilibriumState | None, xi1: float,
                 depth: int = 0) -> tuple[EquilibriumState, dict]:
        """State at xi1, starting from the state st0 at xi0; handles reactant conversions."""
        try:
            st1, reacted = self._solve(xi1, st0)
        except SpeciationError:
            if xi0 is None or depth > 25:
                raise
            mid = 0.5 * (xi0 + xi1)
            st_mid, _ = self._advance(xi0, st0, mid, depth + 1)
            return self._advance(mid, st_mid, xi1, depth + 1)
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
        self._transferred.update(self._saturated(st_b))
        st_conv, reacted = self._solve(boundary, st_b)          # same solution, now with the solid
        if boundary >= xi1:
            return st_conv, reacted
        return self._advance(boundary, st_conv, xi1, depth + 1)

    def run(self, xis: Sequence[float]) -> list[PathPoint]:
        """Solve the system at each reaction-progress value (ascending)."""
        self._transferred = set()
        points: list[PathPoint] = []
        xi_prev, state = None, None
        for xi in xis:
            state, reacted = self._advance(xi_prev, state, xi)
            points.append(PathPoint(xi, self._temperature(xi), state.speciation.pressure_bar,
                                    state, reacted, tuple(sorted(self._transferred))))
            xi_prev = xi
        return points

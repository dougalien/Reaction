"""Run a scenario card on the equilibrium engine and package the numbers for the results window.

A run takes a fluid along a straight line in temperature and pressure and, at every point, finds the
equilibrium of the whole closed system: the aqueous species, plus whichever minerals are stable.  The
equilibrium constants are computed from Gibbs energies at each (T, P), so pressure is as free as
temperature: there is no fixed-pressure database anywhere.  In *static* mode minerals stay and can
redissolve.  In *flow-through* mode minerals are left behind as they form and never redissolve, so only
the fluid goes on; the run reports everything left behind so far.
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import multiprocessing
import os
import re
import threading
import time
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable

import numpy as np

from ..constants import WATER_MOLES_PER_KG
from ..database import load_data0
from ..database.model import Kind
from ..equilibrium import System
from ..kinetics import basalt_glass, rate_law
from ..network import PROTON, WATER
from ..path import PathPoint, ReactionPath, Reactant, _composition_row, mineral_reactant
from ..phasesep import phase_state, salinity_from_chloride
from ..solidsolution import table as solid_solution_table
from ..sorption import Sorption, exchanger_loading
from ..speciation import PH, Couple, Eh, LogFO2, Molal, Pe, Speciation, speciate
from ..thermo import validity
from ..transport import Column
from ..thermo.eosdb import DEFAULT_MASTERS, load_eos_database
from .. import isotopes
from . import analysis, config, rocks
from .fluids import AIR_LOG_FO2, preset
from .names import pretty_mineral, pretty_species
from .scenario import ALWAYS_EXCLUDED, MODE_LABELS, Scenario

Progress = Callable[[float, str], None]

TRACE = 1e-12          # mol per kg of water of a component only a reactant brings in
BAND_PERCENTILES = (5.0, 95.0)            # the shaded range: the middle 90% of the uncertainty runs
BAND_SEED = 20261025                      # fixed, so the same card gives the same range
AMOUNT_GROUPS = ("species", "elements", "minerals", "gas", "sorbed")      # a missing key there means none
AXES = {"temperature_c": ("Temperature", "°C"), "pressure_bar": ("Pressure", "bar")}


class _Evaporation(Reactant):
    """Water removed so that the concentration factor grows by the same ratio at every step: after
    step xi of ``steps``, 1 - factor^(-xi/steps) kg of the starting kilogram is gone."""

    def __init__(self, name, composition, amount, *, factor: float, steps: int):
        super().__init__(name, composition, amount)
        self.factor, self.steps = factor, steps

    def reacted(self, xi: float) -> float:
        return min(self.amount, max(0.0, 1.0 - self.factor ** (-xi / self.steps)))


# starting isotope values when the card gives none: seawater as measured, the rest typical values to replace
ISOTOPE_FLUIDS = {
    "seawater": {"d18O": 0.0, "dD": 0.0, "d13C": 0.0, "d34S": 21.0, "Sr": 7.8, "87Sr/86Sr": 0.70918},
    "vent_fluid": {"d18O": 1.0, "dD": 1.0, "d13C": -4.0, "d34S": 4.0, "Sr": 7.8, "87Sr/86Sr": 0.7035},
    "formation_brine": {"d18O": 2.0, "dD": -20.0, "d13C": -5.0, "d34S": 20.0, "Sr": 600.0, "87Sr/86Sr": 0.7100},
    "rainwater": {"d18O": -8.0, "dD": -54.0, "d13C": -9.0, "d34S": 5.0, "Sr": 0.002, "87Sr/86Sr": 0.7100},
}
ISOTOPE_ROCKS = {"basalt_glass": {"d18O": 5.7, "dD": -70.0, "d13C": -6.0, "d34S": -0.9, "Sr": 115.0,
                                  "87Sr/86Sr": 0.7027}}


class EngineUnavailable(RuntimeError):
    """The thermodynamic data the engine needs could not be found."""


def _sig(x: float | None, digits: int = 6):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return None
    return float(f"{x:.{digits}g}")


class Engine:
    def __init__(self, data_dir: str | Path | None = None):
        self.data_dir = Path(data_dir) if data_dir else config.external_dir()
        self._db = None
        self._load_lock = threading.Lock()
        self.run_lock = threading.Lock()
        self._advice: dict[tuple, list[str]] = {}
        self._prepared = None                    # (card, _prepare's result) of the last uncertainty run
        self._pool, self._pool_size = None, 0

    # ---- data -----------------------------------------------------------------------------------
    @property
    def db(self):
        with self._load_lock:
            if self._db is None:
                obigt, eq36 = self.data_dir / "OBIGT", self.data_dir / "EQ36_DB" / "data0.cmp"
                if not (obigt / "inorganic_aq.csv").exists():
                    raise EngineUnavailable(
                        f"Thermodynamic data not found in {self.data_dir}. Expected OBIGT/ with the "
                        f"CHNOSZ csv files; scripts/fetch_references.ps1 downloads them, or set "
                        f"REACTION_DATA_DIR.")
                ion_sizes = load_data0(eq36) if eq36.exists() else None
                self._db = load_eos_database(obigt, "current", eq36_database=ion_sizes)
            return self._db

    def status(self) -> dict:
        try:
            db = self.db
        except EngineUnavailable as exc:
            return {"ready": False, "message": str(exc)}
        return {"ready": True, "message": f"{db.title}, {len(db)} species"}

    # ---- how far to trust the thermodynamics along a path ---------------------------------------
    def reliability(self, points: list[tuple[float, float]]) -> list[validity.Reliability]:
        db = self.db
        out = []
        for t, p in points:
            try:
                out.append(validity.water_ionization_check(db, t, p))
            except Exception:                          # outside the liquid field of water
                out.append(validity.Reliability(t, p, float("nan"), float("nan"), "poor"))
        return out

    def advisories(self, sc: Scenario) -> list[str]:
        """Warnings about a card's path that do not stop a run, cached per path."""
        if sc.problems():
            return []
        key = (sc.t_start_c, sc.t_end_c, sc.p_start_bar, sc.p_end_bar, sc.mix_with, sc.mix_ratio, sc.mix_t_c)
        if key not in self._advice:
            if len(self._advice) > 64:
                self._advice.clear()
            pts = [sc.tp(i / 40) for i in range(41)]
            self._advice[key] = reliability_notes(pts, self.reliability(pts), sc.varies)
        return list(getattr(self.db, "data_notes", [])) + self._advice[key]

    # ---- one run --------------------------------------------------------------------------------
    def run(self, sc: Scenario, progress: Progress | None = None) -> dict:
        problems = sc.problems()
        if problems:
            raise ValueError("; ".join(problems))
        progress = progress or (lambda f, m: None)
        if sc.transport:
            return self._transport(sc, progress)
        started = time.time()
        db = self.db
        fluid = sc.fluid_spec()
        sp0, fluid_notes, added, mix, system, bulk0 = self._prepare(sc, progress)
        points = sc.points()
        modes = ["static", "flow_through"] if sc.mode == "both" else [sc.mode]
        rels = self.reliability(points)

        notes: list[str] = list(fluid_notes)
        runs: dict[str, dict] = {}
        share = 1.0 / (1 + sc.uncertainty_runs)          # of the progress bar, for the main run
        for k, mode in enumerate(modes):
            def tick(i, n, k=k):
                progress(share * (k + i / n) / len(modes), f"{MODE_LABELS[modes[k]]}: {i} of {n} points")
            runs[mode] = self._march(system, bulk0, points, sc, mode == "flow_through", tick, notes, mix)
        if sc.uncertainty_runs:
            notes.extend(self._uncertainty(sc, modes, runs, lambda f, m: progress(share + (1.0 - share) * f, m)))

        for run in runs.values():
            notes.extend(run.pop("warnings"))
        top = max((v for run in runs.values() for v in run["properties"]["ionic_strength"] if v is not None),
                  default=0.0)
        if sc.activity_model == "pitzer":
            notes.append("Activity coefficients from Pitzer's equations (PHREEQC pitzer.dat parameters); the "
                         "complexes of Na, K, Ca, Mg, Cl, sulfate and carbonate are in those parameters. The "
                         "equilibrium constants are still the database's (HKF), so a solubility can differ "
                         "from a Pitzer-only database by its log K.")
        elif top > 1.0:
            notes.append(f"Ionic strength reaches {top:.1f} mol/kg. The {sc.activity_model.replace('bdot', 'B-dot')} "
                         f"equation is unreliable above about 1 mol/kg; choose Pitzer on the card for brines.")
        for mode, run in runs.items():
            _phase_split(run, points)
            notes.extend(phase_notes(run, points, MODE_LABELS[mode] if len(runs) > 1 else "",
                                     sc.follow_phase))
        if sc.exchange_meq > 0 or sc.hfo_g > 0:
            notes.append("Sorption uses PHREEQC's constants (phreeqc.dat): Gaines-Thomas exchange, the exchange species "
                         "taking their cation's activity coefficient, and the Dzombak and Morel (1990) diffuse-layer model "
                         "for hydrous ferric oxide (89 g/mol Fe, 600 m2/g, 0.2 weak and 0.005 strong sites per Fe). "
                         "The constants are for 25 °C; exchange constants with an enthalpy follow van 't Hoff, the rest "
                         "are held at their 25 °C values.")
        if sc.solid_solutions:
            notes.append("Solid solutions mix ideally on their sites (the EQ3/6 models). Real carbonates and "
                         "feldspars unmix at low temperature, so there the minor endmember is overstated.")
        notes = [n for n in notes if not n.startswith("water density")]       # covered by reliability
        notes.extend(reliability_notes(points, rels, sc.varies))
        start_elements: dict[str, float] = {}
        for name, st in sp0.species.items():
            sp = db.species.get(name)
            for el, k in (sp.elements.items() if sp is not None else ()):
                if el not in ("H", "O"):
                    start_elements[el] = start_elements.get(el, 0.0) + k * st.molality
        names = _catalog(db, runs, sc)
        known = {n.lower() for n in db.species}
        blocked_unknown = [b for b in sc.blocked_minerals if b.lower() not in known
                           and b.lower() not in ("dolomite,ordered", "dolomite,disordered")]  # old defaults
        if blocked_unknown:
            notes.append("Not in the database, so ignored: " + ", ".join(blocked_unknown))
        notes.extend(fluid.caveats)
        if sc.mixing:
            notes.extend(c for c in preset(sc.mix_with).caveats if c not in fluid.caveats)
        seen: list[str] = []
        for n in notes:
            if n not in seen:
                seen.append(n)
        progress(1.0, "Done")

        axes = {key: {"key": key, "label": label, "unit": unit,
                      "values": [_sig(pt[0] if key == "temperature_c" else pt[1]) for pt in points]}
                for key, (label, unit) in AXES.items()}
        varies = sc.varies
        if sc.kinetics:
            last = max(len(points) - 1, 1)
            axes["time"] = {"key": "time", "label": "Time", "unit": sc.duration_unit,
                            "values": [_sig(sc.duration * i / last) for i in range(len(points))]}
            varies = ["time"] + varies
        elif sc.reactants:
            extra = self._reaction_axes(sc, added, len(points))
            axes.update(extra)
            varies = list(extra) + varies
        if sc.evaporating:
            last = max(len(points) - 1, 1)
            axes["evaporation"] = {"key": "evaporation", "label": "Concentration factor", "unit": "times",
                                   "values": [_sig(sc.evaporate ** (i / last)) for i in range(len(points))]}
            varies = ["evaporation"] + varies
        if sc.mixing:
            extra = self._mixing_axes(sc, len(points))
            axes.update(extra)
            varies = list(extra) + varies
        return {
            "id": _dt.datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:4],
            "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "scenario": sc.to_dict(),
            "title": sc.display_title,
            "summary": sc.describe(),
            "database": f"{db.title}; ion sizes from EQ3/6 data0.cmp",
            "fluid": {"key": fluid.key, "label": fluid.label, "description": fluid.description,
                      "measured_t_c": fluid.measured_t_c, "source": fluid.source,
                      "analysis": [e.line() for e in sc.fluid_analysis().entries], "ph": fluid.ph,
                      "totals_molal": {pretty_species(k, db.species[k].charge if k in db.species else None):
                                       _sig(v) for k, v in sp0.totals.items()}},
            "axes": axes,
            "varies": varies,
            "x": axes[varies[0] if varies else "temperature_c"],
            "reliability": {"levels": [r.level for r in rels],
                            "deviation": [_sig(r.deviation, 3) for r in rels],
                            "density": [_sig(r.density, 4) for r in rels]},
            "start_elements": {k: _sig(v) for k, v in start_elements.items()},
            "catalog": names,
            "modes": runs,
            "notes": seen,
            "elapsed_s": round(time.time() - started, 1),
        }

    def _transport(self, sc: Scenario, progress: Progress) -> dict:
        """A column run (see transport.py): the fluid at the outlet against pore volumes, and every cell at the
        end against distance, each as one 'mode' with its own axes."""
        started = time.time()
        db = self.db
        fluid = sc.fluid_spec()
        sp0, fluid_notes, added, inflow, system, bulk0 = self._prepare(sc, progress)
        inflow = inflow or sp0
        notes: list[str] = list(fluid_notes)
        cells = sc.cell_points()
        n = len(cells)
        # the rock in each cell, per kg of pore water, brought to equilibrium with the pore fluid first
        fill, amounts = dict(bulk0), {}
        for r in self._reactants(system, sc, added, 2, notes):
            for comp, c in r.composition.items():
                fill[comp] = fill.get(comp, 0.0) + c * r.amount
            if r.mineral and db.species[r.mineral].kind == Kind.MINERAL:
                amounts[r.mineral] = amounts.get(r.mineral, 0.0) + r.amount
        entering = {k: v for k, v in inflow.component_totals.items() if k in system.basis}
        for b in system.basis:
            entering.setdefault(b, TRACE)
        blocked = {b.lower() for b in sc.blocked_minerals} | {b.lower() for b in ALWAYS_EXCLUDED}
        column = Column(system, cells, fill, entering, dispersion=sc.mixing_fraction, minerals=list(amounts),
                        amounts=amounts, suppress=blocked, guess=sp0)
        progress(0.0, f"Filling {n} cells")
        shifts = sc.shifts
        res = column.run(shifts, snapshots=(0,), progress=lambda s, total: progress(s / total, f"Shift {s} of {total}"))

        def rows(states, where):
            return [None if st is None else self._record(db, PathPoint(0.0, t, p, st, {}))
                    for st, (t, p) in zip(states, where)]
        outlet_at = [cells[-1]] * len(res.outlet)
        outlet = _pack(rows(res.outlet, outlet_at))
        profile = _pack(rows(res.profile, cells))
        start = _pack(rows(res.snapshots[0], cells))
        rel_cells = self.reliability(cells)
        dx = sc.column_m / n
        pv_axis = {"key": "pore_volumes", "label": "Pore volumes", "unit": "of the column",
                   "values": [_sig(v) for v in res.pore_volumes]}
        dist_axis = {"key": "distance", "label": "Distance from the inlet", "unit": "m",
                     "values": [_sig((i + 0.5) * dx) for i in range(n)]}
        prof_axes = {"distance": dist_axis}
        for key, (label, unit) in AXES.items():
            vals = [pt[0] if key == "temperature_c" else pt[1] for pt in cells]
            if max(vals) - min(vals) > 1e-9:
                prof_axes[key] = {"key": key, "label": label, "unit": unit, "values": [_sig(v) for v in vals]}
        outlet.update(axes={"pore_volumes": pv_axis}, varies=["pore_volumes"],
                      reliability={"levels": [rel_cells[-1].level] * len(res.outlet)})
        for run in (start, profile):
            run.update(axes=prof_axes, varies=list(prof_axes), reliability={"levels": [r.level for r in rel_cells]})
        runs = {"outlet": outlet, "start": start, "profile": profile}
        for name, run, where in (("outlet", outlet, outlet_at), ("start", start, cells), ("profile", profile, cells)):
            notes.extend(run.pop("warnings"))
            _phase_split(run, where)
            if name != "start":
                notes.extend(phase_notes(run, where, "the outlet" if name == "outlet" else "the column"))

        k_mix = math.ceil(3.0 * sc.mixing_fraction - 1e-12) if sc.mixing_fraction > 0 else 0
        notes.append(f"Transport as in PHREEQC's TRANSPORT: {n} cells of {dx:.3g} m; each of the {shifts} shifts "
                     f"moves the pore water one cell on, and the dispersivity ({sc.dispersivity_m:g} m) mixes "
                     f"neighbouring cells by {sc.mixing_fraction:.3g} of the difference per shift"
                     + (f" (in {k_mix} mixes, half before and half after the shift)" if k_mix > 1 else "")
                     + ". Minerals, gas and sorbed matter stay in their cell; every cell returns to equilibrium "
                     "at its own temperature and pressure after each move. Kinetics is not used: the column is "
                     "at local equilibrium.")
        if sc.reactants:
            notes.append("The reactants were put into every cell and brought to equilibrium with the pore fluid "
                         "before the flow started (a rock such as basalt glass is replaced by what it alters to).")
        if res.failures:
            first = res.failures[0]
            notes.append(f"No equilibrium was found {len(res.failures)} times (first in cell {first[1] + 1} at "
                         f"shift {first[0]}: {first[2]}); those cells kept their last state, with their matter "
                         f"still counted.")
        if sc.exchange_meq > 0 or sc.hfo_g > 0:
            notes.append("Sorption uses PHREEQC's constants (phreeqc.dat): Gaines-Thomas exchange and the Dzombak "
                         "and Morel (1990) diffuse-layer model for hydrous ferric oxide; the exchanger and the "
                         "surface start in equilibrium with the pore fluid.")
        notes = [x for x in notes if not x.startswith("water density")]
        notes.extend(reliability_notes(cells, rel_cells, sc.varies))
        notes.extend(fluid.caveats)
        seen: list[str] = []
        for x in notes:
            if x not in seen:
                seen.append(x)
        progress(1.0, "Done")
        axes = {"pore_volumes": pv_axis, **prof_axes}
        start_elements: dict[str, float] = {}
        for name, st in inflow.species.items():
            sp = db.species.get(name)
            for el, k in (sp.elements.items() if sp is not None else ()):
                if el not in ("H", "O"):
                    start_elements[el] = start_elements.get(el, 0.0) + k * st.molality
        return {
            "id": _dt.datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:4],
            "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "scenario": sc.to_dict(),
            "title": sc.display_title,
            "summary": sc.describe(),
            "database": f"{db.title}; ion sizes from EQ3/6 data0.cmp",
            "fluid": {"key": fluid.key, "label": fluid.label, "description": fluid.description,
                      "measured_t_c": fluid.measured_t_c, "source": fluid.source,
                      "analysis": [e.line() for e in sc.fluid_analysis().entries], "ph": fluid.ph,
                      "totals_molal": {pretty_species(k, db.species[k].charge if k in db.species else None):
                                       _sig(v) for k, v in inflow.totals.items()}},
            "transport": {"cells": n, "shifts": shifts, "column_m": sc.column_m, "dx_m": dx,
                          "temperatures": [_sig(t) for t, _ in cells], "pressures": [_sig(p) for _, p in cells],
                          "pore_fluid": sc.pore_fluid or sc.fluid, "failures": len(res.failures)},
            "axes": axes,
            "varies": ["pore_volumes"],
            "x": pv_axis,
            "reliability": outlet["reliability"],
            "start_elements": {k: _sig(v) for k, v in start_elements.items()},
            "catalog": _catalog(db, runs, sc),
            "modes": runs,
            "notes": seen,
            "elapsed_s": round(time.time() - started, 1),
        }

    def _prepare(self, sc: Scenario, progress: Progress):
        """The starting fluid, the reactants' names, the second fluid, the system and its bulk."""
        db = self.db
        progress(0.0, f"Speciating {sc.fluid_spec().label.lower()}")
        sp0, fluid_notes = self.starting_fluid(sc)
        inflow = None
        if sc.transport and sc.pore_fluid and sc.pore_fluid != sc.fluid:
            # a column first filled with another fluid: that one is the start, the card's fluid flows in
            inflow = sp0
            sp0, pore_notes = self.starting_fluid(Scenario(fluid=sc.pore_fluid, activity_model=sc.activity_model,
                                                           methane_sulfate=sc.methane_sulfate))
            fluid_notes = fluid_notes + [f"Pore fluid: {n}" for n in pore_notes]
        basis = [b for b in sp0.basis_names if b not in (WATER, PROTON)]
        added = self._reactant_names(sc)
        for name in added.values():                      # components the fluid does not carry yet
            elements = (rocks.elements_per_gram(rocks.ROCKS[name[5:]])[0] if name.startswith("rock:")
                        else db.species[name].elements)
            for el in elements:
                if name.startswith("rock:") and (el == "Fe3" or DEFAULT_MASTERS.get(el) is None):
                    continue                             # ferric iron goes with Fe; Ti is left out
                master = DEFAULT_MASTERS.get(el)
                if el in ("H", "O"):
                    continue
                if master is None:
                    raise ValueError(f"{name} contains {el}, which this calculator cannot track yet")
                if master not in basis:
                    basis.append(master)
            if (sc.hold_methane_apart and name in ("CH4(g)", "CH4") and "CH4" not in basis
                    and "O2(g)" in basis):
                basis.append("CH4")                      # methane added to a fluid without any
        mix = None
        if sc.mixing:
            progress(0.0, f"Speciating {preset(sc.mix_with).label.lower()}")
            mix, mix_notes = self.starting_fluid(Scenario(fluid=sc.mix_with, activity_model=sc.activity_model))
            fluid_notes += [f"Second fluid: {n}" for n in mix_notes]
            basis += [b for b in mix.basis_names if b not in basis and b not in (WATER, PROTON)]
        if inflow is not None:
            mix = inflow
            basis += [b for b in inflow.basis_names if b not in basis and b not in (WATER, PROTON)]
        hold = [b for b in basis if db.species[b].kind != Kind.BASIS]     # CH4, held apart
        blocked = {b.lower() for b in sc.blocked_minerals}
        mixes = [ss for ss in solid_solution_table() if ss.name not in blocked] if sc.solid_solutions else []
        sorb = Sorption(exchange_eq=sc.exchange_meq / 1000.0, hfo_g=sc.hfo_g)
        system = System(db, [b for b in basis if b not in hold], activity_model=sc.activity_model,
                        gas_phase=sc.gas_phase, hold_apart=hold, solid_solutions=mixes,
                        sorption=sorb if sorb.active else None)
        bulk0 = dict(sp0.component_totals)
        for b in basis:
            bulk0.setdefault(b, TRACE)                   # a trace, so the solution has a start
        if sorb.exchange_eq > 0:
            # the exchanger starts loaded with the fluid's cations in equilibrium with it (the fluid unchanged)
            t0 = sp0.temperature_c
            model = system.model(t0, sp0.pressure_bar)
            have = dict(zip(sp0.basis_names, sp0.ln_activities))
            u = np.array([have.get(b, -50.0) for b in model.names])
            gam = {n: st.log_activity * math.log(10.0) - math.log(st.molality)
                   for n, st in sp0.species.items() if st.molality > 0}
            for comp, moles in exchanger_loading(model, u, sorb.exchange_eq, t0, gam).items():
                bulk0[comp] = bulk0.get(comp, 0.0) + moles
        return sp0, fluid_notes, added, mix, system, bulk0

    def _reactant_names(self, sc: Scenario) -> dict[str, str]:
        """The database name of each reactant on the card (matched without regard to case)."""
        db = self.db
        lower: dict[str, list[str]] = {}
        for name, sp in db.species.items():
            lower.setdefault(name.lower(), []).append(name)
        out = {}
        for r in sc.reactants:
            name = r["name"]
            gas_word = re.fullmatch(r"(.+?)\s+gas|(.+?)\s*\(g\)|(methane|carbon dioxide)", name.strip(), re.I)
            if gas_word:
                stem = (gas_word.group(1) or gas_word.group(2) or gas_word.group(3)).strip()
                stem = {"methane": "CH4", "carbon dioxide": "CO2", "hydrogen": "H2", "nitrogen": "N2",
                        "hydrogen sulfide": "H2S", "oxygen": "O2"}.get(stem.lower(), stem)
                if f"{stem}(g)" in db.species:
                    name = f"{stem}(g)"
            if name not in db.species:
                hits = lower.get(name.lower(), [])
                minerals = [h for h in hits if db.species[h].kind == Kind.MINERAL]
                rock = rocks.find(name)
                if not hits and rock is not None:
                    if r["unit"] != "g":
                        raise ValueError(f"give {rock.label.lower()} in grams")
                    out[r["name"]] = f"rock:{rock.key}"
                    continue
                if not hits:
                    raise ValueError(f"'{name}' is not a mineral, rock or species in the database")
                name = (minerals or hits)[0]
            if db.species[name].kind not in (Kind.MINERAL, Kind.AQUEOUS, Kind.BASIS, Kind.GAS):
                raise ValueError(f"{name} cannot be a reactant (only minerals, gases and dissolved species)")
            out[r["name"]] = name
        return out

    def _reactants(self, system: System, sc: Scenario, names: dict[str, str], steps: int,
                   notes: list[str] | None = None) -> list:
        """Each reactant added at a steady rate, all of it by the last step (reaction progress counts
        the steps).  With kinetics, reaction progress is time in seconds: a mineral or glass with a
        measured rate reacts at that rate, anything else is added steadily over the run."""
        out = []
        span = sc.duration_s if sc.kinetics else max(steps - 1, 1)
        area = sc.surface_area_cm2_per_g * 1e-4          # m2 per g
        steady: list[str] = []
        for r in sc.reactants:
            name = names[r["name"]]
            if name.startswith("rock:"):                 # dissolves congruently, never precipitates
                rock = rocks.ROCKS[name[5:]]
                comp, rock_notes = rocks.components_per_gram(rock, self.db, system.basis)
                if notes is not None:
                    notes.extend(n for n in rock_notes if n not in notes)
                law = None
                if sc.kinetics and rock.key == "basalt_glass":
                    law = basalt_glass(rocks.elements_per_gram(rock)[0]["Si"])
                elif sc.kinetics:
                    steady.append(rock.label)
                out.append(Reactant(rock.label, comp, r["amount"], (r["amount"] / span, 0.0, 0.0),
                                    rate_law=law, area=area * r["amount"]))
                continue
            if name not in system.net.phase_names and name not in system.net.dep_names:
                raise ValueError(f"{name} cannot be added: its formula needs the redox state of the "
                                 f"fluid, which is not tracked; set the redox state on the card")
            mw = self.db.molecular_weight(name)
            moles = r["amount"] / mw if r["unit"] == "g" else r["amount"]
            law = rate_law(name) if sc.kinetics and self.db.species[name].kind == Kind.MINERAL else None
            if sc.kinetics and law is None:
                steady.append(pretty_mineral(name) if self.db.species[name].kind == Kind.MINERAL else name)
            out.append(mineral_reactant(system, name, moles, rates=(moles / span, 0.0, 0.0),
                                        rate_law=law, area=area * moles * mw))
        if steady and notes is not None:
            notes.append(f"No measured rate for {', '.join(steady)}: added steadily over the run instead.")
        return out

    def _mixing_axes(self, sc: Scenario, steps: int) -> dict:
        """Second fluid added (kg of its water per kg of the first's) and its fraction of the mixture."""
        kg = [sc.mix_ratio * i / max(steps - 1, 1) for i in range(steps)]
        other = preset(sc.mix_with).label.lower()
        return {"mixed": {"key": "mixed", "label": f"{other.capitalize()} added", "unit": "kg per kg",
                          "values": [_sig(k) for k in kg]},
                "mix_fraction": {"key": "mix_fraction", "label": f"Fraction of {other}", "unit": "by mass",
                                 "values": [_sig(k / (1.0 + k)) for k in kg]}}

    def _reaction_axes(self, sc: Scenario, names: dict[str, str], steps: int) -> dict:
        """Reactant added (g per kg of starting water) and the water/rock mass ratio at each step; for
        reactants given only in moles (a titration), the millimoles added."""
        if all(r["unit"] == "mol" for r in sc.reactants):
            mmol = sum(r["amount"] for r in sc.reactants) * 1000.0
            return {"reacted": {"key": "reacted", "label": "Reactant added", "unit": "mmol per kg water",
                                "values": [_sig(mmol * i / max(steps - 1, 1)) for i in range(steps)]}}
        grams = [0.0] * steps
        for r in sc.reactants:
            g = r["amount"] if r["unit"] == "g" else r["amount"] * self.db.molecular_weight(names[r["name"]])
            for i in range(steps):
                grams[i] += g * i / max(steps - 1, 1)
        return {"reacted": {"key": "reacted", "label": "Reactant added", "unit": "g per kg water",
                            "values": [_sig(g) for g in grams]},
                "water_rock": {"key": "water_rock", "label": "Water/rock", "unit": "by mass",
                               "values": [_sig(1000.0 / g) if g > 0 else None for g in grams]}}

    def starting_fluid(self, sc: Scenario) -> tuple[Speciation, list[str]]:
        """Speciate the card's fluid at its measurement temperature, closing the charge balance as the
        card asks.  Returns the speciation and sentences about what was adjusted."""
        db = self.db
        fluid = sc.fluid_spec()
        read = analysis.parse(fluid.analysis, fluid.units)
        if read.problems:
            raise ValueError("; ".join(read.problems))
        molal, _ = analysis.molalities(read.entries)
        kind, arg = sc.redox_setting()
        # every redox species counts toward the element total of its basis partner; H2 and O2 carry
        # only the redox state.  The couple (if any) first gives the oxygen fugacity on its own.
        apart = sc.hold_methane_apart
        folded: dict[str, float] = {}
        for name, m in molal.items():
            target = analysis.REDOX.get(name, name) if name in analysis.REDOX else name
            if name == "CH4" and apart:
                target = "CH4"                           # its own component: it never meets sulfate
            if target is not None:
                folded[target] = folded.get(target, 0.0) + m
        if apart and kind != "none" and "HCO3-" in folded and "CH4" not in folded and "CH4" in db.species:
            folded["CH4"] = TRACE                        # a column of its own, so a run can add some
        constraints = {name: Molal(m) for name, m in folded.items()}
        constraints[PROTON] = PH(fluid.ph)
        redox = {"air": lambda: LogFO2(AIR_LOG_FO2), "log_fo2": lambda: LogFO2(arg),
                 "eh": lambda: Eh(arg), "pe": lambda: Pe(arg)}.get(kind, lambda: None)()

        def solve(balance_on=None):
            if kind != "couple":
                return speciate(db, fluid.measured_t_c, constraints, redox=redox,
                                pressure_bar=fluid.measured_p_bar, balance_on=balance_on, activity_model=sc.activity_model)
            own = dict(constraints)
            partner = analysis.REDOX[arg]
            if partner is not None:
                rest = folded[partner] - (0.0 if arg == "CH4" and apart else molal[arg])
                own[partner] = Molal(rest)
            own[arg] = Molal(molal[arg])
            first = speciate(db, fluid.measured_t_c, own, redox=Couple(arg),
                             pressure_bar=fluid.measured_p_bar, balance_on=balance_on, activity_model=sc.activity_model)
            return speciate(db, fluid.measured_t_c, constraints, redox=LogFO2(first.log_fo2),
                            pressure_bar=fluid.measured_p_bar, balance_on=balance_on, activity_model=sc.activity_model)
        notes = list(read.notes)
        sp = solve()
        if kind != "none":
            how = {"air": "equilibrium with air", "log_fo2": "the log fO2 given",
                   "eh": "the Eh given", "pe": "the pe given"}.get(
                kind, f"dissolved {analysis.SOLUTES_LABEL(arg) if arg else ''}")
            notes.append(f"Redox state from {how}: log fO2 {sp.log_fo2:.1f}, Eh {sp.eh_volts:+.2f} V "
                         f"at {fluid.measured_t_c:g} C.")
        equivalents = sum(abs(analysis.SOLUTES[n][0]) * m for n, m in molal.items()) or 1.0
        off = sp.charge_imbalance / equivalents * 2.0         # relative to the mean of the two sides
        ion = fluid.charge_balance
        if ion == "auto":
            ion = analysis.balance_ion(molal, sp.charge_imbalance) if abs(off) > 1e-6 else None
        if ion in (None, "none"):
            if abs(off) > 0.05:
                notes.append(f"The analysis is out of charge balance by {100 * abs(off):.0f}% "
                             f"({'cations' if off > 0 else 'anions'} in excess) and was used as given.")
            return sp, notes
        sp = solve(ion)
        before, after = molal[ion], sp.totals[ion]
        change = (after - before) / before
        z = db.species[ion].charge if ion in db.species else None
        notes.append(f"Charge balance closed on {pretty_species(ion, z)}: {before:.4g} to {after:.4g} "
                     f"mol/kg ({100 * change:+.1f}%).")
        if abs(change) > 0.10:
            notes.append(f"That is a large change to {pretty_species(ion, z)}; check the analysis.")
        return sp, notes

    # ---- the march along the path ---------------------------------------------------------------
    def _march(self, system: System, bulk0: dict, points: list[tuple[float, float]], sc: Scenario,
               flow_through: bool, tick, notes: list[str], mix: Speciation | None = None) -> dict:
        """The library's reaction path, with reaction progress counting the steps of the card: the same
        closed-system and flow-through rules as every other path (see ``path.py``)."""
        db = self.db
        n = len(points)
        last = max(n - 1, 1)
        step = sc.duration_s / last if sc.kinetics else 1.0      # reaction progress per step
        blocked = {b.lower() for b in sc.blocked_minerals} | set(ALWAYS_EXCLUDED)
        suppress = [name for name, kind in zip(system.net.phase_names, system.net.phase_kind)
                    if kind == Kind.MINERAL and name.lower() in blocked]
        reactants = self._reactants(system, sc, self._reactant_names(sc), n, notes)
        if mix is not None:          # the second fluid, per kg of its water, as one more reactant
            composition = {k: v for k, v in mix.component_totals.items() if k in system.basis}
            reactants.append(Reactant(f"mix:{sc.mix_with}", composition, sc.mix_ratio,
                                      (sc.mix_ratio / (last * step), 0.0, 0.0)))
        if sc.evaporating:           # water taken away, evenly in the logarithm of the concentration
            reactants.append(_Evaporation("evaporation", {WATER: -WATER_MOLES_PER_KG},
                                          1.0 - 1.0 / sc.evaporate, factor=sc.evaporate, steps=last))
        split: dict[float, tuple] = {}
        separate = self._separator(system, sc, split) if sc.follow_phase != "whole" else None
        path = ReactionPath(system, bulk0, reactants=reactants, separate=separate,
                            temperature=lambda xi: sc.tp(xi / (last * step))[0],
                            pressure=lambda xi: sc.tp(xi / (last * step))[1],
                            suppress_minerals=suppress, flow_through=flow_through)

        def failed(xi, exc):
            t, p = points[int(round(xi / step))]
            notes.append(f"No equilibrium found at {t:.0f} C and {p:.0f} bar ({exc}); that point is "
                         f"left blank.")

        solved = path.run([i * step for i in range(n)], on_failure=failed, progress=tick)
        kinetic = [r for r in reactants if r.kinetic]
        run = _pack([None if pt is None else self._record(db, pt, kinetic) for pt in solved])
        if separate is not None:
            run["split"] = [split.get(round(i * step, 9)) for i in range(n)]
        if sc.isotopes:
            run["isotopes"] = self._isotopes(sc, bulk0, reactants, solved, points, notes)
        return run

    def _isotopes(self, sc: Scenario, bulk0: dict, reactants: list, solved: list, points: list,
                  notes: list[str]) -> dict[str, list]:
        """δ18O, δ13C, δ34S of the fluid and of each mineral, and 87Sr/86Sr of the fluid, at each point (see
        isotopes.py), as {"d18O:water": [...], "d18O:quartz": [...], "87Sr/86Sr:fluid": [...]}."""
        db = self.db
        given = {k.lower(): v for k, v in sc.isotope_inputs.items()}
        used: list[str] = []
        fluid_sr = bulk0.get("Sr+2", 0.0)
        tracker = isotopes.Tracker(db, self._isotope_source(sc.fluid, given.get("fluid", {}), fluid_sr, used,
                                                           "the fluid"), bulk0)
        sources = {}
        for card, r in zip(sc.reactants, reactants):
            sources[r.name] = self._isotope_source(
                card["name"], given.get(card["name"].lower(), given.get(r.name.lower(), {})), 0.0, used, r.name, r)
        for r in reactants[len(sc.reactants):]:
            if r.name.startswith("mix:"):
                sources[r.name] = self._isotope_source(r.name[4:], given.get("mix", {}), r.composition.get("Sr+2", 0.0),
                                                       used, "the second fluid")
        n = len(points)
        out: dict[str, list] = {}
        done = {r.name: 0.0 for r in reactants}
        for i, pt in enumerate(solved):
            if pt is None:
                continue
            t = pt.temperature_c
            for r in reactants:
                now = pt.reacted.get(r.name, 0.0)
                step, done[r.name] = now - done[r.name], now
                if step == 0:
                    continue
                if r.name == "evaporation":
                    tracker.evaporate(-step * r.composition.get(WATER, 0.0), t)
                elif r.name in sources:
                    tracker.add(r.composition, sources[r.name], step)
            deltas = tracker.point(pt.state, t, pt.speciation.component_totals.get("Cl-", 0.0))
            row = {}
            for e, by in deltas.items():
                key = isotopes.SYSTEMS[e]
                for carrier, d in by.items():
                    if carrier in ("CO2", "H2S", "other S", "organic C") or carrier.endswith("(g)"):
                        continue
                    if e == "H" and carrier != "water":   # minerals only echo the water's δD
                        continue
                    row[f"{key}:{carrier}"] = d
            ratio = tracker.sr_fluid_ratio
            if ratio is not None:
                row["87Sr/86Sr:fluid"] = ratio
            for k, v in row.items():
                out.setdefault(k, [None] * n)[i] = round(float(v), 6 if k.startswith("87") else 3)
        said: list[str] = []
        if used:
            said.append("Isotopes: starting values not given on the card were taken as " + "; ".join(used) + ".")
        said.append("Isotopes are in equilibrium among the fluid and the minerals at each step (oxygen: Vho and others "
                     "2020; carbon: Ohmoto and Rye 1979; sulfur: Ohmoto and Rye 1979, Ohmoto and Lasaga 1982); "
                     "δ18O and δ13C of carbonates against VSMOW and VPDB, δ34S against VCDT. Strontium is carried "
                     "with the fluid only (minerals are taken to keep none). δD (VSMOW) is given for the water only: "
                     "it follows evaporation (liquid-vapour fractionation of Horita and Wesolowski 1994), mixing and "
                     "phase separation; hydrous minerals take the water's δD, as mineral-water D/H factors are not in "
                     "yet.")
        for e, names in sorted(tracker.frac.analogues.items()):
            said.append(f"Isotopes ({isotopes.SYSTEMS[e]}): no fractionation of their own for "
                         f"{', '.join(sorted(names))}, so those analogues are used.")
        if tracker.frac.beyond:
            said.append("Isotopes: fractionations used above their calibrated range: "
                        f"{', '.join(sorted(tracker.frac.beyond))}.")
        notes.extend(n for n in said if n not in notes)
        return out

    def _isotope_source(self, name: str, given: dict, sr_moles: float, used: list[str], label: str,
                        reactant: Reactant | None = None) -> isotopes.Source:
        """Starting δ values and Sr of the fluid, the second fluid or a reactant: the card's, or defaults."""
        key = name.lower()
        db = self.db
        rock = rocks.find(name) if reactant is not None else None
        if reactant is None:                          # a fluid: per kg of its water
            base = dict(ISOTOPE_FLUIDS.get(key, ISOTOPE_FLUIDS["seawater"]))
            sr_unit = "mg/kg"
        elif rock is not None:
            base = dict(ISOTOPE_ROCKS.get(rock.key, ISOTOPE_ROCKS["basalt_glass"]))
            sr_unit = "ppm"
        else:
            el = db.species[reactant.name].elements if reactant.name in db.species else {}
            base = ({"d18O": 29.0, "d13C": 0.0} if el.get("C") and el.get("O")
                    else {"d13C": -8.5} if el.get("C")
                    else {"d18O": 9.5, "d34S": 21.0} if el.get("S") and el.get("O", 0) >= 4 * el["S"]
                    else {"d34S": 0.0} if el.get("S") else {"d18O": 6.0})
            if el.get("H"):
                base["dD"] = -60.0                   # an illustrative hydrous mineral
            sr_unit = "ppm"
        defaults = {k: v for k, v in base.items() if k not in given}
        merged = {**base, **given}
        if sr_moles > 0 and "Sr" not in given:        # the fluid's own analysis has Sr
            merged.pop("Sr", None)
            defaults.pop("Sr", None)
        if defaults:
            used.append(f"{label}: " + ", ".join(
                f"{k} {v:g}" + (" " + sr_unit if k == "Sr" else "‰" if k.startswith("d") else "")
                for k, v in defaults.items()))
        if "Sr" in merged:
            if reactant is None:
                sr_moles = merged["Sr"] / 1000.0 / isotopes.SR_GRAMS_PER_MOLE
            else:
                grams = 1.0 if rock is not None else (db.molecular_weight(reactant.name) or 0.0)
                sr_moles = merged["Sr"] * 1e-6 * grams / isotopes.SR_GRAMS_PER_MOLE
        return isotopes.Source({k: v for k, v in merged.items() if k.startswith("d")}, sr_moles,
                               merged.get("87Sr/86Sr"))

    def _separator(self, system: System, sc: Scenario, record: dict):
        """The path's phase-separation hook when the card follows the brine or the vapour: where the fluid
        splits (NaCl-H2O, by its chloride), the other phase leaves.  Water goes by the lever rule; chloride
        and every other solute in proportion to the salt in each phase; dissolved gases (neutral species
        with a gas of the same formula) with the water, as they do near the critical point, where gases
        partition about evenly between the two phases per kg of water.  Minerals stay with the brine and a
        separate gas phase goes with the vapour.  ``record`` gets, for each xi,
        (phase state before the split, kg of fluid that left)."""
        net = system.net
        gases = {n[:-3] for n in net.phase_names if n.endswith("(g)")}
        brine = sc.follow_phase == "brine"
        gone = [0.0]

        def separate(xi, state):
            spec = state.speciation
            x = salinity_from_chloride(spec.component_totals.get("Cl-", 0.0))
            ps = phase_state(spec.temperature_c, spec.pressure_bar, x)
            if not ps.separated:
                record[round(xi, 9)] = (ps, gone[0])
                return None
            w = state.water_kg
            share_w = ps.vapour_mass * (1.0 - ps.x_vap) / (1.0 - x)       # of the water, in the vapour
            share_s = ps.vapour_mass * ps.x_vap / x                         # of the salt, in the vapour
            if not brine:
                share_w, share_s = 1.0 - share_w, 1.0 - share_s
            fluid_kg = w / (1.0 - x)
            gone[0] += fluid_kg * (ps.vapour_mass if brine else 1.0 - ps.vapour_mass)
            record[round(xi, 9)] = (ps, gone[0])
            out = {WATER: share_w * w * WATER_MOLES_PER_KG}
            for name, st in spec.species.items():
                if name == WATER or st.molality <= 0:
                    continue
                volatile = name in gases and self.db.species[name].charge == 0
                moles = (share_w if volatile else share_s) * st.molality * w
                for comp, c in _composition_row(system, name).items():
                    out[comp] = out.get(comp, 0.0) + c * moles
            # the vapour moves on without the minerals; a separate gas phase goes with the vapour
            left = dict(state.gas) if brine else dict(state.minerals)
            if not brine:
                for members in state.solid_solutions.values():
                    left.update(members)
            for name, moles in left.items():
                for comp, c in _composition_row(system, name).items():
                    out[comp] = out.get(comp, 0.0) + c * moles
            return out

        return separate

    def _uncertainty(self, sc: Scenario, modes: list[str], runs: dict[str, dict], progress: Progress
                     ) -> list[str]:
        """Run the card ``sc.uncertainty_runs`` more times with perturbed log K (see ``_sample``), on
        all but one processor, and attach the spread of the results to ``runs`` as ``bands``.  Returns
        notes on the range and on minerals the uncertainty adds or removes."""
        n = sc.uncertainty_runs
        samples: dict[str, list[dict]] = {m: [] for m in modes}
        sigmas: list[float] = []

        def take(out):
            for mode, packed in out["runs"].items():
                if packed is not None:
                    samples[mode].append(packed)
            sigmas.extend(out["sigma"])

        workers = min(n, max((os.cpu_count() or 1) - 1, 1))
        pool = self._sample_pool(workers) if workers > 1 else None
        if pool is None:
            for j in range(n):
                progress(j / n, f"Uncertainty run {j + 1} of {n}")
                take(self._sample(sc, modes, j))
        else:
            card = sc.to_dict()
            futures = [pool.submit(_sample_in_worker, str(self.data_dir), card, modes, j) for j in range(n)]
            for done, fut in enumerate(as_completed(futures), start=1):
                progress(done / n, f"Uncertainty runs: {done} of {n} done")
                take(fut.result())
        notes = []
        for mode, run in runs.items():
            got = samples[mode]
            run["bands"] = _bands(run, got)
            run["band_runs"] = len(got)
            if not got:
                continue
            where = f" ({MODE_LABELS[mode].lower()})" if len(runs) > 1 else ""

            def formed(r, name):
                return any(v and v > 1e-6 for v in r["minerals"].get(name, ()))

            central = [m for m in run["minerals"] if formed(run, m)]
            extra = sorted({m for r in got for m in r["minerals"] if m not in central and formed(r, m)})
            some = [(m, sum(formed(r, m) for r in got)) for m in central]
            if extra:
                notes.append(f"Within the log K uncertainty{where}, these minerals also form in some runs: " +
                             ", ".join(f"{pretty_mineral(m)} ({sum(formed(r, m) for r in got)} of {len(got)})"
                                       for m in extra) + ".")
            shaky = [f"{pretty_mineral(m)} ({k} of {len(got)})" for m, k in some if k < len(got)]
            if shaky:
                notes.append(f"Within the log K uncertainty{where}, these minerals do not form in every run: "
                             + ", ".join(shaky) + ".")
        lost = sum(n - len(v) for v in samples.values())
        sig = sorted(sigmas)
        notes.append(f"Shaded ranges: the middle 90% of {n} runs with every log K shifted at random by "
                     f"{sig[0]:.2f} to {sig[-1]:.2f} log units (one standard deviation, larger where water "
                     f"ionization is reproduced less well)." if sig else "")
        if lost:
            notes.append(f"{lost} of the uncertainty runs did not converge and are left out of the range.")
        return [x for x in notes if x]

    def _sample(self, sc: Scenario, modes: list[str], j: int) -> dict:
        """Uncertainty run ``j`` of a card: every log K (aqueous species, minerals, gases) shifted by its
        own standard-normal draw times ``validity.logk_sigma`` at each point.  A draw is kept along the
        whole path, as a database error would be, and depends only on ``j``, so a card always gives
        the same range.  Returns {"runs": {mode: packed run, or None if it failed}, "sigma": [...]}."""
        key = json.dumps(sc.to_dict(), sort_keys=True)
        if self._prepared is None or self._prepared[0] != key:
            self._prepared = (key, self._prepare(sc, lambda f, m: None))
        _, _, _, mix, system, bulk0 = self._prepared[1]
        db, net, points = self.db, system.net, sc.points()
        cache: dict[tuple, float] = {}

        def sigma(t_c, p_bar):
            k = (round(t_c, 6), round(p_bar, 6))
            if k not in cache:
                cache[k] = validity.logk_sigma(db, t_c, p_bar)
            return cache[k]

        rng = np.random.default_rng([BAND_SEED, j])
        dep = rng.standard_normal(len(net.dep_names))
        # polymorphs (same formula: calcite and aragonite, quartz and its relatives) share one draw: their
        # differences are known far better than their absolute values
        _, group = np.unique(np.round(net.nu_phase, 6), axis=0, return_inverse=True)
        phase = rng.standard_normal(int(group.max()) + 1)[group.ravel()]
        if "H2O(g)" in net.phase_names:          # steam against liquid water: the steam tables, exact
            phase[net.phase_names.index("H2O(g)")] = 0.0
        system.perturb(dep, phase, sigma)
        out: dict[str, dict | None] = {}
        try:
            for mode in modes:
                try:
                    out[mode] = self._march(system, bulk0, points, sc, mode == "flow_through",
                                            lambda i, total: None, [], mix)
                except Exception:                        # a draw the solver cannot follow: left out
                    out[mode] = None
        finally:
            system.perturb()
        return {"runs": out, "sigma": list(cache.values())}

    def _sample_pool(self, workers: int):
        """Worker processes for uncertainty runs, started once (spawned: each loads the data itself)."""
        if self._pool is None or self._pool_size != workers:
            if self._pool is not None:
                self._pool.shutdown(wait=False, cancel_futures=True)
            try:
                self._pool = ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn"))
                self._pool_size = workers
            except (OSError, NotImplementedError, ValueError):
                self._pool = None
        return self._pool

    def _record(self, db, point: PathPoint, kinetic: list | None = None) -> dict:
        state = point.state
        spec = state.speciation
        elements: dict[str, float] = {}
        species = {}
        for name, st in spec.species.items():
            species[name] = st.molality
            sp = db.species.get(name)
            if sp is None:
                continue
            for el, k in sp.elements.items():
                if el not in ("H", "O"):
                    elements[el] = elements.get(el, 0.0) + k * st.molality
        minerals, gas = {}, {}                           # flow-through: everything left behind so far
        for name, moles in point.removed.items():
            target = gas if name.endswith("(g)") else minerals
            target[name] = target.get(name, 0.0) + moles
        for name, moles in state.minerals.items():
            minerals[name] = minerals.get(name, 0.0) + moles
        for name, moles in state.gas.items():
            gas[name] = gas.get(name, 0.0) + moles
        composition = {}                                 # mole fraction of each endmember
        for name, members in state.solid_solutions.items():
            total = sum(members.values())
            minerals[name] = minerals.get(name, 0.0) + total
            for member, moles in members.items():
                composition[f"{name}:{member}"] = moles / total if total > 0 else 0.0
        rates = {}                                       # kinetic reactants: rate and share still there
        for r in kinetic or ():
            left = r.amount - point.reacted.get(r.name, 0.0)
            if r.mineral:
                minerals[r.mineral] = minerals.get(r.mineral, 0.0) + max(left, 0.0)
            rates[f"rate:{r.name}"] = point.rates.get(r.name, 0.0)
            rates[f"left:{r.name}"] = 100.0 * max(left, 0.0) / r.amount
        phases = {"salinity": 100.0 * salinity_from_chloride(spec.component_totals.get("Cl-", 0.0))}
        return {"kinetic": rates, "composition": composition, "phases": phases, "pH": spec.ph, "ionic_strength": spec.ionic_strength, "water_kg": state.water_kg,
                "log_fO2": spec.log_fo2, "Eh": spec.eh_volts,
                "log_aw": spec.log_a_water, "species": species,
                "log_a": {n: spec.species[n].log_activity for n in spec.basis_names
                          if n in spec.species and n not in (WATER, PROTON)}, "elements": elements,
                "minerals": minerals, "gas": gas, "saturation": dict(spec.saturation),
                "warnings": list(spec.warnings), "sorbed": dict(state.sorbed),
                "surface_potential_mV": state.surface_potential_mv}


def reliability_notes(points: list[tuple[float, float]], rels: list[validity.Reliability],
                      varies: list[str]) -> list[str]:
    """One sentence per stretch of the path where the thermodynamics is not 'good'."""
    segments: list[list] = []                      # [level, first index, last index, worst deviation]
    for i, r in enumerate(rels):
        if r.level == "good":
            continue
        dev = abs(r.deviation) if math.isfinite(r.deviation) else float("inf")
        if segments and segments[-1][0] == r.level and segments[-1][2] == i - 1:
            segments[-1][2], segments[-1][3] = i, max(segments[-1][3], dev)
        else:
            segments.append([r.level, i, i, dev])
    notes = []
    for level, i0, i1, worst in segments:
        (t0, p0), (t1, p1) = points[i0], points[i1]
        if i0 == i1:
            where = f"At {t0:.0f} °C and {p0:.0f} bar"
        elif "temperature_c" in varies and "pressure_bar" in varies:
            where = f"From {t0:.0f} °C, {p0:.0f} bar to {t1:.0f} °C, {p1:.0f} bar"
        elif "temperature_c" in varies:
            where = f"From {min(t0, t1):.0f} to {max(t0, t1):.0f} °C at {p0:.0f} bar"
        else:
            where = f"From {min(p0, p1):.0f} to {max(p0, p1):.0f} bar at {t0:.0f} °C"
        worst_text = 99.0 if not math.isfinite(worst) else worst
        notes.append(f"{where}: " + validity.LEVEL_TEXT[level].format(worst=worst_text, tol=validity.GOOD))
    return notes


def _pack(rows: list[dict | None]) -> dict:
    """Turn a list of per-point records into one list per quantity (None where unsolved)."""
    good = [r for r in rows if r is not None]

    def column(getter, keys):
        return {k: [_sig(getter(r, k)) if r is not None else None for r in rows] for k in keys}

    species_keys = sorted({k for r in good for k, v in r["species"].items() if v > 0})
    species_keys = [k for k in species_keys if max(r["species"].get(k, 0.0) for r in good) >= 1e-12]
    element_keys = sorted({k for r in good for k in r["elements"]})
    mineral_keys = sorted({k for r in good for k, v in r["minerals"].items() if v > 0})
    gas_keys = sorted({k for r in good for k, v in r.get("gas", {}).items() if v > 1e-9})
    sat_keys = sorted({k for r in good for k, v in r["saturation"].items()
                       if math.isfinite(v) and v > -3.0})
    return {
        "properties": {
            "pH": [_sig(r["pH"]) if r else None for r in rows],
            "ionic_strength": [_sig(r["ionic_strength"]) if r else None for r in rows],
            "water_kg": [_sig(r["water_kg"]) if r else None for r in rows],
            "log_aw": [_sig(r["log_aw"]) if r else None for r in rows],
            **{k: [_sig(r.get(k)) if r else None for r in rows] for k in ("log_fO2", "Eh")
               if any(r and math.isfinite(r.get(k, math.nan)) for r in rows)},   # only with redox
            **({"surface_potential_mV": [_sig(r.get("surface_potential_mV")) if r else None for r in rows]}
               if any(r and r.get("surface_potential_mV") is not None for r in rows) else {}),
        },
        "sorbed": column(lambda r, k: r.get("sorbed", {}).get(k, 0.0),
                         sorted({k for r in good for k, v in r.get("sorbed", {}).items() if v > 1e-15})),
        "species": column(lambda r, k: r["species"].get(k, 0.0), species_keys),
        "log_a": column(lambda r, k: r.get("log_a", {}).get(k), sorted({k for r in good for k in r.get("log_a", {})})),
        "elements": column(lambda r, k: r["elements"].get(k, 0.0), element_keys),
        "minerals": column(lambda r, k: r["minerals"].get(k, 0.0), mineral_keys),
        "gas": column(lambda r, k: r.get("gas", {}).get(k, 0.0), gas_keys),
        "kinetic": column(lambda r, k: r.get("kinetic", {}).get(k), sorted({k for r in good for k in r.get("kinetic", {})})),
        "phases": column(lambda r, k: r.get("phases", {}).get(k), ["salinity"]),
        "composition": column(lambda r, k: r.get("composition", {}).get(k),
                              sorted({k for r in good for k in r.get("composition", {})})),
        "saturation": column(lambda r, k: r["saturation"].get(k), sat_keys),
        "warnings": sorted({w for r in good for w in r["warnings"]}),
    }


_WORKER: "Engine | None" = None


def _sample_in_worker(data_dir: str, card: dict, modes: list[str], j: int) -> dict:
    """One uncertainty run in a worker process, with that process's own engine."""
    global _WORKER
    if _WORKER is None or str(_WORKER.data_dir) != data_dir:
        _WORKER = Engine(data_dir)
    return _WORKER._sample(Scenario.from_dict(card), modes, j)


def _phase_split(run: dict, points: list[tuple[float, float]]) -> None:
    """Brine and vapour along the path from the fluid's salinity (NaCl-H2O): the phase state of every
    point, and where the fluid splits, the vapour share and the salinity of each phase.  A point the
    chemistry could not solve takes the salinity of the nearest point that was solved, so the split shows
    even where the aqueous model gives up near the critical point."""
    sal = run["phases"]["salinity"]
    known = [i for i, v in enumerate(sal) if v is not None]
    state, share, brine, vapour = [], [], [], []
    split = run.pop("split", None)
    if split is not None:                    # following one phase: the split before it was taken
        run["phases"]["split_off_kg"] = [_sig(s[1]) if s else None for s in split]
    for i, (t, p) in enumerate(points):
        if split is not None and split[i] is not None:
            ps = split[i][0]
            state.append(ps.kind)
            share.append(_sig(100.0 * ps.vapour_mass) if ps.separated else None)
            brine.append(_sig(100.0 * ps.x_liq) if ps.separated else None)
            vapour.append(_sig(100.0 * ps.x_vap) if ps.separated else None)
            continue
        if not known:
            state.append(None)
            share.append(None), brine.append(None), vapour.append(None)
            continue
        x = sal[min(known, key=lambda k: abs(k - i))] / 100.0
        ps = phase_state(t, p, x)
        state.append(ps.kind)
        share.append(_sig(100.0 * ps.vapour_mass) if ps.separated else None)
        brine.append(_sig(100.0 * ps.x_liq) if ps.separated else None)
        vapour.append(_sig(100.0 * ps.x_vap) if ps.separated else None)
    run["phase_state"] = state
    if any(v is not None for v in share):
        run["phases"].update(vapour_share=share, brine_salinity=brine, vapour_salinity=vapour)


def phase_notes(run: dict, points: list[tuple[float, float]], mode: str = "",
                follow: str = "whole") -> list[str]:
    """Where the fluid splits into brine and vapour (NaCl-H2O, Driesner and Heinrich 2007), or meets
    halite, and where the path leaves the table.  ``follow`` is the card's follow_phase."""
    state = run.get("phase_state", [])
    where = f" ({mode.lower()})" if mode else ""
    out = []
    split = [i for i, s in enumerate(state) if s == "liquid + vapour"]
    if split:
        ph = run["phases"]
        i0 = split[0]
        top = max(split, key=lambda i: ph["vapour_share"][i] or 0.0)
        t, p = points[i0]
        x0 = next(v for v in ph["salinity"] if v is not None)
        out.append(f"Phase separation{where}: by its NaCl-equivalent salinity ({x0:.1f} wt%) the fluid "
                   f"splits into brine and vapour from about {t:.0f} °C and {p:.0f} bar (Driesner and Heinrich 2007). "
                   f"At most {ph['vapour_share'][top]:.0f}% of its mass is vapour with "
                   f"{ph['vapour_salinity'][top]:.2g} wt% NaCl, against {ph['brine_salinity'][top]:.3g} wt% in the brine. "
                   + ("The chemistry shown is that of the whole fluid; the Brine and vapour tab has the split."
                      if follow == "whole" else _follow_note(run, follow)))
    halite = [i for i, s in enumerate(state) if s in ("vapour + halite", "liquid + halite")]
    if halite:
        t, p = points[halite[0]]
        out.append(f"By its salinity the fluid would be saturated with halite{where} from about {t:.0f} °C and "
                   f"{p:.0f} bar (NaCl-H2O alone; the calculation's own halite saturation index is the one to trust).")
    outside = [i for i, s in enumerate(state) if s == "outside"]
    if outside and len(outside) < len(state):
        out.append("Phase separation is checked only from 10 to 490 bar (the range of the NaCl-H2O table); "
                   "part of this path is outside it.")
    elif outside and any(t > 300 for t, _ in points):
        out.append("Phase separation is checked only from 10 to 490 bar (the range of the NaCl-H2O table), so it "
                   "is not checked on this path.")
    return out


def _follow_note(run: dict, follow: str) -> str:
    """The phase-separation note's last sentence when the run follows the brine or the vapour."""
    gone = [v for v in run["phases"].get("split_off_kg", []) if v is not None]
    left = f" By the end {gone[-1]:.2g} kg of fluid has left." if gone and gone[-1] > 0 else ""
    if follow == "brine":
        return ("The chemistry shown is that of the brine: at each step the vapour boils off and leaves, taking "
                "the water and dissolved gases with it." + left + " Brines above about 20 wt% NaCl are "
                "beyond the activity model, so points there may be left blank.")
    return ("The chemistry shown is that of the vapour: at each step the brine stays behind with most of the "
            "salt; dissolved gases go with the water." + left)


def _bands(run: dict, samples: list[dict]) -> dict:
    """{group: {key: [low, high]}}: per point, the BAND_PERCENTILES of the samples, widened to take in
    the run's own value, for every line the run has."""
    out: dict[str, dict] = {}
    if len(samples) < 3:
        return out
    for group, lines in run.items():
        if not isinstance(lines, dict) or group == "bands":
            continue
        fill = 0.0 if group in AMOUNT_GROUPS else None
        here = out.setdefault(group, {})
        for key, line in lines.items():
            cols = [s.get(group, {}).get(key) for s in samples]
            lo, hi = [], []
            for i, v in enumerate(line):
                vals = [c[i] if c is not None else fill for c in cols]
                vals = [x for x in vals if x is not None]
                if v is None or len(vals) < 3:
                    lo.append(None)
                    hi.append(None)
                    continue
                a, b = np.percentile(vals, BAND_PERCENTILES)
                lo.append(_sig(min(float(a), v)))
                hi.append(_sig(max(float(b), v)))
            here[key] = [lo, hi]
    return out


def _catalog(db, runs: dict[str, dict], sc: Scenario) -> dict:
    """Display names for everything that can be plotted."""
    species, elements, minerals, gases = {}, {}, {}, {}
    blocked = {b.lower() for b in sc.blocked_minerals}
    for run in runs.values():
        for k in run["species"]:
            sp = db.species.get(k)
            z = sp.charge if sp is not None else None
            species[k] = {"label": pretty_species(k, z), "charge": z, "mw": _mw(db, k)}
        for k in run["elements"]:
            elements[k] = {"label": k, "mw": db.elements.get(k)}
        for k in run.get("gas", {}):
            gases[k] = {"label": pretty_species(k[:-3], 0).replace("(aq)", "") + " gas",
                        "mw": _mw(db, k)}
        for group in ("minerals", "saturation"):
            for k in run[group]:
                minerals[k] = {"label": pretty_mineral(k), "blocked": k.lower() in blocked, "mw": _mw(db, k)}
    kinetic = {}
    for run in runs.values():
        for k in run.get("kinetic", {}):
            what, name = k.split(":", 1)
            shown = pretty_mineral(name) if name in db.species else name
            kinetic[k] = {"label": f"{shown} rate" if what == "rate" else f"{shown} left",
                          "unit": ("mol/s" if name in db.species else "g/s") if what == "rate" else "% of start"}
    composition = {}
    for run in runs.values():
        for k in run.get("composition", {}):
            ss, member = k.split(":", 1)
            composition[k] = {"label": f"{pretty_mineral(member)} in {ss}", "unit": "mole fraction"}
    phases = {"salinity": {"label": "Salinity of the whole fluid", "unit": "wt% NaCl"},
              "vapour_share": {"label": "Vapour share of the fluid", "unit": "% by mass"},
              "brine_salinity": {"label": "Brine salinity", "unit": "wt% NaCl"},
              "vapour_salinity": {"label": "Vapour salinity", "unit": "wt% NaCl"},
              "split_off_kg": {"label": ("Vapour boiled off" if sc.follow_phase == "brine" else "Brine left behind"),
                               "unit": "kg per kg of starting water"}}
    isotope = {}
    words = {"water": "water", "DIC": "dissolved inorganic carbon", "CH4": "dissolved methane",
             "sulfate": "dissolved sulfate", "sulfide": "dissolved sulfide", "fluid": "fluid"}
    for run in runs.values():
        for k in run.get("isotopes", {}):
            what, carrier = k.split(":", 1)
            shown = words.get(carrier) or pretty_mineral(carrier)
            if carrier not in words and not shown[1:2].isupper():
                shown = shown[0].lower() + shown[1:]
            if what == "87Sr/86Sr":
                isotope[k] = {"label": "87Sr/86Sr of the fluid", "unit": "ratio"}
            else:
                isotope[k] = {"label": f"δ{what[1:]} {shown}", "unit": f"‰ {isotopes.STANDARDS[what]}"}
    sorbed = {}
    for run in runs.values():
        for k in run.get("sorbed", {}):
            if k.startswith("Hfo_"):
                site = "strong" if k.startswith("Hfo_s") else "weak"
                sorbed[k] = {"label": f"{k[5:]} on ferric oxide ({site} sites)", "unit": "mol"}
            else:
                ion = re.sub(r"X\d*$", "", k)
                sorbed[k] = {"label": f"{ion} on exchange sites ({k})", "unit": "mol"}
    return {"species": species, "elements": elements, "minerals": minerals, "gas": gases,
            "kinetic": kinetic, "composition": composition, "phases": phases, "isotopes": isotope,
            "sorbed": sorbed}


def _mw(db, name: str) -> float | None:
    """Formula weight (g/mol) for the mass units of the results window."""
    try:
        return _sig(db.molecular_weight(name))
    except (KeyError, TypeError):
        return None

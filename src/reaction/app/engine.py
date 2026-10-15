"""Run a scenario card on the equilibrium engine and package the numbers for the results window.

A heating run takes a fluid, warms it in small temperature steps at constant pressure, and at every
step finds the equilibrium of the whole closed system: the aqueous species, plus whichever minerals
are stable.  In *static* mode minerals stay and can redissolve.  In *flow-through* mode the minerals
formed at each step are left behind, so only the fluid goes on to the next step.
"""
from __future__ import annotations

import datetime as _dt
import math
import threading
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from ..database import load_data0
from ..database.model import Kind
from ..equilibrium import EquilibriumState, System, equilibrate
from ..network import PROTON, WATER
from ..path import _composition_row
from ..speciation import SpeciationError, speciate
from ..thermo.eosdb import load_eos_database
from . import config
from .fluids import preset
from .names import pretty_mineral, pretty_species
from .scenario import ALWAYS_EXCLUDED, MODE_LABELS, Scenario

Progress = Callable[[float, str], None]


class EngineUnavailable(RuntimeError):
    """The thermodynamic data the engine needs could not be found."""


def _sig(x: float | None, digits: int = 6):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return None
    return float(f"{x:.{digits}g}")


@dataclass
class _Context:
    bulk: dict
    state: EquilibriumState | None = None
    t_c: float | None = None
    cumulative: dict | None = None            # flow-through: everything left behind so far


class Engine:
    def __init__(self, data_dir: str | Path | None = None):
        self.data_dir = Path(data_dir) if data_dir else config.external_dir()
        self._db = None
        self._load_lock = threading.Lock()
        self.run_lock = threading.Lock()

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

    # ---- one run --------------------------------------------------------------------------------
    def run(self, sc: Scenario, progress: Progress | None = None) -> dict:
        problems = sc.problems()
        if problems:
            raise ValueError("; ".join(problems))
        progress = progress or (lambda f, m: None)
        started = time.time()
        db = self.db
        fluid = preset(sc.fluid)
        progress(0.0, f"Speciating {fluid.label.lower()}")
        sp0 = speciate(db, fluid.measured_t_c, fluid.constraints, activity_model=sc.activity_model)
        basis = [b for b in sp0.basis_names if b not in (WATER, PROTON)]
        system = System(db, basis, activity_model=sc.activity_model)
        bulk0 = dict(sp0.component_totals)
        temps = [float(t) for t in np.linspace(sc.t_start_c, sc.t_end_c, sc.steps)]
        modes = ["static", "flow_through"] if sc.mode == "both" else [sc.mode]

        notes: list[str] = []
        runs: dict[str, dict] = {}
        for k, mode in enumerate(modes):
            def tick(i, n, k=k):
                progress((k + i / n) / len(modes), f"{MODE_LABELS[modes[k]]}: {i} of {n} temperatures")
            runs[mode] = self._march(system, bulk0, temps, sc, mode == "flow_through", tick, notes)

        for run in runs.values():
            notes.extend(run.pop("warnings"))
        start_elements: dict[str, float] = {}
        for name, st in sp0.species.items():
            sp = db.species.get(name)
            for el, k in (sp.elements.items() if sp is not None else ()):
                if el not in ("H", "O"):
                    start_elements[el] = start_elements.get(el, 0.0) + k * st.molality
        names = _catalog(db, system, runs, sc)
        known = {n.lower() for n in db.species}
        blocked_unknown = [b for b in sc.blocked_minerals if b.lower() not in known]
        if blocked_unknown:
            notes.append("Not in the database, so ignored: " + ", ".join(blocked_unknown))
        notes.extend(fluid.caveats)
        seen = []
        for n in notes:
            if n not in seen:
                seen.append(n)
        progress(1.0, "Done")
        return {
            "id": _dt.datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:4],
            "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "scenario": sc.to_dict(),
            "summary": sc.describe(),
            "database": f"{db.title}; ion sizes from EQ3/6 data0.cmp",
            "fluid": {"key": fluid.key, "label": fluid.label, "description": fluid.description,
                      "measured_t_c": fluid.measured_t_c, "source": fluid.source,
                      "totals_molal": {pretty_species(k, db.species[k].charge if k in db.species else None):
                                       _sig(v) for k, v in sp0.totals.items()}},
            "x": {"key": "temperature_c", "label": "Temperature", "unit": "°C",
                  "values": [_sig(t) for t in temps]},
            "start_elements": {k: _sig(v) for k, v in start_elements.items()},
            "catalog": names,
            "modes": runs,
            "notes": seen,
            "elapsed_s": round(time.time() - started, 1),
        }

    # ---- the temperature march ------------------------------------------------------------------
    def _march(self, system: System, bulk0: dict, temps: list[float], sc: Scenario,
               flow_through: bool, tick, notes: list[str]) -> dict:
        db = self.db
        ctx = _Context(dict(bulk0), cumulative={} if flow_through else None)
        rows: list[dict | None] = []
        n = len(temps)
        for i, t in enumerate(temps):
            tick(i, n)
            try:
                state = self._advance(system, ctx, t, sc)
            except SpeciationError as exc:
                notes.append(f"No equilibrium found at {t:.0f} C ({exc}); that point is left blank.")
                rows.append(None)
                continue
            rows.append(self._record(db, state, ctx))
        tick(n, n)
        return _pack(rows)

    def _advance(self, system: System, ctx: _Context, t_c: float, sc: Scenario, depth: int = 0):
        """Move the system to ``t_c``; on a failed solve, go there in two half steps."""
        model = system.model(t_c, sc.pressure_bar)
        blocked = {b.lower() for b in sc.blocked_minerals} | set(ALWAYS_EXCLUDED)
        eligible = [name for name, kind in zip(model.ph_names, model.ph_kind)
                    if kind == Kind.MINERAL and name.lower() not in blocked]
        start = [] if ctx.state is None or ctx.cumulative is not None else list(ctx.state.minerals)
        try:
            state = equilibrate(system, t_c, ctx.bulk, pressure_bar=sc.pressure_bar, minerals=start,
                                eligible=eligible, guess=ctx.state)
        except SpeciationError:
            if ctx.state is None or ctx.t_c is None or depth >= 6:
                raise
            self._advance(system, ctx, 0.5 * (ctx.t_c + t_c), sc, depth + 1)
            return self._advance(system, ctx, t_c, sc, depth + 1)
        if ctx.cumulative is not None:                       # flow-through: leave the minerals behind
            for name, moles in state.minerals.items():
                ctx.cumulative[name] = ctx.cumulative.get(name, 0.0) + moles
                for comp, c in _composition_row(system, name, t_c).items():
                    ctx.bulk[comp] = ctx.bulk.get(comp, 0.0) - c * moles
        ctx.state, ctx.t_c = state, t_c
        return state

    def _record(self, db, state: EquilibriumState, ctx: _Context) -> dict:
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
        minerals = dict(ctx.cumulative) if ctx.cumulative is not None else dict(state.minerals)
        return {"pH": spec.ph, "ionic_strength": spec.ionic_strength, "water_kg": state.water_kg,
                "log_aw": spec.log_a_water, "species": species, "elements": elements,
                "minerals": minerals, "saturation": dict(spec.saturation),
                "warnings": list(spec.warnings)}


def _pack(rows: list[dict | None]) -> dict:
    """Turn a list of per-temperature records into one list per quantity (None where unsolved)."""
    good = [r for r in rows if r is not None]

    def column(getter, keys):
        return {k: [_sig(getter(r, k)) if r is not None else None for r in rows] for k in keys}

    species_keys = sorted({k for r in good for k, v in r["species"].items() if v > 0})
    species_keys = [k for k in species_keys if max(r["species"].get(k, 0.0) for r in good) >= 1e-12]
    element_keys = sorted({k for r in good for k in r["elements"]})
    mineral_keys = sorted({k for r in good for k, v in r["minerals"].items() if v > 0})
    sat_keys = sorted({k for r in good for k, v in r["saturation"].items()
                       if math.isfinite(v) and v > -3.0})
    return {
        "properties": {
            "pH": [_sig(r["pH"]) if r else None for r in rows],
            "ionic_strength": [_sig(r["ionic_strength"]) if r else None for r in rows],
            "water_kg": [_sig(r["water_kg"]) if r else None for r in rows],
            "log_aw": [_sig(r["log_aw"]) if r else None for r in rows],
        },
        "species": column(lambda r, k: r["species"].get(k, 0.0), species_keys),
        "elements": column(lambda r, k: r["elements"].get(k, 0.0), element_keys),
        "minerals": column(lambda r, k: r["minerals"].get(k, 0.0), mineral_keys),
        "saturation": column(lambda r, k: r["saturation"].get(k), sat_keys),
        "warnings": sorted({w for r in good for w in r["warnings"]}),
    }


def _catalog(db, system: System, runs: dict[str, dict], sc: Scenario) -> dict:
    """Display names for everything that can be plotted."""
    species, elements, minerals = {}, {}, {}
    blocked = {b.lower() for b in sc.blocked_minerals}
    for run in runs.values():
        for k in run["species"]:
            sp = db.species.get(k)
            z = sp.charge if sp is not None else None
            species[k] = {"label": pretty_species(k, z), "charge": z}
        for k in run["elements"]:
            elements[k] = {"label": k}
        for group in ("minerals", "saturation"):
            for k in run[group]:
                minerals[k] = {"label": pretty_mineral(k), "blocked": k.lower() in blocked}
    return {"species": species, "elements": elements, "minerals": minerals}

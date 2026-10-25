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
import math
import threading
import time
import uuid
from pathlib import Path
from typing import Callable

from ..database import load_data0
from ..database.model import Kind
from ..equilibrium import System
from ..network import PROTON, WATER
from ..path import PathPoint, ReactionPath, Reactant, mineral_reactant
from ..speciation import PH, Couple, Eh, LogFO2, Molal, Pe, Speciation, speciate
from ..thermo import validity
from ..thermo.eosdb import DEFAULT_MASTERS, load_eos_database
from . import analysis, config, rocks
from .fluids import AIR_LOG_FO2, preset
from .names import pretty_mineral, pretty_species
from .scenario import ALWAYS_EXCLUDED, MODE_LABELS, Scenario

Progress = Callable[[float, str], None]

TRACE = 1e-12          # mol per kg of water of a component only a reactant brings in
AXES = {"temperature_c": ("Temperature", "°C"), "pressure_bar": ("Pressure", "bar")}


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
        started = time.time()
        db = self.db
        fluid = sc.fluid_spec()
        progress(0.0, f"Speciating {fluid.label.lower()}")
        sp0, fluid_notes = self.starting_fluid(sc)
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
        mix = None
        if sc.mixing:
            progress(0.0, f"Speciating {preset(sc.mix_with).label.lower()}")
            mix, mix_notes = self.starting_fluid(Scenario(fluid=sc.mix_with, activity_model=sc.activity_model))
            fluid_notes += [f"Second fluid: {n}" for n in mix_notes]
            basis += [b for b in mix.basis_names if b not in basis and b not in (WATER, PROTON)]
        system = System(db, basis, activity_model=sc.activity_model)
        bulk0 = dict(sp0.component_totals)
        for b in basis:
            bulk0.setdefault(b, TRACE)                   # a trace, so the solution has a start
        points = sc.points()
        modes = ["static", "flow_through"] if sc.mode == "both" else [sc.mode]
        rels = self.reliability(points)

        notes: list[str] = list(fluid_notes)
        runs: dict[str, dict] = {}
        for k, mode in enumerate(modes):
            def tick(i, n, k=k):
                progress((k + i / n) / len(modes), f"{MODE_LABELS[modes[k]]}: {i} of {n} points")
            runs[mode] = self._march(system, bulk0, points, sc, mode == "flow_through", tick, notes, mix)

        for run in runs.values():
            notes.extend(run.pop("warnings"))
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
        if sc.reactants:
            extra = self._reaction_axes(sc, added, len(points))
            axes.update(extra)
            varies = list(extra) + varies
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

    def _reactant_names(self, sc: Scenario) -> dict[str, str]:
        """The database name of each reactant on the card (matched without regard to case)."""
        db = self.db
        lower: dict[str, list[str]] = {}
        for name, sp in db.species.items():
            lower.setdefault(name.lower(), []).append(name)
        out = {}
        for r in sc.reactants:
            name = r["name"]
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
            if db.species[name].kind not in (Kind.MINERAL, Kind.AQUEOUS, Kind.BASIS):
                raise ValueError(f"{name} cannot be a reactant (only minerals and dissolved species)")
            out[r["name"]] = name
        return out

    def _reactants(self, system: System, sc: Scenario, names: dict[str, str], steps: int,
                   notes: list[str] | None = None) -> list:
        """Each reactant added at a steady rate, all of it by the last step (reaction progress counts
        the steps)."""
        out = []
        for r in sc.reactants:
            name = names[r["name"]]
            if name.startswith("rock:"):                 # dissolves congruently, never precipitates
                rock = rocks.ROCKS[name[5:]]
                comp, rock_notes = rocks.components_per_gram(rock, self.db, system.basis)
                if notes is not None:
                    notes.extend(n for n in rock_notes if n not in notes)
                out.append(Reactant(rock.label, comp, r["amount"],
                                    (r["amount"] / max(steps - 1, 1), 0.0, 0.0)))
                continue
            if name not in system.net.phase_names and name not in system.net.dep_names:
                raise ValueError(f"{name} cannot be added: its formula needs the redox state of the "
                                 f"fluid, which is not tracked; set the redox state on the card")
            moles = r["amount"] / self.db.molecular_weight(name) if r["unit"] == "g" else r["amount"]
            out.append(mineral_reactant(system, name, moles, rates=(moles / max(steps - 1, 1), 0.0, 0.0)))
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
        folded: dict[str, float] = {}
        for name, m in molal.items():
            target = analysis.REDOX.get(name, name) if name in analysis.REDOX else name
            if target is not None:
                folded[target] = folded.get(target, 0.0) + m
        constraints = {name: Molal(m) for name, m in folded.items()}
        constraints[PROTON] = PH(fluid.ph)
        redox = {"air": lambda: LogFO2(AIR_LOG_FO2), "log_fo2": lambda: LogFO2(arg),
                 "eh": lambda: Eh(arg), "pe": lambda: Pe(arg)}.get(kind, lambda: None)()

        def solve(balance_on=None):
            if kind != "couple":
                return speciate(db, fluid.measured_t_c, constraints, redox=redox,
                                balance_on=balance_on, activity_model=sc.activity_model)
            own = dict(constraints)
            partner = analysis.REDOX[arg]
            if partner is not None:
                rest = folded[partner] - molal[arg]
                own[partner] = Molal(rest)
            own[arg] = Molal(molal[arg])
            first = speciate(db, fluid.measured_t_c, own, redox=Couple(arg),
                             balance_on=balance_on, activity_model=sc.activity_model)
            return speciate(db, fluid.measured_t_c, constraints, redox=LogFO2(first.log_fo2),
                            balance_on=balance_on, activity_model=sc.activity_model)
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
        blocked = {b.lower() for b in sc.blocked_minerals} | set(ALWAYS_EXCLUDED)
        suppress = [name for name, kind in zip(system.net.phase_names, system.net.phase_kind)
                    if kind == Kind.MINERAL and name.lower() in blocked]
        reactants = self._reactants(system, sc, self._reactant_names(sc), n, notes)
        if mix is not None:          # the second fluid, per kg of its water, as one more reactant
            composition = {k: v for k, v in mix.component_totals.items() if k in system.basis}
            reactants.append(Reactant(f"mix:{sc.mix_with}", composition, sc.mix_ratio,
                                      (sc.mix_ratio / last, 0.0, 0.0)))
        path = ReactionPath(system, bulk0, reactants=reactants,
                            temperature=lambda xi: sc.tp(xi / last)[0],
                            pressure=lambda xi: sc.tp(xi / last)[1],
                            suppress_minerals=suppress, flow_through=flow_through)

        def failed(xi, exc):
            t, p = points[int(round(xi))]
            notes.append(f"No equilibrium found at {t:.0f} C and {p:.0f} bar ({exc}); that point is "
                         f"left blank.")

        solved = path.run(list(range(n)), on_failure=failed, progress=tick)
        return _pack([None if pt is None else self._record(db, pt) for pt in solved])

    def _record(self, db, point: PathPoint) -> dict:
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
        minerals = dict(point.removed)                    # flow-through: everything left behind so far
        for name, moles in state.minerals.items():
            minerals[name] = minerals.get(name, 0.0) + moles
        return {"pH": spec.ph, "ionic_strength": spec.ionic_strength, "water_kg": state.water_kg,
                "log_fO2": spec.log_fo2, "Eh": spec.eh_volts,
                "log_aw": spec.log_a_water, "species": species,
                "log_a": {n: spec.species[n].log_activity for n in spec.basis_names
                          if n in spec.species and n not in (WATER, PROTON)}, "elements": elements,
                "minerals": minerals, "saturation": dict(spec.saturation),
                "warnings": list(spec.warnings)}


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
        },
        "species": column(lambda r, k: r["species"].get(k, 0.0), species_keys),
        "log_a": column(lambda r, k: r.get("log_a", {}).get(k), sorted({k for r in good for k in r.get("log_a", {})})),
        "elements": column(lambda r, k: r["elements"].get(k, 0.0), element_keys),
        "minerals": column(lambda r, k: r["minerals"].get(k, 0.0), mineral_keys),
        "saturation": column(lambda r, k: r["saturation"].get(k), sat_keys),
        "warnings": sorted({w for r in good for w in r["warnings"]}),
    }


def _catalog(db, runs: dict[str, dict], sc: Scenario) -> dict:
    """Display names for everything that can be plotted."""
    species, elements, minerals = {}, {}, {}
    blocked = {b.lower() for b in sc.blocked_minerals}
    for run in runs.values():
        for k in run["species"]:
            sp = db.species.get(k)
            z = sp.charge if sp is not None else None
            species[k] = {"label": pretty_species(k, z), "charge": z, "mw": _mw(db, k)}
        for k in run["elements"]:
            elements[k] = {"label": k, "mw": db.elements.get(k)}
        for group in ("minerals", "saturation"):
            for k in run[group]:
                minerals[k] = {"label": pretty_mineral(k), "blocked": k.lower() in blocked, "mw": _mw(db, k)}
    return {"species": species, "elements": elements, "minerals": minerals}


def _mw(db, name: str) -> float | None:
    """Formula weight (g/mol) for the mass units of the results window."""
    try:
        return _sig(db.molecular_weight(name))
    except (KeyError, TypeError):
        return None

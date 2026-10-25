"""Run an EQ6 input through Reaction's path engine and compare with EQ6's CSV tables."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from ..database.model import Kind, ThermoDatabase
from ..equilibrium import System
from ..path import PathPoint, Reactant, ReactionPath, _composition_row, mineral_reactant
from .eq6_csv import Eq6Csv, read_eq6_csv
from .eq6_input import Eq6Input, read_eq6_input
from .eq3nr_run import Unsupported


@dataclass
class Eq6Setup:
    system: System
    path: ReactionPath
    xis: list[float]
    notes: list[str] = field(default_factory=list)


def build_case(db: ThermoDatabase, inp: Eq6Input, csv: Eq6Csv) -> Eq6Setup:
    if inp.activity_model not in ("bdot", "davies"):
        raise Unsupported(f"activity model {inp.activity_model}")
    if inp.iopt.get(1, 0) == 1:
        raise Unsupported("titration system model (iopt(1) = 1)")
    if inp.iopt.get(2, 0) == 1:
        raise Unsupported("kinetic (time) mode: milestone M4")
    if inp.iopt.get(4, 0) == 1:
        raise Unsupported("solid solutions permitted (iopt(4) = 1)")
    jt, base, deriv, extras = inp.temperature
    if jt not in (0, 1):
        raise Unsupported(f"temperature option {jt}")
    if inp.pressure[0] != 0:
        raise Unsupported(f"pressure option {inp.pressure[0]}")
    if inp.switches:
        raise Unsupported("special basis switches")
    if inp.other_alter:
        raise Unsupported(f"alter options {inp.other_alter}")

    basis = [c for c in csv["E1"].columns]
    if "H2O" not in basis:
        raise Unsupported("no H2O in the CSV basis")
    bulk = {}
    for name in basis:
        if name not in inp.totals:
            raise Unsupported(f"basis species {name} missing from the mass-balance totals")
        bulk[name] = inp.totals[name][0]

    mineral_names = {s.name for s in db.of_kind(Kind.MINERAL, Kind.LIQUID)}
    species_names = {s.name for s in db.of_kind(Kind.AUX_BASIS, Kind.AQUEOUS)}
    suppress_phases = [n for n in inp.suppress if n in mineral_names or
                       (n in db and db[n].kind in (Kind.GAS,))]
    suppress_species = [n for n in inp.suppress if n in species_names]
    system = System(db, basis, suppress=suppress_species, suppress_phases=suppress_phases,
                    activity_model=inp.activity_model, ph_scale=inp.ph_scale)

    if any(opt not in ("All",) for opt, _ in inp.subset_suppression):
        raise Unsupported(f"mineral subset suppression {inp.subset_suppression}")
    no_minerals = any(opt == "All" for opt, _ in inp.subset_suppression)

    # non-basis species (O2(aq), H2(aq)) are reported apart from the basis totals: add their
    # contribution to the components they are made of
    for sp_name in inp.non_basis:
        amount = inp.totals.get(sp_name, (0.0, 0.0))[0]
        if amount:
            for comp, c in _composition_row(system, sp_name).items():
                bulk[comp] = bulk.get(comp, 0.0) + c * amount
    reactants: list[Reactant] = []
    for r in inp.reactants:
        if r.forward_law != "Relative rate equation":
            raise Unsupported(f"reactant {r.name}: rate law {r.forward_law}")
        if r.kind in ("Pure mineral", "Gas species", "Aqueous species"):
            reactants.append(mineral_reactant(system, r.name, r.remaining, rates=r.rkb))
            if r.kind != "Pure mineral":
                reactants[-1].mineral = None
        elif r.kind == "Special reactant":
            comp = {k: v for k, v in r.reaction.items() if k != r.name}
            missing = [k for k in comp if k not in system.basis]
            if missing:
                raise Unsupported(f"special reactant {r.name}: {missing} not in the basis")
            reactants.append(Reactant(r.name, comp, r.remaining, r.rkb))
        else:
            raise Unsupported(f"reactant type {r.kind}")

    fixed: dict[str, float] = {}
    buffers: list[tuple[str, float]] = []
    if inp.iopt.get(15, 0) == 1:
        fixed["O2(g)"] = inp.log_moles.get("O2(g)", -0.7)
    for gas, _, logf in inp.fixed_fugacity:
        if gas == "O2(g)":
            fixed["O2(g)"] = logf
        else:
            buffers.append((gas, logf))

    temperature = (lambda xi, b=base, d=deriv: b + d * xi) if jt == 1 else base
    path = ReactionPath(system, bulk, reactants=reactants, temperature=temperature,
                        fixed=fixed, buffers=buffers, flow_through=inp.iopt.get(1, 0) == 2,
                        no_precipitation=no_minerals)
    return Eq6Setup(system, path, list(csv.xi))


# --------------------------------------------------------------------------------------------
@dataclass
class Deviation:
    what: str
    xi: float
    ours: float
    reference: float

    @property
    def delta(self) -> float:
        return abs(self.ours - self.reference)


@dataclass
class Eq6Report:
    name: str
    description: str
    n_points: int
    deviations: list[Deviation] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)

    def worst(self, kind: str | None = None, n: int = 4) -> list[Deviation]:
        pool = [d for d in self.deviations if kind is None or d.what.startswith(kind)]
        return sorted(pool, key=lambda d: -d.delta)[:n]


# Tolerances (absolute, relative).  The CSV prints five significant figures, reaction
# progress included, so amounts of a mineral that has only just begun to precipitate carry
# an error of about 1e-3 * xi from the rounding of xi alone.  The solvent mass is printed with three
# decimals (1.000 kg for 1.0004), hence its absolute 5e-4.
TOL = {
    "pH": (3e-4, 0.0), "pmH": (3e-4, 0.0), "log fO2": (5e-3, 5e-4), "aw": (2e-4, 0.0),
    "T": (5e-3, 0.0), "H2O(kg)": (5e-4, 3e-4), "I": (0.0, 2e-3),
    "element": (1e-15, 2e-3), "mineral": (1e-10, 3e-3), "SI": (2e-3, 0.0), "log m": (2e-3, 5e-4),
}


def exceeds(what: str, d: Deviation) -> bool:
    key = next(k for k in TOL if what.startswith(k))
    ab, rel = TOL[key]
    if key == "mineral":
        ab = max(ab, 1e-3 * abs(d.xi))
    return d.delta > ab + rel * abs(d.reference)


def _element_molalities(point: PathPoint, db: ThermoDatabase) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, st in point.speciation.species.items():
        sp = db.species.get(name)
        if sp is None:
            continue
        for el, n in sp.elements.items():
            out[el] = out.get(el, 0.0) + n * st.molality
    return out


def compare_case(db: ThermoDatabase, name: str, description: str, points: list[PathPoint],
                 csv: Eq6Csv, fixed_o2: bool) -> Eq6Report:
    rep = Eq6Report(name, description, len(points))
    b1, b2, c1, e1, pt, q = csv["B1"], csv["B2"], csv["C1"], csv["E1"], csv["P"], csv["Q"]

    def add(what, xi, ours, ref):
        if ref is None or ours is None or (isinstance(ours, float) and math.isnan(ours)):
            return
        rep.deviations.append(Deviation(what, xi, float(ours), float(ref)))

    for i, pnt in enumerate(points):
        xi, spec = pnt.xi, pnt.speciation
        row1, row2 = b1.row_at(i), b2.row_at(i)
        add("T", xi, pnt.temperature_c, row1["Temp(C)"])
        add("pH", xi, spec.ph, row1["pH"])
        add("pmH", xi, spec.pmh, row1["pmH"])
        if not fixed_o2:
            add("log fO2", xi, spec.log_fo2, row1["log fO2"])
        add("aw", xi, 10 ** spec.log_a_water, row1["aw"])
        add("H2O(kg)", xi, pnt.state.water_kg, row2["H2O(kg)"])
        add("I", xi, spec.ionic_strength, row2["I(m)"])
        els = _element_molalities(pnt, db)
        for el, ref in c1.row_at(i).items():
            if el in els and ref is not None:
                add(f"element {el}", xi, els[el], ref)
        for col, ref in e1.row_at(i).items():
            if col == "H2O" or col == "O2(g)" or ref is None:
                continue
            st = spec.species.get(col)
            if st is not None and st.molality > 0:
                add(f"log m {col}", xi, math.log10(st.molality), ref)
        prow = pt.row_at(i)
        names = set(prow) | set(pnt.minerals)
        for mn in names:
            ours = pnt.minerals.get(mn, 0.0)
            add(f"mineral {mn}", xi, ours, prow.get(mn, 0.0) or 0.0)
        for mn, ref in q.row_at(i).items():
            if ref is not None and mn in spec.saturation and mn not in ("H2O", "Gas"):
                add(f"SI {mn}", xi, spec.saturation[mn], ref)
    return rep


def _printed_xis(output_path: Path, csv_xis: list[float]) -> list[float]:
    """The CSV's print points, except that two points the CSV prints alike (1.02047E-04 and 1.02049E-04
    both read 1.0205E-04 there) take the extra digit of the .6o text output.  Only those: at the other
    boundary points EQ6 locates, the extra digit moves the comparison across our own boundary as often
    as it helps, so the CSV values stay."""
    if not output_path.exists():
        return csv_xis
    found = []
    for line in output_path.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] == "Xi=":
            try:
                found.append(float(parts[1]))
            except ValueError:
                return csv_xis
    if len(found) != len(csv_xis):
        return csv_xis
    if any(abs(a - b) > 1e-4 * max(abs(a), abs(b)) for a, b in zip(found, csv_xis)):
        return csv_xis
    n = len(csv_xis)
    tied = [i for i in range(n) if (i > 0 and csv_xis[i] == csv_xis[i - 1])
            or (i + 1 < n and csv_xis[i] == csv_xis[i + 1])]
    return [found[i] if i in tied else x for i, x in enumerate(csv_xis)]


def run_case(db: ThermoDatabase, input_path: str | Path,
             csv_path: str | Path | None = None) -> tuple[Eq6Report, list[PathPoint]]:
    input_path = Path(input_path)
    inp = read_eq6_input(input_path)
    csv_path = Path(csv_path) if csv_path else input_path.parent / "Outputs" / (input_path.stem + ".csv")
    csv = read_eq6_csv(csv_path)
    setup = build_case(db, inp, csv)
    xis = _printed_xis(input_path.parent / "Outputs" / (input_path.stem + ".6o"), setup.xis)
    points = setup.path.run(xis)
    fixed_o2 = "O2(g)" in setup.path.fixed
    return compare_case(db, input_path.stem, inp.description, points, csv, fixed_o2), points

"""A thermodynamic database built from equations of state (OBIGT: HKF + Maier-Kelley).

Unlike a database of tabulated logK values, this one evaluates logK at any
temperature and pressure from the Gibbs energies of the species:

    log K = -(sum_i nu_i G_i) / (R T ln 10)

Formation reactions are generated from the formulas: every species is written
over a set of *master* species (one per element, plus H2O, H+ and O2(g)).

The activity-model parameters (Debye-Hueckel A and B) come from the water model
at the requested T and P; B-dot is tabulated to 300 C (Helgeson, 1969) and held
at its 300 C value above that, which is flagged because it is an extrapolation.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from ..constants import LN10, R
from ..database.model import (GridLogK, Kind, Provenance, Species, ThermoDatabase,
                              ThermoRangeError)
from . import electrostatics, water
from .crystal import CrystalSpecies
from .hkf import HKFSpecies, WaterStateError
from pathlib import Path

from .obigt import (BERMAN_FILES, DEFAULT_FILES, SLOP98_FILES, SUPCRT92_FILES, ObigtData,
                    ObigtSpecies, load_obigt)

WATER, PROTON, OXYGEN = "H2O", "H+", "O2(g)"
CAL = 4.184

# standard atomic weights (g/mol)
ATOMIC_WEIGHTS = {
    "H": 1.00794, "Li": 6.941, "Be": 9.01218, "B": 10.811, "C": 12.011, "N": 14.0067,
    "O": 15.9994, "F": 18.9984, "Na": 22.98977, "Mg": 24.305, "Al": 26.98154,
    "Si": 28.0855, "P": 30.97376, "S": 32.066, "Cl": 35.4527, "K": 39.0983, "Ca": 40.078,
    "Ti": 47.88, "V": 50.9415, "Cr": 51.9961, "Mn": 54.93805, "Fe": 55.847, "Co": 58.9332,
    "Ni": 58.6934, "Cu": 63.546, "Zn": 65.39, "Ga": 69.723, "Ge": 72.61, "As": 74.9216,
    "Se": 78.96, "Br": 79.904, "Rb": 85.4678, "Sr": 87.62, "Zr": 91.224, "Mo": 95.94,
    "Ag": 107.8682, "Cd": 112.411, "Sn": 118.71, "Sb": 121.757, "I": 126.90447,
    "Cs": 132.90543, "Ba": 137.327, "W": 183.84, "Au": 196.96654, "Hg": 200.59, "Pb": 207.2,
    "U": 238.0289, "Ar": 39.948, "He": 4.0026, "Ne": 20.1797, "Kr": 83.8, "Xe": 131.29,
}

# master species by element (OBIGT names); H and O come from H2O, H+ and O2(g)
DEFAULT_MASTERS = {
    "Li": "Li+", "Na": "Na+", "K": "K+", "Rb": "Rb+", "Cs": "Cs+",
    "Be": "Be+2", "Mg": "Mg+2", "Ca": "Ca+2", "Sr": "Sr+2", "Ba": "Ba+2",
    "Cl": "Cl-", "F": "F-", "Br": "Br-", "I": "I-", "S": "SO4-2", "C": "HCO3-",
    "Si": "SiO2", "Al": "Al+3", "Fe": "Fe+2", "Mn": "Mn+2", "Zn": "Zn+2", "Cu": "Cu+2",
    "Pb": "Pb+2", "Cd": "Cd+2", "Ni": "Ni+2", "Co": "Co+2", "Ag": "Ag+", "Hg": "Hg+2",
    "N": "NH3", "P": "HPO4-2", "B": "B(OH)3",
}

# Helgeson (1969) B-dot, kg/mol, on the grid of the EQ3/6 databases
_BDOT_T = np.array([0.01, 25.0, 60.0, 100.0, 150.0, 200.0, 250.0, 300.0])
_BDOT = np.array([0.0394, 0.0410, 0.0438, 0.0460, 0.0470, 0.0470, 0.0340, 0.0000])
# Drummond (1981) salting-out polynomial for nonpolar neutral species
_CCO2 = np.array([-1.0312, 0.0012806, 255.9, 0.4445, -0.001606])
_NONPOLAR = {"CO2", "O2", "H2", "N2", "CH4", "H2S", "Ar", "He", "Ne", "Kr", "Xe", "NH3"}
_DEFAULT_ION_SIZE = 4.0
LOW_DENSITY_G_CM3 = 0.35     # below this the HKF/g-function calibration is poor
_DEFAULT_NEUTRAL_SIZE = 3.0

# the EQ3/6 convention: 1.0132 bar up to 100 C, water saturation pressure above
_ATM = 1.0132


def reference_pressure_bar(t_c: float) -> float:
    """EQ3/6 "reference curve": 1.0132 bar to 100 C, then water saturation."""
    if t_c <= 100.0:
        return _ATM
    t_k = t_c + 273.15
    if t_k >= water.T_CRIT:
        raise ThermoRangeError("no saturation pressure above the critical temperature; "
                               "give an explicit pressure")
    return water.saturation_pressure_mpa(t_k) * water.BAR_PER_MPA


# --- Gibbs-energy providers -----------------------------------------------------------------
class _Water:
    def gibbs(self, t_k, p_bar):
        return water.gibbs(t_k, p_bar, liquid=True)

    def extrapolated(self, t_k):
        return False


class _Aqueous:
    def __init__(self, eos: HKFSpecies):
        self.eos = eos

    def gibbs(self, t_k, p_bar):
        return self.eos.gibbs(t_k, p_bar)

    def extrapolated(self, t_k):
        return False


class _Crystal:
    """One phase, or several polymorphs (cr, cr2, ...) that each hold up to their transition
    temperature, as in CHNOSZ.  Above the limit of the last phase the phase has no data
    (CHNOSZ's ``exceed.Ttr = FALSE``): extrapolating a heat-capacity polynomial far past its
    range made, for example, nesquehonite stable in seawater at 300 C."""

    def __init__(self, eos_list: list[CrystalSpecies]):
        self.eos_list = eos_list

    def _phase(self, t_k):
        for e in self.eos_list:
            if not e.extrapolated(t_k):
                return e
        raise ThermoRangeError(f"{self.eos_list[-1].name} has no data above "
                               f"{self.eos_list[-1].t_max - 273.15:.0f} C")

    def gibbs(self, t_k, p_bar):
        return self._phase(t_k).gibbs(t_k, p_bar)

    def extrapolated(self, t_k):
        return False


class _Electron:
    """Aqueous electron: G = 0 at every temperature (the standard-hydrogen-electrode
    convention; the OBIGT entropy of e- only matters for the H2 couple)."""

    def gibbs(self, t_k, p_bar):
        return 0.0

    def extrapolated(self, t_k):
        return False


class ReactionLogK:
    """log10 K of a reaction from the Gibbs energies of its participants.

    Same call signature as :class:`GridLogK` so the reaction network can use it
    unchanged.  ``terms`` are (coefficient, provider) with the sign convention
    log K = -(sum coefficient * G) / (R T ln 10).
    """

    def __init__(self, terms):
        self.terms = list(terms)

    def __call__(self, t_c: float, p_bar: float | None = None) -> float:
        t_k = t_c + 273.15
        try:
            p = p_bar if p_bar is not None else reference_pressure_bar(t_c)
            total = sum(c * fn.gibbs(t_k, p) for c, fn in self.terms)
        except (WaterStateError, water.WaterError) as exc:
            raise ThermoRangeError(str(exc)) from exc
        return -total / (R * t_k * LN10)

    def extrapolated(self, t_c: float) -> bool:
        t_k = t_c + 273.15
        return any(fn.extrapolated(t_k) for _, fn in self.terms)


class _EhLogK(ReactionLogK):
    """2 H2O = O2(g) + 4 H+ + 4 e-."""

    def __init__(self, providers):
        h2o, o2, hp, e = providers
        super().__init__([(-2.0, h2o), (1.0, o2), (4.0, hp), (4.0, e)])   # products - reactants


# --- formula arithmetic ------------------------------------------------------------------------
def _formation(elements: dict, charge: float, masters: dict):
    """Coefficients of masters (and H2O, H+, O2) whose sum has the composition of a species.

    Returns a dict name -> coefficient, or None if an element has no master.
    """
    nu: dict[str, float] = {}
    h_left = elements.get("H", 0.0)
    o_left = elements.get("O", 0.0)
    z_left = charge
    for el, n in elements.items():
        if el in ("H", "O"):
            continue
        m = masters.get(el)
        if m is None:
            return None
        coef = n / m.elements[el]
        nu[m.name] = nu.get(m.name, 0.0) + coef
        h_left -= coef * m.elements.get("H", 0.0)
        o_left -= coef * m.elements.get("O", 0.0)
        z_left -= coef * m.charge
    nu_h = z_left
    nu_w = (h_left - nu_h) / 2.0
    nu_o2 = (o_left - nu_w) / 2.0
    for name, c in ((PROTON, nu_h), (WATER, nu_w), (OXYGEN, nu_o2)):
        if abs(c) > 1e-12:
            nu[name] = nu.get(name, 0.0) + c
    return nu


@dataclass
class _Master:
    name: str
    elements: dict
    charge: float


class EosDatabase(ThermoDatabase):
    """ThermoDatabase whose logK values come from equations of state."""

    def __init__(self, title: str, source: str, species: dict[str, Species],
                 oxygen: CrystalSpecies):
        providers = (_Water(), _Crystal([oxygen]), _Aqueous(_h_plus()), _Electron())
        super().__init__(
            title=title, source=source, temps_c=_BDOT_T.copy(),
            pressures_bar=np.zeros(8), adh=np.zeros(8), bdh=np.zeros(8), bdot=_BDOT.copy(),
            cco2=_CCO2.copy(), logk_eh=_EhLogK(providers),
            elements=dict(ATOMIC_WEIGHTS), species=species, t_limits=(0.0, 1000.0))
        self._bdot_fit = GridLogK(_BDOT_T, _BDOT, 0.0, 300.0)
        self.data_notes: list[str] = []     # about the data as loaded, shown with every run

    def reference_pressure(self, t_c: float) -> float:
        return reference_pressure_bar(t_c)

    def dh_params(self, t_c: float, p_bar: float | None = None):
        p = p_bar if p_bar is not None else reference_pressure_bar(t_c)
        try:
            a, b = electrostatics.debye_huckel_a_b(t_c + 273.15, p)
        except (WaterStateError, water.WaterError) as exc:
            raise ThermoRangeError(str(exc)) from exc
        bdot = self._bdot_fit(min(t_c, 300.0))
        return a, b, bdot

    def bdot_extrapolated(self, t_c: float) -> bool:
        return t_c > 300.0

    def condition_warnings(self, t_c: float, p_bar: float | None) -> list[str]:
        """Notes about (T, P) regions where the HKF equations are poorly calibrated."""
        p = p_bar if p_bar is not None else reference_pressure_bar(t_c)
        try:
            rho = water.density(t_c + 273.15, p, liquid=True) / 1000.0
        except water.WaterError:
            return []
        if rho < LOW_DENSITY_G_CM3:
            return [f"water density {rho:.2f} g/cm3 is below {LOW_DENSITY_G_CM3}: the HKF "
                    f"equations of state are extrapolated (use a higher pressure)"]
        return []


def _h_plus() -> HKFSpecies:
    return HKFSpecies("H+", 1.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0)


# --- construction ------------------------------------------------------------------------------
def _gas_name(sp: ObigtSpecies) -> str:
    return f"{sp.formula}(g)"


def build_eos_database(data: ObigtData, *, masters: dict[str, str] | None = None,
                       title: str = "OBIGT/HKF", exclude: set[str] | None = None) -> EosDatabase:
    """Assemble an :class:`EosDatabase` from loaded OBIGT species."""
    masters = dict(DEFAULT_MASTERS if masters is None else masters)
    exclude = exclude or set()

    # oxygen gas for the Eh reaction and O2(g) basis species
    oxygen = data.get("oxygen", "gas")
    if oxygen is None:
        raise ValueError("OBIGT data contain no oxygen gas entry")
    providers: dict[str, object] = {WATER: _Water(), OXYGEN: _Crystal([oxygen.eos]),
                                    PROTON: _Aqueous(_h_plus())}
    species: dict[str, Species] = {}
    prov = lambda sp, section: Provenance(  # noqa: E731
        database="OBIGT", file=sp.source, section=section, sources=[s for s in sp.refs if s and s != "NA"])

    def add_basis(name, formula, elements, charge, kind=Kind.BASIS, sp=None, section="basis"):
        species[name] = Species(
            name=name, kind=kind, formula=formula, charge=charge, elements=elements,
            mol_wt=None, a0=_DEFAULT_NEUTRAL_SIZE if charge == 0 else _DEFAULT_ION_SIZE,
            activity_type=0, reaction={}, logk=None, molar_volume=None, ref_state={},
            provenance=prov(sp, section) if sp is not None else Provenance("OBIGT", "built-in", section))

    add_basis(WATER, "H2O", {"H": 2.0, "O": 1.0}, 0.0)
    add_basis(PROTON, "H+", {"H": 1.0}, 1.0)
    species[PROTON].a0 = 9.0
    add_basis(OXYGEN, "O2", {"O": 2.0}, 0.0)

    master_objs: dict[str, _Master] = {}
    for el, mname in masters.items():
        sp = data.get(mname, "aq")
        if sp is None:
            continue
        master_objs[el] = _Master(mname, sp.elements, sp.charge)
        providers[mname] = _Aqueous(sp.eos)
        add_basis(mname, sp.formula, dict(sp.elements), sp.charge, sp=sp, section="master")
        species[mname].activity_type = -1 if mname in _NONPOLAR else 0

    def reaction_for(name: str, sp: ObigtSpecies, provider, kind: Kind, volume=None):
        nu = _formation(sp.elements, sp.charge, master_objs)
        if nu is None or any(k not in species for k in nu):
            return False
        terms = [(-1.0, provider)] + [(c, providers[k]) for k, c in nu.items()]
        reaction = {name: -1.0, **nu}
        if sp.is_aqueous:
            a0 = _DEFAULT_NEUTRAL_SIZE if sp.charge == 0 else _DEFAULT_ION_SIZE
            atype = -1 if sp.charge == 0 and sp.formula in _NONPOLAR else 0
        else:
            a0, atype = None, 0
        species[name] = Species(
            name=name, kind=kind, formula=sp.formula, charge=sp.charge, elements=dict(sp.elements),
            mol_wt=None, a0=a0, activity_type=atype, reaction=reaction,
            logk=ReactionLogK(terms), molar_volume=volume, ref_state={},
            provenance=prov(sp, kind.value))
        return True

    # aqueous species
    for (name, state), sp in sorted(data.species.items()):
        if state != "aq" or name in exclude or name in species or sp.charge is None:
            continue
        if sp.formula.startswith("(Z"):
            continue
        provider = _Aqueous(sp.eos)
        providers[name] = provider
        if not reaction_for(name, sp, provider, Kind.AQUEOUS):
            providers.pop(name, None)

    # gases and liquids
    for (name, state), sp in sorted(data.species.items()):
        if state == "gas" and name != "oxygen":
            gname = _gas_name(sp)
            if gname not in species and gname not in exclude:
                reaction_for(gname, sp, _Crystal([sp.eos]), Kind.GAS)
        elif state == "liq" and name != "water":
            if name not in species and name not in exclude:
                reaction_for(name, sp, _Crystal([sp.eos]), Kind.LIQUID)

    # minerals: all polymorphs of a name share one species (the most stable counts)
    mineral_names = sorted({n for (n, st) in data.species if re.fullmatch(r"cr\d?", st)})
    for name in mineral_names:
        if name in exclude or name in species:
            continue
        poly = data.polymorphs(name)
        if not poly:
            continue
        first = poly[0]
        reaction_for(name, first, _Crystal([p.eos for p in poly]), Kind.MINERAL,
                     volume=getattr(first.eos, "v0", None))

    return EosDatabase(title, "OBIGT", species, oxygen.eos)


def set_ion_sizes_from_eq36(db: EosDatabase, eq36_db: ThermoDatabase) -> int:
    """Copy ion-size parameters and nonpolar-neutral flags from an EQ3/6 database.

    Returns the number of species updated.  Only sizes are used; no logK data.
    """
    n = 0
    for sp36 in eq36_db.species.values():
        if sp36.kind in (Kind.MINERAL, Kind.GAS, Kind.LIQUID) or not sp36.a0:
            continue
        key = eq36_name_to_obigt(sp36.name)
        target = db.species.get(key)
        if target is not None and target.kind in (Kind.BASIS, Kind.AQUEOUS):
            target.a0 = sp36.a0
            target.activity_type = sp36.activity_type
            n += 1
    return n


def eq36_name_to_obigt(name: str) -> str:
    """Translate an EQ3/6 aqueous species name to the OBIGT spelling (Ca++ -> Ca+2)."""
    if name.endswith("(aq)"):
        return name[:-4]
    m = re.fullmatch(r"(.+?)(\+{1,4}|-{1,4})", name)
    if m:
        core, sign = m.groups()
        n = len(sign)
        return core + (sign[0] if n == 1 else f"{sign[0]}{n}")
    return name


PRESETS = {
    # newest curated data win: SLOP98 supplies extra species, SUPCRT92 the classic minerals,
    # the maintained inorganic_* files override both, and the Berman (1988 and later) dataset
    # replaces the Helgeson (1978) silicates, as in OBIGT's default data
    "current": SLOP98_FILES + SUPCRT92_FILES + DEFAULT_FILES + BERMAN_FILES,
    # SUPCRT92-era data win, to reproduce databases built from SUPCRT92 (EQ3/6 data0.com)
    "supcrt92": DEFAULT_FILES + SUPCRT92_FILES + SLOP98_FILES,
}


def load_eos_database(obigt_dir: str | Path, preset: str = "current",
                      eq36_database: ThermoDatabase | None = None,
                      masters: dict[str, str] | None = None) -> EosDatabase:
    """Read OBIGT csv files from ``obigt_dir`` and build an :class:`EosDatabase`.

    ``eq36_database`` optionally supplies ion-size parameters (angstrom) for the
    Debye-Hueckel and B-dot activity models; without it species get 4.0 (ions)
    or 3.0 (neutrals).
    """
    data = load_obigt(obigt_dir, PRESETS[preset])
    db = build_eos_database(data, masters=masters, title=f"OBIGT/HKF ({preset})")
    if "Berman" in data.missing:
        db.data_notes.append(
            "The Berman mineral data are missing from the data folder, so the older Helgeson "
            "(1978) silicates are used; run scripts/fetch_references.ps1 again to add them.")
    if eq36_database is not None:
        set_ion_sizes_from_eq36(db, eq36_database)
    return db

"""Rate laws for minerals and glasses that dissolve (or grow) at a measured rate instead of reaching
equilibrium at once.

A rate law is a sum of mechanisms, each with its own rate constant at 25 C, activation energy and
dependence on the activity of a catalysing species (H+ for the acid and base mechanisms, the fugacity
of CO2 for carbonates, Fe+3 and dissolved O2 for pyrite):

    r = A * sum_m k_m(T) * prod_i a_i^n_i * (1 - Omega)
    k_m(T) = k_m,25 * exp(-E_m / R * (1/T - 1/298.15))

where A is the reactive surface area (m2) and Omega = Q/K is the saturation state of the mineral, so the
rate is zero at equilibrium and changes sign past it.  r is in moles of the mineral per second: positive
dissolves, negative grows.  This is the form of Lasaga (1984) and Palandri and Kharaka (2004), whose
compilation (a public-domain USGS report) supplies the constants in ``data/``.

Basalt glass follows Gislason and Oelkers (2003): r = A * 10^-5.6 exp(-25.5 kJ / RT) (aH+^3/aAl+3)^(1/3)
moles of Si per m2 per second, far from equilibrium (glass is always undersaturated), so there is no
(1 - Omega) term.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .network import PROTON
from .speciation import Speciation

R_KJ = 8.314462618e-3                    # kJ / (mol K)
DATA = Path(__file__).with_name("data") / "rates_palandri_kharaka_2004.tsv"
PK_SOURCE = "Palandri and Kharaka (2004), USGS Open-File Report 2004-1068"


@dataclass(frozen=True)
class Mechanism:
    name: str
    log_k25: float                       # mol m-2 s-1 at 25 C
    e_kj: float                          # activation energy, kJ/mol
    terms: tuple[tuple[str, float], ...] = ()     # (species, exponent); CO2(g) means its fugacity

    def k(self, t_k: float) -> float:
        return 10.0 ** self.log_k25 * math.exp(-self.e_kj / R_KJ * (1.0 / t_k - 1.0 / 298.15))


@dataclass(frozen=True)
class RateLaw:
    mechanisms: tuple[Mechanism, ...]
    source: str = PK_SOURCE
    affinity: bool = True                # multiply by (1 - Omega)
    per_unit: float = 1.0                # reactant units per mole of what the constants count

    def far_from_equilibrium(self, spec: Speciation) -> float:
        """sum_m k_m(T) prod a^n, mol m-2 s-1 (no saturation term).  A mechanism whose catalyst is
        not in the fluid at all (no carbon for the CO2 term) contributes nothing."""
        t_k = spec.temperature_c + 273.15
        total = 0.0
        for m in self.mechanisms:
            log_term = 0.0
            for species, n in m.terms:
                la = log_activity(spec, species)
                if la is None:
                    break
                log_term += n * la
            else:
                total += m.k(t_k) * 10.0 ** max(min(log_term, 300.0), -300.0)
        return total

    def rate(self, spec: Speciation, area_m2: float, saturation: float | None = None) -> float:
        """Reactant units per second (positive dissolves).  ``saturation`` is log Omega."""
        if area_m2 <= 0.0:
            return 0.0
        r = area_m2 * self.far_from_equilibrium(spec) * self.per_unit
        if self.affinity:
            log_omega = 0.0 if saturation is None else saturation
            r *= 1.0 - 10.0 ** min(log_omega, 30.0)
        return r


def log_activity(spec: Speciation, name: str) -> float | None:
    if name == PROTON:
        return -spec.ph
    if name.endswith("(g)"):
        return spec.log_fugacity.get(name)
    st = spec.species.get(name)
    if st is None or st.molality <= 0.0:
        return None
    return st.log_activity


@lru_cache(maxsize=1)
def palandri_kharaka() -> dict[str, RateLaw]:
    """Rate laws by mineral name (as this database names them)."""
    mechs: dict[str, list[Mechanism]] = {}
    for line in DATA.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("mineral\t"):
            continue
        cols = line.split("\t")
        mineral, _pk, mech, log_k, e_kj = cols[:5]
        terms = tuple((t.split(":")[0], float(t.split(":")[1]))
                      for t in (cols[5] if len(cols) > 5 else "").split(";") if t)
        mechs.setdefault(mineral, []).append(Mechanism(mech, float(log_k), float(e_kj), terms))
    return {name: RateLaw(tuple(m)) for name, m in mechs.items()}


def rate_law(mineral: str) -> RateLaw | None:
    return palandri_kharaka().get(mineral)


def basalt_glass(si_moles_per_unit: float) -> RateLaw:
    """Gislason and Oelkers (2003); ``si_moles_per_unit`` converts the moles of Si the constant counts
    into the reactant's own unit (per gram of glass, for a rock reactant)."""
    # their 10^-5.6 is the factor in front of exp(-Ea/RT), not the constant at 25 C
    log_k25 = -5.6 - 25.5 / (R_KJ * 298.15 * math.log(10.0))
    return RateLaw((Mechanism("glass", log_k25, 25.5, ((PROTON, 1.0), ("Al+3", -1.0 / 3.0))),),
                   source="Gislason and Oelkers (2003), Geochim. Cosmochim. Acta 67, 3817-3832",
                   affinity=False, per_unit=1.0 / si_moles_per_unit)


# Calcite and aragonite in seawater, near equilibrium.  Laboratory constants such as Palandri and
# Kharaka's are for pure water far from equilibrium and linear in (1 - Omega); in seawater the rate falls
# much faster as Omega approaches 1, and the measured rates are normalized to the BET (gas-adsorption)
# surface, which for a porous foram or pteropod shell is tens of times the outer surface of a grain.
# r = A * k(T) * (1 - Omega)^n, A the BET surface (m2), k in mol m-2 s-1.  Each regime is a pair
# (log k at 5 C, n, Omega where it starts to apply); the fastest regime that applies sets the rate,
# which keeps the rate continuous where the fits meet.
NAVIAUX_SOURCE = ("Naviaux et al. (2019), Geochim. Cosmochim. Acta 246, 363-384 "
                  "(calcite in seawater, BET-normalized, 5 C fits)")
DONG_SOURCE = ("Dong et al. (2019), Earth Planet. Sci. Lett. 515, 1-12 "
               "(aragonite in seawater, 0.013 g/g/day (1 - Omega)^1.37 at 5 C on 2.46 m2/g)")


@dataclass(frozen=True)
class SeawaterCarbonateLaw:
    regimes: tuple[tuple[float, float, float], ...]      # (log k at 5 C in mol m-2 s-1, n, Omega below which it applies)
    e_kj: float
    source: str
    affinity: bool = True
    per_unit: float = 1.0

    def k_factor(self, t_c: float) -> float:
        return math.exp(-self.e_kj / R_KJ * (1.0 / (t_c + 273.15) - 1.0 / 278.15))

    def far_from_equilibrium(self, spec: Speciation) -> float:
        """The rate per m2 at Omega = 0, mol m-2 s-1."""
        return max(10.0 ** lk * self.k_factor(spec.temperature_c) for lk, _n, _o in self.regimes)

    def rate_per_m2(self, t_c: float, log_omega: float) -> float:
        """mol m-2 s-1, positive dissolves.  The fits measure dissolution only, so a supersaturated
        fluid neither dissolves nor grows the shell."""
        omega = 10.0 ** min(log_omega, 30.0)
        if omega >= 1.0:
            return 0.0
        f = self.k_factor(t_c)
        return max(10.0 ** lk * f * (1.0 - omega) ** n for lk, n, below in self.regimes if omega < below)

    def rate(self, spec: Speciation, area_m2: float, saturation: float | None = None) -> float:
        if area_m2 <= 0.0:
            return 0.0
        return area_m2 * self.rate_per_m2(spec.temperature_c, 0.0 if saturation is None else saturation) \
            * self.per_unit


# Naviaux et al. 2019, Table 2 at 5 C: log k = -13.07 mol cm-2 s-1, n = 0.34 for Omega > 0.75 (step
# retreat and defect-assisted pits) and -10.01, n = 4.81 below (homogeneous etch pits); activation energy
# 25 kJ/mol (their Section 4).  Dong et al. 2019: 0.013 g g-1 day-1 (1 - Omega)^1.37 at 5 C for aragonite
# of 2.46 m2/g, i.e. 6.1e-10 mol m-2 s-1; about twice as fast at 21 C, i.e. about 30 kJ/mol.
# By default aragonite takes the calcite fits on its own saturation state: Dong et al. find the same
# mechanisms and nearly the same critical saturation state (0.9 against calcite's 0.87) and that
# "calcite and aragonite dissolve at similar ratios relative to their abundance" in the water column, so
# what sets them apart is aragonite's higher solubility.  Their own 0.013 g/g/day law ("seawater_dong")
# is far slower per m2 than Naviaux's calcite fits away from equilibrium, which would make a pteropod
# outlast a foram in water more corrosive to it; it is kept as an option.
NAVIAUX_FITS = ((-9.07, 0.34, 1.0), (-6.01, 4.81, 0.75))
SEAWATER_CARBONATE = {
    "calcite": SeawaterCarbonateLaw(NAVIAUX_FITS, 25.0, NAVIAUX_SOURCE),
    "aragonite": SeawaterCarbonateLaw(NAVIAUX_FITS, 25.0,
                                      "Naviaux et al. (2019) calcite fits on aragonite's own saturation state "
                                      "(same mechanisms in seawater: Dong et al. 2019, Earth Planet. Sci. Lett. 515, 1-12)"),
    "aragonite_dong": SeawaterCarbonateLaw(((math.log10(0.013 / 2.46 / 100.087 / 86400.0), 1.37, 1.0),),
                                           30.0, DONG_SOURCE),
}


def seawater_rate_law(mineral: str, dong_aragonite: bool = False) -> SeawaterCarbonateLaw | None:
    if mineral == "aragonite" and dong_aragonite:
        return SEAWATER_CARBONATE["aragonite_dong"]
    return SEAWATER_CARBONATE.get(mineral)

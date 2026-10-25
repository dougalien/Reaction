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

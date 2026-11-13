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


# Olivine turning to serpentine at 250-350 C.  Laboratory dissolution constants (Palandri and Kharaka's
# forsterite) extrapolated to 300 C give rates orders of magnitude faster than olivine actually serpentinizes,
# as Malvoisin et al. (2012) also found with CHESS.  This law takes the measured serpentinization rate as
# the far-from-equilibrium rate and lets the database supply the driving force (1 - Omega).
# Derivation (scripts/malvoisin_rate.py): stage I rates of 38-50 um San Carlos olivine (Fo91) at
# 500 bar, 4.91, 5.96 and 6.86 % per 1000 h at 250, 270 and 300 C (Malvoisin et al. 2012, section 3.4, para. 19), on the
# surface of smooth 44 um spheres of 3.3 g/cm3, divided by 1 - Omega of forsterite in a fluid buffered by
# chrysotile + brucite (this database, a_H2O = 1, 500 bar), and fitted as ln k = a - E/RT.  The 350 C rate
# (0.15 % per 1000 h) was left out of the fit and is predicted within 11% (tests/test_kinetics.py), but 350 C
# is within 0.5 C of the equilibrium temperature, so that check is only as good as the log K.
# The constant at 25 C is only how the fit is stored: the law is measured from 250 to 350 C, at 500 bar, and
# its surface is geometric (about 400 cm2/g for 44 um grains), not BET.
MALVOISIN_SOURCE = ("serpentinization of olivine: rate of Malvoisin et al. (2012), J. Geophys. Res. 117, B04102 "
                    "(San Carlos olivine, 500 bar, 250-300 C, per m2 of geometric surface), times 1 - Omega of the "
                    "olivine from this database; measured from 250 to 350 C, extrapolated outside")
MALVOISIN_LOG_K25 = -11.174              # mol m-2 s-1, ln k = -11.855 - 4136.3 / T
MALVOISIN_E_KJ = 34.39
MALVOISIN_T_RANGE_C = (250.0, 350.0)


@dataclass(frozen=True)
class OlivineRateLaw(RateLaw):
    """A rate law for olivine (Fe,Mg)2SiO4 of fixed composition ``x_fo`` that dissolves congruently: its
    saturation is that of the ideal two-site solid solution of the solid-solution table,
    log Omega = x_fo (SI_fo - 2 log x_fo) + x_fa (SI_fa - 2 log x_fa), from the fluid's saturation
    indices of forsterite and fayalite."""
    x_fo: float = 0.9

    def log_omega(self, spec: Speciation) -> float | None:
        out = 0.0
        for name, x in (("forsterite", self.x_fo), ("fayalite", 1.0 - self.x_fo)):
            if x <= 0.0:
                continue
            si = spec.saturation.get(name)
            if si is None:
                return None
            out += x * (si - 2.0 * math.log10(x))
        return out

    def rate(self, spec: Speciation, area_m2: float, saturation: float | None = None) -> float:
        log_omega = self.log_omega(spec)
        if log_omega is None:
            raise ValueError("olivine rate law: forsterite or fayalite saturation is not available")
        return RateLaw.rate(self, spec, area_m2, log_omega)


def olivine_serpentinization(x_fo: float, mol_per_unit: float) -> OlivineRateLaw:
    """Serpentinization rate law for olivine with forsterite fraction ``x_fo``; ``mol_per_unit`` is moles of
    olivine (one Si each) per reactant unit (per gram, for a rock reactant)."""
    return OlivineRateLaw((Mechanism("serpentinization", MALVOISIN_LOG_K25, MALVOISIN_E_KJ),),
                          source=MALVOISIN_SOURCE, affinity=True, per_unit=1.0 / mol_per_unit, x_fo=x_fo)


# Olivine serpentinization as reaction progress, after Malvoisin et al. (2012), J. Geophys. Res. 117, B04102,
# Appendix A.  Their fit to San Carlos olivine (Fo91, 500 bar, 250-350 C, initial grain sizes 5-150 um) is
#   zeta = 1 - exp(-k f1(T) (t / R^2)^n)                                    (their eq. A2)
#   f1(T) = A exp(-b / T) (1 - exp(-c (1/T - 1/T0)))                        (their eq. A1)
# with t in hours, R the initial grain radius in um, A = 808.3, b = 3640 K, c = 8759 K, T0 = 623.6 K, and
# k, n from their Fig. A2: 0.008514 and 1.459 for 5 um < IGS < 50 um (run 22); for IGS >= 50 um two curves,
# 0.01333 and 1.289 (dissolution) and 8.561e-7 and 5.400 (dissolution enhanced by grain fracturing; run 2).
# Below 5 um the law does not hold (their section 4.1).  The sigmoidal shape is the creation of new surface by
# etch pits and fracturing (their para. 34); brucite armoring limits the measured progress to about 90 %
# (para. 29, 36), which they rescaled to 100 %.
# Adaptations made here (not in the paper): (1) A2 is used in its rate form at the current progress,
#   d zeta / dt = n K^(1/n) (1 - zeta) (-ln(1 - zeta))^((n - 1)/n),  K = k f(T) / R^(2n),
# exact at constant T, so the rate follows changing conditions; with n > 1 it is zero at zeta = 0, so it is
# evaluated at no less than ZETA_SEED.  (2) Their temperature-only equilibrium factor
# 1 - exp(-c (1/T - 1/T0)) is replaced by 1 - Omega of the olivine from this database, so the fluid and the
# pressure set where the reaction stops (T0 = 350.5 C is where this database puts it at 500 bar).  (3) For
# IGS >= 50 um the faster of the two curves at the current progress is used, which is how their Fig. A2d reads.
MALVOISIN_AVRAMI_SOURCE = ("serpentinization of olivine as reaction progress: Avrami law of Malvoisin et al. (2012), "
                           "J. Geophys. Res. 117, B04102, eqs. A1-A2 and Fig. A2 (San Carlos olivine, 500 bar, "
                           "250-350 C, grains 5-150 um), in rate form, with their equilibrium factor replaced by "
                           "1 - Omega of the olivine from this database; they observed at most about 90 % reaction "
                           "(brucite armoring)")
MALVOISIN_A1 = (808.3, 3640.0)          # A, b (K)
MALVOISIN_AVRAMI_SMALL = ((0.008514, 1.459),)                     # 5 um < IGS < 50 um
MALVOISIN_AVRAMI_LARGE = ((0.01333, 1.289), (8.561e-7, 5.400))    # IGS >= 50 um
ZETA_SEED = 1e-4


@dataclass(frozen=True)
class OlivineAvramiLaw(OlivineRateLaw):
    """Malvoisin et al. (2012) eq. A2 in rate form; reads the reaction progress, not the surface."""
    grain_radius_um: float = 22.0
    needs_progress: bool = True

    def curves(self) -> tuple[tuple[float, float], ...]:
        igs = 2.0 * self.grain_radius_um
        if igs < 5.0:
            raise ValueError(f"olivine grains of {igs:.1f} um: Malvoisin et al. (2012) found that their law does "
                             f"not hold below 5 um")
        return MALVOISIN_AVRAMI_SMALL if igs < 50.0 else MALVOISIN_AVRAMI_LARGE

    def rate_progress(self, spec: Speciation, fraction: float, amount: float) -> float:
        """Reactant units per second (positive reacts), at progress ``fraction`` of ``amount``."""
        if amount <= 0.0 or fraction >= 1.0:
            return 0.0
        log_omega = self.log_omega(spec)
        if log_omega is None:
            raise ValueError("olivine rate law: forsterite or fayalite saturation is not available")
        a, b = MALVOISIN_A1
        t_k = spec.temperature_c + 273.15
        f = a * math.exp(-b / t_k) * (1.0 - 10.0 ** min(log_omega, 30.0))
        if f <= 0.0:                     # at or past equilibrium: the law is for the forward reaction only
            return 0.0
        z = min(max(fraction, ZETA_SEED), 1.0 - 1e-12)
        best = 0.0
        for k, n in self.curves():
            big_k = k * f / self.grain_radius_um ** (2.0 * n)
            best = max(best, n * big_k ** (1.0 / n) * (1.0 - z) * (-math.log(1.0 - z)) ** ((n - 1.0) / n))
        return best / 3600.0 * amount


def olivine_avrami(x_fo: float, grain_radius_um: float) -> OlivineAvramiLaw:
    law = OlivineAvramiLaw((), source=MALVOISIN_AVRAMI_SOURCE, affinity=True, x_fo=x_fo,
                           grain_radius_um=grain_radius_um)
    law.curves()                         # refuse grains below 5 um at once
    return law


# Calcite and aragonite in seawater, near equilibrium.  Laboratory constants such as Palandri and
# Kharaka's are for pure water far from equilibrium and linear in (1 - Omega); in seawater the rate falls
# much faster as Omega approaches 1, and the measured rates are normalized to the BET (gas-adsorption)
# surface, which for a porous foram or pteropod shell is tens of times the outer surface of a grain.
# r = A * k(T) * (1 - Omega)^n, A the BET surface (m2), k in mol m-2 s-1.  Each regime is a pair
# (log k at 5 C, n, Omega below which it was fitted); the faster regime sets the rate, through a smooth
# maximum, so the rate and its slope are continuous where the fits meet.
NAVIAUX_SOURCE = ("Naviaux et al. (2019), Mar. Chem. 215, 103684 (calcite in seawater, BET-normalized: the in-situ "
                  "lower bound recommended there for water-column dissolution)")
NAVIAUX_LAB_SOURCE = ("Naviaux et al. (2019), Geochim. Cosmochim. Acta 246, 363-384 "
                      "(calcite in filtered, poisoned seawater in the laboratory, BET-normalized, 5 C fits)")
DONG_SOURCE = ("Dong et al. (2019), Earth Planet. Sci. Lett. 515, 1-12 "
               "(aragonite in seawater, 0.013 g/g/day (1 - Omega)^1.37 at 5 C on 2.46 m2/g)")


SMOOTHNESS = 10.0                       # s in the smooth maximum of the rate regimes


@dataclass(frozen=True)
class SeawaterCarbonateLaw:
    regimes: tuple[tuple[float, float, float], ...]      # (log k at 5 C in mol m-2 s-1, n, Omega below which it was fitted)
    e_kj: float
    source: str
    affinity: bool = True
    per_unit: float = 1.0
    # log Omega from (t_c, p_bar) to use instead of the database's saturation state (the marine
    # carbonate system: see marine.py); None = the database's
    omega: object = None
    mineral: str = ""

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
        # a smooth maximum of the regimes, (sum r_i^s)^(1/s): each fit holds where it is the faster, the
        # corner between them is rounded (a solver sees no kink), and at the crossover the rate is
        # 2^(1/s) times either fit.  The in-situ fits cross at Omega = 0.83, next to the 0.8 at which Naviaux
        # et al. (2019, Mar. Chem.) split their data; the laboratory fits at 0.79, next to 0.75.  (Cala et al.
        # 2024 round the same corner with an error function.)
        log_r = [lk + n * math.log10(1.0 - omega) for lk, n, _below in self.regimes]
        top = max(log_r)
        s = SMOOTHNESS
        return f * 10.0 ** (top + math.log10(sum(10.0 ** (s * (x - top)) for x in log_r)) / s)

    def rate(self, spec: Speciation, area_m2: float, saturation: float | None = None) -> float:
        if area_m2 <= 0.0:
            return 0.0
        if self.omega is not None:
            saturation = self.omega(spec.temperature_c, spec.pressure_bar or 1.01325)
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
NAVIAUX_LAB_FITS = ((-9.07, 0.34, 1.0), (-6.01, 4.81, 0.75))
# The default (Douglas, 2026-11-11): Naviaux et al. (2019), Mar. Chem. 215, 103684, abstract and Section 5: calcite
# dissolves about 4 times slower in situ (N. Pacific, 1.9-4.9 C) than in the laboratory, which they trace to dissolved
# organic carbon on the surface; for particles formed in and sinking out of surface water they recommend the lower
# bound R = 10^-14.3+-0.2 (1 - Omega)^0.11+-0.1 mol cm-2 s-1 for 0.8 < Omega < 1 and 10^-10.8+-0.4 (1 - Omega)^4.7+-0.7
# below 0.8.  Their in-situ Omega is from alkalinity and pH, which they find 5-10% below Omega from alkalinity and
# DIC (what Reaction computes): with Omega(Alk, DIC) these rates set in somewhat deeper than they would.  The
# temperature dependence is still the laboratory 25 kJ/mol, about the 5 C of their fits.
NAVIAUX_FITS = ((-10.3, 0.11, 1.0), (-6.8, 4.7, 0.8))
SEAWATER_CARBONATE = {
    "calcite": SeawaterCarbonateLaw(NAVIAUX_FITS, 25.0, NAVIAUX_SOURCE),
    "aragonite": SeawaterCarbonateLaw(NAVIAUX_FITS, 25.0,
                                      "Naviaux et al. (2019, Mar. Chem. 215) in-situ calcite fits on aragonite's own "
                                      "saturation state (same mechanisms in seawater: Dong et al. 2019, Earth Planet. "
                                      "Sci. Lett. 515, 1-12)"),
    "calcite_lab": SeawaterCarbonateLaw(NAVIAUX_LAB_FITS, 25.0, NAVIAUX_LAB_SOURCE),
    "aragonite_lab": SeawaterCarbonateLaw(NAVIAUX_LAB_FITS, 25.0,
                                          "Naviaux et al. (2019, GCA 246) laboratory calcite fits on aragonite's own "
                                          "saturation state"),
    "aragonite_dong": SeawaterCarbonateLaw(((math.log10(0.013 / 2.46 / 100.087 / 86400.0), 1.37, 1.0),),
                                           30.0, DONG_SOURCE),
}


def seawater_rate_law(mineral: str, dong_aragonite: bool = False,
                      lab: bool = False) -> SeawaterCarbonateLaw | None:
    if mineral == "aragonite" and dong_aragonite:
        return SEAWATER_CARBONATE["aragonite_dong"]
    if lab and mineral in ("calcite", "aragonite"):
        return SEAWATER_CARBONATE[f"{mineral}_lab"]
    return SEAWATER_CARBONATE.get(mineral)

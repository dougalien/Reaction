"""The aqueous model at one temperature and pressure.

Everything that depends only on (T, P, activity model) and not on a particular
set of constraints: the reaction network evaluated at T and P, the
Debye-Hueckel parameters, activity coefficients, water activity and the
molalities implied by a vector of basis-species log activities.  Shared by the
speciation solver (`speciation.py`) and the closed-system equilibrium solver
(`equilibrium.py`).
"""
from __future__ import annotations

import numpy as np

from . import activity, phreeqc_pitzer, pitzer
from .constants import LN10
from .database.model import ThermoDatabase
from .network import Network
from .thermo import bgamma as bgamma_table


MODELS = ("bgamma", "bdot", "davies", "pitzer", "pitzer_hkf", "ideal")
PITZER_MODELS = ("pitzer", "pitzer_hkf")


def _reference_pressure(db, temperature_c: float) -> float:
    """The pressure a database means when none is given (EQ3/6: 1.013 bar to 100 C, then saturation)."""
    p = db.reference_pressure(temperature_c) if hasattr(db, "reference_pressure") else None
    if p is None:
        from .thermo.eosdb import reference_pressure_bar
        p = reference_pressure_bar(temperature_c)
    return p


def _pitzer_names(db) -> dict:
    """Database name -> pitzer.dat name, computed once per database."""
    names = getattr(db, "_pitzer_dat_names", None)
    if names is None:
        names = phreeqc_pitzer.reaction_names(phreeqc_pitzer.read(), db)
        try:
            db._pitzer_dat_names = names
        except AttributeError:
            pass
    return names


class AqueousModel:
    def __init__(self, db: ThermoDatabase, net: Network, temperature_c: float,
                 pressure_bar: float | None, activity_model: str = "bdot",
                 ph_scale: str = "nbs"):
        if activity_model not in MODELS:
            raise ValueError("activity_model must be one of " + ", ".join(repr(m) for m in MODELS))
        if ph_scale not in ("nbs", "internal", "mesmer"):
            raise ValueError("ph_scale must be 'nbs', 'internal' or 'mesmer'")
        self.db, self.net = db, net
        self.temperature_c = temperature_c
        self.pressure_bar = pressure_bar
        self.activity_model, self.ph_scale = activity_model, ph_scale

        ev = net.evaluate(temperature_c, pressure_bar)
        self.ev = ev
        self.names = net.basis
        self.nbf = len(self.names)
        self.ib = {n: i for i, n in enumerate(self.names)}
        dep_index, logk_dep = ev.dep_index, ev.logk_dep
        if activity_model in PITZER_MODELS:
            # the Pitzer system's own complexes (NaCl, CaSO4, NaSO4- ...) are in its parameters
            cols = np.array([b in pitzer.COMPONENTS for b in net.basis])
            keep = [k for k, i in enumerate(dep_index)
                    if net.dep_names[i] in pitzer.KEEP or np.any((net.nu_dep[i] != 0) & ~cols)]
            dep_index, logk_dep = dep_index[keep], logk_dep[keep]
        self.dep_index = dep_index
        self.nu = net.nu_dep[dep_index]
        self.dep_names = [net.dep_names[i] for i in dep_index]
        self.logk_dep = np.array(logk_dep, dtype=float)
        self.z_dep = net.z_dep[dep_index]
        a0_dep = net.a0_dep[dep_index]
        self.type_dep = net.type_dep[dep_index]
        self.nu_ph = net.nu_phase[ev.phase_index]
        self.ph_names = [net.phase_names[i] for i in ev.phase_index]
        self.ph_kind = [net.phase_kind[i] for i in ev.phase_index]
        self.logk_ph = ev.logk_phase

        self.A, self.B, self.bdot = db.dh_params(temperature_c, pressure_bar)
        self.notes: list[str] = []
        self.bgamma = None
        if activity_model == "bgamma":
            p_ref = pressure_bar if pressure_bar is not None else _reference_pressure(db, temperature_c)
            self.bgamma = bgamma_table.bgamma(temperature_c, p_ref)
            why = bgamma_table.extrapolated(temperature_c, p_ref)
            if why:
                self.notes.append(f"b_gamma of the activity model is extrapolated here ({why})")
        elif activity_model == "bdot" and hasattr(db, "bdot_extrapolated") and db.bdot_extrapolated(temperature_c):
            self.notes.append("B-dot is tabulated to 300 C and is held at its 300 C value here "
                              "(activity coefficients are extrapolated)")
        if hasattr(db, "condition_warnings"):
            self.notes += db.condition_warnings(temperature_c, pressure_bar)
        self.logk_eh = db.logk_eh(temperature_c, pressure_bar)

        self.solute = net.is_solute
        self.z_b = net.z_basis
        self.a0_b = np.where(net.a0_basis > 0, net.a0_basis, activity.DEFAULT_ION_SIZE)
        self.a0_d = np.where(a0_dep > 0, a0_dep, activity.DEFAULT_ION_SIZE)
        self.cl_a0 = (db["Cl-"].a0 or 3.0) if "Cl-" in db else 3.0
        self.pitzer = None
        self.pitzer_data = None
        self._vex = None
        if activity_model in PITZER_MODELS:
            solutes = [n for n, s in zip(self.names, self.solute) if s]
            a_phi, p_atm = self.A * LN10 / 3.0, None
            if activity_model == "pitzer":
                a_phi = self._use_pitzer_dat(db, solutes, temperature_c, pressure_bar)
                p_atm = self.pitzer_data.water_state.p_atm
            self.pitzer = pitzer.Pitzer(solutes + self.dep_names,
                                        np.concatenate([self.z_b[self.solute], self.z_dep]),
                                        temperature_c, a_phi, p_atm)
            if temperature_c > pitzer.T_FITTED_C:
                self.notes.append(f"the Pitzer parameters are fitted to about {pitzer.T_FITTED_C:.0f} C "
                                  f"and are extrapolated here")

    @property
    def relaxed(self) -> bool:
        """Whether the solvers should under-relax the activity-coefficient update: Pitzer's coefficients
        depend on every molality, and a negative b_gamma near the critical point feeds back strongly."""
        return self.pitzer is not None or (self.activity_model == "bgamma" and (self.bgamma or 0.0) < 0.0)

    def _use_pitzer_dat(self, db, solutes, temperature_c, pressure_bar) -> float:
        """Take log K, molar volumes and A-phi from PHREEQC's pitzer.dat (``phreeqc_pitzer``) for every
        species and phase it defines, so the Pitzer system is one consistent model; the rest keep the
        database's log K.  Returns A-phi."""
        p_bar = pressure_bar if pressure_bar is not None else _reference_pressure(db, temperature_c)
        data = phreeqc_pitzer.PitzerData.at(temperature_c, p_bar)
        mapping = _pitzer_names(db)
        pot = {}
        for name in set(self.names) | set(self.dep_names) | set(self.ph_names):
            v = phreeqc_pitzer.potential(data, mapping, name)
            if v is not None:
                pot[name] = v

        def relog(row, name):
            if name not in pot:
                return None
            total = -pot[name]
            for j in np.nonzero(row)[0]:
                b = self.names[j]
                if b not in pot:
                    return None
                total += row[j] * pot[b]
            return total

        n_dep = n_ph = 0
        for k, name in enumerate(self.dep_names):
            v = relog(self.nu[k], name)
            if v is not None:
                self.logk_dep[k] = v
                n_dep += 1
        self.logk_ph = np.array(self.logk_ph, dtype=float)
        for k, name in enumerate(self.ph_names):
            v = relog(self.nu_ph[k], name)
            if v is not None:
                self.logk_ph[k] = v
                n_ph += 1
        self.pitzer_logk_counts = (n_dep, n_ph)
        # the ionic-strength part of each species' molar volume acts on its activity coefficient
        names = solutes + self.dep_names
        entries = [data.entries.get(mapping[n][0]) if n in mapping else None for n in names]
        if data.pfac and any(e is not None and e.vm_aq for e in entries):
            self._vex = (data, entries)
        self.pitzer_data = data
        if pressure_bar is not None and pressure_bar / 1.01325 > phreeqc_pitzer.P_MAX_ATM:
            self.notes.append(f"pitzer.dat is fitted to {phreeqc_pitzer.P_MAX_ATM:.0f} atm and is extrapolated "
                              f"here")
        if temperature_c > phreeqc_pitzer.T_CAP_C:
            self.notes.append("PHREEQC's water properties for pitzer.dat stop at 350 C and are held there")
        return data.water_state.a_phi

    def activities(self, m_dep, m_b, ionic: float, sum_m: float):
        """(log gamma of the basis, log gamma of the dependent species, log a(H2O)) for a solution."""
        if self.pitzer is not None:
            m = np.concatenate([m_b[self.solute], m_dep])
            ln_g, ln_aw = self.pitzer.compute(m)
            if self._vex is not None:
                data, entries = self._vex
                for i, e in enumerate(entries):
                    if e is not None and e.vm_aq:
                        ln_g[i] += LN10 * data.pfac * phreeqc_pitzer.vm_excess(e, data.water_state, ionic)
            lg_b = np.zeros(self.nbf)
            ns = int(np.sum(self.solute))
            lg_b[self.solute] = ln_g[:ns] / LN10
            return lg_b, ln_g[ns:] / LN10, ln_aw / LN10
        lg_b, lg_d = self.log_gammas(ionic, sum_m)
        return lg_b, lg_d, self.log_a_water(sum_m, ionic)

    # ---- activity coefficients and water activity ------------------------------------
    def log_gammas(self, ionic: float, sum_m: float | None = None):
        """log10 activity coefficients at ionic strength ``ionic``; ``sum_m`` (all solutes, mol/kg)
        enters only the mole-fraction term of the b_gamma model and defaults to ``ionic``."""
        nbf, t_c = self.nbf, self.temperature_c
        if sum_m is None:
            sum_m = ionic
        lg_b, lg_d = np.zeros(nbf), np.zeros(len(self.dep_names))
        if self.activity_model in ("ideal",) + PITZER_MODELS:  # Pitzer needs the molalities: activities()
            return lg_b, lg_d
        A, B, bdot, model = self.A, self.B, self.bdot, self.activity_model

        bg = self.bgamma

        def one(z, a0, typ):
            if z != 0:
                if model == "davies":
                    return activity.log_gamma_davies(z, ionic, A)
                if model == "bgamma":
                    return activity.log_gamma_bgamma(z, ionic, sum_m, A, B, bg)
                return activity.log_gamma_bdot(z, a0, ionic, A, B, bdot)
            if typ == -1 and model in ("bdot", "bgamma"):
                # Drummond (1981) fitted molal coefficients directly, so no mole-fraction term
                return activity.log_gamma_nonpolar(ionic, t_c, self.db.cco2)
            if model == "bgamma":
                return -activity.log_x_water(sum_m)        # HKF81 Gamma: rational to molal scale
            return 0.0

        for i in range(nbf):
            if self.solute[i]:
                lg_b[i] = one(self.z_b[i], self.a0_b[i], self.net.type_basis[i])
        for k in range(len(self.dep_names)):
            lg_d[k] = one(self.z_dep[k], self.a0_d[k], self.type_dep[k])
        if self.ph_scale == "nbs":
            if model == "davies":
                own = activity.log_gamma_davies(-1.0, ionic, A)
            elif model == "bgamma":
                own = activity.log_gamma_bgamma(-1.0, ionic, sum_m, A, B, bg)
            else:
                own = activity.log_gamma_bdot(-1.0, self.cl_a0, ionic, A, B, bdot)
            f = activity.log_gamma_chloride_bates_guggenheim(ionic, A) - own
            lg_b -= self.z_b * f
            lg_d -= self.z_dep * f
        return lg_b, lg_d

    def log_a_water(self, sum_m: float, ionic: float) -> float:
        if self.activity_model == "bdot":
            return activity.log_activity_water_bdot(sum_m, ionic, self.A, self.B, self.bdot)
        if self.activity_model == "bgamma":
            return activity.log_activity_water_bgamma(sum_m, ionic, self.A, self.B, self.bgamma)
        if self.activity_model == "davies":
            return activity.log_activity_water_davies(sum_m, ionic, self.A)
        return 0.0

    # ---- species molalities from basis log activities ------------------------------------
    def state(self, uvec, lg_b_ln, lg_d_ln):
        """(m_dep, m_basis) for ln activities ``uvec`` and ln activity coefficients."""
        m_dep = np.exp(np.clip(self.nu @ uvec - LN10 * self.logk_dep - lg_d_ln, -600.0, 600.0))
        m_b = np.where(self.solute, np.exp(np.clip(uvec - lg_b_ln, -600.0, 600.0)), 0.0)
        return m_dep, m_b

    def ionic_and_sum(self, m_dep, m_b):
        ionic = 0.5 * (np.sum(m_b * self.z_b ** 2) + np.sum(m_dep * self.z_dep ** 2))
        return float(ionic), float(np.sum(m_b) + np.sum(m_dep))

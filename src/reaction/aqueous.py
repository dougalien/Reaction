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

from . import activity
from .constants import LN10
from .database.model import ThermoDatabase
from .network import Network


class AqueousModel:
    def __init__(self, db: ThermoDatabase, net: Network, temperature_c: float,
                 pressure_bar: float | None, activity_model: str = "bdot",
                 ph_scale: str = "nbs"):
        if activity_model not in ("bdot", "davies", "ideal"):
            raise ValueError("activity_model must be 'bdot', 'davies' or 'ideal'")
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
        self.nu = net.nu_dep[ev.dep_index]
        self.dep_names = [net.dep_names[i] for i in ev.dep_index]
        self.logk_dep = ev.logk_dep
        self.z_dep = net.z_dep[ev.dep_index]
        a0_dep = net.a0_dep[ev.dep_index]
        self.type_dep = net.type_dep[ev.dep_index]
        self.nu_ph = net.nu_phase[ev.phase_index]
        self.ph_names = [net.phase_names[i] for i in ev.phase_index]
        self.ph_kind = [net.phase_kind[i] for i in ev.phase_index]
        self.logk_ph = ev.logk_phase

        self.A, self.B, self.bdot = db.dh_params(temperature_c, pressure_bar)
        self.notes: list[str] = []
        if hasattr(db, "bdot_extrapolated") and db.bdot_extrapolated(temperature_c):
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

    # ---- activity coefficients and water activity ------------------------------------
    def log_gammas(self, ionic: float):
        nbf, t_c = self.nbf, self.temperature_c
        lg_b, lg_d = np.zeros(nbf), np.zeros(len(self.dep_names))
        if self.activity_model == "ideal":
            return lg_b, lg_d
        A, B, bdot, model = self.A, self.B, self.bdot, self.activity_model

        def one(z, a0, typ):
            if z != 0:
                if model == "davies":
                    return activity.log_gamma_davies(z, ionic, A)
                return activity.log_gamma_bdot(z, a0, ionic, A, B, bdot)
            if typ == -1 and model == "bdot":
                return activity.log_gamma_nonpolar(ionic, t_c, self.db.cco2)
            return 0.0

        for i in range(nbf):
            if self.solute[i]:
                lg_b[i] = one(self.z_b[i], self.a0_b[i], self.net.type_basis[i])
        for k in range(len(self.dep_names)):
            lg_d[k] = one(self.z_dep[k], self.a0_d[k], self.type_dep[k])
        if self.ph_scale == "nbs":
            if model == "davies":
                own = activity.log_gamma_davies(-1.0, ionic, A)
            else:
                own = activity.log_gamma_bdot(-1.0, self.cl_a0, ionic, A, B, bdot)
            f = activity.log_gamma_chloride_bates_guggenheim(ionic, A) - own
            lg_b -= self.z_b * f
            lg_d -= self.z_dep * f
        return lg_b, lg_d

    def log_a_water(self, sum_m: float, ionic: float) -> float:
        if self.activity_model == "bdot":
            return activity.log_activity_water_bdot(sum_m, ionic, self.A, self.B, self.bdot)
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

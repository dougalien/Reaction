"""Reaction on Streamlit: an Input window (the scenario card) and a Results window, as in the desktop app.

The desktop app's ``Session`` (``reaction.app.session``) is used unchanged: it holds the card, runs the
engine in a background thread, saves runs and writes their story.  Each user has a private workspace id
in the page address (``?w=...``); the Input window and the Results window opened from it share that id,
so they show the same card and runs, live, in two browser windows.  One ``Engine`` (the thermodynamic
data, loaded once) serves everyone, one run at a time.  Nothing here computes thermodynamics.
"""
from __future__ import annotations

import base64
import csv
import io
import json
import re
import secrets
import shutil
import tempfile
import threading
import time
from pathlib import Path

import plotly.graph_objects as go
import streamlit as st

from .. import __version__
from ..app.engine import Engine
from ..app.session import Session
from ..app.scenario import reactant_text

FEEDBACK_URL = "https://github.com/dougalien/Reaction/issues"
# Where the windows live; a host app (Protolith) that mounts these pages elsewhere changes them, and can
# add panels to the Input window through the hook (a function of the Session).
PATHS = {"input": "./", "results": "./results", "sources": "./sources"}
INPUT_EXTRAS: list = []
WORKSPACE_TTL = 12 * 3600
STATIC = Path(__file__).resolve().parents[1] / "app" / "static"
ICON = str(STATIC / "reaction_512.png")
STRIPE = ("linear-gradient(90deg, #D6B27C 0 22%, #3E5C6B 22% 34%, #D6B27C 34% 62%, #B4532A 62% 70%, "
          "#D6B27C 70% 88%, #3E5C6B 88% 100%)")
PALETTE = ["#3E5C6B", "#B4532A", "#C49A55", "#5E8C6A", "#8A6A9A", "#4F86A0", "#A8644A", "#6B6F45"]
PLOT_BG = "#F3F1EB"
GRID = "#DDD8CC"

# ---- display text: the engine's notes and stories are plain ASCII; show them with proper symbols -------
_SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
_SUP = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
_FORMULA = re.compile(r"(?<![\w.])(fO2|H2O|H2S|H2|O2|CO2|CH4|SO4|SO2|N2|NO3|NH4|HCO3|CO3|SiO2|CaCO3|FeS2|"
                      r"Fe\(OH\)3|Fe2O3|Al2O3|Na2O|K2O)(?![\w])")
_CHARGE = re.compile(r"(?<![\w])((?:[A-Z][a-z]?[₀-₉]*)+)([+-])(\d)(?![\w.])")
_SPELL = {"vapour": "vapor", "Vapour": "Vapor", "aluminium": "aluminum", "Aluminium": "Aluminum",
          "behaviour": "behavior", "colour": "color", "centre": "center", "labelled": "labeled"}
_SPELL_RE = re.compile(r"\b(" + "|".join(_SPELL) + r")")


def tidy(text) -> str:
    """Degree signs, subscripts, charges, δ, minus signs and U.S. spelling in text shown to people."""
    if not isinstance(text, str):
        return text
    t = re.sub(r"(?<=\d) C\b", " °C", text)
    t = _SPELL_RE.sub(lambda m: _SPELL[m.group(1)], t)
    t = re.sub(r"(?<![\w])d(18O|13C|34S|D|7Li|11B|30Si|56Fe|26Mg|44Ca)(?![\w])", r"δ\1", t)
    t = _FORMULA.sub(lambda m: m.group(1).translate(_SUB), t)
    t = _CHARGE.sub(lambda m: m.group(1) + (m.group(3) + m.group(2)).translate(_SUP), t)
    t = re.sub(r"(?:(?<=[\s(=])|^)-(?=\d)", "−", t)
    return t


def sentence(text: str) -> str:
    """Tidy text that should start with a capital letter (mineral names in notes start lowercase)."""
    t = tidy(text)
    first = re.match(r"[a-z]+", t)
    if first and first.group(0) not in ("log", "pH", "pe"):
        t = t[0].upper() + t[1:]
    return t


# ---- workspaces: one Session per user, shared by that user's two windows -----------------------------
@st.cache_resource
def _engine() -> Engine:
    eng = Engine()
    threading.Thread(target=lambda: _try(lambda: eng.db), daemon=True).start()     # load the data early
    return eng


def _try(fn):
    try:
        fn()
    except Exception:
        pass


@st.cache_resource
def _registry() -> dict:
    return {"lock": threading.Lock(), "ws": {}, "root": Path(tempfile.mkdtemp(prefix="reaction-ws-"))}


def workspace() -> tuple[str, Session]:
    reg = _registry()
    wid = st.query_params.get("w")
    if not wid or not re.fullmatch(r"[A-Za-z0-9_-]{8,64}", wid):
        wid = secrets.token_urlsafe(16)
        st.query_params["w"] = wid
    now = time.time()
    with reg["lock"]:
        for k in [k for k, (s, t) in reg["ws"].items() if now - t > WORKSPACE_TTL and k != wid]:
            if reg["ws"][k][0].run_state["state"] != "running":
                del reg["ws"][k]
                shutil.rmtree(reg["root"] / k, ignore_errors=True)
        if wid in reg["ws"]:
            sess = reg["ws"][wid][0]
        else:
            sess = Session(reg["root"] / wid, _engine(), None)
        reg["ws"][wid] = (sess, now)
    return wid, sess


def follow(sess: Session, key: str) -> None:
    """Re-run the page when the session changes (run progress, a new run, an edit in the other window)."""
    seen_key = f"seen-{key}"
    st.session_state.setdefault(seen_key, sess.version)
    running = sess.run_state["state"] == "running"

    @st.fragment(run_every=1.0 if running else 2.0)
    def _watch():
        if sess.version != st.session_state[seen_key]:
            st.session_state[seen_key] = sess.version
            st.rerun()
    _watch()
    st.session_state[seen_key] = sess.version


def _svg(name: str) -> str:
    return "data:image/svg+xml;base64," + base64.b64encode((STATIC / name).read_bytes()).decode()


def header(subtitle: str) -> None:
    st.markdown(
        f"""<div style="display:flex;align-items:center;gap:12px;padding:4px 0 10px">
<img src="{_svg('dougalien_mark.svg')}" width="34" height="34" alt="Dougalien Systems">
<img src="{_svg('reaction_mark.svg')}" width="30" height="30" alt="">
<span style="font-size:1.55rem;font-weight:600;color:#2B3134">Reaction</span>
<a href="https://dougalien.com" target="_blank" style="color:#5b6468;text-decoration:none">by Dougalien Systems</a>
<span style="color:#6f777a;font-size:0.9rem">· {subtitle}</span></div>
<div style="height:4px;border-radius:2px;background:{STRIPE};margin:0 0 14px"></div>""",
        unsafe_allow_html=True)


def footer(wid: str) -> None:
    st.divider()
    st.caption(f"Reaction {__version__} · beta · [Sources]({PATHS['sources']}?w={wid}) · [Report a problem]({FEEDBACK_URL}) · "
               "Every number comes from Reaction's own solver; nothing you enter leaves this app.")


# ---- the Input window: the scenario card ----------------------------------------------------------
NUM_FIELDS = {  # field: (label, min, max, step, help)
    "t_start_c": ("Temperature start (°C)", 0.0, 600.0, 1.0, None),
    "t_end_c": ("Temperature end (°C)", 0.0, 600.0, 1.0, None),
    "p_start_bar": ("Pressure start (bar)", 1.0, 5000.0, 10.0, None),
    "p_end_bar": ("Pressure end (bar)", 1.0, 5000.0, 10.0, None),
    "steps": ("Steps", 2, 200, 1, None),
    "cells": ("Cells", 2, 100, 1, "Number of cells the column is divided into (2 to 100)"),
    "pore_volumes": ("Pore volumes", 0.0, 100.0, 0.5, "How many column volumes of fluid flow through"),
    "column_m": ("Column length (m)", 0.0, None, 0.1, None),
    "dispersivity_m": ("Dispersivity (m)", 0.0, None, 0.01, "Spreads the fronts: D = dispersivity × velocity. 0 gives sharp fronts"),
    "evaporate": ("Evaporate to (times concentrated)", 1.0, 10000.0, 1.0,
                  "Water is taken away until the fluid is this many times more concentrated; 1 means no evaporation"),
    "mix_ratio": ("Up to (kg per kg)", 0.0, 1000.0, 0.1, None),
    "mix_t_c": ("Its temperature (°C)", 0.0, 600.0, 1.0, "Empty: the temperature follows the path"),
    "shell_area_exponent": ("Shell surface exponent", -1.0, 1.0, 0.1, "Forams and pteropods with seawater rates: "
                            "whole surface = start surface x (mass left / mass)^p. 0 keeps it, so the surface "
                            "per gram rises as the shell dissolves; 0.667 is a shrinking solid grain"),
    "surface_area_cm2_per_g": ("Surface area (cm²/g)", 0.0, None, 10.0, "Reactive surface per gram of each reactant; 100 cm²/g is roughly 0.1 mm grains. A reactant can carry its own: 20 ug forams at 2 m2/g (with seawater carbonate rates, forams and pteropods default to 2 m²/g, BET)"),
    "duration": ("Run for", 0.0, None, 1.0, "In a column, the pore volumes pass in this time"),
    "exchange_meq": ("Exchange sites (meq per kg water)", 0.0, None, 1.0,
                     "Cation exchange sites per kg of water (Gaines-Thomas, PHREEQC constants). 0 = none"),
    "hfo_g": ("Ferric oxide surface (g per kg water)", 0.0, None, 0.1,
              "Hydrous ferric oxide with Dzombak and Morel surface complexation. 0 = none"),
    "kaolinite_g": ("Kaolinite", 0.0, None, 1.0, None),
    "illite_g": ("Illite", 0.0, None, 1.0, None),
    "montmorillonite_g": ("Montmorillonite", 0.0, None, 1.0, None),
    "kd_solid_g": ("Kd solid (g per kg water)", 0.0, None, 1.0,
                   "Grams of sediment or soil per kg of water that hold metals by their partition coefficient Kd "
                   "(Allison and Allison 2005, U.S. EPA). 0 = none"),
    "uncertainty_runs": ("Uncertainty range (extra runs, 0 = off)", 0, 200, 1,
                         "Extra runs with every log K shifted at random by its uncertainty; the plots then shade the "
                         "middle 90% of them. Each takes as long as the main run"),
    "ph": ("pH", 0.0, 14.0, 0.1, None),
    "fluid_t_c": ("Measured at (°C)", 0.0, 600.0, 1.0, None),
    "fluid_p_bar": ("Measured at (bar)", 1.0, 5000.0, 1.0, "Empty means 1 bar (or the boiling curve above 100 °C)"),
    "redox_value": ("Value", None, None, 0.1, None),
}
CHOICES = {  # field: (label, [(value, text)], help)
    "activity_model": ("Activity model", [("bdot", "B-dot"), ("davies", "Davies"), ("pitzer", "Pitzer (brines)")], None),
    "units": ("Default unit", [("mg/kg", "mg/kg"), ("mg/l", "mg/L"), ("mmol/kg", "mmol/kg"), ("mol/kg", "mol/kg"),
                               ("mmolal", "mmol/kg water"), ("molal", "mol/kg water"), ("meq/l", "meq/L")], None),
    "mix_spacing": ("Steps even in", [("ratio", "Amount mixed"), ("temperature", "Temperature")], None),
    "kinetics": ("How reactants react", [(False, "All of it, steadily"), (True, "At measured rates")],
                 "Minerals and basalt glass dissolve or grow at measured rates (Palandri and Kharaka 2004; "
                 "Gislason and Oelkers 2003) instead of all at once"),
    "carbonate_rates": ("Calcite and aragonite rates", [("lab", "Laboratory (far from equilibrium)"),
                                                        ("seawater", "Seawater, near saturation"),
                                                        ("seawater_dong", "Seawater, aragonite after Dong et al.")],
                        "Laboratory: Palandri and Kharaka (2004). Seawater: measured near saturation on the BET "
                        "surface, Naviaux et al. (2019) for calcite and for aragonite on its own saturation state; or aragonite "
                        "after Dong et al. (2019)"),
    "depth_axis": ("Plot against", [(False, "Time"), (True, "Depth (sinking shell)")],
                   "Depth reads the pressure as depth in seawater, about 10 m per bar"),
    "duration_unit": ("Time unit", [(u, u) for u in ("seconds", "minutes", "hours", "days", "years")], None),
    "clay_edges": ("Clay edge sites", [(True, "Bind metals (2SPNE)"), (False, "Exchange sites only")], None),
    "kd_medium": ("Kd values for", [("sediment", "Bed sediment"), ("suspended", "Suspended matter"), ("soil", "Soil")], None),
    "solid_solutions": ("Solid solutions", [(False, "Pure minerals only"), (True, "Can form")], None),
    "follow_phase": ("Follow", [("whole", "Whole fluid"), ("brine", "Brine"), ("vapour", "Vapor")],
                     "Where a saline fluid splits into brine and vapor: which one the run follows"),
    "isotopes": ("Isotopes", [(False, "Off"), (True, "On")], None),
    "gas_phase": ("Gas phase", [(True, "Can form"), (False, "Suppressed")], None),
    "methane_sulfate": ("Methane and sulfate", [("apart", "Held apart"), ("react", "Can react")], None),
    "redox_apart": ("Redox couples", [(False, "Reach equilibrium"), (True, "Held apart")], None),
    "minerals_form": ("Minerals", [(True, "Can form"), (False, "Suppressed")], None),
    "energy": ("Catabolic energy", [(False, "Off"), (True, "On")], None),
}
REDOX_LABEL = {"O2": "Dissolved O₂", "H2": "Dissolved H₂", "Fe+3": "Fe(III) / Fe(II)", "H2S": "Sulfide / sulfate",
               "NO3-": "Nitrate / ammonia"}


def _edit(sess: Session, field: str) -> None:
    value = st.session_state[f"f-{field}"]
    if field in ("mix_t_c", "fluid_p_bar", "ph", "redox_value") and value is None:
        value = ""
    st.session_state["problems"] = sess.edit_card({field: value})


def _sync(key: str, value) -> None:
    st.session_state[key] = value


def num(sess: Session, card: dict, field: str, col=st, disabled=False) -> None:
    label, lo, hi, step, help_ = NUM_FIELDS[field]
    v = card.get(field)
    integer = isinstance(step, int)
    if v is not None:
        v = int(v) if integer else float(v)
    _sync(f"f-{field}", v)
    col.number_input(label, key=f"f-{field}", min_value=lo, max_value=hi, step=step, help=help_,
                     disabled=disabled, on_change=_edit, args=(sess, field),
                     format=None if integer else "%g")


def choice(sess: Session, card: dict, field: str, col=st, options=None, label=None) -> None:
    lab, opts, help_ = CHOICES.get(field, (label, options, None))
    opts = options or opts
    values = [o[0] for o in opts]
    text = dict(opts)
    v = card.get(field)
    if v not in values:
        values.append(v)
        text[v] = str(v)
    _sync(f"f-{field}", v)
    col.selectbox(label or lab, values, key=f"f-{field}", format_func=lambda x: text.get(x, str(x)), help=help_,
                  on_change=_edit, args=(sess, field))


def _reset(sess: Session, key: str) -> None:
    sess.reset_card(key)


def _example(sess: Session, key: str) -> None:
    sess.load_example(key)


def _rock(sess: Session, label: str) -> None:
    g = st.session_state.get("rock-g") or 0
    if g > 0:
        st.session_state["problems"] = sess.edit_card({"reactants": f"{g:g} g {label}"})


def _add_reactant(sess: Session) -> None:
    text = (st.session_state.get("add-reactant") or "").strip()
    if not text:
        return
    now = [reactant_text(r) for r in sess.card.to_dict()["reactants"]]
    problems = sess.edit_card({"reactants": "; ".join(now + [text])})
    st.session_state["problems"] = problems
    if not problems:
        st.session_state["add-reactant"] = ""


def _drop_reactant(sess: Session, i: int) -> None:
    rs = [r for k, r in enumerate(sess.card.to_dict()["reactants"]) if k != i]
    sess.edit_card({"reactants": rs})


def _add_blocked(sess: Session) -> None:
    v = (st.session_state.get("add-blocked") or "").strip().lower()
    blocked = sess.card.to_dict()["blocked_minerals"]
    if v and v not in blocked:
        sess.edit_card({"blocked_minerals": blocked + [v]})
    st.session_state["add-blocked"] = ""


def _set(sess: Session, field: str, value) -> None:
    sess.edit_card({field: value})


def _reactant_text(r: dict) -> str:
    return reactant_text(r)


def input_page() -> None:
    st.set_page_config(page_title="Reaction | Input", page_icon=ICON, layout="wide")
    wid, sess = workspace()
    snap = sess.snapshot()
    card = snap["card"]
    col = card.get("kind") == "transport"
    running = snap["run"]["state"] == "running"

    top = st.columns([4, 1.4])
    with top[0]:
        header("water-rock reaction, speciation and reaction paths")
    with top[1]:
        st.link_button("Open results window", f"{PATHS['results']}?w={wid}", width="stretch",
                       help="Opens the Results window in a new browser window; put it on a second screen. "
                            "It follows this window's runs.")

    with st.container(border=True):
        st.markdown("**Scenario card**")
        st.caption("Start from")
        ex_labels = {e["label"] for e in snap["examples"]}
        cards = [c for c in snap["starting_cards"] if c["label"] not in ex_labels]
        for c, s in zip(st.columns(len(cards) + 1), cards + [None]):
            if s is None:
                c.button("Reset", on_click=_reset, args=(sess, "seawater"), help="Back to the default card",
                         width="stretch")
            else:
                c.button(s["label"], key=f"start-{s['key']}", help=tidy(s["description"]), on_click=_reset,
                         args=(sess, s["key"]), width="stretch", disabled=running)
        st.caption("Or load a worked example")
        exs = snap["examples"]
        per_row = 5
        for i in range(0, len(exs), per_row):
            row = exs[i:i + per_row]
            for c, e in zip(st.columns(per_row), row):
                c.button(e["label"], key=f"ex-{e['key']}", help=tidy(f"{e['feature']}: {e['description']}"),
                         on_click=_example, args=(sess, e["key"]), width="stretch", disabled=running,
                         type="primary" if snap["example"] == e["key"] else "secondary")
        ex = next((e for e in exs if e["key"] == snap["example"]), None)
        if ex:
            with st.container(border=True):
                st.markdown(f"**Example: {ex['label']}** · shows {ex['feature'].lower()}")
                st.markdown(tidy(ex["description"]))
                st.caption("How the card sets it up:")
                for _, text in ex["setup"]:
                    st.markdown(f"- {tidy(text)}")
                st.caption("Press Run to see it. Change anything on the card first to see what difference it makes.")
        _sync("f-title", card.get("title") or "")
        st.text_input("Title", key="f-title", placeholder="optional", on_change=_edit, args=(sess, "title"))

    left, right = st.columns(2, gap="large")
    with left:
        _fluid_box(sess, snap, card, col)
        _path_box(sess, card, col)
    with right:
        _reactants_box(sess, snap, card, col, running)
        _surfaces_box(sess, card)
    _options_box(sess, card)

    for p in snap["card_problems"] + list(st.session_state.pop("problems", []) or []):
        st.error(tidy(p))
    for a in snap["advisories"] or []:
        st.warning(tidy(a))
    run_bar(sess, input_page=True)
    for extra in INPUT_EXTRAS:
        extra(sess)
    _history(sess, snap)
    footer(wid)


def _fluid_box(sess, snap, card, col) -> None:
    read = snap["fluid_read"]
    fluids = snap["fluids"]
    with st.container(border=True):
        st.markdown("**Fluid**")
        keys = [f["key"] for f in fluids]
        labels = {f["key"]: f["label"] for f in fluids}
        a, b = st.columns(2)
        choice(sess, card, "fluid", a, options=[(k, labels[k]) for k in keys],
               label="Fluid flowing in" if col else "Fluid")
        if col:
            choice(sess, card, "pore_fluid", b, options=[("", "The same fluid")] +
                   [(k, labels[k]) for k in keys if k != "custom"], label="Column first filled with")
        f = next((x for x in fluids if x["key"] == card["fluid"]), None)
        if f and card["fluid"] != "custom":
            st.caption(tidy(f["description"]) + " Editing it below starts a custom fluid from it.")
        _sync("f-analysis", read.get("analysis") or "")
        st.text_area("Analysis (one solute per line: name, amount, unit if not the default)", key="f-analysis",
                     height=150, placeholder="Na 10768\nCa2+ 10.3 mmol/kg\nSiO2 12", on_change=_edit,
                     args=(sess, "analysis"))
        rd = {**card, "units": read["units"], "ph": read["ph"], "fluid_t_c": read["measured_t_c"],
              "fluid_p_bar": read["measured_p_bar"], "charge_balance": read["charge_balance"],
              "redox": read["redox"], "redox_value": read["redox_value"]}
        a, b = st.columns(2)
        choice(sess, rd, "units", a)
        num(sess, rd, "ph", b)
        num(sess, rd, "fluid_t_c", a)
        num(sess, rd, "fluid_p_bar", b)
        ions = [ln["species"] for ln in read.get("lines", []) if ln.get("charge")]
        choice(sess, rd, "charge_balance", a, label="Charge balance",
               options=[("auto", "Auto"), ("none", "None (use as given)")] + [(i, f"On {i}") for i in ions])
        choice(sess, rd, "redox", b, label="Redox state",
               options=[("auto", "Auto (from the analysis)"), ("none", "Not tracked"), ("air", "Air (O₂ 0.21 bar)"),
                        ("log_fo2", "log fO₂"), ("eh", "Eh (V)"), ("pe", "pe")]
               + [(s, REDOX_LABEL.get(s, s)) for s in read.get("redox_species", [])])
        if read["redox"] in ("log_fo2", "eh", "pe"):
            num(sess, rd, "redox_value", a)
        st.caption("Per kg means per kg of solution; mg/L is read at 1 kg/L. Auto balance adjusts the main ion of the "
                   "side that is short. Redox: Fe(III), sulfide, H₂ or O₂ in the analysis set it; otherwise choose air, "
                   "a log fO₂, Eh or pe.")
        n = len(read.get("lines", []))
        tds = read.get("tds_mg_per_kg") or 0
        tds_text = f"{tds:,.0f}" if tds >= 100 else f"{tds:.3g}"
        summary = (f"As read: {n} solute{'s' if n != 1 else ''}, pH {read['ph'] if read['ph'] is not None else 'not given'} "
                   f"at {read['measured_t_c']:g} °C, {tds_text} mg/kg dissolved solids"
                   if n else "Nothing read yet")
        with st.expander(summary):
            st.dataframe([{"As written": ln["written"], "Value": f"{ln['value']} {ln['unit']}",
                           "mol/kg water": "" if ln["molal"] is None else f"{ln['molal']:.4g}"}
                          for ln in read.get("lines", [])], hide_index=True, width="stretch")
            for t in read.get("notes", []):
                st.caption(tidy(t))
        a, b = st.columns(2)
        choice(sess, card, "activity_model", a)
        num(sess, card, "evaporate", b)


def _path_box(sess, card, col) -> None:
    with st.container(border=True):
        st.markdown("**Path**")
        kind = card.get("kind") or "path"
        _sync("f-kind", kind)
        st.segmented_control("Kind of run", ["path", "transport"], key="f-kind", on_change=_edit, args=(sess, "kind"),
                             format_func=lambda k: {"path": "One fluid along a path", "transport": "Flow through a column"}[k])
        if not col:
            num(sess, card, "steps")
        st.caption("Column: from the inlet to the outlet" if col else "Path: a straight line from start to end")
        a, b = st.columns(2)
        mixing_t = bool(card.get("mix_with")) and card.get("mix_t_c") is not None
        num(sess, card, "t_start_c", a)
        num(sess, card, "t_end_c", b, disabled=mixing_t)
        num(sess, card, "p_start_bar", a)
        num(sess, card, "p_end_bar", b)
        st.caption("Equal start and end hold that variable fixed. Equilibrium constants are computed at every point, so "
                   "any temperature and pressure works (up to 600 °C and 5000 bar).")
        if col:
            a, b = st.columns(2)
            num(sess, card, "cells", a)
            num(sess, card, "pore_volumes", b)
            num(sess, card, "column_m", a)
            num(sess, card, "dispersivity_m", b)
            st.caption("The fluid above flows in at the inlet; minerals and sorbed matter stay in their cell, and every "
                       "cell is at equilibrium at its own temperature and pressure (PHREEQC's TRANSPORT scheme), except "
                       "for minerals reacting at measured rates.")
        else:
            _sync("f-mode", card.get("mode"))
            st.segmented_control("Mode", ["static", "flow_through", "both"], key="f-mode", on_change=_edit,
                                 args=(sess, "mode"),
                                 format_func=lambda m: {"static": "Static", "flow_through": "Flow-through", "both": "Both"}[m])
            st.caption("Static: minerals stay and can redissolve. Flow-through: each step's minerals are left behind.")


def _reactants_box(sess, snap, card, col, running) -> None:
    with st.container(border=True):
        st.markdown("**Reactants and mixing**")
        if not col:
            labels = {f["key"]: f["label"] for f in snap["fluids"] if f["key"] != "custom"}
            a, b = st.columns(2)
            choice(sess, card, "mix_with", a, options=[("", "None")] + list(labels.items()), label="Mix in a second fluid")
            num(sess, card, "mix_ratio", b, disabled=not card.get("mix_with"))
            num(sess, card, "mix_t_c", a, disabled=not card.get("mix_with"))
            choice(sess, card, "mix_spacing", b)
            if card.get("mix_with"):
                st.caption("Mixed in steadily; the temperature follows the path above." if card.get("mix_t_c") is None
                           else "Mixed in steadily; the temperature follows the mixing (mass-weighted, equal heat "
                                "capacities), so the end temperature above is not used.")
        kin = card.get("kinetics")
        hint = ("(in every cell, per kg of pore water)" if col else
                "(what is there at the start, per kg of water; minerals without a measured rate are added steadily)"
                if kin else "(added steadily, all of it by the last step, per kg of water)")
        st.markdown(f"Reactants <span style='color:#6f777a;font-size:0.85rem'>{hint}</span>", unsafe_allow_html=True)
        rocks = snap["rock_cards"]
        st.session_state.setdefault("rock-g", 100.0)
        st.number_input("Rock, g per kg of water (for the rock buttons below)", key="rock-g", min_value=0.0, step=10.0, format="%g",
                        help="Each rock button puts that rock on the card as the reactant, in place of what is there")
        for i in range(0, len(rocks), 3):
            for c, r in zip(st.columns(3), rocks[i:i + 3]):
                c.button(r["label"], key=f"rock-{r['key']}", help=r["description"], on_click=_rock,
                         args=(sess, r["label"]), disabled=running, width="stretch")
        rs = card.get("reactants") or []
        if rs:
            for c, (i, r) in zip(st.columns(min(len(rs), 4) or 1) * 4, enumerate(rs)):
                c.button(f"✕ {_reactant_text(r)}", key=f"drop-r-{i}", on_click=_drop_reactant, args=(sess, i),
                         help=f"Remove {r['name']}")
        else:
            st.caption("None: only temperature and pressure change")
        st.text_input("Add a reactant", key="add-reactant", placeholder="e.g. 10 g albite or 0.05 mol HCl, press Enter",
                      on_change=_add_reactant, args=(sess,), label_visibility="collapsed")
        a, b = st.columns(2)
        choice(sess, card, "kinetics", a)
        if kin:
            num(sess, card, "surface_area_cm2_per_g", b)
            num(sess, card, "duration", a)
            choice(sess, card, "duration_unit", b)
            choice(sess, card, "carbonate_rates", a)
            choice(sess, card, "depth_axis", b)
            if str(card.get("carbonate_rates", "")).startswith("seawater"):
                num(sess, card, "shell_area_exponent", a)
        st.markdown("Blocked minerals <span style='color:#6f777a;font-size:0.85rem'>(supersaturated but too slow to "
                    "form)</span>", unsafe_allow_html=True)
        blocked = card.get("blocked_minerals") or []
        if blocked:
            cols = st.columns(min(len(blocked), 5))
            for i, m in enumerate(blocked):
                cols[i % len(cols)].button(f"✕ {m}", key=f"drop-b-{m}", help=f"Allow {m}", on_click=_set,
                                           args=(sess, "blocked_minerals", [x for x in blocked if x != m]))
        else:
            st.caption("None: every mineral in the database may form")
        st.text_input("Add a blocked mineral", key="add-blocked", placeholder="add a mineral name, press Enter",
                      on_change=_add_blocked, args=(sess,), label_visibility="collapsed")


def _surfaces_box(sess, card) -> None:
    with st.container(border=True):
        st.markdown("**Surfaces and sorption**")
        a, b = st.columns(2)
        num(sess, card, "exchange_meq", a)
        num(sess, card, "hfo_g", b)
        st.caption("Clays (g per kg water, 0 = none). Each carries exchange sites at its cation exchange capacity "
                   "(kaolinite 20, illite 225, montmorillonite 764 meq per kg of clay).")
        a, b, c = st.columns(3)
        num(sess, card, "kaolinite_g", a)
        num(sess, card, "illite_g", b)
        num(sess, card, "montmorillonite_g", c)
        a, b = st.columns(2)
        choice(sess, card, "clay_edges", a)
        num(sess, card, "kd_solid_g", b)
        choice(sess, card, "kd_medium", a)


def _options_box(sess, card) -> None:
    with st.container(border=True):
        st.markdown("**Options**")
        cols = st.columns(4)
        for i, f in enumerate(("solid_solutions", "follow_phase", "isotopes", "gas_phase", "methane_sulfate",
                               "redox_apart", "minerals_form", "energy")):
            choice(sess, card, f, cols[i % 4])
        num(sess, card, "uncertainty_runs", cols[0])


def run_bar(sess: Session, input_page: bool) -> None:
    rs = sess.run_state
    running = rs["state"] == "running"
    cols = st.columns([1, 7])
    if input_page:
        if cols[0].button("Running…" if running else "Run", type="primary", disabled=running, width="stretch"):
            sess.start_run()
            st.rerun()
    if running:
        cols[1].progress(min(max(rs["progress"], 0.0), 1.0), text=rs["message"] or "Calculating")
    elif rs["state"] == "error":
        cols[1].error(rs["error"])
    follow(sess, "input" if input_page else "results")


def _pick_run(sess: Session, key: str) -> None:
    rid = st.session_state.get(key)
    if rid:
        sess.set_current(rid)


def _history(sess: Session, snap: dict) -> None:
    with st.container(border=True):
        st.markdown("**Earlier runs**")
        runs = snap["runs"]
        if not runs:
            st.caption("No runs yet.")
        for r in runs:
            label = (f"**{r['title']}:** " if r["title"] else "") + f"{tidy(r['summary'])} · {r['created'].replace('T', ' ')[:16]}"
            st.button(label, key=f"hist-{r['id']}", on_click=sess.set_current, args=(r["id"],),
                      type="primary" if r["id"] == snap["current"] else "tertiary")


# ---- the Results window -------------------------------------------------------------------------------
MODE_LABEL = {"static": "Static", "flow_through": "Flow-through", "outlet": "At the outlet",
              "start": "Column at the start", "profile": "Column at the end"}
UNIT_SCALE = {"mol": 1, "mmol": 1e3, "umol": 1e6, "mg": 1e3, "ug": 1e6}
UNIT_TEXT = {"mol": "mol", "mmol": "mmol", "umol": "µmol", "mg": "mg", "ug": "µg"}
TABS = [("species", "Species"), ("elements", "Elements"), ("minerals", "Mineral amounts"), ("gas", "Gas phase"),
        ("kinetic", "Rates"), ("energy", "Energy"), ("composition", "Solid solutions"), ("phases", "Brine and vapor"),
        ("isotopes", "Isotopes"), ("sorbed", "Sorbed"), ("partition", "Partition coefficients"),
        ("saturation", "Saturation index"), ("property", "Properties")]
PROPS = {"pH": "pH", "ionic_strength": "Ionic strength", "water_kg": "Solvent mass", "log_aw": "log water activity",
         "log_fO2": "log fO₂", "Eh": "Eh (V)", "surface_potential_mV": "Surface potential (mV)",
         "net_acidity": "Net acidity (mg/kg as CaCO₃)"}
AMOUNT_GROUPS = ("species", "elements", "minerals", "gas", "sorbed")


def _lead(res: dict, view: dict) -> dict:
    mid = next((m for m in view.get("modes", []) if m in res["modes"]), next(iter(res["modes"])))
    return res["modes"][mid]


def _frame(res: dict, view: dict) -> dict:
    r = _lead(res, view)
    return {"axes": r.get("axes") or res["axes"], "varies": r.get("varies") or res.get("varies") or ["temperature_c"],
            "rel": r.get("reliability") or res.get("reliability")}


def label_of(item: dict, res: dict) -> str:
    cat, g, k = res["catalog"], item["group"], item["key"]
    if g == "species":
        return cat["species"][k]["label"]
    if g == "elements":
        return f"{k} (total dissolved)"
    if g == "minerals":
        return cat["minerals"][k]["label"]
    if g == "saturation":
        return "SI " + cat["minerals"][k]["label"]
    if g == "property":
        return PROPS.get(k, k)
    return (cat.get(g) or {}).get(k, {}).get("label", k)


def _values(item: dict, mode: str, view: dict, res: dict):
    run = res["modes"].get(mode)
    if not run:
        return None
    arr = run["properties"].get(item["key"]) if item["group"] == "property" else (run.get(item["group"]) or {}).get(item["key"])
    if arr is None:
        return None
    if item["group"] in AMOUNT_GROUPS:
        k = UNIT_SCALE[view["unit"]]
        if view["unit"] in ("mg", "ug"):
            c = res["catalog"].get("minerals" if item["group"] == "minerals" else item["group"], {})
            mw = (c.get(item["key"]) or {}).get("mw")
            if not mw:
                return None
            k *= mw
        return [None if v is None else v * k for v in arr]
    return arr


def _band(item: dict, mode: str, view: dict, res: dict):
    run = res["modes"].get(mode) or {}
    b = (run.get("bands") or {}).get("properties" if item["group"] == "property" else item["group"])
    pair = b and b.get(item["key"])
    if not pair:
        return None
    line = _values(item, mode, view, res)
    raw = run["properties"][item["key"]] if item["group"] == "property" else run[item["group"]][item["key"]]
    k = 1.0
    for v, w in zip(raw, line or []):
        if v and w is not None:
            k = w / v
            break
    return [[None if v is None else v * k for v in arr] for arr in pair]


def y_label(items: list, view: dict, res: dict) -> str:
    if not items:
        return ""
    g, k = items[0]["group"], items[0]["key"]
    if g in ("species", "elements"):
        return f"{UNIT_TEXT[view['unit']]} per kg of water"
    if g in ("minerals", "gas", "sorbed"):
        return f"{UNIT_TEXT[view['unit']]} per kg of original water"
    if g == "saturation":
        return "log Q/K (0 = at saturation)"
    if g == "partition":
        return "Kd, L/kg"
    if g == "composition":
        return "mole fraction in the solid solution"
    if g == "isotopes":
        c = (res["catalog"].get("isotopes") or {}).get(k, {})
        return "87Sr/86Sr" if k.startswith("87") else f"{k.split(':')[0].replace('d', 'δ')} ({c.get('unit', '‰')})"
    if g == "phases":
        return "% of the fluid's mass in the vapor" if k == "vapour_share" else "wt% NaCl (equivalent)"
    if g == "energy":
        if k.startswith("dG:"):
            return "ΔG, kJ per mole of reaction (below 0: energy released)"
        return "energy available, J per kg of the starting fluid" if k.startswith("Es:") else "energy available, J per kg of fluid"
    if g == "kinetic":
        c = (res["catalog"].get("kinetic") or {}).get(k, {})
        return f"{c.get('unit', 'mol/s')} per kg water (+ dissolves)" if k.startswith("rate:") else "% of starting amount left"
    return PROPS.get(k, k)


def _log_ok(items: list) -> bool:
    if not items:
        return False
    g, k = items[0]["group"], items[0]["key"]
    return (g in ("species", "elements", "minerals", "gas", "sorbed", "partition", "phases")
            or (g == "property" and k == "ionic_strength") or (g == "kinetic" and k.startswith("rate:"))
            or (g == "energy" and not k.startswith("dG:")))


def _candidates(tab: str, res: dict, view: dict) -> list[dict]:
    run = _lead(res, view)
    cat = res["catalog"]
    if tab == "species":
        peak = lambda a: max([v for v in a if v is not None] + [0])
        keys = sorted(run["species"], key=lambda k: -peak(run["species"][k]))
        return [{"group": "species", "key": k, "label": cat["species"][k]["label"]} for k in keys]
    if tab == "saturation":
        return [{"group": "saturation", "key": k, "label": cat["minerals"][k]["label"]
                 + (" (blocked)" if cat["minerals"][k].get("blocked") else "")} for k in run["saturation"]]
    if tab == "property":
        return [{"group": "property", "key": k, "label": v} for k, v in PROPS.items() if k in run["properties"]]
    return [{"group": tab, "key": k, "label": label_of({"group": tab, "key": k}, res)} for k in (run.get(tab) or {})]


def _tabs(res: dict, view: dict) -> list[tuple[str, str]]:
    run0 = _lead(res, view)
    always = {"species", "elements", "minerals", "gas", "saturation", "property"}
    return [(t, lab) for t, lab in TABS if t in always or
            ("vapour_share" in (run0.get("phases") or {}) if t == "phases" else bool(run0.get(t)))]


def _set_view(sess: Session, patch: dict) -> None:
    sess.set_view(patch)


def _fmt(v) -> str:
    return "" if v is None else f"{v:.4g}"


def _series(res, other, view):
    modes = [m for m in view["modes"] if m in res["modes"]]
    out = []
    for i, item in enumerate(view["items"]):
        color = PALETTE[i % len(PALETTE)]
        for m in modes:
            y = _values(item, m, view, res)
            if y is None:
                continue
            out.append({"label": label_of(item, res) + (f" ({MODE_LABEL.get(m, m).lower()})" if len(modes) > 1 else ""),
                        "color": color, "dash": "dash" if m in ("flow_through", "start") else None, "y": y,
                        "band": _band(item, m, view, res), "x": None,
                        "csv": f"{label_of(item, res)}{' ' + MODE_LABEL.get(m, m).lower() if len(modes) > 1 else ''} "
                               f"[{y_label([item], view, res)}]"})
        ox = other and _frame(other, view)["axes"].get(view.get("x"))
        if not ox:
            continue
        shared = [m for m in modes if m in other["modes"]]
        for m in shared or list(other["modes"]):
            oy = _values(item, m, view, other)
            if oy is None:
                continue
            out.append({"label": f"{label_of(item, res)} (compared run)", "color": color, "dash": "dot", "y": oy,
                        "x": ox["values"], "band": None, "compare": True})
    return out


def results_page() -> None:
    st.set_page_config(page_title="Reaction | Results", page_icon=ICON, layout="wide")
    wid, sess = workspace()
    snap = sess.snapshot()
    top = st.columns([4, 1.4])
    with top[0]:
        header("results")
    with top[1]:
        st.link_button("Back to the input window", f"{PATHS['input']}?w={wid}", width="stretch")
    run_bar(sess, input_page=False)
    res = sess.result_with_story(snap["current"]) if snap["current"] else None
    if res is None:
        st.info("No run yet. Set up the scenario card in the Input window and press Run; this window follows it live.")
        if snap["diagram"]:
            _diagram(sess, snap, None)
        footer(wid)
        return
    view = sess.view
    st.subheader(res.get("title") or res["scenario"].get("title") or "Run")
    st.caption(f"{tidy(res['summary'])} · {res['created'].replace('T', ' ')} · {res['elapsed_s']} s")

    runs = snap["runs"]
    names = {r["id"]: f"{r['created'][5:16].replace('T', ' ')} · {tidy(r['summary'])}" for r in runs}
    a, b = st.columns(2)
    st.session_state["pick-run"] = res["id"]
    a.selectbox("Run shown", list(names), key="pick-run", format_func=names.get, on_change=_pick_run,
                args=(sess, "pick-run"))
    cmp_ids = [""] + [r for r in names if r != res["id"]]
    st.session_state["pick-cmp"] = view.get("compare") or ""
    b.selectbox("Compare with", cmp_ids, key="pick-cmp", format_func=lambda r: names.get(r, "Nothing"),
                on_change=lambda: _set_view(sess, {"compare": st.session_state["pick-cmp"] or None}))
    other = sess.result(view.get("compare")) if view.get("compare") else None

    with st.container(border=True):
        st.markdown("**What we found**")
        for p in res.get("story") or []:
            st.markdown(tidy(p))
    _plot_panel(sess, res, other, view)
    _end_lists(res.get("end") or [])
    st.markdown("**Diagrams**")
    cols = st.columns(4)
    for c, (kind, lab) in zip(cols, [("eh_ph", "Eh-pH diagram"), ("fo2_ph", "log fO₂-pH diagram"),
                                      ("activity", "Activity diagram")]):
        c.button(lab, key=f"dia-open-{kind}", on_click=sess.set_diagram, args=({"kind": kind},),
                 width="stretch")
    if snap["diagram"]:
        _diagram(sess, snap, res["id"])
    with st.expander("Assumptions and limits"):
        for n in res.get("notes") or []:
            st.markdown(f"- {tidy(n)}")
        if res.get("database"):
            st.markdown(f"- Thermodynamic data: {res['database']}.")
    footer(wid)


def _plot_panel(sess, res, other, view) -> None:
    fr = _frame(res, view)
    varies = fr["varies"]
    xkey = view.get("x") if view.get("x") in varies else varies[0]
    xax = fr["axes"][xkey] if fr["axes"] else res["x"]
    items = view["items"]
    g = items[0]["group"] if items else ""
    ctrl = st.columns([1.4, 1.4, 1.2, 2.2])
    if len(varies) > 1:
        st.session_state["view-x"] = xkey
        ctrl[0].segmented_control("Horizontal axis", varies, key="view-x", format_func=lambda k: fr["axes"][k]["label"],
                                  on_change=lambda: _set_view(sess, {"x": st.session_state["view-x"] or xkey}))
    mode_ids = list(res["modes"])
    if len(mode_ids) > 1:
        st.session_state["view-modes"] = [m for m in view["modes"] if m in mode_ids]
        ctrl[1].pills("Modes", mode_ids, key="view-modes", selection_mode="multi",
                      format_func=lambda m: MODE_LABEL.get(m, m),
                      on_change=lambda: _set_view(sess, {"modes": st.session_state["view-modes"] or view["modes"]}))
    log_ok = _log_ok(items)
    st.session_state["view-log"] = "Log" if (view.get("log") and log_ok) else "Linear"
    ctrl[2].segmented_control("Axis", ["Log", "Linear"], key="view-log", disabled=not log_ok,
                              on_change=lambda: _set_view(sess, {"log": st.session_state["view-log"] == "Log"}))
    if g in AMOUNT_GROUPS:
        st.session_state["view-unit"] = view["unit"]
        ctrl[3].segmented_control("Units", list(UNIT_TEXT), key="view-unit", format_func=UNIT_TEXT.get,
                                  on_change=lambda: _set_view(sess, {"unit": st.session_state["view-unit"] or view["unit"]}))

    series = _series(res, other, view)
    use_log = bool(view.get("log")) and log_ok
    fig = go.Figure()
    for s in series:
        x = s["x"] or xax["values"]
        if s.get("band"):
            lo, hi = s["band"]
            fig.add_scatter(x=list(x) + list(x)[::-1], y=list(hi) + list(lo)[::-1], fill="toself", mode="lines",
                            line=dict(width=0), fillcolor=s["color"], opacity=0.18, hoverinfo="skip", showlegend=False)
        fig.add_scatter(x=x, y=s["y"], name=s["label"], mode="lines+markers", marker=dict(size=4),
                        line=dict(color=s["color"], dash=s["dash"], width=2.4 if not s.get("compare") else 1.6))
    if g == "saturation" or (g == "energy" and items and items[0]["key"].startswith("dG:")):
        fig.add_hline(y=0, line=dict(color="#8d8b82", width=1, dash="dot"))
    rel = (fr["rel"] or {}).get("levels") or []
    if any(lv != "good" for lv in rel):
        colors = {"good": "rgba(94,140,106,0.10)", "caution": "rgba(196,154,85,0.18)", "poor": "rgba(180,83,42,0.18)"}
        xs = xax["values"]
        for i, lv in enumerate(rel):
            if lv != "good" and i < len(xs):
                x0 = (xs[i - 1] + xs[i]) / 2 if i else xs[i]
                x1 = (xs[i] + xs[i + 1]) / 2 if i + 1 < len(xs) else xs[i]
                fig.add_vrect(x0=x0, x1=x1, fillcolor=colors[lv], line_width=0, layer="below")
    fig.update_layout(height=500, margin=dict(l=10, r=10, t=30, b=10), hovermode="x unified",
                      legend=dict(orientation="h", y=-0.18), paper_bgcolor=PLOT_BG, plot_bgcolor=PLOT_BG,
                      font=dict(color="#2B3134"), title=dict(text=tidy(res["summary"]), font=dict(size=13)),
                      xaxis_title=f"{xax['label']} ({xax['unit']})", yaxis_title=y_label(items, view, res))
    fig.update_xaxes(gridcolor=GRID)
    fig.update_yaxes(gridcolor=GRID, type="log" if use_log else "linear", exponentformat="power")
    st.plotly_chart(fig, width="stretch")
    if any(lv != "good" for lv in rel):
        st.caption("Shaded columns: points where the thermodynamic data are less reliable (amber approximate, rust "
                   "unreliable), from the water-ionization check against Marshall and Franck (1981).")
    banded = next((r for r in res["modes"].values() if r.get("band_runs")), None)
    if banded and any(s.get("band") for s in series):
        st.caption(f"Shaded: the middle 90% of {banded['band_runs']} runs with the log K values shifted at random by "
                   "their uncertainty (see Assumptions and limits).")
    if view.get("log") and not log_ok and items:
        st.caption("Log axis is not used for this quantity.")

    with st.container(border=True):
        st.markdown("**What to plot**")
        tabs = _tabs(res, view)
        tab_ids = [t for t, _ in tabs]
        st.session_state.setdefault("pick-tab", "species")
        if st.session_state["pick-tab"] not in tab_ids:
            st.session_state["pick-tab"] = "species"
        st.segmented_control("Group", tab_ids, key="pick-tab", format_func=dict(tabs).get, label_visibility="collapsed")
        tab = st.session_state["pick-tab"] or "species"
        cands = _candidates(tab, res, view)
        keyed = {f"{c['group']}|{c['key']}": c for c in cands}
        chosen = [f"{i['group']}|{i['key']}" for i in items if f"{i['group']}|{i['key']}" in keyed]
        st.session_state[f"pick-items-{tab}"] = chosen

        def _changed(tab=tab, chosen=chosen):
            now = set(st.session_state[f"pick-items-{tab}"])
            for k in now ^ set(chosen):
                grp, key = k.split("|", 1)
                sess.set_view({"toggle": {"group": grp, "key": key}})
        st.multiselect("Pick what to plot (type to search)", list(keyed), key=f"pick-items-{tab}",
                       format_func=lambda k: keyed[k]["label"], on_change=_changed)
        st.caption("Plotting a quantity from another group replaces the lines shown, since one axis holds one kind "
                   "of quantity.")

    fr_axes = [fr["axes"][k] for k in dict.fromkeys(list(varies) + ["temperature_c", "pressure_bar"]) if k in fr["axes"]]
    main = [s for s in series if not s.get("compare")]
    head = [f"{a['label']} ({a['unit']})" for a in fr_axes]
    for s in main:
        head += [s["csv"], s["csv"] + " 5%", s["csv"] + " 95%"] if s.get("band") else [s["csv"]]
    rows = []
    for i in range(len(fr_axes[0]["values"]) if fr_axes else 0):
        row = [a["values"][i] for a in fr_axes]
        for s in main:
            row.append(s["y"][i])
            if s.get("band"):
                row += [s["band"][0][i], s["band"][1][i]]
        rows.append(row)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(head)
    w.writerows(rows)
    a, b, _ = st.columns([1, 1, 2])
    a.download_button("Download plotted data (CSV)", buf.getvalue(), file_name=f"{res['id']}.csv", mime="text/csv",
                      width="stretch")
    b.download_button("Download the whole run (JSON)", json.dumps(sess.result(res["id"])), file_name=f"{res['id']}.json",
                      mime="application/json", width="stretch")


def _end_lists(ends: list) -> None:
    for e in ends:
        with st.container(border=True):
            st.markdown(f"**{e['title']}**")
            a, b, c = st.columns(3)
            if e["solids"]:
                a.caption("Solids, per kg of starting water")
                a.dataframe([dict(zip(["Mineral", "mmol", "g"], r)) for r in e["solids"]], hide_index=True,
                            width="stretch")
            else:
                a.caption("Solids")
                a.markdown("No minerals are left at this point.")
            b.caption("Fluid")
            b.dataframe([{"Property": tidy(r[0]), "Value": tidy(r[1])} for r in e["fluid"]], hide_index=True,
                        width="stretch")
            c.caption("Dissolved")
            c.dataframe([dict(zip(["Element", "mmol/kg", "mg/kg"], r)) for r in e["dissolved"]], hide_index=True,
                        width="stretch")


def _dia_edit(sess: Session, key: str, field: str) -> None:
    v = st.session_state.get(key)
    sess.set_diagram({field: "" if v is None else v})


def _diagram(sess: Session, snap: dict, run_id: str | None) -> None:
    spec = snap["diagram"]
    with st.container(border=True):
        top = st.columns([5, 1])
        data = sess.diagram_data(run_id)
        top[0].markdown(f"**{tidy((data or {}).get('title') or 'Diagram')}**")
        top[1].button("Close", key="dia-close", on_click=sess.set_diagram, args=({"close": True},),
                      width="stretch")
        act = spec["kind"] == "activity"
        c = st.columns(5)
        if act:
            st.session_state["dia-cat"] = spec.get("cation") or "K+"
            c[0].text_input("Cation", key="dia-cat", on_change=_dia_edit, args=(sess, "dia-cat", "cation"))
        else:
            st.session_state["dia-el"] = spec.get("element") or ""
            c[0].text_input("Element", key="dia-el", on_change=_dia_edit, args=(sess, "dia-el", "element"))
            st.session_state["dia-a"] = float(spec.get("log_a", -6.0))
            c[3].number_input("log activity", key="dia-a", step=1.0, format="%g", on_change=_dia_edit, args=(sess, "dia-a", "log_a"))
        t_c = spec.get("t_c") if spec.get("t_c") is not None else (data or {}).get("t_c")
        p_bar = spec.get("p_bar") if spec.get("p_bar") is not None else (data or {}).get("p_bar")
        st.session_state["dia-t"] = None if t_c is None else float(t_c)
        c[1].number_input("Temperature (°C)", key="dia-t", step=10.0, format="%g", on_change=_dia_edit, args=(sess, "dia-t", "t_c"))
        st.session_state["dia-p"] = None if p_bar is None else float(p_bar)
        c[2].number_input("Pressure (bar)", key="dia-p", step=10.0, format="%g", on_change=_dia_edit, args=(sess, "dia-p", "p_bar"))
        st.session_state["dia-with"] = " ".join(spec.get("with") or [])
        c[4].text_input("With", key="dia-with", help="Other elements to include, e.g. S C",
                        on_change=_dia_edit, args=(sess, "dia-with", "with"))
        if not data or data.get("none"):
            return
        if data.get("error"):
            st.error(tidy(data["error"]))
            return
        _draw_diagram(data)
        rel = data["reliability"]
        level = {"good": "good", "caution": "approximate", "poor": "unreliable"}[rel["level"]]
        st.caption(f"Hatched: the boundary is within {rel['band']} log units, the uncertainty of log K here (water "
                   f"ionization {level} at this temperature and pressure). Dashed: water stability limits. Black: the "
                   "run's path from the open square (start) to the red dot (end).")
        for n in data.get("notes") or []:
            st.caption(f"- {sentence(n)}")


def _draw_diagram(d: dict) -> None:
    nx, ny = d["nx"], d["ny"]
    x0, x1, y0, y1 = d["x"]["lo"], d["x"]["hi"], d["y"]["lo"], d["y"]["hi"]
    xs = [x0 + (x1 - x0) * (i + 0.5) / nx for i in range(nx)]
    ys = [y0 + (y1 - y0) * (j + 0.5) / ny for j in range(ny)]
    fields = {f["index"]: f for f in d["fields"]}
    idx = sorted(fields)
    pos = {k: n for n, k in enumerate(idx)}
    fill = {"aqueous": ["#E7E2D4", "#DCE3E2", "#E9DCD0", "#DDE2D6", "#E4DCE6", "#D9E1E8"],
            "mineral": ["#CDB98F", "#B7C7CB", "#D2AE97", "#BFCBB0", "#CBB9CF", "#AFC3D2"]}
    colors = [fill.get(fields[k]["kind"], fill["mineral"])[n % 6] for n, k in enumerate(idx)]
    z = [[pos.get(v, None) for v in row] for row in d["grid"]]
    n = max(len(idx), 1)
    scale = []
    for i, ccol in enumerate(colors):
        scale += [[i / n, ccol], [(i + 1) / n, ccol]]
    text = [[fields[v]["label"] if v in fields else "" for v in row] for row in d["grid"]]
    fig = go.Figure(go.Heatmap(x=xs, y=ys, z=z, zmin=0, zmax=n, colorscale=scale or [[0, "#ddd"], [1, "#ddd"]],
                               showscale=False, text=text, hovertemplate="%{text}<br>x %{x:.2f}, y %{y:.2f}<extra></extra>"))
    unc = d.get("uncertain") or []
    ux, uy = [], []
    for j, row in enumerate(unc):
        for i, ch in enumerate(row):
            if ch == "1" and (i + j) % 2 == 0:
                ux.append(xs[i])
                uy.append(ys[j])
    if ux:
        fig.add_scatter(x=ux, y=uy, mode="markers", marker=dict(size=2, color="#8a5a00"), hoverinfo="skip",
                        showlegend=False)
    for f in d["fields"]:
        if not f.get("small"):
            fig.add_annotation(x=f["at"][0], y=f["at"][1], text=f["label"] + ("*" if f.get("extrapolated") else ""),
                               showarrow=False, font=dict(size=12, color="#2B3134"))
    for line in (d.get("water") or {}).values():
        fig.add_scatter(x=[p[0] for p in line], y=[p[1] for p in line], mode="lines", showlegend=False,
                        line=dict(color="#3E5C6B", dash="dash", width=1.5), hoverinfo="skip")
    for mode, pts in (d.get("path") or {}).items():
        t0 = d.get("t_c")
        fig.add_scatter(x=[p[0] for p in pts], y=[p[1] for p in pts], mode="lines+markers", name=f"run path ({mode})",
                        line=dict(color="#1b1e20", width=2),
                        marker=dict(size=[7 if t0 is not None and abs(p[2] - t0) <= 10 else 4 for p in pts],
                                    color="#1b1e20"),
                        text=[f"{p[2]:.0f} °C" for p in pts], hovertemplate="%{text}<extra></extra>", showlegend=False)
        if pts:
            fig.add_scatter(x=[pts[0][0]], y=[pts[0][1]], mode="markers", showlegend=False, hoverinfo="skip",
                            marker=dict(symbol="square-open", size=11, color="#1b1e20"))
            fig.add_scatter(x=[pts[-1][0]], y=[pts[-1][1]], mode="markers", showlegend=False, hoverinfo="skip",
                            marker=dict(size=10, color="#B4532A"))
    fig.update_layout(height=520, margin=dict(l=10, r=10, t=10, b=10), paper_bgcolor=PLOT_BG, plot_bgcolor=PLOT_BG,
                      xaxis_title=d["x"]["label"], yaxis_title=d["y"]["label"], font=dict(color="#2B3134"))
    fig.update_xaxes(range=[x0, x1])
    fig.update_yaxes(range=[y0, y1])
    st.plotly_chart(fig, width="stretch")


# ---- the Sources page -------------------------------------------------------------------------------
def sources_page() -> None:
    st.set_page_config(page_title="Reaction | Sources", page_icon=ICON, layout="wide")
    wid = st.query_params.get("w") or ""
    header("sources")
    st.link_button("Back to the input window", f"{PATHS['input']}?w={wid}")
    html = (STATIC / "sources.html").read_text(encoding="utf-8")
    body = re.search(r"<main[^>]*>(.*)</main>", html, re.S)
    content = body.group(1) if body else re.sub(r"(?s).*<body[^>]*>|</body>.*", "", html)
    content = re.sub(r'<div class="bar">.*?</div>', "", content, flags=re.S)
    st.html(f"<div style='max-width:900px;line-height:1.55;color:#2B3134'>{content}</div>")
    footer(wid)

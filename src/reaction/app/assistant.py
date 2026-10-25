"""Turns what the user types into changes to the scenario card and the plot.

The local language model (Ollama) is asked for a fixed JSON shape: a short reply, the card fields
to change, species to show or hide, whether to run, and at most one question.  It never reports
results or does arithmetic.  If the model is not reachable, or answers in the wrong shape, a small
set of plain rules picks out the obvious details (temperatures, pressure, mode, species to plot) so
the app keeps working without it.
"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.request
from dataclasses import dataclass, field

from . import analysis, config
from .fluids import CUSTOM, PRESETS
from .scenario import FIELD_HELP, Scenario

SCHEMA = {
    "type": "object",
    "properties": {
        "reply": {"type": "string"},
        "scenario_changes": {
            "type": "object",
            "properties": {
                "title": {"type": "string"},
                "fluid": {"type": "string", "enum": [*PRESETS, CUSTOM]},
                "ph": {"type": "number"},
                "fluid_t_c": {"type": "number"},
                "t_start_c": {"type": "number"},
                "t_end_c": {"type": "number"},
                "steps": {"type": "integer"},
                "p_start_bar": {"type": "number"},
                "p_end_bar": {"type": "number"},
                "mode": {"type": "string", "enum": ["static", "flow_through", "both"]},
                "blocked_minerals": {"type": "array", "items": {"type": "string"}},
            },
            "additionalProperties": False,
        },
        "show": {"type": "array", "items": {"type": "string"}},
        "hide": {"type": "array", "items": {"type": "string"}},
        "axis": {"type": "string", "enum": ["log", "linear", "unchanged"]},
        "x_axis": {"type": "string", "enum": ["temperature", "pressure", "unchanged"]},
        "run": {"type": "boolean"},
        "question": {"type": "string"},
    },
    "required": ["reply", "scenario_changes", "show", "hide", "axis", "x_axis", "run", "question"],
}

SYSTEM = """You are the assistant inside a personal geochemical calculator. The user describes a \
model in plain language. You do two jobs: fill in the scenario card, and control what the plot shows.

Rules:
- Change only the card fields the user actually asked for or clearly implied. Never invent values.
- You never state or guess any chemical results (concentrations, pH, minerals). The calculator does that.
- Temperatures are degrees C. Pressure is in bar (1 MPa = 10 bar, 1 kbar = 1000 bar).
- A run follows a straight line from a start point to an end point in temperature and pressure. To hold
  pressure fixed at X, set both p_start_bar and p_end_bar to X; likewise for temperature. Either or both
  may change (heating, cooling, compression, decompression).
- A pasted water analysis is read by the calculator itself; do not copy it into the card.
- Only these path runs exist so far. If the user asks for something else, say so briefly in reply.
- "show" lists species, minerals, or properties to plot exactly as the user wrote them (for example
  "CaCl+", "anhydrite", "pH", "SO4--"). "hide" lists ones to remove. "axis" is log, linear or unchanged.
  "x_axis" is temperature, pressure or unchanged.
- Set run to true only when the user asks to run, calculate, simulate, react, or heat something.
- If something essential is unclear, put ONE short question in "question"; otherwise leave it empty.
- reply: one or two short sentences saying what you changed. No chemistry claims.

Card fields:
{fields}

Current card:
{card}

Currently plotted: {plotted}

Examples (user message -> JSON):
"heat seawater to 300 C at 400 bar, static" ->
{{"reply":"Set seawater heated to 300 C at 400 bar, static.","scenario_changes":{{"t_end_c":300,"p_start_bar":400,"p_end_bar":400,"mode":"static"}},"show":[],"hide":[],"axis":"unchanged","x_axis":"unchanged","run":true,"question":""}}
"decompress from 1000 to 200 bar at 300 C" ->
{{"reply":"Set 300 C held fixed, pressure from 1000 down to 200 bar.","scenario_changes":{{"t_start_c":300,"t_end_c":300,"p_start_bar":1000,"p_end_bar":200}},"show":[],"hide":[],"axis":"unchanged","x_axis":"pressure","run":true,"question":""}}
"plot CaCl+ instead of Ca++" ->
{{"reply":"Plotting CaCl+ in place of Ca++.","scenario_changes":{{}},"show":["CaCl+"],"hide":["Ca++"],"axis":"unchanged","x_axis":"unchanged","run":false,"question":""}}
"""


@dataclass
class Interpretation:
    reply: str = ""
    scenario_changes: dict = field(default_factory=dict)
    show: list[str] = field(default_factory=list)
    hide: list[str] = field(default_factory=list)
    axis: str = "unchanged"
    x_axis: str = "unchanged"              # "temperature_c", "pressure_bar" or "unchanged"
    run: bool = False
    question: str = ""
    source: str = "rules"                  # "model" or "rules"
    mode_mentioned: bool = False


class Ollama:
    def __init__(self, url: str = config.OLLAMA_URL, model: str = config.OLLAMA_MODEL):
        self.url, self.model = url.rstrip("/"), model

    def available(self, timeout: float = 1.5) -> bool:
        try:
            with urllib.request.urlopen(f"{self.url}/api/tags", timeout=timeout) as r:
                names = [m.get("name", "") for m in json.load(r).get("models", [])]
            return any(n == self.model or n.split(":")[0] == self.model.split(":")[0] for n in names)
        except (OSError, ValueError):
            return False

    def chat(self, messages: list[dict], timeout: float = 120.0) -> str:
        body = json.dumps({
            "model": self.model, "messages": messages, "stream": False, "format": SCHEMA,
            "keep_alive": "30m", "options": {"temperature": 0.1, "num_ctx": 4096},
        }).encode()
        req = urllib.request.Request(f"{self.url}/api/chat", data=body,
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.load(r)["message"]["content"]


def interpret(text: str, card: Scenario, plotted: list[str], history: list[dict],
              model: Ollama | None) -> Interpretation:
    """The model's reading of ``text`` if it works, otherwise the rules' reading."""
    rules = interpret_rules(text, card)
    if model is None or "analysis" in rules.scenario_changes:     # a pasted analysis is read exactly
        return rules
    system = SYSTEM.format(
        fields="\n".join(f"- {k}: {v}" for k, v in FIELD_HELP.items()),
        card=json.dumps(card.to_dict()), plotted=", ".join(plotted) or "nothing yet")
    messages = [{"role": "system", "content": system}]
    for turn in history[-6:]:
        messages.append({"role": "user" if turn["role"] == "user" else "assistant",
                         "content": turn["text"]})
    messages.append({"role": "user", "content": text})
    try:
        data = json.loads(model.chat(messages))
        got = Interpretation(
            reply=str(data.get("reply", "")).strip(),
            scenario_changes={k: v for k, v in dict(data.get("scenario_changes") or {}).items()
                              if v not in (None, "")},
            show=[str(s) for s in data.get("show", []) if str(s).strip()],
            hide=[str(s) for s in data.get("hide", []) if str(s).strip()],
            axis=data.get("axis", "unchanged") if data.get("axis") in ("log", "linear") else "unchanged",
            x_axis={"temperature": "temperature_c", "pressure": "pressure_bar"}.get(
                data.get("x_axis"), "unchanged"),
            run=bool(data.get("run", False)), question=str(data.get("question", "")).strip(),
            source="model")
    except (OSError, ValueError, KeyError, TypeError):
        rules.reply = ("The language model did not answer, so I read your message with simple rules "
                       "and may have missed details. " + rules.reply).strip()
        return rules
    # A small model sometimes invents values.  Numbers and modes are kept only if the user's own
    # words contain them, and the rules' reading of those fields wins over the model's.
    got.scenario_changes = _guard(got.scenario_changes, text)
    for group in (T_KEYS, P_KEYS):                    # the rules' reading of T or P replaces the model's
        if any(k in rules.scenario_changes for k in group):
            for k in group:
                got.scenario_changes.pop(k, None)
    for key in (*T_KEYS, *P_KEYS, "mode", "fluid", "steps", "analysis", "units", "ph", "reactants",
                "mix_with", "mix_ratio", "mix_t_c", "redox", "redox_value"):
        if key in rules.scenario_changes:
            got.scenario_changes[key] = rules.scenario_changes[key]
    if got.x_axis == "unchanged":
        got.x_axis = rules.x_axis
    got.mode_mentioned = rules.mode_mentioned
    plot_request = bool(rules.show or rules.hide)
    got.run = rules.run if plot_request else (got.run or rules.run)
    got.show = got.show or rules.show
    got.hide = got.hide or rules.hide
    if got.axis == "unchanged":
        got.axis = rules.axis
    return got


T_KEYS = ("t_start_c", "t_end_c", "temperature_c")
P_KEYS = ("p_start_bar", "p_end_bar", "pressure_bar")
_P_UNITS = {"bar": 1.0, "bars": 1.0, "mpa": 10.0, "kbar": 1000.0, "gpa": 10000.0, "atm": 1.01325}
_BLOCK_WORDS = re.compile(r"block|suppress|prevent|exclude|forbid|not allowed|can'?t precipitate|cannot|"
                          r"without|\bno\b", re.I)


def _guard(changes: dict, text: str) -> dict:
    nums = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", text)]

    def said(value, scales=(1.0,)) -> bool:
        return any(abs(float(value) - n * k) <= 1e-9 * max(1.0, abs(float(value)))
                   for n in nums for k in scales)

    out = {}
    for key, value in changes.items():
        try:
            if key in (*T_KEYS, "steps") and not said(value):
                continue
            if key in P_KEYS and not said(value, (1.0, 10.0, 1000.0, 10000.0, 1.01325)):
                continue
            if key in ("ph", "fluid_t_c") and not said(value):
                continue
        except (TypeError, ValueError):
            continue
        if key in ("mode", "fluid"):
            continue                                   # taken from the rules only
        if key == "blocked_minerals" and not _BLOCK_WORDS.search(text):
            continue
        out[key] = value
    return out


# ---- plain rules ------------------------------------------------------------------------------------
_NUM = r"(-?\d+(?:\.\d+)?)"
_DEG = r"\s*(?:°|deg(?:rees?)?)?\s*c(?![a-z])"
RUN_WORDS = re.compile(r"\b(run|go|start|simulate|calculate|compute|react|titrate|mix|add|heat|warm|cool|decompress|evaporate|"
                       r"compress|depressuri[sz]e|pressuri[sz]e|do it|what happens)\b", re.I)


_FLUID_WORDS = r"seawater|sea water|vent fluid|vent|black smoker fluid|rainwater|rain|formation brine|brine"
_MIX = re.compile(rf"\bmix(?:es|ed|ing)?\s+(?:(?:the\s+)?(?P<a>{_FLUID_WORDS})\b[^.;]*?\s+)?(?:with|into|and)\s+"
                  rf"(?:(?:up\s+to\s+)?\d*\.?\d+\s*(?:kg|parts?)\s+(?:of\s+)?)?(?:the\s+)?(?P<b>{_FLUID_WORDS})\b")


def _fluid_key(words: str) -> str | None:
    w = words.strip()
    if not w:
        return None
    if "vent" in w or "smoker" in w:
        return "vent_fluid"
    if "brine" in w:
        return "formation_brine"
    if "rain" in w:
        return "rainwater"
    return "seawater"


_MORE_REACTANT = re.compile(
    r"\s*(?:,\s*(?:and\s+)?|\s+and\s+|\s*\+\s*)(\d*\.?\d+)\s*(g|grams?|mg|kg|mol|moles?|mmol)\s+(?:of\s+)?"
    r"([A-Za-z][\w()+\-]*(?:,(?:ordered|disordered|low|high|alpha|beta|maximum)|\s+glass|\s+gas)?)", re.I)
_REACTANT = re.compile(
    r"\b(?:react(?:ing)?|add(?:ing)?|dissolv(?:e|ing)|titrat(?:e|ing)(?:\s+\w+)?\s+with|with)\s+"
    r"(?:up\s+to\s+)?(\d*\.?\d+)\s*(g|grams?|mg|kg|mol|moles?|mmol)\s+(?:of\s+)?"
    r"([A-Za-z][\w()+\-]*(?:,(?:ordered|disordered|low|high|alpha|beta|maximum)|\s+glass|\s+gas)?)", re.I)


def interpret_rules(text: str, card: Scenario) -> Interpretation:
    low = text.lower()
    out = Interpretation(source="rules")
    ch = out.scenario_changes

    unit_re = r"(bars?|mpa|kbar|gpa|atm)\b"
    t_range = re.search(rf"from\s+{_NUM}(?:{_DEG})?\s*(?:to|-|until)\s*{_NUM}{_DEG}", low)
    p_range = re.search(rf"(?:from\s+)?{_NUM}(?:\s*{unit_re})?\s*(?:to|-|until)\s*{_NUM}\s*{unit_re}", low)
    pressure_verb = bool(re.search(r"compress|pressuri|(lower|drop|raise|reduce|increase) the pressure", low))
    isothermal = bool(re.search(r"isotherm|constant temperature|fixed temperature|hold the temperature", low))

    t_single = None
    if t_range:
        ch["t_start_c"], ch["t_end_c"] = float(t_range.group(1)), float(t_range.group(2))
    else:
        m = re.search(rf"(?:to|up to|until|at)\s*{_NUM}{_DEG}", low) or re.search(rf"{_NUM}{_DEG}", low)
        if m:
            t_single = float(m.group(1))
    if t_single is not None:
        # a pressure sweep at "300 C" holds 300 C; "heat to 300 C" ends there
        if isothermal or p_range or pressure_verb:
            ch["temperature_c"] = t_single
        else:
            ch["t_end_c"] = t_single
    if p_range:
        k_end = _P_UNITS[p_range.group(4)]
        k_start = _P_UNITS[p_range.group(2)] if p_range.group(2) else k_end     # "from 2 kbar to 500 bar"
        ch["p_start_bar"], ch["p_end_bar"] = float(p_range.group(1)) * k_start, float(p_range.group(3)) * k_end
    else:
        m = re.search(rf"{_NUM}\s*{unit_re}", low)
        if m:
            value = float(m.group(1)) * _P_UNITS[m.group(2)]
            if pressure_verb:
                ch["p_end_bar"] = value                    # "decompress to 200 bar"
            else:
                ch["pressure_bar"] = value                 # "at 400 bar" holds it
    m = re.search(r"\b(?:against|versus|vs\.?|as a function of|over)\s+(?:the\s+)?(pressure|temperature)\b", low)
    if m:
        out.x_axis = "pressure_bar" if m.group(1) == "pressure" else "temperature_c"
    m = re.search(r"(\d+)\s*(?:steps|points)", low)
    if m:
        ch["steps"] = int(m.group(1))
    if re.search(r"flow[\s-]*through", low) and re.search(r"static", low) or re.search(r"\bboth\b|\bcompare\b", low):
        ch["mode"] = "both"
    elif re.search(r"flow[\s-]*through", low):
        ch["mode"] = "flow_through"
    elif re.search(r"\bstatic\b|\bclosed\b|\bbatch\b", low):
        ch["mode"] = "static"
    out.mode_mentioned = "mode" in ch
    for key, p in PRESETS.items():
        if key in low or p.label.lower() in low or (key == "vent_fluid" and re.search(r"\bvent\b|black smoker", low)) \
                or (key == "formation_brine" and re.search(r"\bbrine\b", low)) or (key == "rainwater" and re.search(r"\brain\b", low)):
            ch["fluid"] = key
    mix = _MIX.search(low)
    if mix:
        first, second = _fluid_key(mix.group("a") or ""), _fluid_key(mix.group("b"))
        if second:
            ch["mix_with"] = second
            if first:
                ch["fluid"] = first
            elif ch.get("fluid") == second:
                ch.pop("fluid")                        # "mix seawater ..." named only the second fluid
            ratio = re.search(r"(?:up to|ratio(?: of)?|to)\s*(\d*\.?\d+)\s*(?:kg|:\s*1\b|parts?)", low[mix.end():]) \
                or re.search(r"(\d*\.?\d+)\s*(?:kg|parts?)\s+(?:of\s+)?" + re.escape(mix.group("b")), low)
            if ratio:
                ch["mix_ratio"] = float(ratio.group(1))
            hot = re.search(rf"{re.escape(mix.group('b'))}\s+(?:at\s+)?{_NUM}{_DEG}", low)
            if hot:
                ch["mix_t_c"] = float(hot.group(1))
                for k in ("t_end_c", "temperature_c"):
                    if ch.get(k) == ch["mix_t_c"]:
                        ch.pop(k)
            first_hot = mix.group("a") and re.search(rf"{re.escape(mix.group('a'))}\s+(?:at\s+)?{_NUM}{_DEG}", low)
            if first_hot:
                ch["t_start_c"] = float(first_hot.group(1))
                for k in ("t_end_c", "temperature_c"):
                    if ch.get(k) == ch["t_start_c"]:
                        ch.pop(k)
    elif re.search(r"\b(no|stop|without)\s+mixing\b", low):
        ch["mix_with"] = ""
    pasted = analysis.parse(text, card.units, strict=False)
    if len(pasted.entries) >= 3:                       # a pasted water analysis
        ch["fluid"] = CUSTOM
        ch["analysis"] = "\n".join(e.line() for e in pasted.entries)
        if pasted.unit:
            ch["units"] = pasted.unit
        ch["ph"] = pasted.ph                           # None asks for it rather than keeping an old pH
        names = ", ".join(e.written for e in pasted.entries)
        out.reply = f"Read a custom fluid with {len(pasted.entries)} solutes ({names})."
        if pasted.ph is None:
            out.question = "What is its pH?"
        out.run = False
    elif pasted.ph is not None and not re.search(r"\b(plot|show|hide)\b", low):
        ch["ph"] = pasted.ph
    ch.update(_redox(low))
    ch.update(_gas_options(low))
    ch.update(_kinetic_options(low))
    ch.update(_brine_options(low))
    held = re.search(rf"\bat\s*{_NUM}{_DEG}", low)
    if ch.get("evaporate", 1.0) > 1.0 and held and "t_end_c" in ch and "t_start_c" not in ch \
            and float(held.group(1)) == ch["t_end_c"]:
        ch["temperature_c"] = ch.pop("t_end_c")              # "evaporate ... at 25 C" holds 25 C
    found = []
    for m in _REACTANT.finditer(text):                 # "with 10 g albite, 1 g calcite and 5 g quartz"
        found.append(m.groups())
        rest = text[m.end():]
        while (more := _MORE_REACTANT.match(rest)):
            found.append(more.groups())
            rest = rest[more.end():]
    found = [f for f in found if not re.fullmatch(_FLUID_WORDS, f[2].lower())]   # "with 10 kg of seawater" is mixing
    if found:
        ch["reactants"] = "; ".join(f"{amount} {unit} {name}" for amount, unit, name in found)
        held = re.search(rf"\bat\s*{_NUM}{_DEG}", low)
        if held and "t_end_c" in ch and "t_start_c" not in ch and float(held.group(1)) == ch["t_end_c"]:
            ch["temperature_c"] = ch.pop("t_end_c")         # "react ... at 25 C" holds 25 C
    elif re.search(r"\b(no|remove the|without|clear the)\s+reactants?\b", low):
        ch["reactants"] = []
    out.run = bool(RUN_WORDS.search(low)) and not re.search(r"\b(don't|do not|dont)\s+run\b", low)

    # plot requests: "plot A instead of B", "show A", "hide B", "remove B"
    m = re.search(r"(?:plot|show|display|use|add)\s+(.+?)\s+(?:instead of|rather than|in place of)\s+(.+)", text, re.I)
    if m:
        out.show, out.hide = [_clean(m.group(1))], [_clean(m.group(2))]
        out.run = False
    else:
        m = re.search(r"\b(?:plot|show|display|add)\s+(.+)", text, re.I)
        if m and not ch:
            out.show = [_clean(s) for s in re.split(r",|\band\b", m.group(1)) if _clean(s)]
            out.run = False
        m = re.search(r"\b(?:hide|remove|drop)\s+(.+)", text, re.I)
        if m:
            out.hide = [_clean(s) for s in re.split(r",|\band\b", m.group(1)) if _clean(s)]
    if re.search(r"\blog(arithmic)?\b", low):
        out.axis = "log"
    elif re.search(r"\blinear\b", low):
        out.axis = "linear"
    if ch:
        bits = []
        if "t_start_c" in ch and "t_end_c" in ch:
            bits.append(f"temperature {ch['t_start_c']:g} to {ch['t_end_c']:g} C")
        elif "t_end_c" in ch:
            bits.append(f"end temperature {ch['t_end_c']:g} C")
        if "temperature_c" in ch:
            bits.append(f"temperature held at {ch['temperature_c']:g} C")
        if "p_start_bar" in ch and "p_end_bar" in ch:
            bits.append(f"pressure {ch['p_start_bar']:g} to {ch['p_end_bar']:g} bar")
        elif "p_end_bar" in ch:
            bits.append(f"end pressure {ch['p_end_bar']:g} bar")
        if "pressure_bar" in ch:
            bits.append(f"pressure held at {ch['pressure_bar']:g} bar")
        if "mode" in ch:
            bits.append(ch["mode"].replace("_", "-"))
        if ch.get("mix_with"):
            bits.append(f"mixing in up to {ch.get('mix_ratio', card.mix_ratio):g} kg of "
                        f"{ch['mix_with'].replace('_', ' ')} per kg"
                        + (f" at {ch['mix_t_c']:g} C" if "mix_t_c" in ch else ""))
        if ch.get("reactants"):
            verb = "reacting with " if ch.get("kinetics", card.kinetics) else "adding "
            bits.append(verb + ch["reactants"].replace(";", ","))
        elif "reactants" in ch:
            bits.append("no reactants")
        if ch.get("kinetics") and "duration" in ch:
            bits.append(f"measured rates over {ch['duration']:g} {ch['duration_unit']}")
        elif ch.get("kinetics"):
            bits.append("measured rates")
        elif "kinetics" in ch:
            bits.append("no kinetics (reactants react completely)")
        if "evaporate" in ch:
            bits.append(f"evaporated to {ch['evaporate']:g} times" if ch["evaporate"] > 1 else "no evaporation")
        if "activity_model" in ch:
            bits.append({"pitzer": "Pitzer activity coefficients", "bdot": "B-dot activity coefficients",
                         "davies": "Davies activity coefficients"}[ch["activity_model"]])
        if "surface_area_cm2_per_g" in ch:
            bits.append(f"surface area {ch['surface_area_cm2_per_g']:g} cm2/g")
        if "ph" in ch and ch["ph"] is not None and "analysis" not in ch:
            bits.append(f"pH {ch['ph']:g}")
        if "redox" in ch:
            bits.append({"air": "redox from equilibrium with air", "none": "redox not tracked",
                         "log_fo2": f"log fO2 {ch.get('redox_value', 0):g}",
                         "eh": f"Eh {ch.get('redox_value', 0):g} V",
                         "pe": f"pe {ch.get('redox_value', 0):g}"}[ch["redox"]])
        if "analysis" not in ch:
            out.reply = "Set " + ", ".join(bits) + "." if bits else ""
        elif bits:
            out.reply += " Set " + ", ".join(bits) + "."
    return out


_TIME = re.compile(r"\b(?:for|over|during|after|lasting)\s+(?:a\s+period\s+of\s+)?(\d*\.?\d+(?:e\d+)?)\s*"
                   r"(seconds?|secs?|minutes?|mins?|hours?|hrs?|days?|weeks?|months?|years?|yrs?)\b")
_TIME_SCALE = {"week": ("days", 7.0), "month": ("years", 1.0 / 12.0)}


def _kinetic_options(low: str) -> dict:
    """Kinetics from "for 10 years", "at measured rates" or "kinetically"; off with "instantly" or
    "to equilibrium"; and the reactive surface area."""
    out: dict = {}
    m = _TIME.search(low)
    if m:
        unit = re.sub(r"s$", "", m.group(2))
        unit = {"sec": "second", "min": "minute", "hr": "hour", "yr": "year"}.get(unit, unit)
        name, k = _TIME_SCALE.get(unit, (unit + "s", 1.0))
        out["kinetics"] = True
        out["duration"] = float(m.group(1)) * k
        out["duration_unit"] = name
    if re.search(r"kinetic(s|ally)?|measured rates?|rate laws?", low):
        out["kinetics"] = not re.search(r"(no|without|turn off|switch off)\s+kinetics", low)
    elif re.search(r"(instantly|instantaneous(ly)?|all at once)", low):
        out["kinetics"] = False
    m = re.search(r"(?:surface area|area)\s*(?:of\s+)?(?:to\s+)?(\d*\.?\d+)\s*(cm2|cm\^2|m2|m\^2)\s*(?:/|per)\s*g", low)
    if m:
        out["surface_area_cm2_per_g"] = float(m.group(1)) * (1e4 if m.group(2).startswith("m") else 1.0)
    return out


def _brine_options(low: str) -> dict:
    """Pitzer activity coefficients, and evaporation ("evaporate seawater 50 times", "to 20x")."""
    out: dict = {}
    if re.search(r"\bpitzer\b", low):
        out["activity_model"] = "pitzer"
    elif re.search(r"\bb[- ]?dot\b", low):
        out["activity_model"] = "bdot"
    elif re.search(r"\bdavies\b", low):
        out["activity_model"] = "davies"
    if re.search(r"\b(evaporat|concentrat)\w*", low):
        m = re.search(r"(\d*\.?\d+)\s*(?:x\b|times|×|-?fold)", low)
        if m:
            out["evaporate"] = float(m.group(1))
        elif re.search(r"\bno evaporation\b|\bwithout evaporat", low):
            out["evaporate"] = 1.0
    return out


def _gas_options(low: str) -> dict:
    """The gas phase switch and whether methane may reduce sulfate."""
    out: dict = {}
    if re.search(r"\b(no|without|suppress\w*|disable\w*|turn off)\b[^.]*\bgas phase\b|\bgas phase\b[^.]*\b(off|suppressed)\b", low):
        out["gas_phase"] = False
    elif re.search(r"\b(allow|with|enable|turn on)\b[^.]*\bgas phase\b", low):
        out["gas_phase"] = True
    if re.search(r"methane", low) and re.search(r"sulfate|sulphate", low):
        if re.search(r"\b(apart|separate|not react|don'?t react|no reaction|block)", low):
            out["methane_sulfate"] = "apart"
        elif re.search(r"\b(react|equilibrat|reduc)", low):
            out["methane_sulfate"] = "react"
    return out


def _redox(low: str) -> dict:
    """The redox state named in a message: air, log fO2, Eh, pe, or not tracked."""
    m = re.search(r"\b(?:(?:log\s*)?f\s?o2|oxygen fugacity)\s*(?:of|=|:|is|at|to)?\s*([-+]?\d+(?:\.\d+)?)", low)
    if m:
        return {"redox": "log_fo2", "redox_value": float(m.group(1))}
    m = re.search(r"\beh\s*(?:of|=|:|is|at|to)?\s*([-+]?\d*\.?\d+)\s*(mv|v)?\b", low)
    if m:
        volts = float(m.group(1)) / (1000.0 if m.group(2) == "mv" else 1.0)
        return {"redox": "eh", "redox_value": volts}
    m = re.search(r"\bpe\s*(?:of|=|:|is|at|to)?\s*([-+]?\d*\.?\d+)\b", low)
    if m:
        return {"redox": "pe", "redox_value": float(m.group(1))}
    if re.search(r"equilibri\w* with (?:the )?(?:air|atmosphere)|air[- ]saturated|oxygenated", low):
        return {"redox": "air"}
    if re.search(r"\b(?:no|ignore|without|don'?t track)\s+(?:the\s+)?redox\b|redox\s+(?:off|not tracked)", low):
        return {"redox": "none"}
    return {}


def _clean(s: str) -> str:
    s = re.sub(r"\b(?:against|versus|vs\.?|as a function of|over)\s+(?:the\s+)?(?:pressure|temperature)\b",
               "", s, flags=re.I)
    s = re.sub(r"\b(on|using|with|in)\s+(?:a\s+|an\s+|the\s+)?(?:linear|log(?:arithmic)?)(?:\s+(?:axis|scale))?",
               "", s, flags=re.I)
    s = re.sub(r"\b(please|the|on the plot|in the plot|from the plot|on the graph)\b", "", s, flags=re.I)
    return re.sub(r"\s+", " ", s).strip(" .,:;")

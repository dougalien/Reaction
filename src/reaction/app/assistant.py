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

from . import config
from .fluids import PRESETS
from .scenario import FIELD_HELP, Scenario

SCHEMA = {
    "type": "object",
    "properties": {
        "reply": {"type": "string"},
        "scenario_changes": {
            "type": "object",
            "properties": {
                "title": {"type": "string"},
                "fluid": {"type": "string", "enum": list(PRESETS)},
                "t_start_c": {"type": "number"},
                "t_end_c": {"type": "number"},
                "steps": {"type": "integer"},
                "pressure_bar": {"type": "number"},
                "mode": {"type": "string", "enum": ["static", "flow_through", "both"]},
                "blocked_minerals": {"type": "array", "items": {"type": "string"}},
            },
            "additionalProperties": False,
        },
        "show": {"type": "array", "items": {"type": "string"}},
        "hide": {"type": "array", "items": {"type": "string"}},
        "axis": {"type": "string", "enum": ["log", "linear", "unchanged"]},
        "run": {"type": "boolean"},
        "question": {"type": "string"},
    },
    "required": ["reply", "scenario_changes", "show", "hide", "axis", "run", "question"],
}

SYSTEM = """You are the assistant inside a personal geochemical calculator. The user describes a \
model in plain language. You do two jobs: fill in the scenario card, and control what the plot shows.

Rules:
- Change only the card fields the user actually asked for or clearly implied. Never invent values.
- You never state or guess any chemical results (concentrations, pH, minerals). The calculator does that.
- Temperatures are degrees C. Pressure is in bar (1 MPa = 10 bar, 1 kbar = 1000 bar).
- Only heating runs exist so far. If the user asks for something else, say so briefly in reply.
- "show" lists species, minerals, or properties to plot exactly as the user wrote them (for example
  "CaCl+", "anhydrite", "pH", "SO4--"). "hide" lists ones to remove. "axis" is log, linear or unchanged.
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
{{"reply":"Set seawater heated to 300 C at 400 bar, static.","scenario_changes":{{"t_end_c":300,"pressure_bar":400,"mode":"static"}},"show":[],"hide":[],"axis":"unchanged","run":true,"question":""}}
"plot CaCl+ instead of Ca++" ->
{{"reply":"Plotting CaCl+ in place of Ca++.","scenario_changes":{{}},"show":["CaCl+"],"hide":["Ca++"],"axis":"unchanged","run":false,"question":""}}
"""


@dataclass
class Interpretation:
    reply: str = ""
    scenario_changes: dict = field(default_factory=dict)
    show: list[str] = field(default_factory=list)
    hide: list[str] = field(default_factory=list)
    axis: str = "unchanged"
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
    if model is None:
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
            run=bool(data.get("run", False)), question=str(data.get("question", "")).strip(),
            source="model")
    except (OSError, ValueError, KeyError, TypeError):
        rules.reply = ("The language model did not answer, so I read your message with simple rules "
                       "and may have missed details. " + rules.reply).strip()
        return rules
    # A small model sometimes invents values.  Numbers and modes are kept only if the user's own
    # words contain them, and the rules' reading of those fields wins over the model's.
    got.scenario_changes = _guard(got.scenario_changes, text)
    for key in ("t_start_c", "t_end_c", "pressure_bar", "mode", "fluid", "steps"):
        if key in rules.scenario_changes:
            got.scenario_changes[key] = rules.scenario_changes[key]
    got.mode_mentioned = rules.mode_mentioned
    plot_request = bool(rules.show or rules.hide)
    got.run = rules.run if plot_request else (got.run or rules.run)
    got.show = got.show or rules.show
    got.hide = got.hide or rules.hide
    if got.axis == "unchanged":
        got.axis = rules.axis
    return got


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
            if key in ("t_start_c", "t_end_c", "steps") and not said(value):
                continue
            if key == "pressure_bar" and not said(value, (1.0, 10.0, 1000.0, 1.01325)):
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
RUN_WORDS = re.compile(r"\b(run|go|start|simulate|calculate|compute|react|heat|warm|do it|what happens)\b", re.I)


def interpret_rules(text: str, card: Scenario) -> Interpretation:
    low = text.lower()
    out = Interpretation(source="rules")
    ch = out.scenario_changes

    m = re.search(rf"from\s+{_NUM}(?:{_DEG})?\s*(?:to|-|until)\s*{_NUM}{_DEG}", low)
    if m:
        ch["t_start_c"], ch["t_end_c"] = float(m.group(1)), float(m.group(2))
    else:
        m = re.search(rf"(?:to|up to|until|at)\s*{_NUM}{_DEG}", low) or re.search(rf"{_NUM}{_DEG}", low)
        if m:
            ch["t_end_c"] = float(m.group(1))
    m = re.search(rf"{_NUM}\s*(bars?|mpa|kbar|atm)\b", low)
    if m:
        v, unit = float(m.group(1)), m.group(2)
        ch["pressure_bar"] = v * (10.0 if unit == "mpa" else 1000.0 if unit == "kbar" else
                                  1.01325 if unit == "atm" else 1.0)
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
    for key in PRESETS:
        if key in low:
            ch["fluid"] = key
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
        if "t_end_c" in ch:
            bits.append(f"end temperature {ch['t_end_c']:g} C")
        if "pressure_bar" in ch:
            bits.append(f"pressure {ch['pressure_bar']:g} bar")
        if "mode" in ch:
            bits.append(ch["mode"].replace("_", "-"))
        out.reply = "Set " + ", ".join(bits) + "." if bits else ""
    return out


def _clean(s: str) -> str:
    s = re.sub(r"\b(on|using|with|in)\s+(?:a\s+|an\s+|the\s+)?(?:linear|log(?:arithmic)?)(?:\s+(?:axis|scale))?",
               "", s, flags=re.I)
    s = re.sub(r"\b(please|the|on the plot|in the plot|from the plot|on the graph)\b", "", s, flags=re.I)
    return re.sub(r"\s+", " ", s).strip(" .,:;")

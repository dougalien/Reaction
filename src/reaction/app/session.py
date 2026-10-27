"""The state the two windows share: the card, the chat, the running job, saved runs, and the view."""
from __future__ import annotations

import json
import queue
import re
import threading
import time
import traceback
from pathlib import Path

from .assistant import Interpretation, Ollama, change_bits, interpret
from . import analysis, carryover, diagrams
from .engine import Engine, EngineUnavailable
from .fluids import CUSTOM, CUSTOM_LABEL, PRESETS
from .narrate import story
from .scenario import EXAMPLES, STARTING_CARDS, Scenario, example_card, reactant_text, starting_card
from .view import apply_patch, default_view, resolve_item, sanitize

KEEP_CHAT = 80


# what a previous card added to its fluid, cleared when a message starts a new experiment from a named fluid
PROCESS_FIELDS = {"reactants": "reactants", "mix_with": "mixing", "evaporate": "evaporation", "kinetics": "kinetics",
                  "kind": "the column", "pore_fluid": "the column's pore fluid", "exchange_meq": "exchange sites",
                  "hfo_g": "ferric oxide surface", "kaolinite_g": "kaolinite", "illite_g": "illite",
                  "montmorillonite_g": "montmorillonite", "kd_solid_g": "partition-coefficient sediment",
                  "follow_phase": "phase following"}


GREETING = ("Tell me what to model, for example: \"heat seawater to 350 C at 500 bar\". "
            "You can also edit the card by hand.")


class Session:
    def __init__(self, workspace: Path, engine: Engine, model: Ollama | None):
        self.workspace = Path(workspace)
        (self.workspace / "runs").mkdir(parents=True, exist_ok=True)
        self.engine, self.model = engine, model
        self.lock = threading.RLock()
        self.subscribers: list[queue.Queue] = []
        self.card = Scenario()
        self.chat: list[dict] = []
        self.view: dict = {"items": [], "log": True, "modes": [], "unit": "mmol"}
        self.second: dict = {"current": None, "view": {}}      # the second results window
        self.diagram: dict | None = None                       # the diagram both windows show
        self.example: str | None = None                        # the example the card was loaded from
        self.example_pending = False                           # its first run opens on the example's plot
        self.current: str | None = None
        self.run_state = {"state": "idle", "progress": 0.0, "message": "", "error": ""}
        self.mode_confirmed = False
        self.rerun = False
        self.pending_run = False
        self.pending_x: str | None = None
        self.version = 0
        self._result_cache: dict[str, dict] = {}
        self._llm_ok = False
        self._load()
        if not self.chat:
            self.say("assistant", GREETING, log=False)

    # ---- persistence --------------------------------------------------------------------------
    @property
    def _file(self) -> Path:
        return self.workspace / "session.json"

    def _load(self) -> None:
        try:
            data = json.loads(self._file.read_text(encoding="utf-8"))
            self.card = Scenario.from_dict(data.get("card", {}))
            self.chat = data.get("chat", [])[-KEEP_CHAT:]
            self.current = data.get("current")
            self.view = data.get("view", self.view)
            self.mode_confirmed = bool(data.get("mode_confirmed", False))
            self.diagram = data.get("diagram") if isinstance(data.get("diagram"), dict) else None
            self.example = data.get("example") if data.get("example") in EXAMPLES else None
            second = data.get("second") or {}
            if isinstance(second, dict) and self.result(second.get("current")) is not None:
                self.second = {"current": second["current"],
                               "view": sanitize(second.get("view", {}), self.result(second["current"]))}
        except (OSError, ValueError):
            pass
        if self.current and not (self.workspace / "runs" / f"{self.current}.json").exists():
            self.current = None
        if self.current:
            self.view = sanitize(self.view, self.result(self.current))

    def _save(self) -> None:
        data = {"card": self.card.to_dict(), "chat": self.chat[-KEEP_CHAT:], "current": self.current,
                "view": self.view, "mode_confirmed": self.mode_confirmed, "second": self.second,
                "diagram": self.diagram, "example": self.example}
        tmp = self._file.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=1), encoding="utf-8")
        tmp.replace(self._file)

    # ---- events -------------------------------------------------------------------------------
    def subscribe(self) -> queue.Queue:
        q: queue.Queue = queue.Queue(maxsize=50)
        with self.lock:
            self.subscribers.append(q)
        q.put(self.snapshot())
        return q

    def unsubscribe(self, q: queue.Queue) -> None:
        with self.lock:
            if q in self.subscribers:
                self.subscribers.remove(q)

    def _publish(self) -> None:
        with self.lock:
            self.version += 1
            snap = self.snapshot()
            for q in list(self.subscribers):
                try:
                    q.put_nowait(snap)
                except queue.Full:
                    self.subscribers.remove(q)

    def snapshot(self) -> dict:
        with self.lock:
            advice: list[str] = []
            advisories = getattr(self.engine, "advisories", None)
            if advisories is not None and self.engine._db is not None:
                try:
                    advice = advisories(self.card)
                except Exception:                      # advice is a courtesy, never a reason to fail
                    advice = []
            spec = self.card.fluid_spec()
            return {"version": self.version, "card": self.card.to_dict(),
                    "fluids": [{"key": p.key, "label": p.label, "description": p.description,
                                "source": p.source} for p in PRESETS.values()]
                    + [{"key": CUSTOM, "label": CUSTOM_LABEL,
                        "description": "Type or paste an analysis, one solute per line."}],
                    "fluid_read": {**analysis.preview(spec.analysis, spec.units),
                                   "analysis": spec.analysis, "units": spec.units,
                                   "ph": spec.ph if spec.ph == spec.ph else None,   # NaN is not JSON
                                   "measured_t_c": spec.measured_t_c, "measured_p_bar": spec.measured_p_bar,
                                   "charge_balance": spec.charge_balance,
                                   "redox": spec.redox, "redox_value": spec.redox_value,
                                   "redox_species": _couples(analysis.parse(spec.analysis, spec.units).species)},
                    "starting_cards": [{"key": c.key, "label": c.label, "description": c.description}
                                       for c in STARTING_CARDS.values()],
                    "examples": [{"key": e.key, "label": e.label, "feature": e.feature,
                                  "description": e.description, "setup": [list(x) for x in e.setup]}
                                 for e in EXAMPLES.values()],
                    "example": self.example,
                    "card_problems": self.card_problems(), "advisories": advice,
                    "chat": self.chat[-KEEP_CHAT:],
                    "run": dict(self.run_state), "runs": self.list_runs(), "current": self.current,
                    "view": self.view, "windows": {"2": self.window("2")}, "diagram": self.diagram, "llm": {"model": getattr(self.model, "model", None),
                                               "ready": self._llm_ok},
                    "engine": self.engine.status() if self.engine._db is not None else
                    {"ready": None, "message": "Thermodynamic data load on the first run."}}

    def card_problems(self) -> list[str]:
        """The card's own checks, plus reactant names once the database is loaded to check them."""
        out = self.card.problems()
        problem = self._reactant_problem(self.card)
        if problem:
            out.append(problem)
        return out

    def _reactant_problem(self, card: Scenario) -> str:
        """Why the card's reactants cannot be used ("'xyz' is not a mineral ..."), or "".  Loads the database
        when a reactant needs checking, so a run never starts on a name it cannot find."""
        if not card.reactants or not hasattr(self.engine, "_reactant_names"):
            return ""
        try:
            self.engine.db
            self.engine._reactant_names(card)
        except EngineUnavailable:
            return ""                                   # reported when a run is attempted
        except ValueError as exc:
            return str(exc)
        return ""

    # ---- chat ---------------------------------------------------------------------------------
    def say(self, role: str, text: str, log: bool = True) -> None:
        with self.lock:
            self.chat.append({"role": role, "text": text, "t": time.strftime("%H:%M")})
            if log:
                self._save()

    def clear_chat(self) -> None:
        """Empty the chat window (and the history the language model reads); the card and runs stay."""
        with self.lock:
            self.chat = []
            self.say("assistant", GREETING)
        self._publish()

    def check_model(self) -> bool:
        ok = bool(self.model and self.model.available())
        changed = ok != self._llm_ok
        self._llm_ok = ok
        if changed:
            self._publish()
        return ok

    def preload_engine(self) -> None:
        """Load the thermodynamic data in the background so path advice is ready before the first run."""
        def load():
            try:
                self.engine.db
            except Exception:
                pass                                   # reported when a run is attempted
            self._publish()
        if getattr(self.engine, "_db", None) is None:
            threading.Thread(target=load, daemon=True).start()

    def watch_model(self, every: float = 20.0) -> None:
        """Keep the "language model ready" indicator current in a background thread."""
        def loop():
            while True:
                self.check_model()
                time.sleep(every)
        threading.Thread(target=loop, daemon=True).start()

    def handle_chat(self, text: str) -> None:
        text = text.strip()
        if not text:
            return
        with self.lock:
            self.say("user", text)
            card, history = self.card, self.chat[:-1]
            result = self.result(self.current) if self.current else None
            plotted = [i["key"] for i in self.view.get("items", [])]
        self._publish()
        if re.fullmatch(r"\s*(new|reset|start over|new run|clear)\s*[.!]?\s*", text, re.I):
            self.reset_card()
            self.say("assistant", "Card reset to the defaults.")
            self._publish()
            return
        start = _starting_card_named(text)
        if start is not None:
            self.reset_card(start)
            card = STARTING_CARDS[start]
            self.say("assistant", f"{card.label} card loaded. {card.description} Change anything on the card, "
                                  "then run it.")
            self._publish()
            return
        request = diagrams.parse_request(text)
        if request is not None:
            self.say("assistant", self.set_diagram(request, publish=False))
            self._publish()
            return
        model = self.model if self.check_model() else None
        got = interpret(text, card, plotted, history, model)
        with self.lock:
            self._apply(got, result)
        self._publish()

    def _fresh(self, got: Interpretation) -> tuple[dict, str]:
        """A message that names its starting fluid ("react the hydrothermal fluid with olivine") describes a new
        experiment: what the last card added to the fluid (mixing, reactants, evaporation, a column, sorbents)
        and the message does not mention is cleared.  Returns those changes and a sentence naming them."""
        if "fluid" not in got.scenario_changes or got.from_result:
            return {}, ""
        default, card = Scenario(), self.card
        out, said = {}, []
        for key, label in PROCESS_FIELDS.items():
            if key in got.scenario_changes or getattr(card, key) == getattr(default, key):
                continue
            if key == "mix_with":
                label = f"mixing with {card.mix_ratio:g} kg of {card.mix_with.replace('_', ' ')} per kg"
            elif key == "reactants":
                label = ", ".join(reactant_text(r) for r in card.reactants)
            out[key] = getattr(default, key)
            said.append(label)
        return out, ("Cleared from the last card: " + "; ".join(said) + ".") if said else ""

    def _apply(self, got: Interpretation, result: dict | None) -> None:
        notes: list[str] = []
        fresh, cleared = self._fresh(got)
        got.scenario_changes = {**fresh, **got.scenario_changes}
        changes, carried, blocked = self._carry(got, result)
        new, problems = self.card.changed(changes)
        notes += problems
        if changes and not problems and new.to_dict() == self.card.to_dict() and not carried:
            got.reply = f"The card already says that: {new.describe()}."
        if "reactants" in changes and new.reactants:
            bad = self._reactant_problem(new)
            if bad:                                     # keep the old reactants rather than a name that fails
                new.reactants = self.card.reactants
                had = ", ".join(reactant_text(r) for r in self.card.reactants) or "none"
                blocked = blocked or f"I can't use that: {bad}. The reactants stay as they were ({had})."
                rest = {k: v for k, v in got.scenario_changes.items() if k != "reactants"}
                got.reply = ("Set " + ", ".join(change_bits(rest, self.card)) + ".") if change_bits(rest, self.card) else " "
        if changes:
            self.card = new
        if got.mode_mentioned or "mode" in got.scenario_changes:
            self.mode_confirmed = True
        patch: dict = {}
        if result is not None and (got.show or got.hide):
            adds, removes, missing = [], [], []
            for s in got.show:
                item = resolve_item(s, result)
                (adds if item else missing).append(item or s)
            for s in got.hide:
                item = resolve_item(s, result)
                if item:
                    removes.append(item)
            patch.update(add=adds, remove=removes)
            if missing:
                notes.append("I could not find " + ", ".join(f"'{m}'" for m in missing)
                             + " in this result.")
        elif got.show or got.hide:
            notes.append("There is no result to plot yet.")
        if got.axis in ("log", "linear"):
            patch["log"] = got.axis == "log"
        if got.column_view and result is not None:
            part = result["modes"].get(got.column_view)
            if part is None:
                notes.append("That view is for a column run; this result is not one.")
            else:
                patch.update(modes=[got.column_view], x=part["varies"][0])
        if patch and result is not None:
            self.view = apply_patch(self.view, patch, result)
        reply = got.reply or ("Done." if (patch or got.scenario_changes) else "")
        problems_now = self.card_problems()
        # answering the mode question counts as the go-ahead for the run that prompted it
        want_run = got.run or (self.pending_run and (got.mode_mentioned or "mode" in got.scenario_changes))
        self.pending_run = False
        run_now = want_run and not got.question and not blocked
        if run_now and result is not None and not got.again and self.card.to_dict() == result.get("scenario"):
            run_now = False                             # the chat understood nothing new: say so, don't repeat
            notes.append("Nothing on the card changed from the last run, so I have not run it again. "
                         "Say what to change, or \"run it again\".")
        if blocked:
            notes.append(blocked)
        if run_now and problems_now:
            reply = (reply + " I can't run that yet: " + "; ".join(problems_now) + ".").strip()
            run_now = False
        elif run_now and not self.mode_confirmed and not self.card.transport:
            got.question = "Static, flow-through, or both?"
            run_now = False
            self.pending_run = True
        x_req = got.x_axis if got.x_axis in ("temperature_c", "pressure_bar") else None
        if x_req:
            name = "pressure" if x_req == "pressure_bar" else "temperature"
            if run_now:
                self.pending_x = x_req                  # applies to the result of the run about to start
            elif result is None:
                notes.append("There is no result to plot yet.")
            elif x_req in _varies(result, self.view):
                self.view = apply_patch(self.view, {"x": x_req}, result)
                reply = reply or f"Plotting against {name}."
            else:
                notes.append(f"{name.capitalize()} does not change in this run, so it cannot be the "
                             f"horizontal axis.")
        if carried:
            reply = (carried + " " + reply).strip()
        if cleared:
            reply = (reply + " " + cleared).strip()
        if got.question:
            reply = (reply + " " + got.question).strip()
        if notes:
            reply = (reply + " " + " ".join(notes)).strip()
        if run_now and self.run_state["state"] == "running":
            self.rerun = True
            reply = (reply + " A run is still going; I will start this one as soon as it finishes.").strip()
            run_now = False
        elif run_now:
            reply = (reply + f" Running: {self.card.describe()}.").strip()     # exactly what will run
        if not reply:
            reply = "I did not find anything to change. Try something like \"heat seawater to 300 C\"."
        self.say("assistant", reply)
        if run_now:
            self.start_run()

    def _carry(self, got: Interpretation, result: dict | None) -> tuple[dict, str, str]:
        """The card changes to make, with the fluid at the end of the last run as the starting fluid when the
        user asked for it ("take that fluid and cool it to 25 C"), or when they asked to cool and only the
        last run's fluid is hotter than where they want to go.  Returns (changes, what was carried, why the
        run must wait)."""
        changes = dict(got.scenario_changes)
        trial, _ = self.card.changed(changes)
        t0, t1 = trial.t_start_c, trial.t_end_c
        wrong_way = ((got.direction == "cool" and t1 > t0) or (got.direction == "heat" and t1 < t0)) \
            and "t_start_c" not in changes
        end = carryover.end_point(result) if result is not None else None
        carry = got.from_result or (wrong_way and end is not None and "fluid" not in changes and (
            end[2] > t1 if got.direction == "cool" else end[2] < t1))
        if not carry and wrong_way and "fluid" in changes:
            start = Scenario().t_start_c                # a newly named fluid starts where a fresh card does
            if (start < t1) == (got.direction == "heat"):
                return {**changes, "t_start_c": start}, f"Starting the {got.direction}ing at {start:g} °C.", ""
        if not carry:
            if wrong_way:
                verb = "cool" if got.direction == "cool" else "heat"
                return changes, "", (f"That would {'heat' if verb == 'cool' else 'cool'} the fluid from {t0:g} to "
                                     f"{t1:g} °C, not {verb} it, so I have not run it. Give the start temperature, "
                                     f"or say \"take the fluid from the last run\".")
            return changes, "", ""
        if result is None:
            return changes, "", "There is no run yet to take the fluid from, so I have not run anything."
        try:
            fluid, note = carryover.fluid_from_result(result)
        except ValueError as exc:
            return changes, "", f"I could not take the fluid from the last run: {exc}."
        for single, pair in (("temperature_c", ("t_start_c", "t_end_c")), ("pressure_bar", ("p_start_bar", "p_end_bar"))):
            if single in changes:                       # the user's held value replaces both ends
                for k in pair:
                    fluid.pop(k)
        if "p_end_bar" in changes and "p_start_bar" not in changes:
            fluid.pop("p_end_bar")
        return {**fluid, **changes}, note, ""

    # ---- diagrams ---------------------------------------------------------------------------------
    def set_diagram(self, request: dict, publish: bool = True) -> str:
        """Open, change or close the diagram; returns a sentence for the chat."""
        with self.lock:
            if request.get("close"):
                self.diagram = None
                reply = "Diagram closed."
            else:
                self.diagram = diagrams.merge(self.diagram, request)
                d = self.diagram
                what = (f"an activity diagram for {d['cation']} with aluminium and silica"
                        if d["kind"] == "activity" else
                        f"{'an' if d['kind'] == 'eh_ph' else 'a'} {diagrams.KINDS[d['kind']]} diagram for "
                        f"{d['element']}")
                where = (f"at {d['t_c']:g} C" if d.get("t_c") is not None else
                         "at the last point of the run" if self.current else "at 25 C")
                reply = (f"Drawing {what} {where} in the results window"
                         + (", with the run's path on it." if self.current else "."))
            self._save()
        if publish:
            self._publish()
        return reply

    def diagram_data(self, run_id: str | None) -> dict | None:
        if self.diagram is None:
            return None
        return diagrams.compute(self.engine.db, self.diagram, self.result(run_id) if run_id else None)

    # ---- card, view, runs -------------------------------------------------------------------------
    def edit_card(self, changes: dict) -> list[str]:
        with self.lock:
            new, problems = self.card.changed(changes)
            self.card = new
            if "mode" in changes:
                self.mode_confirmed = True
            self._save()
        self._publish()
        return problems

    def reset_card(self, start: str = "seawater") -> None:
        """Back to a fresh card: the default one, or one of the starting cards."""
        with self.lock:
            self.card = starting_card(start)
            self.mode_confirmed = False
            self.example = None
            self._save()
        self._publish()

    def load_example(self, key: str) -> bool:
        """Load one of the examples onto the card and run it.  False if a run is going (nothing changes)."""
        with self.lock:
            if self.run_state["state"] == "running":
                return False
            self.card = example_card(key)
            self.mode_confirmed = True
            self.example = key
            self.example_pending = True
            self._save()
        started = self.start_run()
        if not started:
            with self.lock:
                self.example_pending = False
        return started

    def window(self, w: str = "1") -> dict:
        """What results window ``w`` shows: {"current": run id, "view": ...}.  Window 1 is the one the
        chat drives; window 2 starts as a copy of it and goes its own way once something changes it."""
        if w == "1":
            return {"current": self.current, "view": self.view}
        if self.second.get("current") is None:
            return {"current": self.current, "view": self.view}
        return self.second

    def set_view(self, patch: dict, window: str = "1") -> None:
        with self.lock:
            win = self.window(window)
            result = self.result(win["current"]) if win["current"] else None
            view = apply_patch(win["view"], patch, result)
            if view.get("compare") and (view["compare"] == win["current"] or self.result(view["compare"]) is None):
                view["compare"] = None
            if window == "1":
                self.view = view
            else:
                self.second = {"current": win["current"], "view": view}
            self._save()
        self._publish()

    def set_current(self, run_id: str, window: str = "1") -> bool:
        with self.lock:
            result = self.result(run_id)
            if result is None:
                return False
            old = self.window(window)["view"]
            view = sanitize(old, result) if old.get("items") else default_view(result)
            if not view["items"]:
                view = {**default_view(result), "compare": view.get("compare")}
            if view.get("compare") == run_id:
                view["compare"] = None
            if window == "1":
                self.current, self.view = run_id, view
            else:
                self.second = {"current": run_id, "view": view}
            self._save()
        self._publish()
        return True

    def result(self, run_id: str | None) -> dict | None:
        if not run_id:
            return None
        with self.lock:
            if run_id in self._result_cache:
                return self._result_cache[run_id]
            path = self.workspace / "runs" / f"{run_id}.json"
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                return None
            self._result_cache[run_id] = data
            if len(self._result_cache) > 12:
                self._result_cache.pop(next(iter(self._result_cache)))
            return data

    def result_with_story(self, run_id: str) -> dict | None:
        data = self.result(run_id)
        if data is None:
            return None
        return {**data, "story": story(data)}

    def list_runs(self) -> list[dict]:
        out = []
        for path in sorted((self.workspace / "runs").glob("*.json"), reverse=True)[:40]:
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
                out.append({"id": data["id"], "created": data["created"], "summary": data["summary"],
                            "title": data["scenario"].get("title", "")})
            except (OSError, ValueError, KeyError):
                continue
        return out

    # ---- running ------------------------------------------------------------------------------
    def start_run(self) -> bool:
        with self.lock:
            if self.run_state["state"] == "running":
                return False
            problems = self.card.problems()
            if problems:
                self.run_state = {"state": "error", "progress": 0.0, "message": "",
                                  "error": "; ".join(problems)}
                self._publish()
                return False
            card = Scenario.from_dict(self.card.to_dict())
            self.run_state = {"state": "running", "progress": 0.0, "message": "Starting", "error": ""}
        self._publish()
        threading.Thread(target=self._run_thread, args=(card,), daemon=True).start()
        return True

    def _run_thread(self, card: Scenario) -> None:
        last = [0.0]

        def progress(frac: float, message: str) -> None:
            now = time.time()
            if now - last[0] < 0.25 and frac < 1.0:
                return
            last[0] = now
            with self.lock:
                self.run_state.update(progress=round(frac, 3), message=message)
            self._publish()

        try:
            with self.engine.run_lock:
                result = self.engine.run(card, progress)
            path = self.workspace / "runs" / f"{result['id']}.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            with self.lock:
                self._result_cache[result["id"]] = result
                self.current = result["id"]
                self.view = default_view(result)
                ex = EXAMPLES.get(self.example) if self.example_pending else None
                self.example_pending = False
                if ex is not None and ex.plot:          # an example opens on what it is there to show
                    items = [i for i in (resolve_item(t, result) for t in ex.plot) if i]
                    if items:
                        self.view = sanitize({**self.view, "items": items, "log": ex.log,
                                              "unit": ex.unit or self.view["unit"]}, result)
                if self.pending_x in _varies(result, self.view):
                    self.view["x"] = self.pending_x
                elif self.pending_x:
                    asked = "pressure" if self.pending_x == "pressure_bar" else "temperature"
                    self.say("assistant", f"{asked.capitalize()} does not change in this run, so the "
                                          f"horizontal axis is {result['x']['label'].lower()}.")
                self.pending_x = None
                self.say("assistant", f"Done in {result['elapsed_s']} s. The results are in the other window.")
                self._save()
                # idle only once saved: anything that waits for the run then finds it on disk too
                self.run_state = {"state": "idle", "progress": 1.0, "message": "", "error": ""}
        except EngineUnavailable as exc:
            self._fail(str(exc))
        except Exception as exc:                       # keep the app alive whatever the engine does
            traceback.print_exc()
            self._fail(f"{type(exc).__name__}: {exc}")
        self._publish()
        with self.lock:
            again, self.rerun = self.rerun, False
        if again:
            self.start_run()

    def _fail(self, message: str) -> None:
        with self.lock:
            self.run_state = {"state": "error", "progress": 0.0, "message": "", "error": message}
            self.say("assistant", f"The run failed: {message}")


_STARTING_WORDS = {"seawater": r"sea\s*water", "freshwater": r"fresh\s*water|river\s*water",
                   "brine": r"brine|(?:great\s+)?salt\s+lake", "hydrothermal": r"hydrothermal(?:\s+fluid)?",
                   "stream": r"(?:contaminated\s+)?stream|contaminated\s+site"}


def _starting_card_named(text: str) -> str | None:
    """The starting card a message asks for on its own: "brine card", "load the fresh water scenario"."""
    for key, words in _STARTING_WORDS.items():
        if re.fullmatch(rf"\s*(?:(?:load|use|open|start\s+(?:from|with)|give\s+me|new)\s+)?(?:the\s+|a\s+)?"
                        rf"(?:{words})\s+(?:scenario\s+)?(?:card|scenario)\s*[.!]?\s*", text, re.I):
            return key
    return None


def _couples(species: list[str]) -> list[str]:
    """Redox species of an analysis that can set its redox state (their partner is there too)."""
    return [s for s in species if s in analysis.REDOX
            and (analysis.REDOX[s] is None or analysis.REDOX[s] in species)]


def _varies(result: dict, view: dict) -> list[str]:
    """What can be the horizontal axis of the view: a column run's modes have axes of their own."""
    modes = [m for m in view.get("modes", []) if m in result["modes"]] or list(result["modes"])
    return result["modes"][modes[0]].get("varies", result.get("varies", []))

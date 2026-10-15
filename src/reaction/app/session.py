"""The state the two windows share: the card, the chat, the running job, saved runs, and the view."""
from __future__ import annotations

import json
import queue
import re
import threading
import time
import traceback
from pathlib import Path

from .assistant import Interpretation, Ollama, interpret
from .engine import Engine, EngineUnavailable
from .narrate import story
from .scenario import Scenario
from .view import apply_patch, default_view, resolve_item, sanitize

KEEP_CHAT = 80


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
        self.current: str | None = None
        self.run_state = {"state": "idle", "progress": 0.0, "message": "", "error": ""}
        self.mode_confirmed = False
        self.rerun = False
        self.pending_run = False
        self.version = 0
        self._result_cache: dict[str, dict] = {}
        self._llm_ok = False
        self._load()
        if not self.chat:
            self.say("assistant", "Tell me what to model, for example: \"heat seawater to 350 C at "
                                  "500 bar\". You can also edit the card by hand.", log=False)

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
        except (OSError, ValueError):
            pass
        if self.current and not (self.workspace / "runs" / f"{self.current}.json").exists():
            self.current = None
        if self.current:
            self.view = sanitize(self.view, self.result(self.current))

    def _save(self) -> None:
        data = {"card": self.card.to_dict(), "chat": self.chat[-KEEP_CHAT:], "current": self.current,
                "view": self.view, "mode_confirmed": self.mode_confirmed}
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
            return {"version": self.version, "card": self.card.to_dict(),
                    "card_problems": self.card.problems(), "chat": self.chat[-KEEP_CHAT:],
                    "run": dict(self.run_state), "runs": self.list_runs(), "current": self.current,
                    "view": self.view, "llm": {"model": getattr(self.model, "model", None),
                                               "ready": self._llm_ok},
                    "engine": self.engine.status() if self.engine._db is not None else
                    {"ready": None, "message": "Thermodynamic data load on the first run."}}

    # ---- chat ---------------------------------------------------------------------------------
    def say(self, role: str, text: str, log: bool = True) -> None:
        with self.lock:
            self.chat.append({"role": role, "text": text, "t": time.strftime("%H:%M")})
            if log:
                self._save()

    def check_model(self) -> bool:
        ok = bool(self.model and self.model.available())
        changed = ok != self._llm_ok
        self._llm_ok = ok
        if changed:
            self._publish()
        return ok

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
        model = self.model if self.check_model() else None
        got = interpret(text, card, plotted, history, model)
        with self.lock:
            self._apply(got, result)
        self._publish()

    def _apply(self, got: Interpretation, result: dict | None) -> None:
        notes: list[str] = []
        new, problems = self.card.changed(got.scenario_changes)
        notes += problems
        if got.scenario_changes:
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
        if patch and result is not None:
            self.view = apply_patch(self.view, patch, result)
        reply = got.reply or ("Done." if (patch or got.scenario_changes) else "")
        problems_now = self.card.problems()
        # answering the mode question counts as the go-ahead for the run that prompted it
        want_run = got.run or (self.pending_run and (got.mode_mentioned or "mode" in got.scenario_changes))
        self.pending_run = False
        run_now = want_run and not got.question
        if run_now and problems_now:
            reply = (reply + " I can't run that yet: " + "; ".join(problems_now) + ".").strip()
            run_now = False
        elif run_now and not self.mode_confirmed:
            got.question = "Static, flow-through, or both?"
            run_now = False
            self.pending_run = True
        if got.question:
            reply = (reply + " " + got.question).strip()
        if notes:
            reply = (reply + " " + " ".join(notes)).strip()
        if run_now and self.run_state["state"] == "running":
            self.rerun = True
            reply = (reply + " A run is still going; I will start this one as soon as it finishes.").strip()
            run_now = False
        elif run_now:
            reply = (reply + " Running now.").strip()
        if not reply:
            reply = "I did not find anything to change. Try something like \"heat seawater to 300 C\"."
        self.say("assistant", reply)
        if run_now:
            self.start_run()

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

    def reset_card(self) -> None:
        with self.lock:
            self.card = Scenario()
            self.mode_confirmed = False
            self._save()
        self._publish()

    def set_view(self, patch: dict) -> None:
        with self.lock:
            result = self.result(self.current) if self.current else None
            self.view = apply_patch(self.view, patch, result)
            self._save()
        self._publish()

    def set_current(self, run_id: str) -> bool:
        with self.lock:
            result = self.result(run_id)
            if result is None:
                return False
            self.current = run_id
            self.view = sanitize(self.view, result) if self.view.get("items") else default_view(result)
            if not self.view["items"]:
                self.view = default_view(result)
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
                self.run_state = {"state": "idle", "progress": 1.0, "message": "", "error": ""}
                self.say("assistant", f"Done in {result['elapsed_s']} s. The results are in the other window.")
                self._save()
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

"""A small local web server for the two windows (chat and card, results).

Standard library only.  The windows stay in step through a server-sent-event stream: every change
to the shared state is pushed to both.  The server listens on the loopback address only and
rejects requests whose Host or Origin is not the loopback address, so a web page you happen to
have open cannot drive it.
"""
from __future__ import annotations

import json
import mimetypes
import queue
import threading
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from . import config
from .assistant import Ollama
from .engine import Engine, EngineUnavailable
from .scenario import STARTING_CARDS
from .session import Session

STATIC = Path(__file__).resolve().parent / "static"
PAGES = {"/": "control.html", "/control": "control.html", "/results": "results.html"}
MAX_BODY = 1_000_000


class Handler(BaseHTTPRequestHandler):
    server_version = "Reaction"
    protocol_version = "HTTP/1.1"

    @property
    def session(self) -> Session:
        return self.server.session            # type: ignore[attr-defined]

    def log_message(self, fmt, *args):          # quiet
        pass

    # ---- plumbing ----------------------------------------------------------------------------
    def _allowed(self) -> bool:
        port = self.server.server_address[1]
        ok_hosts = {f"127.0.0.1:{port}", f"localhost:{port}", f"[::1]:{port}"}
        if self.headers.get("Host", "") not in ok_hosts:
            return False
        origin = self.headers.get("Origin")
        return origin is None or origin in {f"http://{h}" for h in ok_hosts}

    def _send(self, status: int, body: bytes, ctype: str, extra: dict | None = None) -> None:
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, status: int = 200) -> None:
        self._send(status, json.dumps(obj).encode(), "application/json; charset=utf-8")

    def _error(self, status: int, message: str) -> None:
        self._json({"error": message}, status)

    # ---- GET ---------------------------------------------------------------------------------
    def do_GET(self) -> None:
        if not self._allowed():
            return self._error(HTTPStatus.FORBIDDEN, "forbidden")
        path = urlparse(self.path).path
        if path in PAGES:
            return self._file(STATIC / PAGES[path])
        if path.startswith("/static/"):
            return self._file(STATIC / path[len("/static/"):])
        if path == "/api/state":
            return self._json(self.session.snapshot())
        if path.startswith("/api/run/"):
            data = self.session.result_with_story(path[len("/api/run/"):])
            return self._json(data) if data else self._error(404, "no such run")
        if path == "/api/diagram":
            run = (parse_qs(urlparse(self.path).query).get("run") or [None])[0]
            try:
                data = self.session.diagram_data(run)
            except (ValueError, KeyError) as exc:
                return self._error(400, str(exc).strip("'\""))
            except EngineUnavailable as exc:
                return self._error(503, str(exc))
            return self._json(data or {"none": True})
        if path == "/events":
            return self._events()
        if path == "/favicon.ico":
            return self._send(204, b"", "image/x-icon")
        self._error(404, "not found")

    def _file(self, path: Path) -> None:
        try:
            path = path.resolve()
            path.relative_to(STATIC.resolve())
            body = path.read_bytes()
        except (OSError, ValueError):
            return self._error(404, "not found")
        ctype = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        if ctype.startswith("text/") or ctype.endswith("javascript"):
            ctype += "; charset=utf-8"
        self._send(200, body, ctype)

    def _events(self) -> None:
        q = self.session.subscribe()
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Connection", "keep-alive")
        self.end_headers()
        try:
            while True:
                try:
                    snap = q.get(timeout=15)
                    chunk = f"data: {json.dumps(snap)}\n\n"
                except queue.Empty:
                    chunk = ": keepalive\n\n"
                self.wfile.write(chunk.encode())
                self.wfile.flush()
        except (OSError, ValueError):
            pass
        finally:
            self.session.unsubscribe(q)

    # ---- POST --------------------------------------------------------------------------------
    def do_POST(self) -> None:
        if not self._allowed():
            return self._error(HTTPStatus.FORBIDDEN, "forbidden")
        if not self.headers.get("Content-Type", "").startswith("application/json"):
            return self._error(HTTPStatus.UNSUPPORTED_MEDIA_TYPE, "send application/json")
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > MAX_BODY:
            return self._error(HTTPStatus.REQUEST_ENTITY_TOO_LARGE, "too large")
        try:
            body = json.loads(self.rfile.read(length) or b"{}")
        except ValueError:
            return self._error(400, "invalid json")
        s, path = self.session, urlparse(self.path).path
        if path == "/api/chat":
            threading.Thread(target=s.handle_chat, args=(str(body.get("text", "")),),
                             daemon=True).start()
            return self._json({"ok": True})
        if path == "/api/clear-chat":
            s.clear_chat()
            return self._json({"ok": True})
        if path == "/api/card":
            return self._json({"problems": s.edit_card(dict(body.get("changes", {})))})
        if path == "/api/reset":
            start = str(body.get("start") or "seawater")
            if start not in STARTING_CARDS:
                return self._error(400, f"unknown starting card {start!r}")
            s.reset_card(start)
            return self._json({"ok": True})
        if path == "/api/run":
            return self._json({"started": s.start_run()})
        if path == "/api/current":
            return self._json({"ok": s.set_current(str(body.get("id", "")), _window(body))})
        if path == "/api/diagram":
            request = body.get("request") or {}
            if body.get("close"):
                request = {"close": True}
            return self._json({"reply": s.set_diagram(dict(request))})
        if path == "/api/view":
            s.set_view(dict(body.get("patch", {})), _window(body))
            return self._json({"ok": True})
        self._error(404, "not found")


def _window(body: dict) -> str:
    return "2" if str(body.get("window", "1")) == "2" else "1"


class _Server(ThreadingHTTPServer):
    daemon_threads = True
    # On Windows the default lets a second server share a busy port and silently steal or lose
    # requests; refuse instead so the caller can pick another port.
    allow_reuse_address = False

    def handle_error(self, request, client_address):
        import sys
        if isinstance(sys.exc_info()[1], (ConnectionError, TimeoutError)):
            return                       # a window was closed or reloaded mid-request
        super().handle_error(request, client_address)


DEFAULT_PORT = 8790


def make_server(port: int | None = DEFAULT_PORT, workspace: Path | None = None,
                engine: Engine | None = None, model: Ollama | None = None,
                host: str = "127.0.0.1") -> ThreadingHTTPServer:
    """Bind to ``port``; if it is taken, try the next ones.  ``port=0`` asks the system for any."""
    candidates = [0] if port == 0 else range(port or DEFAULT_PORT, (port or DEFAULT_PORT) + 20)
    last: OSError | None = None
    for p in candidates:
        try:
            httpd = _Server((host, p), Handler)
            break
        except OSError as exc:
            last = exc
    else:
        raise OSError(f"no free port from {port} to {port + 19}: {last}")
    httpd.session = Session(workspace or config.workspace_dir(), engine or Engine(),   # type: ignore[attr-defined]
                            model if model is not None else Ollama())
    httpd.session.watch_model()          # type: ignore[attr-defined]
    httpd.session.preload_engine()       # type: ignore[attr-defined]
    return httpd

"""``python -m reaction.app`` starts the local server and opens the control window."""
from __future__ import annotations

import argparse
import multiprocessing
import threading
import webbrowser

from .server import DEFAULT_PORT, make_server


def main() -> None:
    multiprocessing.freeze_support()          # Reaction.exe: a spawned worker runs its task, not the app
    ap = argparse.ArgumentParser(description="Reaction: a personal geochemical calculator")
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--no-browser", action="store_true", help="do not open the windows")
    args = ap.parse_args()
    httpd = make_server(args.port)
    url = f"http://127.0.0.1:{httpd.server_address[1]}"
    print(f"Reaction is running at {url}/control  (results window: {url}/results)")
    print("Press Ctrl+C to stop.")
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(f"{url}/control")).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()


if __name__ == "__main__":
    main()

"""Where the app finds thermodynamic data and keeps its runs."""
from __future__ import annotations

import os
from pathlib import Path


def external_dir() -> Path:
    """Folder holding the third-party data the app reads but never ships.

    Expected contents: ``OBIGT/`` (the CHNOSZ OBIGT csv files) and ``EQ36_DB/data0.cmp``
    (used here only for ion-size parameters).  ``scripts/fetch_references.ps1`` fills it.
    """
    env = os.environ.get("REACTION_DATA_DIR")
    if env:
        return Path(env)
    base = os.environ.get("LOCALAPPDATA")
    if base:
        return Path(base) / "Reaction" / "external"
    return Path.home() / ".local" / "share" / "reaction" / "external"


def workspace_dir() -> Path:
    """Folder where saved runs and the chat session live (``workspace/`` in the repository)."""
    env = os.environ.get("REACTION_WORKSPACE")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[3] / "workspace"


OLLAMA_URL = os.environ.get("REACTION_OLLAMA_URL", "http://127.0.0.1:11434")
OLLAMA_MODEL = os.environ.get("REACTION_OLLAMA_MODEL", "qwen2.5:7b")

"""Where the app finds thermodynamic data and keeps its runs."""
from __future__ import annotations

import os
import sys
from pathlib import Path


BUNDLED_DIR = Path(__file__).resolve().parents[1] / "data" / "reference"


def local_dir() -> Path:
    """The per-user data folder: ``%LOCALAPPDATA%\\Reaction\\external`` (``~/.local/share/reaction/external``
    elsewhere), where ``scripts/fetch_references.ps1`` puts the full EQ3/6 and OBIGT downloads."""
    base = os.environ.get("LOCALAPPDATA")
    if base:
        return Path(base) / "Reaction" / "external"
    return Path.home() / ".local" / "share" / "reaction" / "external"


def external_dir() -> Path:
    """Folder holding the third-party thermodynamic data: ``OBIGT/`` (the CHNOSZ OBIGT csv files) and
    ``EQ36_DB/data0.cmp`` (used here only for ion-size parameters).

    ``REACTION_DATA_DIR`` wins; then a filled per-user folder (``local_dir()``, which also holds the
    EQ3/6 test library the validation tests read); otherwise the copy bundled with the app in
    ``reaction/data/reference`` (see its NOTICE.md), so the app always has its data.
    """
    env = os.environ.get("REACTION_DATA_DIR")
    if env:
        return Path(env)
    local = local_dir()
    if (local / "OBIGT" / "inorganic_aq.csv").exists() or not (BUNDLED_DIR / "OBIGT").is_dir():
        return local
    return BUNDLED_DIR


def workspace_dir() -> Path:
    """Folder where saved runs and the chat session live (``workspace/`` in the repository, or
    ``%LOCALAPPDATA%\\Reaction\\workspace`` for the packaged Reaction.exe, which has no repository)."""
    env = os.environ.get("REACTION_WORKSPACE")
    if env:
        return Path(env)
    if getattr(sys, "frozen", False):
        return local_dir().parent / "workspace"
    return Path(__file__).resolve().parents[3] / "workspace"


OLLAMA_URL = os.environ.get("REACTION_OLLAMA_URL", "http://127.0.0.1:11434")
OLLAMA_MODEL = os.environ.get("REACTION_OLLAMA_MODEL", "qwen2.5:7b")

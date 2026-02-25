"""REPL session persistence for ChemCoworker CLI."""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


SESSIONS_DIR = Path.home() / ".chemcoworker" / "sessions"
_SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9._-]+")


def _slugify(name: str) -> str:
    base = _SAFE_NAME_RE.sub("-", name.strip()).strip("-.")
    return base or "session"


def session_path(name: str) -> Path:
    """Return canonical session file path for a given user-facing name/id."""
    stem = _slugify(name)
    return SESSIONS_DIR / f"{stem}.json"


def list_sessions(limit: int = 50) -> List[Path]:
    """Return recent session files (newest first)."""
    if not SESSIONS_DIR.exists():
        return []
    files = sorted(SESSIONS_DIR.glob("*.json"), key=lambda p: p.stat().st_mtime, reverse=True)
    return files[:limit]


def save_session(
    *,
    name: str | None,
    history: List[Dict[str, Any]],
    model: str,
    provider: str,
    verbose: bool,
    plan_mode: bool,
) -> Path:
    """Persist REPL session state to JSON and return file path."""
    SESSIONS_DIR.mkdir(parents=True, exist_ok=True)
    if not name:
        ts = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
        name = f"session-{ts}"
    path = session_path(name)
    payload = {
        "id": path.stem,
        "saved_at": datetime.now(timezone.utc).isoformat(),
        "model": model,
        "provider": provider,
        "verbose": bool(verbose),
        "plan_mode": bool(plan_mode),
        "history": history,
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def load_session(name: str) -> Dict[str, Any]:
    """Load persisted REPL session JSON by id/name."""
    path = session_path(name)
    if not path.exists():
        raise FileNotFoundError(f"Session not found: {name}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or "history" not in data:
        raise ValueError("Invalid session file format")
    if not isinstance(data.get("history"), list):
        raise ValueError("Invalid session history")
    return data

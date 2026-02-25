"""CLI config persistence."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict


CONFIG_PATH = Path.home() / ".chemcoworker" / "config.json"
DEFAULT_MODEL: Dict[str, str] = {"name": "o4-mini", "provider": "openai"}


def load_config() -> Dict[str, str]:
    """Load saved model config or return defaults."""
    try:
        data = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        if "name" in data and "provider" in data:
            return {"name": str(data["name"]), "provider": str(data["provider"])}
    except Exception:
        pass
    return DEFAULT_MODEL.copy()


def save_config(model: str, provider: str) -> None:
    """Persist model choice."""
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CONFIG_PATH.write_text(
        json.dumps({"name": model, "provider": provider}, indent=2),
        encoding="utf-8",
    )

"""Loading and validation helpers for reaction grammar definitions."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List


_PATH = Path(__file__).with_name("definitions") / "reaction_grammars.v1.json"


@lru_cache(maxsize=1)
def load_reaction_grammars() -> List[Dict[str, Any]]:
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return sorted(payload.get("grammars") or [], key=lambda item: (-int(item.get("priority", 0)), str(item.get("id", ""))))


__all__ = ["load_reaction_grammars"]

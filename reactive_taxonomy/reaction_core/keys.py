"""Deterministic identity helpers for reaction-core projections.

Display labels must never participate in these identities.  Payloads supplied
to :func:`digest` are normalized chemistry fields assembled by the builder.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any


def canonical_json(value: Any) -> str:
    """Serialize an identity payload deterministically."""
    return json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    )


def digest(prefix: str, value: Any, *, length: int = 32) -> str:
    """Return a namespaced, truncated SHA-256 identity for a payload."""
    encoded = canonical_json(value).encode("utf-8")
    return f"{prefix}:" + hashlib.sha256(encoded).hexdigest()[:length]


__all__ = ["canonical_json", "digest"]

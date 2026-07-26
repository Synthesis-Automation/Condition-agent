"""Shared interpretable feature extraction from serialized reaction signatures."""

from __future__ import annotations

from typing import Any, Mapping, Tuple


def environment_tokens(signature: Mapping[str, Any]) -> Tuple[str, ...]:
    """Return stable steric, electronic, and nearby-group tokens."""
    tokens = []
    for partner in signature.get("partners") or ():
        for category in ("steric", "electronic"):
            value = (partner.get(category) or {}).get("class")
            if value:
                tokens.append(f"{category}:{value}")
        for group in partner.get("nearby_groups") or ():
            group_id = group.get("group_id")
            if group_id:
                tokens.append(f"nearby:{group_id}")
    return tuple(sorted(tokens))


__all__ = ["environment_tokens"]

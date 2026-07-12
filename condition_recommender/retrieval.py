"""Hierarchical chemistry retrieval and transparent structured similarity."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from .indexing import IndexedReaction

_RULES_PATH = Path(__file__).with_name("definitions") / "retrieval.v1.json"


@lru_cache(maxsize=1)
def load_retrieval_rules() -> Dict[str, Any]:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        return dict(json.load(handle))


def _value(partner: Dict[str, Any], key: str) -> str:
    if key == "steric":
        return str((partner.get("steric") or {}).get("class") or "")
    if key == "electronic":
        return str((partner.get("electronic") or {}).get("class") or "")
    return str(partner.get(key) or "")


def _jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    a, b = set(left), set(right)
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b) if a | b else 1.0


def structured_similarity(query: Dict[str, Any], row: IndexedReaction) -> Tuple[float, Dict[str, float]]:
    weights = load_retrieval_rules()["similarity_weights"]
    qe, qt = query["electrophile"], query["transfer_partner"]
    pairs = {
        "electrophile_handle": (_value(qe, "handle_token"), _value(row.electrophile, "handle_token")),
        "electrophile_context": (_value(qe, "anchor_context"), _value(row.electrophile, "anchor_context")),
        "transfer_handle": (_value(qt, "handle_token"), _value(row.transfer_partner, "handle_token")),
        "transfer_context": (_value(qt, "anchor_context"), _value(row.transfer_partner, "anchor_context")),
        "electrophile_steric": (_value(qe, "steric"), _value(row.electrophile, "steric")),
        "transfer_steric": (_value(qt, "steric"), _value(row.transfer_partner, "steric")),
        "electrophile_electronic": (_value(qe, "electronic"), _value(row.electrophile, "electronic")),
        "transfer_electronic": (_value(qt, "electronic"), _value(row.transfer_partner, "electronic")),
    }
    components = {name: float(left == right) for name, (left, right) in pairs.items()}
    components["spectators"] = _jaccard(query.get("spectator_group_ids", ()), row.spectator_group_ids)
    components["flags"] = _jaccard(query.get("family_flags", ()), row.family_flags)
    score = sum(float(weights[name]) * components[name] for name in weights)
    return round(score, 6), components


def retrieve_pool(query: Dict[str, Any], rows: Tuple[IndexedReaction, ...]) -> Tuple[str, Tuple[IndexedReaction, ...]]:
    """Choose the narrowest chemistry pool meeting configured support."""
    rules = load_retrieval_rules()
    minimum = int(rules.get("minimum_pool_size", 12))
    qe, qt = query["electrophile"], query["transfer_partner"]
    exact = tuple(row for row in rows if
        _value(row.electrophile, "handle_token") == _value(qe, "handle_token") and
        _value(row.electrophile, "anchor_context") == _value(qe, "anchor_context") and
        _value(row.transfer_partner, "handle_token") == _value(qt, "handle_token") and
        _value(row.transfer_partner, "anchor_context") == _value(qt, "anchor_context"))
    if len(exact) >= minimum:
        return "exact_partner_signature", exact
    compatible_handles = set((rules.get("transfer_handle_compatibility") or {}).get(_value(qt, "handle_token"), [_value(qt, "handle_token")]))
    context_relaxed = tuple(row for row in rows if
        _value(row.electrophile, "handle_token") == _value(qe, "handle_token") and
        _value(row.electrophile, "anchor_context") == _value(qe, "anchor_context") and
        _value(row.transfer_partner, "handle_token") in compatible_handles)
    if len(context_relaxed) >= minimum:
        return "transfer_context_relaxed", context_relaxed
    handle_family = tuple(row for row in rows if
        _value(row.electrophile, "handle_token") == _value(qe, "handle_token") and
        _value(row.transfer_partner, "handle_token") in compatible_handles)
    if len(handle_family) >= minimum:
        return "context_relaxed", handle_family
    return "family_prior", rows


__all__ = ["load_retrieval_rules", "retrieve_pool", "structured_similarity"]

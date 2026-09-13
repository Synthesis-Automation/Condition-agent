"""Ordinal precedent distance, separate from recipe scores and probabilities."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from .conversion.identities import canonical_reaction_identity
from .generic_indexing import GenericIndexedReaction


@lru_cache(maxsize=1)
def load_match_level_rules() -> dict[str, Any]:
    """Load the ordinal evidence contract independently of ranking weights."""
    path = Path(__file__).with_name("definitions") / "precedent_match_levels.v1.json"
    rules = json.loads(path.read_text(encoding="utf-8"))
    names = ("same_reaction", "close_precedent", "related_handle", "broader_analogue")
    levels = rules.get("levels") or {}
    if (rules.get("definition_id") != "precedent_match_levels.v1"
            or rules.get("schema_version") != "1.0"
            or rules.get("is_probability") is not False
            or set(levels) != set(names)
            or any(levels[name].get("level") != i or not levels[name].get("label")
                   for i, name in enumerate(names, 1))
            or set(rules.get("close_retrieval_levels") or ()) != {
                "reaction_facet_exact", "reaction_facet_attachment_relaxed",
                "exact_signature", "handle_signature",
            }):
        raise ValueError("invalid precedent match-level definition")
    return rules


def precedent_match_level(
    query_identity: str, row: GenericIndexedReaction, retrieval_level: str,
) -> tuple[int, str, tuple[str, ...]]:
    """Describe the strongest demonstrated relationship to an observed row."""
    rules = load_match_level_rules()

    def result(name: str, detail: str) -> tuple[int, str, tuple[str, ...]]:
        level = rules["levels"][name]
        return level["level"], level["label"], (detail,)

    identity = row.canonical_reaction_id
    if not identity.startswith("CRX1:") and query_identity:
        canonical = canonical_reaction_identity(row.reaction_smiles)
        identity = canonical.reaction_id if canonical else ""
    if query_identity and identity == query_identity:
        return result("same_reaction", "Same normalized reactants and products")
    if retrieval_level.startswith("related_handle_"):
        source, target = retrieval_level.removeprefix("related_handle_").split("_to_")
        return result("related_handle", f"Query Ar-{source}; precedent Ar-{target}")
    if retrieval_level in rules["close_retrieval_levels"]:
        return result("close_precedent", "Matched reaction features; molecules may differ")
    return result("broader_analogue", "Broader structural evidence; inspect the reaction differences")

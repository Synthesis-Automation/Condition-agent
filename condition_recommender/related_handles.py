"""Conservative query-time lookup of explicitly permitted handle analogues."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal, Mapping

from reactive_taxonomy import departing_fragment_tokens, featurize_reaction
from reactive_taxonomy.handle_variants import aromatic_leaving_handle_variant

from .generic_indexing import GenericReactionIndex
from .reaction_facets import reaction_facet_keys

SearchScope = Literal["same_handle", "automatic", "broad"]


def validate_search_scope(value: str) -> None:
    """Validate search breadth independently of ranking preferences."""
    if value not in {"same_handle", "automatic", "broad"}:
        raise ValueError(f"unknown condition search scope: {value}")


@lru_cache(maxsize=1)
def load_related_handle_rules() -> dict[str, Any]:
    """Validate the bounded, versioned search hypothesis definition."""
    path = Path(__file__).with_name("definitions") / "related_handle_retrieval.v1.json"
    rules = json.loads(path.read_text(encoding="utf-8"))
    if (rules.get("schema_version") != "1.0"
            or rules.get("definition_id") != "related_handle_retrieval.v1"
            or not rules.get("calibration_status")
            or rules.get("directions") != [["I", "Br"], ["Br", "I"]]
            or rules.get("formed_partner_elements") != ["C", "N", "O", "S"]
            or rules.get("allowed_lookup_levels") != ["reaction_facet_exact"]
            or rules.get("require_known_departing_fragments") is not True
            or rules.get("require_equal_departing_fragment_multiset") is not True
            or rules.get("minimum_independent_support") != 2
            or rules.get("match_level") != 3
            or rules.get("match_label") != "Related handle"):
        raise ValueError("invalid related-handle retrieval definition")
    return rules


def related_handle_positions(
    reaction_smiles: str,
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
) -> tuple[tuple[str, set[int]], ...]:
    """Find analogue rows through existing indices without changing the query.

    A virtual Br/I substitution is used only to obtain existing exact facet
    keys. Every departing fragment must match that variant, including the
    partner's leaving fragment. No whole-index scan or broadened chemistry
    gate is used, and no hypothetical signature is reported as observed.
    """
    rules = load_related_handle_rules()
    results = []
    for source, target in rules["directions"]:
        variant = aromatic_leaving_handle_variant(
            reaction_smiles, signature, query_element=source,
            precedent_element=target,
            formed_partner_elements=tuple(rules["formed_partner_elements"]),
        )
        if variant is None:
            continue
        analysis = featurize_reaction(variant.reaction_smiles)
        if not analysis.valid or analysis.reaction_signature is None or analysis.reaction_core is None:
            continue
        variant_signature = asdict(analysis.reaction_signature)
        keys = reaction_facet_keys(
            variant_signature, asdict(analysis.reaction_core),
            asdict(analysis.fallback_descriptor) if analysis.fallback_descriptor else None,
        )
        departing = Counter(departing_fragment_tokens(variant.reaction_smiles, variant_signature))
        if not departing:
            continue
        positions = index.facet_exact.get(keys.get("reaction_facet_exact", ""), ())
        matching = {
            position for position, row in zip(positions, index.select(positions))
            if Counter(departing_fragment_tokens(row.reaction_smiles, row.signature)) == departing
        }
        results.append((f"related_handle_{source}_to_{target}", matching))
    return tuple(results)

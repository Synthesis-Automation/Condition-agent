"""Inspect selected fragment evidence without requiring executable operators."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
from typing import Any

from reactive_taxonomy.fragment_search import compile_fragment_query, validate_fragment_target
from reactive_taxonomy.molecule_comparison import compare_molecules


def inspect_fragment_precedent(
    search: dict[str, Any], observation_id: str, target_smiles: str,
) -> dict[str, Any]:
    """Preserve one actual returned observation, scoped source text and differences.

Callers supply a trusted saved or freshly executed search, never client-authored
source records. Partial search evidence remains inspectable, but cannot authorize
transfer. No new retrieval, condition inference or atom mapping takes place.
"""
    description = search["query"]
    compiled = compile_fragment_query(description["expression"], description["query_format"], description["topology"])
    validation = validate_fragment_target(compiled, target_smiles)
    supplied = search.get("target_validation") or {}
    if (not validation.matches_target or supplied.get("matches_target") is not True
            or supplied.get("query_id") != compiled.query_id
            or supplied.get("target_smiles") != validation.target_smiles):
        raise ValueError("Saved query evidence does not validate this target")
    hits = [hit for hit in search["hits"] if hit["observation_id"] == observation_id]
    if len(hits) != 1:
        raise ValueError("Select one observation from the returned search sample")
    hit = deepcopy(hits[0])
    comparison: dict[str, Any] = {"status": "unavailable", "reason": "no_product_structure"}
    if hit.get("product_smiles"):
        try:
            comparison = compare_molecules(hit["product_smiles"], validation.target_smiles, timeout_seconds=2).to_dict()
        except ValueError as exc:
            comparison = {"status": "unavailable", "reason": str(exc)}
    return {
        "schema_version": "fragment_precedent_inspection.v1",
        "target_smiles": validation.target_smiles, "target_validation": asdict(validation),
        "query": deepcopy(description), "source": hit, "comparison": comparison,
        "search_scope": {key: deepcopy(search.get(key)) for key in (
            "index_id", "search_status", "stop_reason", "source_scope", "source_coverage_complete",
            "ranking_scope", "output_truncated", "counts", "relationship_groups")},
        "limitations": ["Observed conditions and procedures describe the source substrate, not the target.",
                        "Structural alignments are hypotheses, not observed reaction atom correspondence.",
                        "A source remains inspectable when operator compilation or transfer fails."],
    }

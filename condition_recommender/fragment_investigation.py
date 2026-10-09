"""Inspect selected fragment evidence without requiring executable operators."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
from typing import Any

from reactive_taxonomy.fragment_search import compile_fragment_query, validate_fragment_target
from reactive_taxonomy.molecule_comparison import compare_molecules


def selected_fragment_search(
    search: dict[str, Any], observation_id: str, target_smiles: str,
) -> dict[str, Any]:
    """Project a discovery hit through its recorded query without new retrieval.

    Automatic discovery merges several queries. Only the selected hit's original
    query and attempt can authorize inspection; partial scope is never upgraded.
    Ordinary fragment search results retain their existing validation path.
    """
    if search.get("schema_version") != "synthesis_precedent_discovery.v1":
        return search
    hits = [hit for hit in search.get("hits", []) if hit.get("observation_id") == observation_id]
    if len(hits) != 1:
        raise ValueError("Select one observation from the returned discovery sample")
    discovery = hits[0].get("discovery", {})
    query = compile_fragment_query(discovery.get("query"), discovery.get("query_format"),
                                   discovery.get("topology"))
    validation = validate_fragment_target(query, target_smiles)
    if not validation.matches_target or validation.target_smiles != search.get("target_smiles"):
        raise ValueError("Saved discovery query does not validate this target")
    attempts = [attempt for attempt in search.get("attempts", [])
                if attempt.get("query_id") == query.query_id
                and all(attempt.get(key) == discovery.get(key) for key in ("core_id", "level", "query", "topology"))
                and observation_id in attempt.get("returned_observation_ids", [])]
    if len(attempts) != 1:
        raise ValueError("Discovery hit requires one matching recorded query attempt")
    attempt = attempts[0]
    complete = search.get("search_status") == attempt.get("search_status") == "complete"
    return {
        **deepcopy(search), "query": query.describe(), "target_validation": asdict(validation),
        "hits": deepcopy(hits), "counts": deepcopy(attempt["counts"]),
        "search_status": "complete" if complete else "partial",
        "stop_reason": (search.get("stop_reason") if search.get("search_status") != "complete"
                        else attempt.get("stop_reason")),
    }


def inspect_fragment_precedent(
    search: dict[str, Any], observation_id: str, target_smiles: str,
) -> dict[str, Any]:
    """Preserve one actual returned observation, scoped source text and differences.

Callers supply a trusted saved or freshly executed search, never client-authored
source records. Partial search evidence remains inspectable, but cannot authorize
transfer. No new retrieval, condition inference or atom mapping takes place.
"""
    search = selected_fragment_search(search, observation_id, target_smiles)
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

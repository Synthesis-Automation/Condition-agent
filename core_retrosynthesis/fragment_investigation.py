"""Bounded selected-source investigation using canonical transfer gates."""

from __future__ import annotations

from typing import Any

from condition_recommender.fragment_investigation import inspect_fragment_precedent

from .fragment_guidance import DEFAULT_POLICY, evaluate_transfers, load_policy, project_selected_query_guidance
from .generic_library import build_generic_library


def investigate_fragment_precedent(
    search: dict[str, Any], observation_id: str, target_smiles: str,
) -> dict[str, Any]:
    """Inspect one source and test only its admitted transformations in memory.

Production libraries and admission policies are unchanged. Source compilation
rejections, incomplete evidence and ambiguous alignments retain the inspection.
"""
    inspection = inspect_fragment_precedent(search, observation_id, target_smiles)
    policy = load_policy(DEFAULT_POLICY)
    policy.update(max_focus_bonds=3, top_k=3)
    if search.get("search_side", "product") != "product":
        return {**inspection, "transfer": {"status": "product_search_required"}}
    if search.get("search_status") != "complete":
        return {**inspection, "transfer": {"status": "incomplete_search", "stop_reason": search.get("stop_reason")}}
    projected = project_selected_query_guidance(target_smiles, search, [observation_id], policy)
    guidance = projected["guidance"]
    # Reuse the canonical bounded experiment with no corpus-library search.
    transfers = evaluate_transfers(
        {"guidance": guidance, "policy": policy, "include_baseline": False},
        library=build_generic_library([], admission_mode="data_driven"),
        repeat=False, source_library_path=None,
    )
    arms = [{"target_atom_ids": arm["target_atom_ids"], **arm["direct_source_transfer"]}
            for arm in transfers["guided"]]
    if not guidance["focus_bonds"]:
        status = "no_resolved_construction"
    elif not transfers["compiled_source_template_count"]:
        status = "source_compilation_rejected"
    elif any(arm["candidates"] for arm in arms):
        status = "graph_validated_proposals"
    else:
        status = "no_verified_transfer_within_budget"
    return {
        **inspection, "transfer": {
            "status": status, "core_admission_policy": "pass_only", "policy": policy,
            "compiled_source_template_count": transfers["compiled_source_template_count"],
            "source_admissions": transfers["source_admissions"], "arms": arms,
            "guidance": guidance,
            "limitations": ["Only this selected source is tested; no production library is loaded or modified.",
                            "Graph validation is not evidence of experimental feasibility.",
                            "An empty bounded search does not establish chemical impossibility."],
        },
    }

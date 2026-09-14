"""Bounded application composition for condition-advisory graph planning."""

from __future__ import annotations

from typing import Any, Callable

from condition_recommender.reaction_context import (
    classify_planned_reaction,
    context_core,
    load_reaction_context_policy,
    reaction_sides,
)
from core_retrosynthesis import GenericTemplateLibrary, disconnect_strategies_detailed
from forward_synthesis import ForwardOperatorLibrary, assess_proposed_step


def explore_reaction_context(
    reaction: str,
    *,
    forward_library: Callable[[], ForwardOperatorLibrary],
    retro_library: Callable[[], GenericTemplateLibrary],
    recommend: Callable[..., Any],
) -> dict[str, Any]:
    """Run independent bounded channels; retain partial failures and provenance.

    No recipe or operator agreement is promoted into evidence for the original
    request. Alternative conditions are queried only for their own proposals.
    """
    inputs, product = reaction_sides(reaction)
    policy = load_reaction_context_policy()
    result: dict[str, Any] = {
        "schema_version": "1.0",
        "definition_id": policy["definition_id"],
        "query_reaction_smiles": reaction,
        "advisory_only": True,
        "search_limits": policy,
        "warnings": [
            "Generated structures are hypotheses, not observed reactions.",
            "Search is bounded; missing products or routes are inconclusive.",
            "Operator agreement is structural consistency, not independent evidence "
            "that conditions will work. No recipe-specific selectivity was assessed.",
        ],
    }
    if "." in product:
        unavailable = {"status": "out_of_scope", "error": "SINGLE_PRODUCT_REQUIRED"}
        return {
            **result,
            "reactant_analysis": unavailable,
            "product_analysis": {**unavailable, "alternatives": []},
        }
    try:
        forward = assess_proposed_step(
            inputs,
            product,
            forward_library(),
            top_k=policy["forward_products"],
            max_operators_to_apply=policy["forward_operators"],
            max_assignments_per_operator=policy["forward_assignments"],
            max_outcomes_per_operator=policy["forward_outcomes"],
        )
        result["reactant_analysis"] = {"status": "complete", **forward.to_dict()}
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        result["reactant_analysis"] = {"status": "unavailable", "error": str(exc)}
    try:
        retro = disconnect_strategies_detailed(
            product,
            retro_library(),
            top_k_strategies=policy["retro_strategies"],
            max_realizations_per_strategy=policy["retro_realizations"],
            max_templates_to_apply=policy["retro_templates_per_tier"],
            max_candidates_to_validate=policy["retro_validations_per_tier"],
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        result["product_analysis"] = {
            "status": "unavailable",
            "error": str(exc),
            "alternatives": [],
        }
        return result
    query_core = context_core(reaction)
    alternatives = []
    # Representatives first: one prolific strategy cannot consume every slot.
    candidates = [s.representative for s in retro.strategies]
    candidates += [c for s in retro.strategies for c in s.alternate_realizations]
    grouped: dict[tuple[str, str], list[Any]] = {}
    excluded_unverified = 0
    for candidate in candidates:
        if candidate.forward_validation_status != "verified_signature":
            excluded_unverified += 1
            continue
        proposal = (
            candidate.condition_query_reaction_smiles
            or candidate.proposed_reaction_smiles
        )
        identity = reaction_sides(proposal)
        grouped.setdefault(identity, []).append(candidate)
    for group in list(grouped.values())[: policy["alternatives"]]:
        candidate = group[0]
        proposal = (
            candidate.condition_query_reaction_smiles
            or candidate.proposed_reaction_smiles
        )
        relation = classify_planned_reaction(reaction, proposal, query_core=query_core)
        item: dict[str, Any] = {
            "reaction_smiles": proposal,
            "relation": relation.to_dict(),
            "operator_id": candidate.operator_id,
            "strategy_id": candidate.strategy_id,
            "precedent_reaction_ids": sorted(
                {rid for c in group for rid in c.precedent_reaction_ids}
            ),
            "supporting_operator_ids": sorted({c.operator_id for c in group}),
            "selectivity_warnings": candidate.to_dict()["selectivity_warnings"],
            "forward_validation_status": candidate.forward_validation_status,
            "evidence_kind": "generated_hypothesis",
            "conditions": None,
        }
        if (
            candidate.forward_validation_status == "verified_signature"
            and relation.kind in {"precursor_alternative", "route_alternative"}
        ):
            try:
                conditions = recommend(
                    proposal,
                    top_k=policy["recipes"],
                    preferred_reaction_ids=tuple(item["precedent_reaction_ids"]),
                )
                item["conditions"] = conditions.to_dict()
            except (FileNotFoundError, RuntimeError, ValueError) as exc:
                item["condition_error"] = str(exc)
        alternatives.append(item)
    result["product_analysis"] = {
        "status": "complete",
        "alternatives": alternatives,
        "search_diagnostics": retro.to_dict()["search_diagnostics"],
        "unique_proposals_considered": len(grouped),
        "excluded_unverified_count": excluded_unverified,
        "display_limit_reached": len(grouped) > len(alternatives),
    }
    return result

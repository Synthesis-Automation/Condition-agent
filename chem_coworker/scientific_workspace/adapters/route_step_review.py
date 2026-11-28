"""Compose bounded step evidence from the pinned corpus and existing structural checks."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from sqlite3 import DatabaseError
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from condition_recommender.recipe_assessment import ReactionRecipeAssessment
    from core_retrosynthesis.external_proposal_assessment import ExternalRetrosynthesisAssessment, ExternalRetrosynthesisProposal
    from .operations import ScientificOperations


def compact_precedent(match: Mapping[str, Any]) -> dict[str, Any]:
    """Keep the source identity, measurement and substrate differences for one observation."""
    row, context = match["precedent"], match["context"]
    result = {key: row[key] for key in ("reference_id", "reaction_id", "observation_id", "reaction_smiles", "level")}
    result.update(source_yield_pct=row["yield_pct"], outcome=match["outcome"],
                  contributes_to_support=match["contributes_to_support"],
                  differences=row["differences"][:3], difference_count=len(row["differences"]))
    result["context"] = {key: context[key] for key in ("status", "reverse_status", "minimum_partner_similarity")}
    for key in ("missing", "unresolved", "source_only", "source_unresolved"):
        result["context"][key] = context[key][:3]
        if len(context[key]) > 3:
            result["context"][f"{key}_count"] = len(context[key])
    return result


def _corpus(operations: ScientificOperations) -> tuple[Any, dict[str, Any]]:
    artifacts = operations.store.manifest["baseline"].get("artifacts", {})
    missing = [name for name in ("condition_index", "shared_core_index")
               if artifacts.get(name, {}).get("status") != "present"]
    if missing:
        return None, {"status": "not_configured" if len(missing) == 2 else "incomplete", "missing_indexes": missing}
    try:
        return operations._conditions(), {}
    except (FileNotFoundError, ValueError, DatabaseError) as exc:
        return None, {"status": "unavailable", "reason": str(exc)}


def _precedent(engine: Any, reaction: str) -> dict[str, Any]:
    from condition_recommender import assess_step_plausibility

    try:
        evidence = assess_step_plausibility(reaction, engine.index, engine.shared_core_index,
                                            candidate_limit=128, match_limit=1)
    except (FileNotFoundError, ValueError, DatabaseError) as exc:
        return {"status": "unavailable", "reason": str(exc)}
    return {
        "status": "assessed", "grade": evidence.grade, "evidence_rank": evidence.evidence_rank,
        "measured_evidence_status": evidence.measured_evidence_status,
        "candidate_count": evidence.candidate_count, "candidate_limit": evidence.candidate_limit,
        "candidate_truncated": evidence.candidate_truncated, "qualified_count": evidence.qualified_count,
        "support_counts": asdict(evidence.support_counts), "warnings": list(evidence.warnings),
        "definition_id": evidence.policy.definition_id, "definition_version": evidence.policy.definition_version,
        "closest_match": compact_precedent(asdict(evidence.closest_match)) if evidence.closest_match else None,
    }


def review_steps(
    operations: ScientificOperations,
    steps: Sequence[tuple[str | None, ExternalRetrosynthesisProposal, ExternalRetrosynthesisAssessment, ReactionRecipeAssessment | None]],
) -> list[dict[str, Any]]:
    """Audit each physical step and query each distinct validated mapping once in this call."""
    corpus = None
    corpus_by_reaction: dict[str, dict[str, Any]] = {}
    reviews = []
    for step_id, proposal, assessment, recipe in steps:
        gates = {g.gate_id: g for g in assessment.gates}
        missing = [name for name in ("input_structure", "reaction_side_consistency", "atom_correspondence",
                                     "reaction_completeness") if gates[name].status != "pass"]
        mapped = assessment.normalized_mapped_reaction_smiles
        if missing or not mapped:
            precedent = {"status": "unresolved_structure", "unverified_gates": missing}
        else:
            # Use the same chemistry input as an explicit saved-step plausibility call.
            # Generated maps can omit departing atoms and change the observed-core comparison.
            reaction = mapped if proposal.mapped_reaction_smiles else f"{proposal.precursor_smiles}>>{proposal.target_smiles}"
            if corpus is None:
                corpus = _corpus(operations)
            engine, unavailable = corpus
            if reaction not in corpus_by_reaction:
                corpus_by_reaction[reaction] = unavailable if engine is None else _precedent(engine, reaction)
            precedent = corpus_by_reaction[reaction]
            if precedent["status"] == "assessed":
                precedent = {**precedent, "query_reaction_smiles": reaction,
                             "query_basis": "validated_supplied_mapping" if proposal.mapped_reaction_smiles else "proposal_graphs"}
        checks = {name: {"status": gates[name].status, "summary": gates[name].summary}
                  for name in ("precursor_compatibility", "reaction_compatibility", "functional_group_selectivity")}
        compatibility = assessment.compatibility_evidence
        alerts = []
        if compatibility:
            for group in (compatibility.precursor, compatibility.reaction):
                alerts.extend({"rule_id": issue.rule_id, "message": issue.message,
                               "severity": issue.intrinsic_severity} for issue in group.assessments)
            if competition := compatibility.functional_group_competition:
                alerts.append({"code": competition.code, "message": competition.message})
            if geometry := compatibility.migration_geometry:
                checks["migration_geometry"] = asdict(geometry)
                if geometry.requires_caution:
                    alerts.append({"code": geometry.audit_id, "message": geometry.summary})
        reviews.append({
            "step_id": step_id, "structural_status": assessment.status,
            "compatibility": {"checks": checks, "alerts": alerts},
            "recipe_status": recipe.status if recipe is not None else "not_supplied",
            "recipe_conflicts": list(recipe.hard_conflicts) if recipe is not None else [],
            "precedent": precedent,
        })
    return reviews


def weakest_steps(reviews: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Identify blocked steps or ties for least documented support and expose unevaluated steps."""
    blocked = [r["step_id"] for r in reviews if r["structural_status"] in ("invalid", "rejected_by_compatibility")
               or r["recipe_conflicts"]]
    unresolved = [r["step_id"] for r in reviews if r["structural_status"] == "ambiguous"]
    missing = [r["step_id"] for r in reviews if r["precedent"].get("evidence_rank") is None]
    alerts = [r["step_id"] for r in reviews if r["compatibility"]["alerts"]]
    graded = [r for r in reviews if r["precedent"].get("evidence_rank") is not None]
    duplicates = [key for key, count in Counter(r["step_id"] for r in reviews).items() if count > 1]
    if duplicates:
        basis, ids = "ambiguous_step_ids", []
    elif blocked:
        basis, ids = "structural_or_recipe_conflict", blocked
    elif unresolved:
        basis, ids = "unresolved_structure", unresolved
    elif graded:
        rank = min(r["precedent"]["evidence_rank"] for r in graded)
        basis = "lowest_assessed_precedent_grade"
        ids = [r["step_id"] for r in graded if r["precedent"]["evidence_rank"] == rank]
    else:
        basis, ids = "insufficient_evidence", []
    return {"step_ids": ids, "basis": basis, "ambiguous_step_ids": duplicates, "unassessed_precedent_step_ids": missing,
            "compatibility_alert_step_ids": alerts}

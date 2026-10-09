"""Recorded composition for concrete retro validity and saved forward evidence."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from core_retrosynthesis.chemistry import canonical_smiles
from core_retrosynthesis.retro_validity import RetroForwardEvidence, assess_retro_validity

from .route_investigation import _evidence, _step, _strings
from .step_selection import SOURCE_OPERATIONS, _call, _selection, resolve_saved_candidate

if TYPE_CHECKING:
    from .operations import ScientificOperations


def assess_validity(
    operations: ScientificOperations, proposal: dict[str, Any] | None,
    source_ref: str | None, step_id: str | None, realization_id: str | None,
    forward_ref: str | None, candidate_limit: int, match_limit: int,
    evidence_refs: list[str] | None,
    strategy_id: str | None = None, precursor_smiles: str | None = None,
) -> dict[str, Any]:
    """Resolve saved inputs and pinned artifacts; delegate all chemistry to packages."""
    if (proposal is None) == (source_ref is None):
        raise ValueError("Supply exactly one of proposal or source_ref")
    mapping_evidence = None
    if source_ref:
        payload = _call(operations.store, source_ref, SOURCE_OPERATIONS)
        selected, _ = _selection(payload, step_id, realization_id, strategy_id, precursor_smiles)
        proposal = {key: selected[key] for key in (
            "target_smiles", "precursor_smiles", "mapped_reaction_smiles", "proposed_conditions",
        ) if key in selected}
        if payload["operation"] == "disconnect_target" and selected.get("condition_query_reaction_smiles"):
            proposal["saved_candidate"] = {"source_ref": source_ref, "realization_id": realization_id,
                                           "strategy_id": strategy_id}
            proposal, mapping_evidence = resolve_saved_candidate(operations.store, proposal)
    elif any(value is not None for value in (step_id, realization_id, strategy_id, precursor_smiles)):
        raise ValueError("step_id and realization_id require a saved source_ref")
    step = _step(proposal)
    refs = [*_strings(evidence_refs, "evidence_refs"), *([source_ref] if source_ref else []),
            *([forward_ref] if forward_ref else [])]
    context = _evidence(operations, refs)
    forward = None
    forward_execution = "not_run"
    forward_record = None
    if forward_ref:
        # Partial executions are valid evidence of an unfinished check, not of
        # chemical failure. They must still bind to the selected reaction.
        if not any(event.artifact_ref == forward_ref and event.kind == "call"
                   for event in operations.store.events()):
            raise ValueError("forward_ref must identify a recorded forward call")
        saved = operations.store.read_artifact(forward_ref)
        if saved.get("operation") != "assess_route_step_forward":
            raise ValueError("forward_ref must identify assess_route_step_forward")
        forward_record = saved.get("result")
        if not isinstance(forward_record, dict):
            raise ValueError("Forward call has no recorded result")
        from .forward_check import _selected_step

        audited_proposal, audited_assessment, _ = _selected_step(
            operations, forward_record["source_ref"], forward_record["step_id"],
        )
        if (canonical_smiles(step.target_smiles) != audited_assessment["canonical_target_smiles"]
                or canonical_smiles(step.precursor_smiles) != audited_assessment["canonical_precursor_smiles"]
                or step.proposed_conditions != audited_proposal.get("proposed_conditions")):
            raise ValueError("Saved forward check does not match this realization and recipe")
        forward_execution = forward_record["execution_status"]
        if saved.get("execution_status") != forward_execution:
            raise ValueError("Saved forward execution statuses disagree")
        if forward_execution == "completed":
            audit = forward_record["assessment"]
            forward = RetroForwardEvidence(
                starting_materials=audit["starting_materials"], intended_product=audit["intended_product"],
                validity=audit["validity"], targeted_replay_status=audit["targeted_replay_status"],
                intended_match=audit["intended_match"], best_competitor_product=audit["best_competitor_product"],
                warnings=tuple(audit["warnings"]), audited_recipe=step.proposed_conditions,
                checks=tuple(audit.get("checks", ())),
                intended_product_rank=audit.get("intended_product_rank"),
                score_margin=audit.get("score_margin"),
                ranking_definition_id=audit.get("blind_prediction", {}).get("ranking_definition_id"),
                audit_schema_version=audit.get("schema_version"),
            )
    artifacts = operations.store.manifest["baseline"]["artifacts"]
    available = [artifacts.get(name, {}).get("status") == "present"
                 for name in ("condition_index", "shared_core_index")]
    arguments = {}
    if all(available):
        engine = operations._conditions()
        arguments = {"condition_index": engine.index, "shared_core_index": engine.shared_core_index}
    result = assess_retro_validity(
        step, operations._external_route_library(), **arguments,
        forward_evidence=forward, forward_execution_status=forward_execution,
        candidate_limit=candidate_limit, match_limit=match_limit,
    )
    return {
        "schema_version": "retro_validity_investigation.v1", "proposal": step.to_dict(),
        "validity": result.to_dict(), "source_ref": source_ref,
        "selection": {"step_id": step_id, "realization_id": realization_id, "strategy_id": strategy_id},
        "mapping_evidence": mapping_evidence,
        "forward_ref": forward_ref,
        "artifact_warnings": (["CORPUS_ARTIFACT_PAIR_INCOMPLETE"]
                              if any(available) and not all(available) else []),
        "forward_execution": ({"execution_status": forward_execution, "error": forward_record.get("error")}
                              if forward_record else None),
        "experimental_feasibility": "not_established", **context,
    }

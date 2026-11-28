"""Saved step dossiers composed from existing assessments and attributed questions."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from reactive_taxonomy.chemistry.reaction_smiles import split_reaction_fields

from ..adapters.step_selection import SOURCE_OPERATIONS, _call, _selection
from ..core.store import InvestigationEvent
from .step_assessments import answer_step_assessments

if TYPE_CHECKING:
    from ..workspace import ScientificWorkspace


def step_plausibility_summary(result: Mapping[str, Any]) -> dict[str, Any]:
    """Show the recorded grade, search coverage and first source returned by the assessor."""
    assessment = result["assessment"]
    summary = {key: assessment[key] for key in (
        "query_reaction_smiles", "grade", "evidence_rank", "measured_evidence_status", "qualified_count", "candidate_count",
        "candidate_limit", "candidate_truncated", "matches_truncated", "warnings",
    )}
    summary.update(source_ref=result["source_ref"], step_id=result["step_id"], first_match=None)
    summary["support_counts"] = {key: assessment["support_counts"][key] for key in (
        "supporting_source_group_count", "comparable_observation_count", "low_yield_count",
        "zero_yield_count", "unreported_yield_count", "unknown_yield_count",
    )}
    if assessment["matches"]:
        match = assessment["matches"][0]
        row, context = match["precedent"], match["context"]
        first = {key: row[key] for key in ("reference_id", "reaction_id", "observation_id", "reaction_smiles", "level")}
        first.update(source_yield_pct=row["yield_pct"], outcome=match["outcome"],
                     contributes_to_support=match["contributes_to_support"],
                     differences=row["differences"][:3], difference_count=len(row["differences"]))
        first["context"] = {key: context[key] for key in ("status", "reverse_status", "minimum_partner_similarity")}
        for key in ("missing", "unresolved", "source_only", "source_unresolved"):
            first["context"][key] = context[key][:3]
            if len(context[key]) > 3:
                first["context"][f"{key}_count"] = len(context[key])
        summary["first_match"] = first
    return summary


def step_dossier(
    workspace: ScientificWorkspace, source_ref: str, *, step_id: str | None = None,
    evidence_refs: tuple[str, ...] = (),
) -> InvestigationEvent:
    """Save an inspectable view, retaining exact scope and independent assessment dimensions."""
    store = workspace.store
    if len(evidence_refs) > 10:
        raise ValueError("Attach at most 10 relevant saved evidence references")
    condition_operations = {"recommend_conditions", "assess_recipe", "assess_proposed_recipe", "analyze_reaction"}
    payload = _call(store, source_ref, (SOURCE_OPERATIONS - {"disconnect_target"}) | condition_operations)
    if payload["operation"] in condition_operations:
        if step_id is not None:
            raise ValueError("A single reaction call does not require step_id")
        sides = split_reaction_fields(payload["arguments"]["reaction_smiles"])
        if len(sides) != 3 or not sides[0] or not sides[2]:
            raise ValueError("Dossier requires a product-specified reaction")
        proposal = {"precursor_smiles": sides[0], "target_smiles": sides[2]}
    else:
        proposal, _ = _selection(payload, step_id, None)
    refs = tuple(dict.fromkeys((source_ref, *evidence_refs)))
    # The established projection owns exact graph matching and recipe attribution.
    # Linked summaries below are context only, never inferred support for this step.
    answer = {
        "molecules": [{"id": "reactants", "smiles": proposal["precursor_smiles"]},
                      {"id": "products", "smiles": proposal["target_smiles"]}],
        "steps": [{"id": "step", "reactant_ids": ["reactants"], "product_ids": ["products"]}],
        "evidence_refs": list(refs),
    }
    source_summary = workspace.call_summary(source_ref)
    linked = workspace.batch_summary(list(evidence_refs)) if evidence_refs else None
    assessments = answer_step_assessments(store, answer)["step"]
    assessments["step_plausibility"] = []
    if evidence_refs:
        cited = set(evidence_refs)
        for event in store.events():
            if event.kind != "call" or event.artifact_ref not in cited:
                continue
            call = store.read_artifact(event.artifact_ref)
            if call.get("operation") != "assess_step_plausibility" or call.get("execution_status") != "completed":
                continue
            result = call["result"]
            if (result["source_ref"] == source_ref and result["step_id"] == step_id
                    and result["proposal"] == proposal):
                assessments["step_plausibility"].append({
                    "artifact_ref": event.artifact_ref, **step_plausibility_summary(result),
                })
    questions = workspace.list_questions(source_ref=source_ref, step_id=step_id)
    assessment_refs = tuple(item["artifact_ref"] for name in ("structural_assessments", "recipe_assessments")
                            for item in assessments[name])
    question_refs = tuple(item["revision_ref"] for item in questions["items"])
    return store.append("step_dossier", {
        "schema_version": "scientific_step_dossier.v1", "source_ref": source_ref, "step_id": step_id,
        "proposal": proposal, "source_summary": source_summary, "assessments": assessments,
        "linked_evidence": linked, "questions": questions,
        "experimental_feasibility": "not_established", "review_status": "unreviewed",
        "limitations": [
            "This saved view runs no chemistry and does not change admission or ranking.",
            "Linked evidence is agent-selected context, not verified semantic support.",
            "Recipe checks apply only to their saved recipes; changed conditions require reassessment.",
            "Questions apply only to this exact source and step; other revisions require their own investigation.",
            "Missing checks and empty searches do not establish feasibility or impossibility.",
            "Cite underlying calls and sources, not this dossier or agent question resolutions.",
        ],
    }, evidence_refs=tuple(dict.fromkeys((*refs, *assessment_refs, *question_refs))))

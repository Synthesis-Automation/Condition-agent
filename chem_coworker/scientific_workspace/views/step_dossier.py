"""Saved step dossiers composed from existing assessments and attributed questions."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ..adapters.step_selection import SOURCE_OPERATIONS, _call, _selection
from ..core.store import InvestigationEvent
from .step_assessments import answer_step_assessments

if TYPE_CHECKING:
    from ..workspace import ScientificWorkspace


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
        sides = payload["arguments"]["reaction_smiles"].split(">")
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

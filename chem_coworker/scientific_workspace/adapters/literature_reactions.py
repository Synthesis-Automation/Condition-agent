"""Prepare source-linked drawing inputs using recorded structures and exact passages.

This adapter composes taxonomy graph audits with immutable source evidence. It
does not translate names, recognize scheme images, admit precedents or establish
experimental feasibility. Source acquisition and structure origin remain separate.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
import re
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from condition_registry.quantity_audit import audit_source_quantities
from reactive_taxonomy.structure_audit import audit_structure

from ..core.store import InvestigationStore
from .literature import load_source_passage


class ParticipantInput(BaseModel):
    """An explicitly identified participant and its immutable source passage."""

    model_config = ConfigDict(extra="forbid", strict=True)
    name: str = Field(min_length=1, max_length=300)
    compound_id: str = Field(min_length=1, max_length=300)
    evidence_ref: str = Field(pattern=r"^sha256:[0-9a-f]{64}$")
    material_form: str = Field(min_length=1, max_length=500)
    smiles: str | None = Field(default=None, min_length=1, max_length=4000)
    reported_formula: str | None = Field(default=None, min_length=1, max_length=100)


class SourceConflict(BaseModel):
    """Agent-identified conflicting passages, preserved without silently choosing one."""

    model_config = ConfigDict(extra="forbid", strict=True)
    description: str = Field(min_length=1, max_length=1000)
    excerpt_refs: list[str] = Field(min_length=2, max_length=5)


def _indexed_record(store: InvestigationStore, reference: str, reaction_id: str) -> dict[str, Any]:
    kinds = {event.artifact_ref: event.kind for event in store.events()}
    payload = store.read_artifact(reference)
    if kinds.get(reference) != "call" or payload.get("execution_status") != "completed":
        raise ValueError("Indexed structures require a completed recorded corpus lookup")
    result = payload.get("result", {})
    operation = payload.get("operation")
    if operation == "search_fragment_precedents":
        records = [hit.get("record", {}) for hit in result.get("hits", [])]
    elif operation == "get_precedents":
        records = result.get("records", [])
    else:
        raise ValueError("Reuse indexed graphs from search_fragment_precedents or get_precedents")
    matches = [record for record in records if record.get("reaction_id") == reaction_id]
    if len(matches) != 1:
        raise ValueError("Select one unambiguous saved reaction record; use an observation-specific lookup")
    record = matches[0]
    if not record.get("reference_id") or not record.get("reaction_smiles"):
        raise ValueError("Indexed graph reuse requires reaction structures and an explicit publication identity")
    return record


def prepare_literature_reaction(
    store: InvestigationStore, *, source_ref: str, source_id: str, title: str, locator: str,
    structure_evidence: str, reactants: list[dict[str, Any]], products: list[dict[str, Any]],
    structure_origin: Literal["source_explicit", "reconstructed_from_description"] = "reconstructed_from_description",
    indexed_ref: str | None = None, reaction_id: str | None = None,
    conditions: list[dict[str, Any]] | None = None, yield_info: dict[str, Any] | None = None,
    source_conflicts: list[dict[str, Any]] | None = None,
    limitations: list[str] | None = None,
    scheme_refs: list[str] | None = None,
) -> dict[str, Any]:
    """Prepare a literature_reaction.v2 block; every participant needs an exact excerpt.

    Participants specify name, compound_id, evidence_ref, material_form and either
    explicit smiles or indexed_ref/reaction_id for the whole source reaction. For
    indexed reuse omit all participant SMILES; sides/counts come from that record.
    The captured root's reported_reference_id must match the indexed publication.
    Optional reported_formula must occur literally in the participant passage;
    disagreements are retained as conflicts, never silently repaired. Conditions
    and yields remain explicit reported claims, not inferred from proposed steps.
    """
    from ..answers.answer_contracts import LiteratureReaction

    source, text, root_ref = load_source_passage(store, source_ref)
    if not isinstance(structure_evidence, str) or structure_evidence not in text:
        raise ValueError("structure_evidence must quote the captured source exactly")
    if not 1 <= len(reactants) <= 20 or not 1 <= len(products) <= 20:
        raise ValueError("Supply 1 to 20 participants on each reaction side")
    if structure_origin not in {"source_explicit", "reconstructed_from_description"}:
        raise ValueError("Declare source_explicit or reconstructed_from_description")
    inputs = {side: [ParticipantInput.model_validate(item) for item in items]
              for side, items in (("reactants", reactants), ("products", products))}
    record = None
    structures = {}
    if indexed_ref is not None or reaction_id is not None:
        if not indexed_ref or not reaction_id:
            raise ValueError("indexed_ref and reaction_id must be supplied together")
        record = _indexed_record(store, indexed_ref, reaction_id)
        if source.get("reported_reference_id") != record["reference_id"]:
            raise ValueError("Captured publication attribution must match the indexed reference_id")
        parts = record["reaction_smiles"].split(">")
        if len(parts) != 3 or not parts[0] or not parts[2]:
            raise ValueError("Indexed source reaction must contain explicit reactant/product graphs")
        structures = {"reactants": parts[0].split("."), "products": parts[2].split(".")}
        if any(item.smiles is not None for items in inputs.values() for item in items):
            raise ValueError("Omit participant SMILES when reusing indexed structures")
        if any(len(inputs[side]) != len(structures[side]) for side in inputs):
            raise ValueError("Participant counts must match the saved source reaction; salts are not stripped")
        structure_origin = "indexed_record"
    checks, quantity_checks, participants = [], [], {"reactants": [], "products": []}
    refs = [source_ref]
    if not isinstance(scheme_refs or [], list) or len(scheme_refs or []) > 5:
        raise ValueError("scheme_refs must contain at most five captured source images")
    from .source_images import load_source_image

    for reference in scheme_refs or []:
        load_source_image(store, reference, root_ref)
        refs.append(reference)
    notes = list(limitations or [])
    notes.append("Graph checks do not verify the paper's structure assignment, tautomer, or experimental feasibility.")
    for side, items in inputs.items():
        for index, item in enumerate(items):
            _, passage, passage_root = load_source_passage(store, item.evidence_ref)
            if passage_root != root_ref:
                raise ValueError("Participant passages must belong to this captured source")
            if not any(event.kind == "literature_excerpt" and event.artifact_ref == item.evidence_ref
                       for event in store.events()):
                raise ValueError("Each participant requires an exact recorded literature_excerpt")
            if not re.search(r"(?<!\w)" + re.escape(item.compound_id) + r"(?!\w)", passage):
                raise ValueError("Participant compound_id must occur in its captured passage")
            if item.reported_formula is not None and item.reported_formula not in passage:
                raise ValueError("reported_formula must occur literally in the participant passage")
            smiles = structures[side][index] if record else item.smiles
            if smiles is None:
                raise ValueError("Supply participant SMILES or reuse a saved indexed reaction")
            if structure_origin == "source_explicit" and smiles not in passage:
                raise ValueError("Source-explicit SMILES must occur in that participant's passage")
            audit = audit_structure(smiles).to_dict()
            quantities = audit_source_quantities(smiles, passage, (item.name, item.compound_id))
            for quantity in quantities:
                quantity_checks.append({"side": side, "component_index": index,
                                        "evidence_ref": item.evidence_ref, **asdict(quantity)})
                if quantity.status == "conflicting":
                    notes.append(
                        f"{side}[{index}] source mass/amount conflict in {quantity.source_text!r}; "
                        "neither reported value is selected or corrected."
                    )
            status = "invalid" if not audit["valid"] else "graph_checked_assignment_unverified"
            if item.reported_formula and audit["formula"] != item.reported_formula:
                if audit["valid"]:
                    status = "conflicting"
                notes.append(f"{side}[{index}] formula conflict: graph {audit['formula']}; source {item.reported_formula}.")
            checks.append({"side": side, "component_index": index, "status": status,
                           "reported_formula": item.reported_formula, "audit": audit})
            participants[side].append({"name": item.name, "smiles": smiles,
                                       "compound_id": item.compound_id, "evidence_ref": item.evidence_ref,
                                       "material_form": item.material_form, "graph_status": status,
                                       "formula": audit["formula"]})
            refs.append(item.evidence_ref)
    conflicts = [SourceConflict.model_validate(item).model_dump() for item in source_conflicts or []]
    for conflict in conflicts:
        if len(set(conflict["excerpt_refs"])) < 2:
            raise ValueError("A source discrepancy requires at least two distinct exact passages")
        for reference in conflict["excerpt_refs"]:
            _, _, conflict_root = load_source_passage(store, reference)
            if conflict_root != root_ref or not any(event.kind == "literature_excerpt" and event.artifact_ref == reference
                                                    for event in store.events()):
                raise ValueError("Conflicts require exact passages from the same source")
            refs.append(reference)
        notes.append("Unresolved source discrepancy: " + conflict["description"])
    if indexed_ref:
        refs.append(indexed_ref)
        notes.append("Indexed graph identity is reused; association with this captured publication remains agent-attributed.")
    provenance = {"acquisition": source["acquisition"], "retrieval_status": source["retrieval_status"],
                  "extraction_status": source["extraction"]["status"], "claim_support": "not_assessed",
                  "limitations": source["extraction"].get("limitations", [])}
    block = LiteratureReaction.model_validate({
        "schema_version": "literature_reaction.v2", "title": title, "source_id": source_id,
        "locator": locator, "structure_origin": structure_origin, "structure_evidence": structure_evidence,
        **participants,
        "conditions": [{"limitations": [], **claim} for claim in conditions or []],
        "yield_info": {"limitations": [], **yield_info} if yield_info is not None else None,
        "limitations": notes,
        "source_provenance": provenance,
    }).model_dump()
    for claim in [*block["conditions"], *([block["yield_info"]] if block["yield_info"] else [])]:
        if claim["basis"] != "reported" or source_id not in claim["source_ids"]:
            raise ValueError("Source conditions/yield require reported attribution to source_id")
        if claim["text"] not in text:
            raise ValueError("Prepared source conditions/yields must quote the captured text; adaptations belong on the proposed step")
    return {"schema_version": "literature_reaction_preparation.v2", "literature_reaction": block,
            "source_ref": source_ref, "source_provenance": provenance,
            "indexed_record": {"artifact_ref": indexed_ref, "reaction_id": reaction_id,
                                  "reference_id": record["reference_id"]} if record else None,
            "structure_origin": structure_origin, "participant_count": len(checks),
            "structure_checks": checks, "source_conflicts": conflicts, "limitations": notes,
            "quantity_checks": quantity_checks,
            "quantity_check_scope": "explicit_adjacent_mass_amount_pairs_only",
            "scheme_refs": scheme_refs or [], "assignment_verification": "not_performed",
            "evidence_refs": list(dict.fromkeys(refs))}


def prepared_literature_reaction(store: InvestigationStore, reference: str) -> dict[str, Any]:
    """Return the exact recorded block and bind it to its preparation artifact."""
    payload = store.read_artifact(reference)
    if (not any(event.kind == "call" and event.artifact_ref == reference for event in store.events())
            or payload.get("operation") != "prepare_literature_reaction"
            or payload.get("execution_status") != "completed"):
        raise ValueError("A completed prepare_literature_reaction call is required")
    block = deepcopy(payload["result"]["literature_reaction"])
    block["preparation_ref"] = reference
    return block

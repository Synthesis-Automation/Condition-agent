"""Inspect source precedents through canonical domain lookups and pinned catalogs."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from core_retrosynthesis.chemistry import canonical_smiles
from core_retrosynthesis.step_precedents import lookup_reaction_precedents

from .source_catalogs import reference_records

if TYPE_CHECKING:
    from .operations import ScientificOperations

from .step_selection import SCHEMA_VERSION, SOURCE_OPERATIONS, _call, _selection


def inspect_step_precedents(
    operations: ScientificOperations, source_ref: str, step_id: str | None = None,
    realization_id: str | None = None, offset: int = 0, limit: int = 3,
) -> dict[str, Any]:
    """Inspect at most five reactions on a page from one saved step's support.

    Disconnections use the canonical top-20 template lookup; assessments retain
    their saved bounded match set. Counts never imply whole-corpus coverage.
    """
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 5:
        raise ValueError("offset must be nonnegative and limit must be between one and five")
    payload = _call(operations.store, source_ref, SOURCE_OPERATIONS)
    selected, assessment = _selection(payload, step_id, realization_id)
    target = canonical_smiles(selected["target_smiles"])
    precursors = canonical_smiles(selected["precursor_smiles"])
    if not target or not precursors:
        raise ValueError("Selected step requires valid reactant and product structures")
    library = operations._external_route_library()
    if realization_id:
        lookup = lookup_reaction_precedents(
            step_id=realization_id, template_id=selected["template_id"], operator_id=selected["operator_id"],
            product_smiles=target, precursor_smiles=precursors, library=library, limit=20,
        )
        matches = [item.to_dict() for item in lookup.matches]
        available = lookup.available_precedent_count
        scope = "selected_template_top_20"
    else:
        matches = assessment.get("precedent_matches", [])
        available = None
        scope = "saved_assessment_matches"
    # A repeated record through several operators is one supporting reaction here.
    unique = {}
    for match in matches:
        key = (match["reaction_id"], match["reference_id"], match["mapped_reaction_smiles"])
        unique.setdefault(key, match)
    matches = list(unique.values())
    if offset > len(matches):
        raise ValueError("offset exceeds the saved precedent selection")
    page = matches[offset:offset + limit]
    references, reference_status = reference_records(operations, {item["reference_id"] for item in page})
    ids = list(dict.fromkeys(item["reaction_id"] for item in page))
    observations, procedures = {}, {}
    if ids:
        try:
            observations = operations.get_precedents(ids, limit=20)
        except FileNotFoundError:
            observations = {"records": [], "availability": "index_unavailable"}
        try:
            procedures = operations.get_procedures(ids)
        except FileNotFoundError:
            procedures = {"records": [], "availability": "catalog_unavailable"}
    templates = {item.template_id: item for item in library.templates}
    records = []
    for match in page:
        template = templates.get(match["template_id"])
        try:
            comparison = operations.compare_molecules(target, match["product_smiles"], timeout_seconds=1).to_dict()
        except ValueError as exc:
            comparison = {"status": "unavailable", "warnings": [str(exc)]}
        rows = [{key: row.get(key) for key in (
            "reaction_id", "observation_id", "reference_id", "reaction_smiles", "resolved_recipe", "yield_pct",
            "condition_uncertain", "chemistry_status", "condition_status",
        )} for row in observations.get("records", []) if row["reaction_id"] == match["reaction_id"]]
        source_procedures = [row for row in procedures.get("records", [])
                             if row["reaction_id"] == match["reaction_id"]]
        records.append({
            **match, "reaction_smiles": f"{match['precursor_smiles']}>>{match['product_smiles']}",
            "support_kind": "template_precedent", "reference_record": references.get(match["reference_id"]),
            "same_recorded_product": canonical_smiles(match["product_smiles"]) == target,
            "same_recorded_precursors": canonical_smiles(match["precursor_smiles"]) == precursors,
            "product_comparison": comparison,
            "template_context": {"edit_tokens": list(template.edit_tokens), "handle_signature": template.handle_signature,
                                 "stereo_policy": template.stereo_policy} if template else None,
            "observations": rows, "procedures": source_procedures,
            "experimental_link_scope": "reaction_id_only_template_has_no_observation_id",
            "limitations": ["Template membership and similarity do not establish transfer or experimental feasibility.",
                            "Associated conditions and yields belong to their source observations, not the proposed step.",
                            "Substrate, functional-group and stereochemical transfer require inspection."],
        })
    return {
        "schema_version": SCHEMA_VERSION, "source_ref": source_ref,
        "selection": {"step_id": step_id, "realization_id": realization_id,
                      "target_smiles": target, "precursor_smiles": precursors},
        "status": "precedents_available" if matches else "no_precedents_retrieved",
        "assessment_status": assessment.get("status"), "assessment_warnings": assessment.get("warnings", []),
        "scope": scope, "available_template_records": available, "saved_match_count": len(matches),
        "retrieval_truncated": available > 20 if available is not None else None,
        "page": {"offset": offset, "limit": limit, "returned": len(page),
                 "next_offset": offset + limit if offset + limit < len(matches) else None},
        "distinct_references_on_page": len({item["reference_id"] for item in page if item["reference_id"]}),
        "reference_catalog_status": reference_status,
        "observation_page": {key: value for key, value in observations.items() if key != "records"},
        "procedure_catalog_status": procedures.get("availability", "catalog_available" if ids else "not_requested"),
        "precedents": records, "experimental_feasibility": "not_established",
        "limitations": ["Counts describe the selected template or saved assessment, not an exhaustive literature search.",
                        "An empty result does not establish that no experimental precedent exists.",
                        "This inspection does not change route admission or validate proposed conditions."],
    }

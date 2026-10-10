"""Detailed operation-specific views of saved calls; never replacement scientific contracts."""

from __future__ import annotations

from collections.abc import Mapping
from itertools import islice
from typing import Any

from .projection import (
    _COMMON,
    _COMPATIBILITY,
    _PRECEDENT,
    _RECOMMENDATION,
    _Projection,
)


def _result_summary(operation: str, value: Mapping[str, Any], view: _Projection) -> dict[str, Any]:
    path = "$.result"
    summary = view.pick(value, _COMMON, path)
    if operation == "inspect_step_precedents":
        summary.update(view.pick(value, ("assessment_proposal",), path))
        summary.update(view.pick(value, ("selection", "scope", "saved_match_count", "available_template_records",
                                         "retrieval_truncated", "page", "distinct_references_on_page",
                                         "assessment_status", "assessment_warnings"), path))
        view.add_list(summary, value, "precedents", ("match_id", "reaction_id", "reference_id", "reaction_smiles",
                      "support_kind", "product_similarity", "precursor_similarity", "same_recorded_product",
                      "same_recorded_precursors", "product_comparison", "limitations"), path)
    elif operation == "analyze_molecule":
        summary.update(view.pick(value, (
            "canonical_smiles", "component_count", "atom_count", "heavy_atom_count", "formal_charge",
        ), path))
        for key, fields in (
            ("stereocenters", ("atom_index", "assignment", "assigned")),
            ("double_bond_stereo", ("bond_index", "assignment")),
            ("motifs", ("motif_id", "chemist_label")),
            ("reactive_sites", ("site_type", "chemist_label", "availability", "warnings")),
        ):
            view.add_list(summary, value, key, fields, path)
    elif operation == "assess_starting_material":
        summary.update(view.pick(value, (
            "canonical_smiles", "stop_expansion", "stop_reason", "availability",
            "molecular_weight", "mw_threshold", "definition_version", "policy",
            "unavailable_starting_materials",
        ), path))
        view.add_nested(summary, value, "registry", (
            "status", "source", "candidates", "matched_identifier", "warnings", "reason",
            "error_type", "message",
        ), path)
        view.add_nested(summary, value, "literature", (
            "status", "index_id", "match_scope", "source_scope", "source_coverage_complete",
            "precedents", "has_more", "warnings", "reason", "error_type", "message",
        ), path)
    elif operation == "compare_molecules":
        summary.update(view.pick(value, (
            "same_constitution", "stereo_relationship", "core_method", "core_atom_count",
            "left_coverage", "right_coverage", "alignment_count_observed", "alignment_ambiguous",
            "embeddings_truncated", "alignments_truncated", "search_timed_out", "definition_version",
        ), path))
        for side in ("left", "right"):
            view.add_nested(summary, value, side, ("canonical_smiles", "warnings", "atom_id_scope"), path)
        view.add_list(summary, value, "alignments", (
            "left_only_atom_ids", "right_only_atom_ids", "left_boundaries", "right_boundaries",
        ), path, limit=1)
    elif operation == "inspect_reactive_sites":
        summary.update(view.pick(value, ("selected_atom_ids", "radius", "definition_version"), path))
        view.add_nested(summary, value, "molecule", ("canonical_smiles", "warnings", "atom_id_scope"), path)
        for key in ("sites", "other_sites"):
            summary[f"{key}_count"] = len(value.get(key, []))
            view.add_list(summary, value, key, (
                "hypothesis_id", "site_type", "atom_indices", "chemist_label", "availability", "warnings",
            ), path)
        view.add_list(summary, value, "motifs", ("motif_id", "chemist_label", "atom_indices"), path)
        def site_profile(item: Any, child: str) -> Any:
            compact = view.pick(item, ("hypothesis_id", "center_atom_index"), child)
            profile = item.get("reactivity_profile", {}) if isinstance(item, Mapping) else {}
            location = f"{child}.reactivity_profile"
            compact.update(view.pick(profile, ("status", "context_kind", "flags"), location))
            view.add_nested(compact, profile, "steric", ("accessibility_class", "evidence"), location)
            view.add_nested(compact, profile, "electronic", ("activation_class", "evidence"), location)
            return compact
        summary["site_profiles"] = view.preview(
            value.get("environments", []), f"{path}.environments", site_profile, limit=3,
        )
    elif operation == "analyze_reaction":
        summary.update(view.pick(value, (
            "evidence_quality", "edit_archetype", "transformation_class", "named_family", "compatible_named_families",
        ), path))
        for key, fields in (
            ("reaction_signature", ("signature_id", "evidence_quality", "transformation_class", "schema_version")),
            ("observation", ("valid", "evidence_quality", "evidence_confidence", "warnings", "error")),
            ("interpretation", ("evidence_quality", "named_family", "compatible_named_families", "warnings")),
            ("reaction_completeness", ("status", "warnings", "evidence", "suspected_missing_reactant",
                                       "suspected_insufficient_reactant_multiplicity", "product_heavy_atom_coverage")),
        ):
            view.add_nested(summary, value, key, fields, path)
    elif operation == "generate_weak_label_screening_array":
        summary.update(view.pick(value, (
            "recommendation_mode", "reaction_type_id", "reaction_type_hint_id",
            "source_reaction_type_candidates", "candidate_count", "compatible_candidate_count",
            "excluded_candidate_count", "recipe_count",
        ), path))
        summary["recommendation_count"] = len(value.get("recommendations", []))
        view.add_list(summary, value, "query_participants", (
            "role", "signature", "display_label",
        ), path)
        view.add_list(summary, value, "recommendations", (
            "rank", "recipe_id", "support", "source_reaction_types", "source_row_numbers",
            "source_matches", "compatibility_score", "compatibility_evidence",
            "historical_yield_pct", "cautions", "explanation",
        ), path, limit=3)
    elif operation == "recommend_conditions":
        summary.update(view.pick(value, (
            "recommendation_mode", "retrieval_level", "search_scope", "candidate_count",
            "independent_candidate_count", "compatible_candidate_count",
            "independent_compatible_candidate_count", "excluded_candidate_count", "external_mapping_status",
        ), path))
        view.add_list(summary, value, "recommendations", _RECOMMENDATION, path, limit=3)
        view.add_list(summary, value, "retrieval_trace", (
            "level", "status", "candidate_count", "compatible_candidate_count", "excluded_candidate_count",
            "independent_compatible_candidate_count", "broadening_reason",
        ), path)
    elif operation == "assess_recipe":
        summary.update(view.pick(value, _COMPATIBILITY, path))
    elif operation == "assess_proposed_recipe":
        summary.update(view.pick(value, ("reaction_smiles", "process_coverage"), path))
        view.add_nested(summary, value, "compatibility", _COMPATIBILITY, path)
    elif operation == "prepare_literature_reaction":
        summary.update(view.pick(value, ("structure_origin", "participant_count", "quantity_check_scope"), path))
        summary["quantity_conflict_count"] = sum(item.get("status") == "conflicting" for item in value.get("quantity_checks", []))
        view.add_list(summary, value, "quantity_checks", (
            "side", "component_index", "status", "source_text", "evidence_ref", "expected_mass_g", "reason",
        ), path)
    elif operation == "assess_route_step_forward":
        summary.update(view.pick(value, (
            "execution_status", "source_ref", "step_id", "question",
        ), path))
        view.add_nested(summary, value, "assessment", (
            "targeted_replay_status", "intended_match", "intended_product_rank", "best_competitor_product",
            "score_margin", "disposition", "validity", "advisory_only", "warnings",
        ), path)
        execution = value.get("execution")
        if isinstance(execution, Mapping):
            compact: dict[str, Any] = {}
            location = f"{path}.execution"
            view.add_nested(compact, execution, "timings", ("elapsed_seconds", "timeout_seconds"), location)
            view.add_list(compact, execution, "stages", (
                "stage", "status", "elapsed_seconds", "duration_seconds",
            ), location, limit=16)
            view.add_nested(compact, execution, "diagnostics", ("stages.jsonl", "stderr.log"), location)
            summary["execution"] = compact
        view.add_nested(summary, value, "provenance", (
            "source_assessment_id", "selected_operator_match_id", "operator_id",
            "input_hashes_verified", "library_source_verified",
        ), path)
    elif operation == "disconnect_composite":
        summary.update(view.pick(value, ("target_smiles", "catalog_id", "diagnostics"), path))
        def composite(item: Any, child: str) -> Any:
            result = view.pick(item, (
                "action_id", "strategy_id", "intermediate_smiles", "terminal_precursor_smiles",
                "physical_step_count", "physical_step_cost", "condition_compatibility_status",
                "one_pot_status", "score",
            ), child)
            if isinstance(item, Mapping):
                view.add_nested(result, item, "dependency", (
                    "admitted", "status", "relationship_class", "dependency_class", "lineage_status", "warnings",
                ), child)
                view.add_list(result, item, "physical_steps", (
                    "forward_step_number", "reaction_smiles", "operator_id", "template_id", "precedent_reaction_ids",
                    "forward_validation_status", "precursor_compatibility_disposition", "reaction_compatibility_disposition",
                    "condition_status", "evidence_kind", "selectivity_warnings",
                ), child, limit=2)
            return result
        summary["actions"] = view.preview(value.get("actions", []), f"{path}.actions", composite, limit=3)
        view.add_list(summary, value, "dependency_reviews", (
            "strategy_id", "intermediate_smiles", "dependency",
        ), path, limit=3)
        view.add_list(summary, value, "one_step_fallbacks", (
            "precursor_smiles", "proposed_reaction_smiles", "forward_validation_status",
        ), path, limit=3)
    elif operation == "disconnect_target":
        summary.update(view.pick(value, ("bond_focus", "search_diagnostics"), path))
        view.add_nested(summary, value, "request", (
            "target_smiles", "top_k", "max_realizations_per_strategy", "max_templates_to_apply",
            "max_candidates_to_validate", "use_context", "include_l0", "include_conditions",
            "required_disconnection_bond", "focus_target_smiles",
        ), path)
        realization_fields = (
            "realization_id", "target_smiles", "precursor_smiles", "proposed_reaction_smiles",
            "forward_validation_status", "precedent_reaction_ids", "selectivity_warnings",
            "precursor_compatibility_disposition", "reaction_compatibility_disposition",
            "bond_focus_check",
            "strategic_class", "strategic_complexity_score", "strategic_candidate",
        )

        def strategy(item: Any, child: str) -> Any:
            result = view.pick(item, (
                "strategy_id", "strategy_rank", "operator_id", "target_smiles",
                "independent_reference_support", "precedent_reaction_ids",
                "returned_realization_count", "total_realization_count",
            ), child)
            if isinstance(item, Mapping):
                view.add_nested(result, item, "representative", realization_fields, child)
                representative = item.get("representative")
                if isinstance(representative, Mapping):
                    view.add_nested(result["representative"], representative, "strategic_complexity", (
                        "definition_id", "evidence", "graph_complexity_reduction_fraction",
                        "extra_precursor_heavy_atom_count", "warnings",
                    ), f"{child}.representative")
                view.add_list(result, item, "alternate_realizations", realization_fields, child, limit=2)
            return result

        if "strategies" in value:
            summary["strategies"] = view.preview(value["strategies"], f"{path}.strategies", strategy, limit=3)
        view.add_list(summary, value, "condition_evidence", (
            "strategy_id", "evidence", "condition_selectivity_assessment",
        ), path, limit=3)
    elif operation == "assess_retro_validity":
        summary.update(view.pick(value, ("source_ref", "selection", "forward_ref", "artifact_warnings"), path))
        view.add_nested(summary, value, "validity", (
            "status", "structural_status", "precedent_grade", "evidence_rank", "suggested_action",
            "forward_status", "forward_execution_status", "cautions", "unresolved_checks", "warnings",
            "definition_id", "definition_version", "score_semantics", "ranking_influence",
        ), path)
        validity = value.get("validity", {})
        if isinstance(validity, Mapping):
            for name in ("operator_precedent_support", "corpus_precedent_support"):
                view.add_nested(summary["validity"], validity, name, (
                    "status", "strongest_level", "level_counts", "distinct_reference_count",
                    "candidate_count", "qualified_count", "candidate_truncated", "matches_truncated", "scope",
                ), f"{path}.validity")
    elif operation in {"assess_route_step", "inspect_route_step", "assess_route_proposal",
                       "revise_route_branch"}:
        route = operation in {"assess_route_proposal", "revise_route_branch"}
        summary.update(view.pick(value, (
            "source_ref", "step_id", "upstream_step_ids", "downstream_step_ids", "assessment_options", "evidence_refs", "evidence_warnings",
        ), path))
        view.add_list(summary, value, "evidence_provenance", (
            "artifact_ref", "kind", "role", "source_url", "acquisition", "retrieval_status",
            "extraction_status", "claim_support", "warnings",
        ), path)
        if "assessment" in value:
            summary["assessment"] = view.assessment(value["assessment"], f"{path}.assessment", route=route)
        view.add_nested(summary, value, "material_constraints", (
            "status", "stock_availability", "blocked_leaf_smiles", "unavailable_starting_materials",
        ), path)
        view.add_nested(summary, value, "proposed_recipe_assessment", _COMPATIBILITY, path)
        recipes = value.get("proposed_recipe_assessments")
        if isinstance(recipes, Mapping):
            selected = list(islice(recipes, 5))
            summary["proposed_recipe_assessments"] = {
                str(key)[:80]: view.pick(recipes[key], _COMPATIBILITY, f"{path}.proposed_recipe_assessments.{str(key)[:80]}")
                for key in selected
            }
            if len(recipes) > len(selected):
                view.truncated(f"{path}.proposed_recipe_assessments", "mapping_preview", total_fields=len(recipes), shown_fields=len(selected))
        view.add_nested(summary, value, "revision", (
            "improvement_status", "reassessment_scope", "reason", "risks", "assumptions",
            "removed_step_ids", "added_step_ids", "replaced_step_ids", "preserved_step_ids",
        ), path)
    elif operation == "compare_route_proposals":
        def alternative(item: Any, child: str) -> Any:
            result = view.pick(item, ("status", "source_ref", "route_id", "step_count", "unresolved_step_ids", "limitations"), child)
            if isinstance(item, Mapping):
                view.add_nested(result, item, "material_constraints", ("status", "stock_availability", "blocked_leaf_smiles"), child)
                view.gates(result, item, "topology_gates", child)
                view.add_list(result, item, "steps", ("step_id", "status", "strongest_evidence_tier", "warnings"), child)
            return result
        if "alternatives" in value:
            summary["alternatives"] = view.preview(value["alternatives"], f"{path}.alternatives", alternative)
    elif operation == "propose_condition_adaptation":
        summary.update(view.pick(value, ("transfer_status", "source_ref", "observation_id", "assumptions", "risks"), path))
        view.add_nested(summary, value, "compatibility", _COMPATIBILITY, path)
        view.add_list(summary, value, "changes", ("field", "basis", "reason", "evidence_refs"), path)
    elif operation == "suggest_search_fragments":
        summary.update(view.pick(value, ("target_id", "target_smiles", "definition_version", "query_compiler_version",
                                         "selection_mode", "generated_count", "rejected_count",
                                         "generation_truncated", "output_truncated", "limitations"), path))
        view.add_list(summary, value, "candidates", ("candidate_id", "kind", "query", "query_format",
                      "topology", "target_atom_ids", "features", "reasons", "cautions", "matches_target"), path)
    elif operation == "propose_fragment_queries":
        summary.update(view.pick(value, ("target_smiles", "definition_version", "parent_query",
                                         "query_atoms", "aromatic_atom_ids", "output_truncated", "limitations"), path))
        view.add_list(summary, value, "variants", ("variant_id", "query", "query_format", "topology",
                      "relaxations", "reason", "alignment_ambiguous", "target_alignments_truncated"), path)
    elif operation == "investigate_fragment_precedent":
        summary.update(view.pick(value, ("source_ref", "target_smiles", "search_scope", "limitations"), path))
        view.add_nested(summary, value, "source", ("observation_id", "reaction_id", "reference_id",
                        "relationships", "procedure_availability", "warnings"), path)
        view.add_nested(summary, value, "comparison", ("status", "core_atom_count", "left_coverage",
                        "right_coverage", "alignment_ambiguous", "warnings"), path)
        view.add_nested(summary, value, "transfer", ("status", "core_admission_policy",
                        "compiled_source_template_count", "source_admissions", "arms", "limitations"), path)
    elif operation in {"search_fragment_precedents", "find_synthesis_precedents"}:
        summary.update(view.pick(value, (
            "search_status", "stop_reason", "query", "index_id", "source_scope", "source_coverage_complete",
            "counts", "relationship_groups", "group_count_scope", "ranking_scope", "returned_count",
            "refinement_hints", "output_truncated", "target_validation", "target_smiles", "definition_version", "limitations",
        ), path))
        view.add_list(summary, value, "hits", (
            "hit_id", "observation_id", "reaction_id", "reference_id", "relationships",
            "matched_side", "matched_sides", "match_extents", "matched_molecule_smiles",
            "product_smiles", "relationship_summary",
            "citation_availability",
            "warnings", "admission_tier", "admission_reasons", "procedure_availability",
            "procedure_match_scope", "inspect_paths", "discovery",
        ), path)
    elif operation == "inspect_condition_precedents":
        summary.update(view.pick(value, (
            "observation_count", "distinct_reference_count", "missing_reference_observation_ids", "page",
        ), path))
        view.add_nested(summary, value, "query_analysis", ("valid", "error", "warnings", "evidence_quality"), path)
        view.add_nested(summary, value, "procedure_catalog", ("availability", "missing_reaction_ids"), path)

        def precedent(item: Any, child: str) -> Any:
            result = view.pick(item, ("missing_operating_fields", "procedure_link_scope"), child)
            if isinstance(item, Mapping):
                view.add_nested(result, item, "observation", _PRECEDENT, child)
                view.add_nested(result, item, "compatibility", _COMPATIBILITY, child)
            return result

        if "precedents" in value:
            summary["precedents"] = view.preview(value["precedents"], f"{path}.precedents", precedent)
    else:
        summary.update(view.pick(value, (
            "retrieval_level", "candidate_count", "total", "offset", "next_offset", "page",
            "record_scope", "index_scope", "missing_reaction_ids",
        ), path))
        view.add_list(summary, value, "records", _PRECEDENT, path)
    return summary

def summarize_call(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Project a serialized call into a compact, evidence-preserving inspection view.

    The caller adds the event/artifact reference. ``inspection`` names JSON paths
    within that artifact and reports preview bounds. Execution completion is kept
    separate from scientific validity; no status is inferred from missing fields.
    """
    view = _Projection()
    summary = {"summary_schema_version": "scientific_call_summary.v1",
               "operation": view.value(payload.get("operation"), "$.operation"),
               "execution_status": view.value(payload.get("execution_status"), "$.execution_status"),
               "error": view.value(payload.get("error"), "$.error")}
    for key in ("duration_seconds", "result_bytes"):
        if key in payload:
            summary[key] = view.value(payload[key], f"$.{key}")
    result = payload.get("result")
    if isinstance(result, Mapping):
        summary["result_summary"] = _result_summary(str(payload.get("operation", "")), result, view)
    else:
        summary["result_summary"] = {}
    keys = list(islice(result, 25)) if isinstance(result, Mapping) else []
    summary["inspection"] = {
        "result_path": "$.result", "projection_only": True,
        "result_present": "result" in payload, "result_type": type(result).__name__,
        "result_fields": [str(key)[:80] for key in keys],
        "result_field_count": len(result) if isinstance(result, Mapping) else None,
        "collections": view.collections, "truncations": view.truncations,
        "truncation_count": view.truncation_count,
        "hint": "Inspect the saved call artifact for full evidence, omitted fields and unabridged warnings. "
                "Counts describe saved result collections; use source pagination for remaining records.",
    }
    return summary

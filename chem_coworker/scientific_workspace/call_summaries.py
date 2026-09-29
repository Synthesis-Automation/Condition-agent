"""Bounded views of recorded scientific results, without chemistry interpretation.

Summaries are projections, not replacement contracts. Original status strings,
uncertainties and evidence identifiers are retained; full results remain in the
call artifact. Collection counts refer to that artifact, not an entire corpus.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from itertools import islice
import json
from typing import Any


_COMMON = (
    "status", "valid", "error", "warnings", "schema_version", "review_status",
    "origin", "experimental_feasibility", "limitations", "ranking",
)
_COMPATIBILITY = (
    "status", "analysis_status", "compatible", "hard_conflicts",
    "unresolved_requirements", "checked_requirements", "evidence",
    "analysis_warnings", "penalty_ids", "schema_version",
)
_STEP = (
    "status", "strongest_evidence_tier", "actionable", "admission_eligible",
    "warnings", "canonical_target_smiles", "canonical_precursor_smiles",
    "assessment_id", "proposal_id", "ranking_influence",
)
_ROUTE = (
    "status", "actionable", "admission_eligible", "warnings", "route_id",
    "canonical_target_smiles", "unresolved_step_ids", "disconnected_step_ids",
    "leaf_smiles", "ranking_influence",
)
_RECOMMENDATION = (
    "rank", "recipe_id", "compatibility_status", "match_label", "match_level",
    "evidence_relation", "retrieval_level", "cautions", "compatibility_evidence",
    "match_details", "explanation", "support", "reference_support",
    "observation_support", "precedent_reaction_ids", "precedent_reference_ids",
    "historical_yield_pct",
)
_PRECEDENT = (
    "reaction_id", "reference_id", "observation_id", "source_uri", "source_url",
    "match_level", "match_label", "evidence_relation", "warnings", "status",
    "admission_status", "product_similarity", "precursor_similarity",
)


class _Projection:
    """Per-call bounds and inspection metadata; never retains mutable globals."""

    def __init__(self) -> None:
        self.collections: list[dict[str, Any]] = []
        self.truncations: list[dict[str, Any]] = []
        self.truncation_count = 0
        self.remaining_nodes = 450
        self.remaining_text = 12000
        self.text_limit = 400

    def truncated(self, path: str, reason: str, **counts: Any) -> None:
        self.truncation_count += 1
        if len(self.truncations) < 30:
            self.truncations.append({"path": path, "reason": reason, **counts})

    def value(self, value: Any, path: str, depth: int = 0) -> Any:
        """Bound a selected JSON value and disclose every omitted subtree."""
        # Preserve short status/validity values even if surrounding evidence used
        # the preview budget. These keys are selected by the fixed views below.
        key = path.rsplit(".", 1)[-1]
        if (key.endswith("status") or key in {"valid", "compatible", "actionable", "admission_eligible"}) and (
            value is None or isinstance(value, bool) or isinstance(value, str) and len(value) <= 80
        ):
            return value
        self.remaining_nodes -= 1
        if self.remaining_nodes < 0:
            self.truncated(path, "summary_node_budget")
            return {"summary_omitted": True}
        if isinstance(value, str):
            size = min(self.text_limit, self.remaining_text)
            shown = value[:size]
            self.remaining_text -= len(shown)
            if len(value) > size:
                self.truncated(path, "text_preview", total_characters=len(value), shown_characters=size)
                return shown + "…"
            return value
        if value is None or isinstance(value, (bool, int, float)):
            return value
        if depth >= 3:
            self.truncated(path, "nested_detail")
            return {"summary_omitted": True, "value_type": type(value).__name__}
        if isinstance(value, Mapping):
            keys = list(islice(value, 12))
            if len(value) > len(keys):
                self.truncated(path, "mapping_preview", total_fields=len(value), shown_fields=len(keys))
            projected = {}
            for key in keys:
                if not isinstance(key, str) or len(key) > 80:
                    self.truncated(path, "unexpected_or_long_mapping_key")
                    continue
                projected[key] = self.value(value[key], f"{path}.{key}", depth + 1)
            return projected
        if isinstance(value, (list, tuple)):
            return self.preview(value, path, lambda item, child: self.value(item, child, depth + 1),
                                describe_collection=False)
        self.truncated(path, "unexpected_value_type")
        return {"summary_omitted": True, "value_type": type(value).__name__}

    def pick(self, value: Any, fields: tuple[str, ...], path: str) -> Any:
        if not isinstance(value, Mapping):
            return self.value(value, path)
        return {key: self.value(value[key], f"{path}.{key}") for key in fields if key in value}

    def preview(
        self, value: Any, path: str,
        project: Callable[[Any, str], Any] | None = None, *, limit: int = 5,
        describe_collection: bool = True,
    ) -> Any:
        if not isinstance(value, (list, tuple)):
            return self.value(value, path)
        shown = min(limit, len(value))
        if (describe_collection or len(value) > shown) and len(self.collections) < 40:
            self.collections.append({"path": path, "total": len(value), "shown": shown,
                                     "omitted": len(value) - shown, "count_scope": "saved_result"})
        if len(value) > shown:
            self.truncated(path, "collection_preview", total=len(value), shown=shown)
        project = project or self.value
        return [project(item, f"{path}[{index}]") for index, item in enumerate(value[:shown])]

    def add_nested(
        self, target: dict[str, Any], source: Mapping[str, Any], key: str,
        fields: tuple[str, ...], path: str,
    ) -> None:
        if key in source:
            target[key] = self.pick(source[key], fields, f"{path}.{key}")

    def add_list(
        self, target: dict[str, Any], source: Mapping[str, Any], key: str,
        fields: tuple[str, ...], path: str, *, limit: int = 5,
    ) -> None:
        if key in source:
            target[key] = self.preview(source[key], f"{path}.{key}",
                                       lambda item, child: self.pick(item, fields, child), limit=limit)

    def gates(self, target: dict[str, Any], source: Mapping[str, Any], key: str, path: str) -> None:
        self.add_list(target, source, key, ("gate_id", "status", "summary", "warnings", "evidence_ids"),
                      path, limit=16)

    def assessment(self, value: Any, path: str, *, route: bool = False) -> Any:
        summary = self.pick(value, _ROUTE if route else _STEP, path)
        if not isinstance(value, Mapping):
            return summary
        self.gates(summary, value, "topology_gates" if route else "gates", path)
        if route:
            def step(item: Any, child: str) -> Any:
                result = self.pick(item, ("external_step_id",), child)
                if isinstance(item, Mapping):
                    self.add_nested(result, item, "assessment", (
                        "status", "strongest_evidence_tier", "actionable", "admission_eligible", "warnings",
                    ), child)
                return result
            if "step_assessments" in value:
                summary["step_assessments"] = self.preview(
                    value["step_assessments"], f"{path}.step_assessments", step, limit=10,
                )
        else:
            self.add_list(summary, value, "precedent_matches", _PRECEDENT, path)
            self.add_nested(summary, value, "forward_assessment", (
                "targeted_replay_status", "intended_match", "disposition", "validity", "advisory_only", "warnings",
            ), path)
            self.add_nested(summary, value, "condition_evidence", (
                "status", "recommender_valid", "retrieval_level", "recommendation_mode",
                "candidate_count", "compatible_candidate_count", "warnings", "error",
            ), path)
            if "compatibility_evidence" in value:
                evidence = value["compatibility_evidence"]
                if isinstance(evidence, Mapping):
                    compact: dict[str, Any] = {}
                    for key in ("precursor", "reaction"):
                        self.add_nested(compact, evidence, key, ("disposition", "warning_strength"), f"{path}.compatibility_evidence")
                    self.add_nested(compact, evidence, "functional_group_competition", (
                        "code", "message", "assessment_mode", "conditions_evaluated", "ranking_impact",
                    ), f"{path}.compatibility_evidence")
                    summary["compatibility_evidence"] = compact
                else:
                    summary["compatibility_evidence"] = self.value(evidence, f"{path}.compatibility_evidence")
        return summary


def _result_summary(operation: str, value: Mapping[str, Any], view: _Projection) -> dict[str, Any]:
    path = "$.result"
    summary = view.pick(value, _COMMON, path)
    if operation == "inspect_step_precedents":
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
    elif operation == "disconnect_target":
        view.add_nested(summary, value, "request", (
            "target_smiles", "top_k", "max_realizations_per_strategy", "max_templates_to_apply",
            "max_candidates_to_validate", "use_context", "include_l0", "include_conditions",
        ), path)
        realization_fields = (
            "realization_id", "target_smiles", "precursor_smiles", "proposed_reaction_smiles",
            "forward_validation_status", "precedent_reaction_ids", "selectivity_warnings",
            "precursor_compatibility_disposition", "reaction_compatibility_disposition",
        )

        def strategy(item: Any, child: str) -> Any:
            result = view.pick(item, (
                "strategy_id", "strategy_rank", "operator_id", "target_smiles",
                "independent_reference_support", "precedent_reaction_ids",
                "returned_realization_count", "total_realization_count",
            ), child)
            if isinstance(item, Mapping):
                view.add_nested(result, item, "representative", realization_fields, child)
                view.add_list(result, item, "alternate_realizations", realization_fields, child, limit=2)
            return result

        if "strategies" in value:
            summary["strategies"] = view.preview(value["strategies"], f"{path}.strategies", strategy, limit=3)
        view.add_list(summary, value, "condition_evidence", (
            "strategy_id", "evidence", "condition_selectivity_assessment",
        ), path, limit=3)
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
    elif operation == "search_fragment_precedents":
        summary.update(view.pick(value, (
            "search_status", "stop_reason", "query", "index_id", "source_scope", "source_coverage_complete",
            "counts", "relationship_groups", "group_count_scope", "ranking_scope", "returned_count",
            "refinement_hints", "output_truncated", "target_validation",
        ), path))
        view.add_list(summary, value, "hits", (
            "hit_id", "observation_id", "reaction_id", "reference_id", "relationships",
            "product_smiles", "relationship_summary",
            "citation_availability",
            "warnings", "admission_tier", "admission_reasons", "procedure_availability",
            "procedure_match_scope", "inspect_paths",
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


def summarize_call_brief(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Project saved results directly; never count or select from an earlier preview."""
    operation = payload.get("operation")
    source = payload.get("result")
    source = source if isinstance(source, Mapping) else {}
    result = {"summary_schema_version": "scientific_call_brief.v2"}
    result.update(_brief_fields(payload, (
        "operation", "execution_status", "error", "duration_seconds", "result_bytes",
    )))
    if not isinstance(payload.get("result"), Mapping):
        result["result_summary"] = {}
        result["inspection"] = {"result_path": "$.result", "projection_only": True,
                                "detail_available": True, "result_present": "result" in payload}
        return result
    overview = _brief_fields(source, (
        "status", "valid", "error", "warnings", "limitations", "review_status",
        "experimental_feasibility", "ranking", "search_status", "stop_reason",
        "canonical_smiles", "canonical_target_smiles", "retrieval_level",
        "candidate_count", "compatible_candidate_count", "returned_count",
        "definition_version", "selection_mode", "generated_count", "rejected_count",
        "query", "target_smiles", "source_ref", "step_id",
        "source_coverage_complete", "group_count_scope", "ranking_scope",
    ))
    if "warnings" in source and source["warnings"] in ([], ()):
        # Keep the existing root-level warning contract for workspace clients.
        overview["warnings"] = []
    if operation == "inspect_step_precedents":
        overview.update(_brief_fields(source, ("scope", "saved_match_count", "available_template_records",
                                              "retrieval_truncated", "page", "distinct_references_on_page",
                                              "observation_page", "reference_catalog_status", "procedure_catalog_status",
                                              "assessment_status", "assessment_warnings")))
        overview["selection"] = _brief_fields(source.get("selection", {}), (
            "step_id", "realization_id", "target_smiles", "precursor_smiles",
        ))
        overview["precedents"] = [_brief_fields(item, (
            "match_id", "reaction_id", "reference_id", "reaction_smiles", "support_kind",
            "product_similarity", "precursor_similarity", "same_recorded_product", "same_recorded_precursors",
            "limitations",
        )) for item in source.get("precedents", [])[:3]]
        for brief, item in zip(overview["precedents"], source.get("precedents", [])):
            brief["product_comparison"] = _brief_fields(item.get("product_comparison", {}), (
                "status", "same_constitution", "stereo_relationship", "left_coverage", "right_coverage",
                "alignment_ambiguous", "search_timed_out", "warnings",
            ))
            brief["observation_count"] = len(item.get("observations", []))
            brief["procedure_count"] = len(item.get("procedures", []))
    elif operation == "analyze_molecule":
        overview.update(_brief_fields(source, ("formal_charge", "stereocenters")))
        for key, fields in (("motifs", ("motif_id", "chemist_label")),
                            ("reactive_sites", ("chemist_label", "availability", "warnings"))):
            overview[f"{key}_count"] = len(source.get(key, []))
            overview[key] = [_brief_fields(item, fields) for item in source.get(key, [])[:3]]
    elif operation == "compare_molecules":
        overview.update(_brief_fields(source, (
            "same_constitution", "stereo_relationship", "core_method", "core_atom_count",
            "left_coverage", "right_coverage", "alignment_count_observed", "alignment_ambiguous",
            "embeddings_truncated", "alignments_truncated", "search_timed_out",
        )))
        for side in ("left", "right"):
            overview[side] = _brief_fields(source.get(side, {}), ("canonical_smiles", "warnings"))
    elif operation == "inspect_reactive_sites":
        overview.update(_brief_fields(source, ("selected_atom_ids", "radius")))
        overview["sites_count"] = len(source.get("sites", []))
        overview["other_sites_count"] = len(source.get("other_sites", []))
        overview["molecule"] = _brief_fields(source.get("molecule", {}), ("canonical_smiles", "warnings"))
        overview["sites"] = [
            _brief_fields(item, ("hypothesis_id", "chemist_label", "atom_indices", "availability", "warnings"))
            for item in source.get("sites", [])[:3] if isinstance(item, Mapping)
        ]
        overview["other_sites"] = [
            _brief_fields(item, ("hypothesis_id", "chemist_label", "atom_indices", "warnings"))
            for item in source.get("other_sites", [])[:3] if isinstance(item, Mapping)
        ]
        overview["site_profiles"] = []
        for item in source.get("environments", [])[:2]:
            profile = item.get("reactivity_profile", {})
            compact = _brief_fields(item, ("hypothesis_id", "center_atom_index"))
            compact.update(_brief_fields(profile, ("status", "context_kind", "flags")))
            for key, field in (("steric", "accessibility_class"), ("electronic", "activation_class")):
                descriptor = profile.get(key, {})
                compact[key] = _brief_fields(descriptor, (field,))
                if isinstance(descriptor.get("evidence"), Mapping):
                    compact[key]["evidence"] = _brief_fields(descriptor["evidence"], (
                        "source", "method", "confidence", "warnings",
                    ))
            overview["site_profiles"].append(compact)
    elif operation == "analyze_reaction":
        overview.update(_brief_fields(source, ("evidence_quality", "transformation_class", "named_family")))
        for key in ("observation", "interpretation", "reaction_completeness"):
            if isinstance(source.get(key), Mapping):
                overview[key] = _brief_fields(source[key], ("status", "valid", "evidence_quality", "warnings"))
        if isinstance(source.get("reaction_signature"), Mapping):
            signature = source["reaction_signature"]
            overview["bond_changes"] = _brief_fields(signature, ("formed_bond_types", "broken_bond_types"))
            overview["bond_changes"].update({
                f"{key}_count": len(signature.get(key, []))
                for key in ("order_changes", "hydrogen_changes", "stereo_changes")
                if key in signature
            })
    elif operation == "disconnect_target":
        strategies = source.get("strategies", [])
        overview["strategy_count"] = len(strategies)
        overview["strategies"] = []
        for strategy in strategies[:3]:
            item = _brief_fields(strategy, ("strategy_id", "strategy_rank", "independent_reference_support"))
            if isinstance(strategy.get("representative"), Mapping):
                item["representative"] = _brief_fields(strategy["representative"], (
                    "realization_id", "precursor_smiles", "forward_validation_status",
                    "precursor_compatibility_disposition", "reaction_compatibility_disposition",
                    "template_id", "operator_id", "precedent_reaction_ids",
                ))
                warnings = strategy["representative"].get("selectivity_warnings", [])
                if warnings:
                    item["representative"]["selectivity_warnings"] = [
                        _brief_fields(warning, ("code", "message", "conditions_evaluated"))
                        for warning in warnings[:2] if isinstance(warning, Mapping)
                    ]
                    item["representative"]["selectivity_warning_count"] = len(warnings)
            overview["strategies"].append(item)
    elif operation in {"assess_route_step", "inspect_route_step", "assess_route_proposal", "revise_route_branch"}:
        assessment = source.get("assessment")
        if isinstance(assessment, Mapping):
            overview["assessment"] = _brief_fields(assessment, (
                "status", "actionable", "admission_eligible", "strongest_evidence_tier",
                "warnings", "unresolved_step_ids", "disconnected_step_ids",
            ))
            gate_key = "topology_gates" if operation in {"assess_route_proposal", "revise_route_branch"} else "gates"
            gates = assessment.get(gate_key, [])
            non_pass_gates = [
                gate for gate in gates
                if isinstance(gate, Mapping) and gate.get("status") != "pass"
            ]
            non_pass_gates.sort(key=lambda gate: 0 if gate.get("status") in {
                "failed", "invalid", "conflicting", "unresolved", "unknown",
            } else 1)
            overview["assessment"]["non_pass_gates"] = [
                _brief_fields(gate, ("gate_id", "status", "summary", "warnings"))
                for gate in non_pass_gates[:4]
            ]
            overview["assessment"]["non_pass_gate_count"] = len(non_pass_gates)
            if "step_assessments" in assessment:
                steps = assessment["step_assessments"]
                overview["assessment"]["step_count"] = len(steps)
                if len(steps) > 10:
                    steps = sorted(steps, key=lambda step: (
                        (step.get("assessment") or {}).get("actionable") is not False
                        and not (step.get("assessment") or {}).get("warnings")
                    ))
                    overview["assessment"]["steps_preview_order"] = "non_actionable_or_warned_first"
                overview["assessment"]["steps"] = [
                    {"step_id": step.get("external_step_id"),
                     **(_brief_fields(step["assessment"], ("status", "actionable", "warnings"))
                        if isinstance(step.get("assessment"), Mapping) else {})}
                    for step in steps[:10] if isinstance(step, Mapping)
                ]
        if isinstance(source.get("material_constraints"), Mapping):
            overview["material_constraints"] = _brief_fields(source["material_constraints"], (
                "status", "stock_availability", "unavailable_starting_materials",
            ))
    elif operation == "recommend_conditions":
        recommendations = source.get("recommendations", [])
        overview["recommendation_count"] = len(recommendations)
        overview["recommendations"] = []
        for item in recommendations[:3]:
            if not isinstance(item, Mapping):
                continue
            compact = _brief_fields(item, (
                "rank", "recipe_id", "match_label", "retrieval_level", "compatibility_status",
                "reference_support", "observation_support", "precedent_reaction_ids", "cautions",
            ))
            if isinstance(item.get("resolved_recipe"), Mapping):
                compact["recipe"] = _brief_recipe(item["resolved_recipe"])
            overview["recommendations"].append(compact)
    elif operation == "resolve_recipe":
        overview.update(_brief_fields(source, ("recipe_id",)))
        overview["recipe"] = _brief_recipe(source)
    elif operation in {"get_precedents", "get_procedures"}:
        overview.update(_brief_fields(source, (
            "total", "offset", "next_offset", "missing_reaction_ids", "record_scope", "index_scope",
        )))
        records = source.get("records", [])
        overview["saved_record_count"] = len(records)
        overview["records"] = []
        for index, item in enumerate(records[:3]):
            compact = _brief_fields(item, (
                "reaction_id", "observation_id", "reference_id", "source_url", "source_uri",
                "chemistry_status", "condition_status", "condition_uncertain", "warnings",
            ))
            compact["saved_index"] = index
            if operation == "get_procedures":
                text = item.get("procedure_text")
                compact["procedure_characters"] = len(text) if isinstance(text, str) else 0
                compact["procedure_text_present"] = isinstance(text, str) and bool(text.strip())
            else:
                compact.update(_brief_fields(item, ("reaction_smiles", "recipe_id", "yield_pct")))
            overview["records"].append(compact)
    elif operation == "assess_recipe":
        overview.update(_brief_fields(source, ("compatible", "hard_conflicts", "unresolved_requirements")))
    elif operation == "assess_route_step_forward":
        if isinstance(source.get("assessment"), Mapping):
            overview["assessment"] = _brief_fields(source["assessment"], (
                "targeted_replay_status", "intended_match", "disposition", "validity",
                "advisory_only", "warnings",
            ))
        if isinstance(source.get("execution"), Mapping):
            execution = source["execution"]
            overview["execution"] = _brief_fields(execution, ("timings",))
            stages = execution.get("stages", [])
            overview["execution"]["stages"] = [
                _brief_fields(item, ("stage", "status", "duration_seconds")) for item in stages[-2:]
            ]
            overview["execution"]["stage_count"] = len(stages)
    elif operation == "compare_route_proposals":
        alternatives = source.get("alternatives", [])
        overview["alternative_count"] = len(alternatives)
        overview["alternatives"] = [
            _brief_fields(item, ("status", "source_ref", "route_id", "step_count",
                                 "unresolved_step_ids", "limitations"))
            for item in alternatives[:3] if isinstance(item, Mapping)
        ]
    elif operation == "propose_condition_adaptation":
        overview.update(_brief_fields(source, ("transfer_status", "assumptions", "risks")))
        if isinstance(source.get("compatibility"), Mapping):
            overview["compatibility"] = _brief_fields(source["compatibility"], (
                "status", "compatible", "hard_conflicts", "unresolved_requirements",
            ))
        overview["changes"] = [
            _brief_fields(item, ("field", "basis", "reason"))
            for item in source.get("changes", [])[:3] if isinstance(item, Mapping)
        ]
    elif operation == "search_fragment_precedents":
        overview.update(_brief_fields(source, ("counts", "refinement_hints", "output_truncated")))
        if isinstance(source.get("target_validation"), Mapping):
            overview["target_validation"] = _brief_fields(source["target_validation"], (
                "matches_target", "target_smiles", "query_id", "schema_version",
                "definition_version", "compiler_version",
            ))
        if isinstance(source.get("query"), Mapping):
            overview["query"] = _brief_fields(source["query"], ("expression", "query_format", "topology", "query_id"))
        overview["hits"] = [
            _brief_fields(item, ("hit_id", "product_smiles", "relationship_summary",
                                 "procedure_availability", "warnings"))
            for item in source.get("hits", [])[:3] if isinstance(item, Mapping)
        ]
    elif operation == "suggest_search_fragments":
        overview.update(_brief_fields(source, ("generation_truncated", "output_truncated")))
        overview["candidate_count"] = len(source.get("candidates", []))
        overview["candidates"] = [
            _brief_fields(item, ("candidate_id", "kind", "query", "reasons", "cautions"))
            for item in source.get("candidates", [])[:3] if isinstance(item, Mapping)
        ]
    elif operation == "inspect_condition_precedents":
        overview.update(_brief_fields(source, (
            "observation_count", "distinct_reference_count", "page", "missing_reference_observation_ids",
        )))
        if isinstance(source.get("procedure_catalog"), Mapping):
            overview["procedure_catalog"] = _brief_fields(source["procedure_catalog"], ("availability",))
        overview["precedents"] = []
        for item in source.get("precedents", [])[:2]:
            if not isinstance(item, Mapping):
                continue
            precedent = _brief_fields(item, ("missing_operating_fields", "procedure_link_scope"))
            for key, fields in (
                ("observation", ("observation_id", "reaction_id", "reference_id", "chemistry_status", "condition_status")),
                ("compatibility", ("status", "compatible", "hard_conflicts")),
            ):
                if isinstance(item.get(key), Mapping):
                    precedent[key] = _brief_fields(item[key], fields)
            observation = item.get("observation", {})
            if isinstance(observation.get("resolved_recipe"), Mapping):
                precedent["recipe"] = _brief_recipe(observation["resolved_recipe"])
            overview["precedents"].append(precedent)
    else:
        overview.update(_brief_fields(source, (
            "compatible", "hard_conflicts", "unresolved_requirements", "total", "next_offset",
            "observation_count", "distinct_reference_count", "transfer_status", "question",
        )))
    result["result_summary"] = overview
    result["inspection"] = {
        "result_path": "$.result", "projection_only": True,
        "detail_available": True,
        "count_scope": "saved_result",
    }
    if not overview and source:
        result["inspection"]["result_fields"] = [str(key)[:80] for key in islice(source, 20)]
    return result


def _brief_fields(
    source: Mapping[str, Any], keys: tuple[str, ...], *, _depth: int = 0,
) -> dict[str, Any]:
    """Bound nested values, keeping counts from the saved source and literal IDs."""
    selected: dict[str, Any] = {}
    for key in keys:
        if key not in source:
            continue
        value = source[key]
        if key == "error" and value is None:
            continue
        if key in {"warnings", "cautions", "limitations"} and value in ([], (), None):
            continue
        if isinstance(value, str):
            # Structure strings and evidence locators must stay usable. If even
            # this larger bound is exceeded, explicitly label the partial text.
            literal = key.endswith(("smiles", "_id", "_ref", "_url", "_uri")) or key in {"expression", "query"}
            limit = 2000 if literal else 180
            selected[key] = value[:limit] + ("…" if len(value) > limit else "")
            if len(value) > limit:
                selected[f"{key}_truncated"] = True
        elif isinstance(value, (list, tuple)):
            if key.endswith(("smiles", "_ids", "_refs")):
                selected[key] = [item[:2000] + ("…" if len(item) > 2000 else "")
                                 if isinstance(item, str) else _brief_value(item, _depth + 1)
                                 for item in value[:3]]
                if any(isinstance(item, str) and len(item) > 2000 for item in value[:3]):
                    selected[f"{key}_truncated"] = True
            else:
                selected[key] = [_brief_value(item, _depth + 1) for item in value[:3]]
            if len(value) > 3:
                selected[f"{key}_count"] = len(value)
        elif isinstance(value, Mapping):
            selected[key] = _brief_value(value, _depth + 1)
            if len(value) > 5:
                selected[f"{key}_field_count"] = len(value)
        else:
            selected[key] = value
    return selected


def _brief_value(value: Any, depth: int) -> Any:
    """Bound arbitrary nested payloads inside explicitly selected brief fields."""
    if isinstance(value, str):
        return value[:180] + ("…" if len(value) > 180 else "")
    if isinstance(value, (Mapping, list, tuple)) and depth >= 3:
        return {"summary_omitted": True, "item_count": len(value)}
    if isinstance(value, Mapping):
        keys = tuple(key for key in islice(value, 5) if isinstance(key, str) and len(key) <= 80)
        result = _brief_fields(value, keys, _depth=depth)
        if len(value) > len(keys):
            result["omitted_field_count"] = len(value) - len(keys)
        return result
    if isinstance(value, (list, tuple)):
        return {"preview": [_brief_value(item, depth + 1) for item in value[:3]], "item_count": len(value)}
    return value


def _brief_recipe(recipe: Mapping[str, Any]) -> dict[str, Any]:
    """Show ingredients and reported operating fields, leaving provenance on disk."""
    groups = ("catalysts", "ligands", "bases", "acids", "oxidants", "reductants",
              "condensation_agents", "additives", "solvents", "other_components")
    components = [(group, index, item) for group in groups
                  for index, item in enumerate(recipe.get(group, [])) if isinstance(item, Mapping)]
    # An unresolved ingredient late in a long recipe must not disappear behind
    # six ordinary components. This orders the preview, not the recipe itself.
    components.sort(key=lambda entry: entry[2].get("identity_status") == "resolved" and not entry[2].get("warnings"))
    result = _brief_fields(recipe, (
        "temperature_c", "time_h", "concentration_m", "atmosphere", "declared_absences", "warnings",
    ))
    result = {key: value for key, value in result.items() if value is not None and value not in ([], ())}
    result["component_count"] = len(components)
    result["components"] = []
    for group, index, component in components[:6]:
        compact = _brief_fields(component, (
            "primary_role", "identity_status", "role_status", "warnings",
        ))
        compact.update(_brief_fields({"name": component.get("canonical_name") or component.get("raw_identifier")}, ("name",)))
        compact["saved_path"] = [group, index]
        if component.get("amount") is not None:
            compact.update(_brief_fields(component, ("amount", "amount_unit", "quantity_status")))
        result["components"].append(compact)
    if len(components) > 6:
        result["omitted_component_count"] = len(components) - 6
    if recipe.get("stages"):
        result["stage_count"] = len(recipe["stages"])
        result["stage_details_omitted"] = True
    return result


def _context_fields(view: _Projection, item: Mapping[str, Any], fields: tuple[str, ...], path: str) -> dict[str, Any]:
    """Keep ancestor cautions visible without repeating surrounding result trees."""
    def shallow(value: Any, location: str) -> Any:
        if isinstance(value, (Mapping, list, tuple)):
            view.truncated(location, "ancestor_context_detail", total=len(value))
            return {"summary_omitted": True, "value_type": type(value).__name__, "total": len(value)}
        return view.value(value, location)

    result = {}
    for key in fields:
        if key not in item:
            continue
        value, location = item[key], f"{path}.{key}"
        if isinstance(value, (list, tuple)):
            result[key] = view.preview(value, location, shallow, limit=2, describe_collection=False)
        elif isinstance(value, Mapping) and key == "error":
            result[key] = view.pick(value, ("type", "message"), location)
            if set(value) - {"type", "message"}:
                view.truncated(location, "ancestor_context_detail", total_fields=len(value))
        else:
            result[key] = shallow(value, location)
    return result


def inspect_artifact_payload(
    payload: Any, *, artifact_ref: str, path: tuple[str | int, ...] | list[str | int] = (),
    offset: int = 0, limit: int = 5,
) -> dict[str, Any]:
    """Project one literal JSON path, retaining bounded context and pagination.

    Paths are key/index sequences, never expressions. Collection counts refer to
    this saved artifact; nested previews and context can themselves be truncated.
    """
    if (not isinstance(path, (tuple, list)) or len(path) > 12 or any(
            not (isinstance(part, str) and len(part) <= 80 or type(part) is int and part >= 0)
            for part in path)):
        raise ValueError("path must contain at most 12 literal keys (up to 80 characters) or nonnegative indices")
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 20:
        raise ValueError("offset must be nonnegative and limit must be between 1 and 20")
    value = payload
    selected_path = "$"
    ancestors = [(selected_path, value)]
    for part in path:
        if isinstance(value, Mapping) and isinstance(part, str):
            value = value[part]
        elif isinstance(value, list) and type(part) is int:
            value = value[part]
        else:
            raise ValueError(f"Path component has the wrong container type at {selected_path}")
        selected_path += f"[{json.dumps(part, ensure_ascii=False)}]"
        ancestors.append((selected_path, value))
    if not isinstance(value, (Mapping, list)) and offset:
        raise ValueError("offset is supported only for a selected list or mapping")

    view = _Projection()
    view.remaining_nodes = 120
    view.remaining_text = 4000
    context_view = _Projection()
    context_view.remaining_nodes = 80
    context_view.remaining_text = 1200
    context_view.text_limit = 200
    fields = tuple(dict.fromkeys((*_COMMON, "operation", "execution_status", "errors", "uncertainties",
                                 "evidence_refs", "source_ref", "retrieval_status", "acquisition", "claim_support",
                                 "compatible", "actionable", "admission_eligible", "cautions", "risks",
                                 "unresolved_requirements", "analysis_warnings", "evidence_warnings")))
    context = []
    for index, (location, item) in enumerate(ancestors):
        # A requested warning/error field belongs in the preview, not twice in
        # the context budget as well. Other ancestor cautions stay visible.
        selected_fields = tuple(key for key in fields if index >= len(path) or key != path[index])
        if isinstance(item, Mapping) and any(key in item for key in selected_fields):
            context.append({"path": location, "fields": _context_fields(context_view, item, selected_fields, location)})
    total = len(value) if isinstance(value, (Mapping, list)) else 1
    shown = min(limit, max(0, total - offset)) if isinstance(value, (Mapping, list)) else 1
    page = {"offset": offset, "limit": limit, "total": total, "shown": shown,
            "next_offset": offset + shown if offset + shown < total else None,
            "omitted_before": min(offset, total), "omitted_after": max(0, total - offset - shown),
            "count_scope": "saved_artifact", "value_type": type(value).__name__}
    if isinstance(value, Mapping):
        preview = {key: view.value(value[key], f"{selected_path}[{json.dumps(key, ensure_ascii=False)}]")
                   for key in islice(value, offset, offset + shown)}
    elif isinstance(value, list):
        preview = [view.value(item, f"{selected_path}[{index}]")
                   for index, item in enumerate(value[offset:offset + shown], offset)]
    else:
        preview = view.value(value, selected_path)
    if page["omitted_before"] or page["omitted_after"]:
        view.truncated(selected_path, "collection_page", **page)
    # Selection and ancestor context have separate budgets, so large surrounding
    # warnings cannot consume the requested field's preview allowance.
    view.collections.extend(context_view.collections)
    view.truncations.extend(context_view.truncations[:max(0, 30 - len(view.truncations))])
    view.truncation_count += context_view.truncation_count
    result = {
        "inspection_schema_version": "scientific_artifact_inspection.v1", "artifact_ref": artifact_ref,
        "path": list(path), "json_path": selected_path, "preview": preview, "page": page, "context": context,
        "inspection": {"projection_only": True, "collections": view.collections,
                       "truncations": view.truncations, "truncation_count": view.truncation_count,
                       "text_budget_characters": 4000, "serialized_byte_limit": 24000,
                       "context_text_budget_characters": 1200,
                       "hint": "Inspect a more specific path or the saved artifact for complete evidence and warnings. "
                               "Page counts describe saved collections, not the complete source dataset."},
    }
    # Long JSON keys and deeply nested paths can outweigh the text budget. Report
    # an omitted preview rather than silently returning an oversized console dump.
    size = lambda: len(json.dumps(result, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
    if size() > 23500:
        result["preview"] = {"summary_omitted": True}
        page.update({"selected": shown, "shown": 0, "next_offset": offset if offset < total else None,
                     "omitted_after": max(0, total - offset)})
        result["inspection"].update({
            "collections": [], "truncations": [{"path": selected_path, "reason": "serialized_preview_budget"}],
            "truncation_count": view.truncation_count + 1, "metadata_omitted": True,
        })
        while size() > 23500 and len(context) > 1:
            context.pop()
            result["inspection"]["context_omitted"] = True
        if size() > 23500 and context:
            context[0]["fields"] = {key: item for key, item in context[0]["fields"].items()
                                    if item is None or isinstance(item, (bool, int, float, str))}
            result["inspection"]["context_omitted"] = True
    result["inspection"]["truncated"] = bool(result["inspection"]["truncation_count"])
    return result

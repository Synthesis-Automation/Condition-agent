"""Bounded views of recorded scientific results, without chemistry interpretation.

Summaries are projections, not replacement contracts. Original status strings,
uncertainties and evidence identifiers are retained; full results remain in the
call artifact. Collection counts refer to that artifact, not an entire corpus.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from itertools import islice
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
            size = min(400, self.remaining_text)
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
            return self.preview(value, path, lambda item, child: self.value(item, child, depth + 1))
        self.truncated(path, "unexpected_value_type")
        return {"summary_omitted": True, "value_type": type(value).__name__}

    def pick(self, value: Any, fields: tuple[str, ...], path: str) -> Any:
        if not isinstance(value, Mapping):
            return self.value(value, path)
        return {key: self.value(value[key], f"{path}.{key}") for key in fields if key in value}

    def preview(
        self, value: Any, path: str,
        project: Callable[[Any, str], Any] | None = None, *, limit: int = 5,
    ) -> Any:
        if not isinstance(value, (list, tuple)):
            return self.value(value, path)
        shown = min(limit, len(value))
        if len(self.collections) < 40:
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
                    self.add_nested(result, item, "assessment", _STEP, child)
                return result
            if "step_assessments" in value:
                summary["step_assessments"] = self.preview(value["step_assessments"], f"{path}.step_assessments", step)
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
    if operation == "analyze_molecule":
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
    elif operation in {"assess_route_step", "inspect_route_step", "assess_route_proposal",
                       "prepare_route_proposal", "revise_route_branch"}:
        route = operation in {"assess_route_proposal", "prepare_route_proposal", "revise_route_branch"}
        summary.update(view.pick(value, (
            "source_ref", "step_id", "upstream_step_ids", "downstream_step_ids", "assessment_options", "evidence_refs",
        ), path))
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

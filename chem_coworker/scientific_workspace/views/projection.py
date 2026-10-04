"""Bounded result projections preserving statuses, evidence and truncation metadata."""

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
    "coverage", "reaction_completeness",
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

"""Read-only evidence projections from authoritative typed domain results."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Iterable, Mapping, Tuple

from condition_recommender.models import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
)
from core_retrosynthesis import collect_route_refinement_issues
from chem_coworker.contracts import (
    MultistepRetrosynthesisResponse,
    RetrosynthesisResponse,
)

from .contracts import EvidenceItem, stable_assistance_id


def _without_none(value: Mapping[str, Any]) -> Dict[str, Any]:
    return {key: item for key, item in value.items() if item is not None}


def _schema_versions(result: GenericRecommendationResult) -> Dict[str, str]:
    versions = {"recommendation_result": result.schema_version}
    if result.retrieval_definition_version:
        versions["retrieval_definition"] = result.retrieval_definition_version
    return versions


def _item(
    *,
    evidence_id: str,
    source_id: str,
    payload_type: str,
    payload: Mapping[str, Any],
    result: GenericRecommendationResult,
    uncertainty: str | None = None,
    layer: str = "recommendation",
) -> EvidenceItem:
    return EvidenceItem(
        evidence_id=evidence_id,
        layer=layer,  # type: ignore[arg-type]
        source_id=source_id,
        payload_type=payload_type,
        payload=payload,
        provenance="deterministic_inference",
        schema_versions=_schema_versions(result),
        uncertainty=uncertainty,
    )


def _candidate_items(
    result: GenericRecommendationResult,
    result_ref: str,
    alias: str,
    candidate: GenericConditionRecommendation,
) -> Tuple[EvidenceItem, ...]:
    base = f"{alias}"
    score_trace = asdict(candidate.score_trace)
    recipe = {
        "candidate_alias": alias,
        "recipe_id": candidate.recipe_id,
        "recipe_core_id": candidate.recipe_core_id,
        "recipe_variant_ids": list(candidate.recipe_variant_ids),
        "resolved_recipe": candidate.resolved_recipe,
        "synthesis_protocol": candidate.synthesis_protocol,
    }
    ranking = {
        "candidate_alias": alias,
        "rank": candidate.rank,
        "score": candidate.score,
        "similarity_score": candidate.similarity_score,
        "compatibility_score": candidate.compatibility_score,
        "expected_yield_pct": candidate.expected_yield_pct,
        "default_rank": candidate.default_rank,
        "default_score": candidate.default_score,
        "rank_change": candidate.rank_change,
        "retrieval_level": candidate.retrieval_level,
        "score_trace": score_trace,
        "factor_evidence": candidate.factor_evidence,
    }
    support = {
        "candidate_alias": alias,
        "support": candidate.support,
        "observation_support": candidate.observation_support,
        "reference_support": candidate.reference_support,
        "condition_series_support": candidate.condition_series_support,
        "dataset_support": candidate.dataset_support,
        "precedent_reaction_ids": list(candidate.precedent_reaction_ids),
        "precedent_reference_ids": list(candidate.precedent_reference_ids),
        "precedent_reaction_smiles": list(candidate.precedent_reaction_smiles),
        "precedent_reaction_contexts": list(candidate.precedent_reaction_contexts),
    }
    compatibility = {
        "candidate_alias": alias,
        "compatibility_score": candidate.compatibility_score,
        "compatibility_evidence": list(candidate.compatibility_evidence),
    }
    explanation = {
        "candidate_alias": alias,
        "explanation": list(candidate.explanation),
        "cautions": list(candidate.cautions),
    }
    return (
        _item(
            evidence_id=f"{base}.recipe",
            source_id=result_ref,
            payload_type="canonical_condition_recipe",
            payload=recipe,
            result=result,
        ),
        _item(
            evidence_id=f"{base}.ranking",
            source_id=result_ref,
            payload_type="recommendation_score_trace",
            payload=ranking,
            result=result,
        ),
        _item(
            evidence_id=f"{base}.support",
            source_id=result_ref,
            payload_type="precedent_support",
            payload=support,
            result=result,
        ),
        _item(
            evidence_id=f"{base}.compatibility",
            source_id=result_ref,
            payload_type="compatibility_evidence",
            payload=compatibility,
            result=result,
            uncertainty=(
                "Review listed compatibility cautions."
                if candidate.cautions
                else None
            ),
        ),
        _item(
            evidence_id=f"{base}.explanation",
            source_id=result_ref,
            payload_type="deterministic_explanation",
            payload=explanation,
            result=result,
            uncertainty="; ".join(candidate.cautions) or None,
        ),
    )


@dataclass(frozen=True)
class ConditionEvidenceProjection:
    """A compact initial packet plus lossless on-demand candidate evidence."""

    result_ref: str
    result_schema_version: str
    evidence: Tuple[EvidenceItem, ...]
    candidate_aliases: Tuple[Tuple[str, str], ...]

    @property
    def allowed_evidence_ids(self) -> frozenset[str]:
        return frozenset(item.evidence_id for item in self.evidence)

    def resolve_candidate(self, alias: str) -> str:
        """Resolve a request-local alias without accepting internal IDs."""

        matches = dict(self.candidate_aliases)
        if alias not in matches:
            raise ValueError(f"unknown condition candidate alias: {alias}")
        return matches[alias]

    def candidate_evidence(self, alias: str) -> Tuple[EvidenceItem, ...]:
        """Return the detailed evidence slice for one known candidate alias."""

        self.resolve_candidate(alias)
        prefix = f"{alias}."
        return tuple(item for item in self.evidence if item.evidence_id.startswith(prefix))

    def initial_packet(self) -> Dict[str, Any]:
        """Return bounded query and candidate-summary evidence for a first turn."""

        selected = []
        for item in self.evidence:
            if item.evidence_id.startswith("query.") or item.evidence_id.endswith(
                ".summary"
            ):
                selected.append(_evidence_dict(item))
        return {
            "result_ref": self.result_ref,
            "result_schema_version": self.result_schema_version,
            "candidate_aliases": [alias for alias, _ in self.candidate_aliases],
            "evidence": selected,
        }

    def inspection_packet(self, alias: str) -> Dict[str, Any]:
        """Return every projected field for one candidate."""

        return {
            "result_ref": self.result_ref,
            "candidate_alias": alias,
            "evidence": [
                _evidence_dict(item) for item in self.candidate_evidence(alias)
            ],
        }

    def comparison_packet(self, aliases: Iterable[str]) -> Dict[str, Any]:
        """Return existing evidence for two or more candidates without rescoring."""

        normalized = tuple(dict.fromkeys(aliases))
        if len(normalized) < 2:
            raise ValueError("comparison requires at least two distinct candidates")
        for alias in normalized:
            self.resolve_candidate(alias)
        evidence = tuple(
            item
            for alias in normalized
            for item in self.candidate_evidence(alias)
        )
        return {
            "result_ref": self.result_ref,
            "candidate_aliases": list(normalized),
            "evidence": [_evidence_dict(item) for item in evidence],
        }


def _evidence_dict(item: EvidenceItem) -> Dict[str, Any]:
    return {
        "evidence_id": item.evidence_id,
        "layer": item.layer,
        "payload_type": item.payload_type,
        "payload": dict(item.payload),
        "provenance": item.provenance,
        "uncertainty": item.uncertainty,
    }


def project_condition_result(
    result: GenericRecommendationResult,
    *,
    max_candidates: int = 20,
) -> ConditionEvidenceProjection:
    """Project a recommendation result without changing or reinterpreting it."""

    if max_candidates < 1:
        raise ValueError("max_candidates must be at least one")
    result_ref = stable_assistance_id("CONDRESULT", result.to_dict())
    query_items = (
        _item(
            evidence_id="query.reaction",
            source_id=result_ref,
            payload_type="reaction_query",
            payload=_without_none(
                {
                    "query_reaction_smiles": result.query_reaction_smiles,
                    "effective_query_reaction_smiles": (
                        result.effective_query_reaction_smiles
                    ),
                    "valid": result.valid,
                    "error": result.error,
                }
            ),
            result=result,
            uncertainty=result.error,
            layer="observation",
        ),
        _item(
            evidence_id="query.identity",
            source_id=result_ref,
            payload_type="reaction_identity",
            payload=_without_none(
                {
                    "query_signature_id": result.query_signature_id,
                    "query_reaction_core_id": result.query_reaction_core_id,
                    "query_fallback_descriptor_id": (
                        result.query_fallback_descriptor_id
                    ),
                    "query_edit_hypothesis_ids": list(
                        result.query_edit_hypothesis_ids
                    ),
                    "recommendation_mode": result.recommendation_mode,
                }
            ),
            result=result,
            layer="observation",
        ),
        _item(
            evidence_id="query.mapping",
            source_id=result_ref,
            payload_type="atom_mapping_status",
            payload=_without_none(
                {
                    "status": result.external_mapping_status,
                    "provider": result.external_mapping_provider,
                    "confidence": result.external_mapping_confidence,
                }
            ),
            result=result,
            uncertainty=(
                "External atom mapping was not available or was not accepted."
                if result.external_mapping_status not in {None, "accepted", "valid"}
                else None
            ),
            layer="observation",
        ),
        _item(
            evidence_id="query.interpretation",
            source_id=result_ref,
            payload_type="reaction_interpretation",
            payload=_without_none(
                {
                    "reaction_label": result.reaction_label,
                    "named_family": result.named_family,
                    "transformation_class": result.transformation_class,
                    "spectator_groups": list(result.spectator_groups),
                    "reaction_partners": list(result.reaction_partners),
                    "completion_proposal": result.completion_proposal,
                    "completion_selections": list(result.completion_selections),
                }
            ),
            result=result,
            layer="interpretation",
        ),
        _item(
            evidence_id="query.retrieval",
            source_id=result_ref,
            payload_type="retrieval_trace",
            payload={
                "retrieval_definition_version": result.retrieval_definition_version,
                "retrieval_strategy": result.retrieval_strategy,
                "retrieval_level": result.retrieval_level,
                "candidate_count": result.candidate_count,
                "independent_candidate_count": result.independent_candidate_count,
                "compatible_candidate_count": result.compatible_candidate_count,
                "independent_compatible_candidate_count": (
                    result.independent_compatible_candidate_count
                ),
                "excluded_candidate_count": result.excluded_candidate_count,
                "retrieval_trace": [asdict(item) for item in result.retrieval_trace],
                "ranking_preferences": result.ranking_preferences,
            },
            result=result,
        ),
        _item(
            evidence_id="query.warnings",
            source_id=result_ref,
            payload_type="domain_warnings",
            payload={"warnings": list(result.warnings)},
            result=result,
            uncertainty="; ".join(result.warnings) or None,
        ),
        _item(
            evidence_id="query.condition_constraints",
            source_id=result_ref,
            payload_type="condition_constraint_trace",
            payload={
                "condition_constraints": result.condition_constraints,
                "excluded_recommendations": list(result.condition_constraint_trace),
            },
            result=result,
        ),
    )
    candidate_items = []
    aliases = []
    for index, candidate in enumerate(result.recommendations[:max_candidates], start=1):
        alias = f"candidate-{index}"
        aliases.append((alias, candidate.recipe_id))
        candidate_items.append(
            _item(
                evidence_id=f"{alias}.summary",
                source_id=result_ref,
                payload_type="condition_candidate_summary",
                payload={
                    "candidate_alias": alias,
                    "rank": candidate.rank,
                    "score": candidate.score,
                    "similarity_score": candidate.similarity_score,
                    "compatibility_score": candidate.compatibility_score,
                    "expected_yield_pct": candidate.expected_yield_pct,
                    "support": candidate.support,
                    "reference_support": candidate.reference_support,
                    "retrieval_level": candidate.retrieval_level,
                    "cautions": list(candidate.cautions),
                },
                result=result,
                uncertainty="; ".join(candidate.cautions) or None,
            )
        )
        candidate_items.extend(_candidate_items(result, result_ref, alias, candidate))
    return ConditionEvidenceProjection(
        result_ref=result_ref,
        result_schema_version=result.schema_version,
        evidence=query_items + tuple(candidate_items),
        candidate_aliases=tuple(aliases),
    )


@dataclass(frozen=True)
class RetrosynthesisEvidenceProjection:
    """Compact and expandable evidence for validated one-step strategies."""

    result_ref: str
    evidence: Tuple[EvidenceItem, ...]
    strategy_aliases: Tuple[Tuple[str, str], ...]

    def strategy_evidence(self, alias: str) -> Tuple[EvidenceItem, ...]:
        aliases = dict(self.strategy_aliases)
        if alias not in aliases:
            raise ValueError(f"unknown retrosynthesis strategy alias: {alias}")
        return tuple(
            item for item in self.evidence if item.evidence_id.startswith(f"{alias}.")
        )

    def initial_packet(self) -> Dict[str, Any]:
        selected = tuple(
            item
            for item in self.evidence
            if item.evidence_id.startswith("query.")
            or item.evidence_id.endswith(".summary")
        )
        return _projection_packet(
            self.result_ref,
            "strategy_aliases",
            (alias for alias, _ in self.strategy_aliases),
            selected,
        )

    def inspection_packet(self, alias: str) -> Dict[str, Any]:
        return _projection_packet(
            self.result_ref,
            "strategy_aliases",
            (alias,),
            self.strategy_evidence(alias),
        )

    def comparison_packet(self, aliases: Iterable[str]) -> Dict[str, Any]:
        values = tuple(dict.fromkeys(aliases))
        if len(values) < 2:
            raise ValueError("comparison requires at least two distinct strategies")
        evidence = tuple(
            item for alias in values for item in self.strategy_evidence(alias)
        )
        return _projection_packet(
            self.result_ref,
            "strategy_aliases",
            values,
            evidence,
        )


@dataclass(frozen=True)
class MultistepEvidenceProjection:
    """Compact and expandable evidence that keeps partial routes explicit."""

    result_ref: str
    evidence: Tuple[EvidenceItem, ...]
    route_aliases: Tuple[Tuple[str, str], ...]

    def route_evidence(self, alias: str) -> Tuple[EvidenceItem, ...]:
        aliases = dict(self.route_aliases)
        if alias not in aliases:
            raise ValueError(f"unknown multistep route alias: {alias}")
        return tuple(
            item for item in self.evidence if item.evidence_id.startswith(f"{alias}.")
        )

    def initial_packet(self) -> Dict[str, Any]:
        selected = tuple(
            item
            for item in self.evidence
            if item.evidence_id.startswith("query.")
            or item.evidence_id.endswith(".summary")
        )
        return _projection_packet(
            self.result_ref,
            "route_aliases",
            (alias for alias, _ in self.route_aliases),
            selected,
        )

    def inspection_packet(self, alias: str, step_index: int | None = None) -> Dict[str, Any]:
        evidence = self.route_evidence(alias)
        if step_index is not None:
            evidence = tuple(
                item
                for item in evidence
                if item.evidence_id == f"{alias}.step-{step_index}"
                or (
                    item.payload_type == "route_refinement_issue"
                    and item.payload.get("step_index") == step_index
                )
            )
            if not evidence:
                raise ValueError(f"unknown route step index: {step_index}")
        return _projection_packet(
            self.result_ref,
            "route_aliases",
            (alias,),
            evidence,
        )

    def comparison_packet(self, aliases: Iterable[str]) -> Dict[str, Any]:
        values = tuple(dict.fromkeys(aliases))
        if len(values) < 2:
            raise ValueError("comparison requires at least two distinct routes")
        evidence = tuple(item for alias in values for item in self.route_evidence(alias))
        return _projection_packet(
            self.result_ref,
            "route_aliases",
            values,
            evidence,
        )


def _projection_packet(
    result_ref: str,
    alias_key: str,
    aliases: Iterable[str],
    evidence: Iterable[EvidenceItem],
) -> Dict[str, Any]:
    return {
        "result_ref": result_ref,
        alias_key: list(aliases),
        "evidence": [_evidence_dict(item) for item in evidence],
    }


def _application_evidence(
    evidence_id: str,
    result_ref: str,
    payload_type: str,
    payload: Mapping[str, Any],
    *,
    layer: str = "route",
    uncertainty: str | None = None,
) -> EvidenceItem:
    return EvidenceItem(
        evidence_id=evidence_id,
        layer=layer,  # type: ignore[arg-type]
        source_id=result_ref,
        payload_type=payload_type,
        payload=payload,
        provenance="deterministic_inference",
        schema_versions={
            "source": str(payload.get("schema_version") or "application.v1")
        },
        uncertainty=uncertainty,
    )


def project_retrosynthesis_response(
    response: RetrosynthesisResponse,
    *,
    max_strategies: int = 20,
) -> RetrosynthesisEvidenceProjection:
    """Project every returned strategy field without inventing chemistry."""

    snapshot = response.to_dict()
    result_ref = stable_assistance_id("RETRORESULT", snapshot)
    evidence = [
        _application_evidence(
            "query.retrosynthesis",
            result_ref,
            "retrosynthesis_result",
            {
                "target_smiles": response.request.target_smiles,
                "valid": response.valid,
                "error": response.error,
                "warnings": list(response.warnings),
                "library_path": response.library_path,
            },
            uncertainty="; ".join(response.warnings) or response.error,
        )
    ]
    condition_by_id = {
        item.strategy_id: item.evidence.to_dict()
        for item in response.condition_evidence
    }
    aliases = []
    for index, strategy in enumerate(response.strategies[:max_strategies], start=1):
        alias = f"strategy-{index}"
        aliases.append((alias, strategy.strategy_id))
        payload = strategy.to_dict()
        evidence.extend(
            (
                _application_evidence(
                    f"{alias}.summary",
                    result_ref,
                    "retrosynthesis_strategy_summary",
                    {
                        "strategy_alias": alias,
                        "strategy_rank": strategy.strategy_rank,
                        "representative_score": strategy.representative_score,
                        "independent_reference_support": (
                            strategy.independent_reference_support
                        ),
                        "total_realization_count": strategy.total_realization_count,
                        "forward_validation_status": (
                            strategy.representative.forward_validation_status
                        ),
                    },
                ),
                _application_evidence(
                    f"{alias}.strategy",
                    result_ref,
                    "retrosynthesis_strategy",
                    {"strategy_alias": alias, **payload},
                ),
                _application_evidence(
                    f"{alias}.conditions",
                    result_ref,
                    "retrosynthesis_condition_evidence",
                    {
                        "strategy_alias": alias,
                        "condition_evidence": condition_by_id.get(strategy.strategy_id),
                    },
                    uncertainty=(
                        None
                        if strategy.strategy_id in condition_by_id
                        else "No deterministic condition evidence is attached."
                    ),
                ),
            )
        )
    return RetrosynthesisEvidenceProjection(
        result_ref=result_ref,
        evidence=tuple(evidence),
        strategy_aliases=tuple(aliases),
    )


def project_multistep_response(
    response: MultistepRetrosynthesisResponse,
    *,
    max_routes: int = 10,
) -> MultistepEvidenceProjection:
    """Project solved routes or explicitly partial routes and diagnostics."""

    snapshot = response.to_dict()
    result_ref = stable_assistance_id("ROUTERESULT", snapshot)
    result = response.result
    route_kind = "none"
    routes = ()
    diagnostics: Mapping[str, Any] = {}
    if result is not None:
        if result.routes:
            route_kind = "solved"
            routes = result.routes
        else:
            route_kind = "partial"
            routes = result.partial_routes
        diagnostics = result.diagnostics.to_dict()
    evidence = [
        _application_evidence(
            "query.multistep",
            result_ref,
            "multistep_search_result",
            {
                "target_smiles": response.request.target_smiles,
                "valid": response.valid,
                "error": response.error,
                "route_kind": route_kind,
                "warnings": list(response.warnings),
                "search_diagnostics": diagnostics,
                "result_schema_version": (
                    result.schema_version if result is not None else None
                ),
            },
            uncertainty=(
                "No fully solved route was returned; all exposed routes are partial."
                if route_kind == "partial"
                else "; ".join(response.warnings) or response.error
            ),
        )
    ]
    aliases = []
    for index, route in enumerate(routes[:max_routes], start=1):
        alias = f"route-{index}"
        aliases.append((alias, route.route_id))
        route_issues = collect_route_refinement_issues(route)
        evidence.extend(
            (
                _application_evidence(
                    f"{alias}.summary",
                    result_ref,
                    "multistep_route_summary",
                    {
                        "route_alias": alias,
                        "solved": route.solved,
                        "route_cost": route.route_cost,
                        "reaction_count": route.reaction_count,
                        "maximum_depth": route.maximum_depth,
                        "evidence_summary": route.evidence_summary.to_dict(),
                        "refinement_issue_count": len(route_issues),
                        "strong_refinement_issue_count": sum(
                            item.severity == "strong" for item in route_issues
                        ),
                        "warnings": list(route.warnings),
                    },
                    uncertainty=(
                        "This is a partial route."
                        if not route.solved
                        else "; ".join(route.warnings) or None
                    ),
                ),
                _application_evidence(
                    f"{alias}.route",
                    result_ref,
                    "multistep_route",
                    {"route_alias": alias, **route.to_dict()},
                    uncertainty="This route remains partial." if not route.solved else None,
                ),
            )
        )
        for step_index, step in enumerate(route.steps, start=1):
            evidence.append(
                _application_evidence(
                    f"{alias}.step-{step_index}",
                    result_ref,
                    "multistep_route_step",
                    {
                        "route_alias": alias,
                        "step_index": step_index,
                        "step": step.to_dict(),
                        "route_solved": route.solved,
                    },
                    uncertainty="Parent route remains partial." if not route.solved else None,
                )
            )
        for issue_index, issue in enumerate(route_issues, start=1):
            evidence.append(
                _application_evidence(
                    f"{alias}.issue-{issue_index}",
                    result_ref,
                    "route_refinement_issue",
                    {
                        "route_alias": alias,
                        **issue.to_dict(),
                    },
                    uncertainty=(
                        "This issue is a deterministic search or evidence signal; "
                        "resolving it does not establish experimental feasibility."
                    ),
                )
            )
    return MultistepEvidenceProjection(
        result_ref=result_ref,
        evidence=tuple(evidence),
        route_aliases=tuple(aliases),
    )

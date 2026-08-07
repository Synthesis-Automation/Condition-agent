"""Public type-agnostic condition recommendation API."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

from reactive_taxonomy import (
    AtomMappingProvider,
    ExternalMappingAssessment,
    RxnMapperProvider,
    analyze_reaction_with_external_mapping,
    featurize_reaction,
)

from .generic_indexing import (
    GenericReactionIndex,
    load_generic_index,
)
from .core_retrieval import (
    load_reaction_core_retrieval_rules,
    retrieve_core_pool_with_trace,
)
from .generic_retrieval import (
    RetrievalStrategy,
    load_generic_retrieval_rules,
    retrieve_compatible_generic_pool_with_trace,
)
from .fallback_retrieval import retrieve_fallback_pool_with_trace
from .fallback_similarity import (
    assess_fallback_similarity,
    compatibility_signature_from_fallback,
    load_fallback_retrieval_rules,
)
from .hypothesis_retrieval import (
    assess_hypothesis_consensus_similarity,
    build_hypothesis_retrieval_query,
    load_edit_hypothesis_retrieval_rules,
    retrieve_hypothesis_consensus_pool_with_trace,
)
from .models import (
    ChemistRankingPreferences,
    FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION,
    GenericRecommendationResult,
    PrecedentIndexScope,
    ReactionCompletionSelection,
)
from .reaction_completion import (
    build_completed_reaction_smiles,
    build_reaction_completion_proposal,
    validate_completion_selections,
)
from .ranking_preferences import resolve_ranking_preferences
from .recipe_ranking import rank_condition_recipes


def _fragment_source_artifact_is_current(source: Path) -> bool | None:
    manifest_path = (
        source
        if source.name.casefold() == "shard_manifest.json"
        else source.parent / "shard_manifest.json"
    )
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    definition_contract = manifest.get("definition_contract") or {}
    return str(
        definition_contract.get(
            "fragment_source_capability_definition_version"
        )
        or ""
    ) == FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION


def _trusted_review_reuse_row_count(source: Path) -> int | None:
    """Return a verified trusted-row count when no review index is needed."""
    conversion_path = source.parent / "conversion_report.json"
    if conversion_path.is_file():
        try:
            conversion = json.loads(conversion_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            conversion = {}
        integrity = conversion.get("integrity") or {}
        tiers = conversion.get("precedent_tier_counts") or {}
        trusted_count = tiers.get("trusted")
        verified_rows = integrity.get("verified_row_count")
        output_rows = conversion.get("output_row_count")
        conversion_supports_reuse = (
            bool(integrity.get("valid"))
            and int(conversion.get("failed_shard_count") or 0) == 0
            and trusted_count is not None
            and int(tiers.get("review_core") or 0) == 0
            and verified_rows is not None
            and output_rows is not None
            and int(verified_rows) == int(output_rows)
        )
        if conversion_supports_reuse:
            return int(trusted_count)

    report_path = source.parent / "recommendation_artifacts_report.json"
    if not report_path.is_file():
        return None
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    trusted_count = report.get("trusted_precedent_count")
    unrestricted_count = report.get("unrestricted_precedent_count")
    counts_match = (
        trusted_count is not None
        and unrestricted_count is not None
        and int(trusted_count) == int(unrestricted_count)
    )
    fast_index_path = str(
        ((report.get("artifacts") or {}).get("fast_index") or {}).get("path")
        or ""
    ).replace("\\", "/")
    if (
        report.get("review_index_reuses_trusted")
        and int(report.get("review_core_precedent_count") or 0) == 0
        and counts_match
        and fast_index_path.endswith(source.name)
    ):
        return int(trusted_count)
    return None


def _reaction_label_payload(analysis: Any) -> Dict[str, Any]:
    """Serialize the one canonical rendered reaction label."""
    return (
        asdict(analysis.reaction_label)
        if analysis.reaction_label is not None
        else {}
    )


def recommend_indexed_signature(
    signature: Dict[str, Any],
    index: GenericReactionIndex,
    *,
    reaction_core: Mapping[str, Any] | None = None,
    query_reaction_smiles: str = "",
    reaction_label: Mapping[str, Any] | None = None,
    top_k: int = 5,
    minimum_pool_size: int | None = None,
    retrieval_strategy: RetrievalStrategy = "hybrid",
    ranking_weights: Mapping[str, float] | None = None,
    ranking_preferences: ChemistRankingPreferences | None = None,
) -> GenericRecommendationResult:
    """Recommend from an existing signature and index without re-featurization."""
    resolved_preferences = resolve_ranking_preferences(ranking_preferences)
    query_context = {
        "reaction_label": dict(reaction_label or {}),
        "spectator_groups": tuple(
            dict(group)
            for group in (signature.get("spectator_groups") or ())
            if isinstance(group, Mapping)
        ),
        "reaction_partners": tuple(
            dict(partner)
            for partner in (signature.get("partners") or ())
            if isinstance(partner, Mapping)
        ),
        "ranking_preferences": resolved_preferences.to_dict(),
    }
    retrieval_definition_version = str(load_generic_retrieval_rules()["schema_version"])
    if top_k < 1:
        return GenericRecommendationResult(
            query_reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE"
        )
    if not index.rows:
        return GenericRecommendationResult(
            query_reaction_smiles, False, error="EMPTY_GENERIC_INDEX"
        )
    if str(signature.get("schema_version") or "") != (
        index.reaction_signature_schema_version
    ):
        return GenericRecommendationResult(
            query_reaction_smiles,
            False,
            error="INCOMPATIBLE_REACTION_SIGNATURE_SCHEMA",
        )
    query_definitions = signature.get("definition_versions") or {}
    if (
        not isinstance(query_definitions, dict)
        or tuple(
            sorted((str(key), str(value)) for key, value in query_definitions.items())
        )
        != index.taxonomy_definition_versions
    ):
        return GenericRecommendationResult(
            query_reaction_smiles,
            False,
            error="INCOMPATIBLE_REACTION_TAXONOMY_DEFINITIONS",
        )
    retrieval = retrieve_compatible_generic_pool_with_trace(
        signature,
        index,
        minimum_pool_size=minimum_pool_size,
        reaction_core=reaction_core,
        strategy=retrieval_strategy,
    )
    level = retrieval.level
    compatible_pool = retrieval.pool
    if not compatible_pool:
        compatibility_failure = level == "no_compatible_condition_precedent"
        return GenericRecommendationResult(
            query_reaction_smiles=query_reaction_smiles,
            valid=False,
            query_signature_id=str(signature.get("signature_id") or ""),
            query_reaction_core_id=str(
                (reaction_core or {}).get("core_id") or ""
            )
            or None,
            named_family=signature.get("named_family"),
            transformation_class=signature.get("transformation_class"),
            **query_context,
            retrieval_definition_version=retrieval_definition_version,
            retrieval_strategy=retrieval_strategy,
            retrieval_level=level,
            candidate_count=retrieval.candidate_count,
            independent_candidate_count=retrieval.independent_candidate_count,
            excluded_candidate_count=retrieval.excluded_candidate_count,
            retrieval_trace=retrieval.trace,
            warnings=("ALL_RETRIEVED_RECIPES_FAILED_COMPATIBILITY",)
            if compatibility_failure
            else (),
            error="NO_COMPATIBLE_CONDITION_PRECEDENT"
            if compatibility_failure
            else "NO_CHEMICALLY_COMPATIBLE_PRECEDENT",
        )
    warnings = []
    if level.endswith("limited_support"):
        warnings.append("LIMITED_PRECEDENT_SUPPORT")
    if level not in {"exact_signature", "handle_signature", "named_family"}:
        warnings.append("TYPE_AGNOSTIC_FALLBACK_USED")
    if level.startswith("reaction_core_"):
        warnings.append("REACTION_CORE_RETRIEVAL_USED")
    if retrieval.excluded_candidate_count:
        warnings.append(
            f"INCOMPATIBLE_PRECEDENTS_EXCLUDED:{retrieval.excluded_candidate_count}"
        )
    recommendations = rank_condition_recipes(
        signature,
        compatible_pool,
        retrieval_level=level,
        top_k=top_k,
        ranking_profile=(
            "calibration_candidate"
            if ranking_weights is not None
            else retrieval_strategy
            if retrieval_strategy in {"transformation_prior", "legacy_pilot"}
            else "default"
        ),
        ranking_weights=ranking_weights,
        ranking_preferences=resolved_preferences,
        query_reaction_core=reaction_core,
    )
    if any(
        caution.startswith("Reaction-scope mismatch:")
        for recommendation in recommendations
        for caution in recommendation.cautions
    ):
        warnings.append("REACTION_TOPOLOGY_FALLBACK_USED")
    return GenericRecommendationResult(
        query_reaction_smiles=query_reaction_smiles,
        valid=True,
        query_signature_id=str(signature.get("signature_id") or ""),
        query_reaction_core_id=str(
            (reaction_core or {}).get("core_id") or ""
        )
        or None,
        named_family=signature.get("named_family"),
        transformation_class=signature.get("transformation_class"),
        **query_context,
        retrieval_definition_version=retrieval_definition_version,
        retrieval_strategy=retrieval_strategy,
        retrieval_level=level,
        candidate_count=retrieval.candidate_count,
        independent_candidate_count=retrieval.independent_candidate_count,
        compatible_candidate_count=len(compatible_pool),
        independent_compatible_candidate_count=(
            retrieval.independent_compatible_candidate_count
        ),
        excluded_candidate_count=retrieval.excluded_candidate_count,
        retrieval_trace=retrieval.trace,
        recommendations=recommendations,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True)
class GenericConditionRecommender:
    """Reusable recommender backed by one preloaded, validated index."""

    index: GenericReactionIndex
    source_path: str = ""
    mapping_provider: AtomMappingProvider | None = None
    includes_review_precedents: bool = False
    review_index_reuses_trusted: bool = False
    fragment_source_artifact_current: bool | None = None

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        *,
        mapping_provider: AtomMappingProvider | None = None,
        include_review: bool = False,
    ) -> "GenericConditionRecommender":
        source = Path(path)
        index_source = source
        paired_review_names = {
            "generic_index.sqlite": "generic_review_index.sqlite",
        }
        paired_review_name = paired_review_names.get(source.name.casefold())
        review_index_reuses_trusted = False
        review_reuse_row_count: int | None = None
        if include_review and paired_review_name:
            index_source = source.with_name(paired_review_name)
            if not index_source.is_file():
                review_reuse_row_count = _trusted_review_reuse_row_count(source)
                review_index_reuses_trusted = review_reuse_row_count is not None
                if review_reuse_row_count is not None:
                    index_source = source
                else:
                    raise FileNotFoundError(
                        "Review-core index is unavailable. Rebuild recommendation "
                        f"artifacts to create {paired_review_name} or record safe "
                        "trusted-index reuse."
                    )
        index = (
            load_generic_index(index_source, include_review=True)
            if include_review
            else load_generic_index(index_source)
        )
        expected_scope = (
            PrecedentIndexScope.TRUSTED_AND_REVIEW_CORE
            if include_review and not review_index_reuses_trusted
            else PrecedentIndexScope.TRUSTED
        )
        if index.precedent_scope != expected_scope:
            raise ValueError(
                f"Index scope is {index.precedent_scope.value!r}, expected "
                f"{expected_scope.value!r}; select the matching mode or rebuild "
                "recommendation artifacts."
            )
        if (
            review_reuse_row_count is not None
            and len(index.rows) != review_reuse_row_count
        ):
            raise ValueError(
                "Trusted-index review reuse row count does not match the "
                "artifact report; rebuild recommendation artifacts"
            )
        return cls(
            index=index,
            source_path=str(source),
            mapping_provider=mapping_provider,
            includes_review_precedents=include_review,
            review_index_reuses_trusted=review_index_reuses_trusted,
            fragment_source_artifact_current=(
                _fragment_source_artifact_is_current(source)
            ),
        )

    def recommend(
        self,
        reaction_smiles: str,
        *,
        top_k: int = 5,
        minimum_pool_size: int | None = None,
        unrestricted_fallback: bool = False,
        ranking_preferences: ChemistRankingPreferences | None = None,
        completion_selections: Tuple[ReactionCompletionSelection, ...] = (),
    ) -> GenericRecommendationResult:
        """Featurize a query and recommend without reloading the index."""
        resolved_preferences = resolve_ranking_preferences(ranking_preferences)
        result = _recommend_with_index(
            reaction_smiles,
            self.index,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
            unrestricted_fallback=unrestricted_fallback,
            mapping_provider=self.mapping_provider,
            ranking_preferences=resolved_preferences,
            completion_selections=completion_selections,
        )
        result = replace(
            result,
            ranking_preferences=resolved_preferences.to_dict(),
        )
        if (
            completion_selections
            and self.fragment_source_artifact_current is False
        ):
            result = replace(
                result,
                error=(
                    "RECOMMENDATION_INDEX_REBUILD_REQUIRED_FOR_COMPLETION"
                    if not result.valid
                    and result.error
                    in {
                        "NO_SAFE_FALLBACK_PRECEDENT",
                        "INDEX_HAS_NO_FALLBACK_DESCRIPTORS",
                    }
                    else result.error
                ),
                warnings=tuple(
                    dict.fromkeys(
                        (
                            *result.warnings,
                            "RECOMMENDATION_ARTIFACT_PREDATES_"
                            "FRAGMENT_SOURCE_COMPLETION",
                        )
                    )
                ),
            )
        if self.includes_review_precedents:
            review_warning = (
                "UNRESTRICTED_MODE_REUSES_TRUSTED_INDEX"
                if self.review_index_reuses_trusted
                else "UNRESTRICTED_CORE_REVIEW_INDEX_ENABLED"
            )
            result = replace(
                result,
                warnings=tuple(
                    dict.fromkeys(
                        (
                            *result.warnings,
                            review_warning,
                        )
                    )
                ),
            )
        return result


def _recommend_with_index(
    reaction_smiles: str,
    index: GenericReactionIndex,
    *,
    top_k: int,
    minimum_pool_size: int | None,
    unrestricted_fallback: bool,
    mapping_provider: AtomMappingProvider | None = None,
    ranking_preferences: ChemistRankingPreferences | None = None,
    completion_selections: Tuple[ReactionCompletionSelection, ...] = (),
) -> GenericRecommendationResult:
    if top_k < 1:
        return GenericRecommendationResult(
            reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE"
        )
    base_analysis = featurize_reaction(reaction_smiles)
    completion_proposal = (
        build_reaction_completion_proposal(base_analysis)
        if base_analysis.valid
        else None
    )
    if completion_proposal is not None:
        validate_completion_selections(
            completion_proposal,
            completion_selections,
        )
    effective_reaction_smiles, completion_build_warnings = (
        build_completed_reaction_smiles(
            reaction_smiles,
            completion_selections,
        )
    )
    analysis_reaction_smiles = effective_reaction_smiles or reaction_smiles
    query_analysis = (
        featurize_reaction(analysis_reaction_smiles)
        if effective_reaction_smiles is not None
        else base_analysis
    )
    mapping_skip_warning = None
    fallback_descriptor = query_analysis.fallback_descriptor
    if (
        mapping_provider is not None
        and fallback_descriptor is not None
        and query_analysis.partial_product_transformation is not None
        and fallback_descriptor.retrieval_eligible
    ):
        mapping_skip_warning = (
            "EXTERNAL_MAPPING_SKIPPED:PARTIAL_PRODUCT_TRANSFORMATION_AVAILABLE"
        )
    elif (
        mapping_provider is not None
        and fallback_descriptor is not None
        and "incomplete_product_atom_provenance"
        in fallback_descriptor.ineligibility_reasons
    ):
        # Atom correspondence cannot supply elements that are absent from the
        # input. Avoid starting an external mapping model for a query that the
        # chemistry gates must reject regardless of the proposed mapping.
        mapping_skip_warning = (
            "EXTERNAL_MAPPING_SKIPPED:INCOMPLETE_PRODUCT_ATOM_PROVENANCE"
        )
    assessment = (
        analyze_reaction_with_external_mapping(
            analysis_reaction_smiles,
            mapping_provider,
            base_analysis=query_analysis,
        )
        if (
            mapping_provider is not None
            and query_analysis.valid
            and mapping_skip_warning is None
        )
        else None
    )

    def finalize(result: GenericRecommendationResult) -> GenericRecommendationResult:
        attached = _attach_external_mapping_assessment(result, assessment)
        warnings = list(attached.warnings)
        warnings.extend(completion_build_warnings)
        if mapping_skip_warning is not None:
            warnings.append(mapping_skip_warning)
        if completion_proposal is not None and completion_proposal.requirements:
            warnings.append("QUERY_PRODUCT_SOURCE_COMPONENT_MISSING")
            if not completion_selections:
                warnings.append(
                    "SYSTEM_PROPOSED_SOURCES_REQUIRE_USER_CONFIRMATION"
                )
        selection_cautions = []
        for selection in completion_selections:
            if selection.selection_kind == "compatible_source_class":
                warnings.append(
                    "QUERY_SOURCE_CLASS_USER_CONFIRMED:"
                    f"{selection.capability_id}"
                )
            elif selection.selection_kind == "registered_substance":
                warnings.append(
                    "QUERY_SOURCE_SUBSTANCE_USER_CONFIRMED:"
                    f"{selection.substance_id}"
                )
            elif selection.selection_kind == "custom_identifier":
                warnings.append("QUERY_SOURCE_EDIT_UNRESOLVED")
            else:
                warnings.append("QUERY_SOURCE_LEFT_UNRESOLVED")
        if completion_selections:
            warnings.append("REACTION_COMPLETION_SELECTION_RECORDED")
            if effective_reaction_smiles is not None:
                warnings.append("REACTION_COMPLETION_SHADOW_QUERY_USED")
            selection_cautions.extend(
                (
                    "A missing condition source was confirmed or edited by the "
                    "user; it was not observed in the submitted reaction",
                    "The submitted reaction SMILES was retained unchanged",
                )
            )
        recommendations = tuple(
            replace(
                recommendation,
                cautions=tuple(
                    dict.fromkeys(
                        (*selection_cautions, *recommendation.cautions)
                    )
                ),
            )
            for recommendation in attached.recommendations
        )
        return replace(
            attached,
            query_reaction_smiles=reaction_smiles,
            effective_query_reaction_smiles=effective_reaction_smiles,
            completion_proposal=(
                completion_proposal.to_dict()
                if completion_proposal is not None
                and completion_proposal.requirements
                else None
            ),
            completion_selections=tuple(
                selection.to_dict() for selection in completion_selections
            ),
            recommendations=recommendations,
            warnings=tuple(dict.fromkeys(warnings)),
        )

    analysis = assessment.analysis if assessment is not None else query_analysis
    if not analysis.valid:
        return GenericRecommendationResult(
            reaction_smiles,
            False,
            error=analysis.error or "INVALID_REACTION",
        )
    if (
        assessment is not None
        and assessment.status
        in {"external_mapping_internal_consensus", "external_mapping_only"}
        and analysis.reaction_core is not None
    ):
        mapped_core_result = _recommend_core_with_index(
            analysis,
            index,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
            ranking_preferences=ranking_preferences,
        )
        if mapped_core_result.valid:
            return finalize(
                mapped_core_result,
            )
    if analysis.reaction_signature is None:
        core_attempt = None
        if analysis.edit_hypotheses:
            hypothesis_result = _recommend_hypotheses_with_index(
                analysis,
                index,
                top_k=top_k,
                minimum_pool_size=minimum_pool_size,
                ranking_preferences=ranking_preferences,
            )
            if (
                hypothesis_result.valid
                or hypothesis_result.error
                != "QUERY_EDIT_HYPOTHESES_NOT_RETRIEVABLE"
            ):
                return finalize(
                    hypothesis_result,
                )
        if (
            analysis.reaction_core is not None
            and analysis.partial_product_transformation is None
        ):
            core_attempt = _recommend_core_with_index(
                analysis,
                index,
                top_k=top_k,
                minimum_pool_size=minimum_pool_size,
                ranking_preferences=ranking_preferences,
            )
            if core_attempt.valid:
                return finalize(
                    core_attempt,
                )
        result = _recommend_fallback_with_index(
            analysis,
            index,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
            unrestricted=unrestricted_fallback,
            ranking_preferences=ranking_preferences,
            completion_selections=completion_selections,
        )
        if core_attempt is not None:
            result = replace(
                result,
                retrieval_trace=(
                    *core_attempt.retrieval_trace,
                    *result.retrieval_trace,
                ),
                warnings=tuple(
                    dict.fromkeys(
                        (
                            *result.warnings,
                            "REACTION_CORE_RETRIEVAL_ATTEMPTED",
                            f"REACTION_CORE_RETRIEVAL_RESULT:"
                            f"{core_attempt.retrieval_level}",
                        )
                    )
                ),
            )
    else:
        result = recommend_indexed_signature(
            asdict(analysis.reaction_signature),
            index,
            reaction_core=(
                asdict(analysis.reaction_core)
                if analysis.reaction_core is not None
                else None
            ),
            query_reaction_smiles=reaction_smiles,
            reaction_label=_reaction_label_payload(analysis),
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
            ranking_preferences=ranking_preferences,
        )
        result = replace(
            result,
            spectator_groups=tuple(
                asdict(group) for group in analysis.spectator_groups
            ),
        )
    return finalize(result)


def _recommend_core_with_index(
    analysis: Any,
    index: GenericReactionIndex,
    *,
    top_k: int,
    minimum_pool_size: int | None,
    ranking_preferences: ChemistRankingPreferences | None = None,
) -> GenericRecommendationResult:
    """Recommend from a review-qualified core against verified precedents."""
    core_model = analysis.reaction_core
    if core_model is None:
        return GenericRecommendationResult(
            query_reaction_smiles=str(analysis.input_reaction_smiles),
            valid=False,
            recommendation_mode="abstained",
            error="QUERY_REACTION_CORE_NOT_RETRIEVABLE",
        )
    core = asdict(core_model)
    fallback_model = analysis.fallback_descriptor
    fallback = asdict(fallback_model) if fallback_model is not None else {}
    compatibility_signature = compatibility_signature_from_fallback(fallback)
    retrieval = retrieve_core_pool_with_trace(
        core,
        compatibility_signature,
        index,
        minimum_pool_size=minimum_pool_size,
    )
    rules = load_reaction_core_retrieval_rules()
    retrieval_definition_version = (
        f"{rules['definition_id']}@{rules['schema_version']}"
    )
    if not retrieval.pool:
        return GenericRecommendationResult(
            query_reaction_smiles=str(analysis.input_reaction_smiles),
            valid=False,
            query_reaction_core_id=str(core.get("core_id") or "") or None,
            recommendation_mode="reaction_core_review",
            reaction_label=_reaction_label_payload(analysis),
            transformation_class=analysis.transformation_class,
            retrieval_definition_version=retrieval_definition_version,
            retrieval_strategy="reaction_core_ladder",
            retrieval_level=retrieval.level,
            candidate_count=retrieval.candidate_count,
            independent_candidate_count=retrieval.independent_candidate_count,
            excluded_candidate_count=retrieval.excluded_candidate_count,
            retrieval_trace=retrieval.trace,
            warnings=(
                "QUERY_TRANSFORMATION_NOT_VERIFIED",
                "REACTION_CORE_RETRIEVAL_ATTEMPTED",
            ),
            error="QUERY_REACTION_CORE_NOT_RETRIEVABLE",
        )

    query = {
        "reaction_core": core,
        "fallback_descriptor": fallback,
        "spectator_groups": tuple(
            asdict(group) for group in analysis.spectator_groups
        ),
    }

    def core_similarity(
        query_value: Mapping[str, Any],
        row: Any,
    ) -> Any:
        return assess_fallback_similarity(
            query_value["fallback_descriptor"],
            row.fallback_descriptor,
        )

    recommendations = rank_condition_recipes(
        query,
        retrieval.pool,
        retrieval_level=retrieval.level,
        top_k=top_k,
        similarity_assessor=core_similarity,
        query_reaction_core=core,
        ranking_preferences=ranking_preferences,
    )
    cautions = (
        "Query retrieval used minimized reaction-core evidence, not a verified "
        "query reaction signature",
        "Only independently admitted trusted or explicit review-core "
        "precedents supplied conditions",
        "The center-transition key alone was not used for retrieval",
        "Reaction-core recommendations require expert review",
    )
    recommendations = tuple(
        replace(
            recommendation,
            cautions=tuple(
                dict.fromkeys((*cautions, *recommendation.cautions))
            ),
        )
        for recommendation in recommendations
    )
    warnings = [
        "QUERY_TRANSFORMATION_NOT_VERIFIED",
        "REACTION_CORE_RETRIEVAL_USED",
        "RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW",
    ]
    if retrieval.level.endswith("limited_support"):
        warnings.append("LIMITED_PRECEDENT_SUPPORT")
    if retrieval.excluded_candidate_count:
        warnings.append(
            f"INCOMPATIBLE_PRECEDENTS_EXCLUDED:"
            f"{retrieval.excluded_candidate_count}"
        )
    return GenericRecommendationResult(
        query_reaction_smiles=str(analysis.input_reaction_smiles),
        valid=True,
        query_reaction_core_id=str(core.get("core_id") or "") or None,
        recommendation_mode="reaction_core_review",
        reaction_label=_reaction_label_payload(analysis),
        transformation_class=analysis.transformation_class,
        spectator_groups=tuple(
            asdict(group) for group in analysis.spectator_groups
        ),
        retrieval_definition_version=retrieval_definition_version,
        retrieval_strategy="reaction_core_ladder",
        retrieval_level=retrieval.level,
        candidate_count=retrieval.candidate_count,
        independent_candidate_count=retrieval.independent_candidate_count,
        compatible_candidate_count=len(retrieval.pool),
        independent_compatible_candidate_count=(
            retrieval.independent_compatible_candidate_count
        ),
        excluded_candidate_count=retrieval.excluded_candidate_count,
        retrieval_trace=retrieval.trace,
        recommendations=recommendations,
        warnings=tuple(warnings),
    )


def _attach_external_mapping_assessment(
    result: GenericRecommendationResult,
    assessment: ExternalMappingAssessment | None,
) -> GenericRecommendationResult:
    """Expose mapper provenance and mandatory review cautions on query results."""
    if assessment is None or assessment.mapping_result is None:
        return result
    mapping = assessment.mapping_result
    warnings = list(result.warnings)
    warnings.extend(assessment.warnings)
    warnings.append(f"EXTERNAL_MAPPING_STATUS:{assessment.status}")
    if assessment.status in {
        "external_mapping_internal_consensus",
        "external_mapping_only",
    }:
        warnings.extend(
            (
                "QUERY_SIGNATURE_USES_EXTERNAL_MAPPING",
                "RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW",
            )
        )
    cautions = (
        "Query atom correspondence was generated by an external mapper, not "
        "supplied with the source reaction",
        (
            "RXNMapper selected one internally enumerated edit hypothesis"
            if assessment.status == "external_mapping_internal_consensus"
            else "External mapping has no unique internal correspondence consensus"
        ),
        "Mapper-supported recommendations require expert review",
    )
    recommendations = (
        tuple(
            replace(
                recommendation,
                cautions=tuple(
                    dict.fromkeys((*cautions, *recommendation.cautions))
                ),
            )
            for recommendation in result.recommendations
        )
        if assessment.status
        in {"external_mapping_internal_consensus", "external_mapping_only"}
        else result.recommendations
    )
    return replace(
        result,
        external_mapping_status=assessment.status,
        external_mapping_provider=mapping.metadata.provider_id,
        external_mapping_confidence=mapping.mapper_confidence,
        recommendation_mode=(
            assessment.status
            if assessment.status
            in {"external_mapping_internal_consensus", "external_mapping_only"}
            else result.recommendation_mode
        ),
        recommendations=recommendations,
        warnings=tuple(dict.fromkeys(warnings)),
    )


def _recommend_hypotheses_with_index(
    analysis: Any,
    index: GenericReactionIndex,
    *,
    top_k: int,
    minimum_pool_size: int | None,
    ranking_preferences: ChemistRankingPreferences | None = None,
) -> GenericRecommendationResult:
    """Recommend only from precedents robust across every retained edit hypothesis."""
    reaction_smiles = str(analysis.input_reaction_smiles)
    hypotheses = tuple(asdict(value) for value in analysis.edit_hypotheses)
    fallback_model = analysis.fallback_descriptor
    fallback = asdict(fallback_model) if fallback_model is not None else {}
    query = build_hypothesis_retrieval_query(hypotheses, fallback)
    hypothesis_ids = tuple(
        str(value.get("hypothesis_id") or "") for value in hypotheses
    )
    if query is None:
        return GenericRecommendationResult(
            query_reaction_smiles=reaction_smiles,
            valid=False,
            query_edit_hypothesis_ids=hypothesis_ids,
            recommendation_mode="abstained",
            warnings=("AMBIGUOUS_EDIT_HYPOTHESES_RETAINED",),
            error="QUERY_EDIT_HYPOTHESES_NOT_RETRIEVABLE",
        )
    rules = load_edit_hypothesis_retrieval_rules()
    retrieval = retrieve_hypothesis_consensus_pool_with_trace(
        query,
        index,
        minimum_pool_size=minimum_pool_size,
    )
    base_warnings = [
        "QUERY_TRANSFORMATION_NOT_VERIFIED",
        "AMBIGUOUS_EDIT_HYPOTHESES_RETAINED",
        "EDIT_HYPOTHESIS_CONSENSUS_REQUIRED",
    ]
    if not retrieval.pool:
        compatibility_failure = (
            retrieval.level == "no_compatible_condition_precedent"
        )
        return GenericRecommendationResult(
            query_reaction_smiles=reaction_smiles,
            valid=False,
            query_edit_hypothesis_ids=hypothesis_ids,
            recommendation_mode="ambiguous_edit_hypotheses",
            reaction_label=_reaction_label_payload(analysis),
            transformation_class=analysis.transformation_class,
            retrieval_definition_version=str(rules["schema_version"]),
            retrieval_strategy="edit_hypothesis_consensus",
            retrieval_level=retrieval.level,
            candidate_count=retrieval.candidate_count,
            independent_candidate_count=retrieval.independent_candidate_count,
            excluded_candidate_count=retrieval.excluded_candidate_count,
            retrieval_trace=retrieval.trace,
            warnings=tuple(
                base_warnings
                + (
                    ["ALL_RETRIEVED_RECIPES_FAILED_COMPATIBILITY"]
                    if compatibility_failure
                    else []
                )
            ),
            error=(
                "NO_COMPATIBLE_CONDITION_PRECEDENT"
                if compatibility_failure
                else "NO_ROBUST_EDIT_HYPOTHESIS_PRECEDENT"
            ),
        )
    query_mapping = {
        **query.to_mapping(),
        "spectator_groups": tuple(
            asdict(group) for group in analysis.spectator_groups
        ),
    }
    recommendations = rank_condition_recipes(
        query_mapping,
        retrieval.pool,
        retrieval_level=retrieval.level,
        top_k=top_k,
        similarity_assessor=assess_hypothesis_consensus_similarity,
        query_reaction_core=(
            asdict(analysis.reaction_core)
            if analysis.reaction_core is not None
            else None
        ),
        ranking_preferences=ranking_preferences,
    )
    hypothesis_cautions = (
        "Query atom correspondence is ambiguous; no hypothesis is presented "
        "as the verified reaction center",
        "Every returned precedent passed the anonymous edit-graph threshold "
        "for every retained query hypothesis",
        "Recipe ranking uses the worst-case edit similarity across hypotheses",
        "Ambiguous query records remain excluded from the verified precedent index",
    )
    recommendations = tuple(
        replace(
            recommendation,
            cautions=tuple(
                dict.fromkeys(
                    (*hypothesis_cautions, *recommendation.cautions)
                )
            ),
        )
        for recommendation in recommendations
    )
    warnings = [
        *base_warnings,
        "AMBIGUOUS_EDIT_HYPOTHESIS_RETRIEVAL_USED",
        "RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW",
    ]
    if retrieval.excluded_candidate_count:
        warnings.append(
            f"INCOMPATIBLE_PRECEDENTS_EXCLUDED:{retrieval.excluded_candidate_count}"
        )
    return GenericRecommendationResult(
        query_reaction_smiles=reaction_smiles,
        valid=True,
        query_edit_hypothesis_ids=hypothesis_ids,
        recommendation_mode="ambiguous_edit_hypotheses",
        reaction_label=_reaction_label_payload(analysis),
        transformation_class=analysis.transformation_class,
        retrieval_definition_version=str(rules["schema_version"]),
        retrieval_strategy="edit_hypothesis_consensus",
        retrieval_level=retrieval.level,
        candidate_count=retrieval.candidate_count,
        independent_candidate_count=retrieval.independent_candidate_count,
        compatible_candidate_count=len(retrieval.pool),
        independent_compatible_candidate_count=(
            retrieval.independent_compatible_candidate_count
        ),
        excluded_candidate_count=retrieval.excluded_candidate_count,
        retrieval_trace=retrieval.trace,
        recommendations=recommendations,
        warnings=tuple(warnings),
    )


def _recommend_fallback_with_index(
    analysis: Any,
    index: GenericReactionIndex,
    *,
    top_k: int,
    minimum_pool_size: int | None,
    unrestricted: bool = False,
    ranking_preferences: ChemistRankingPreferences | None = None,
    completion_selections: Tuple[ReactionCompletionSelection, ...] = (),
) -> GenericRecommendationResult:
    """Recommend through the gated or explicit unrestricted fallback route."""
    reaction_smiles = str(analysis.input_reaction_smiles)
    descriptor_model = analysis.fallback_descriptor
    if descriptor_model is None:
        return GenericRecommendationResult(
            reaction_smiles,
            False,
            recommendation_mode="abstained",
            error="QUERY_HAS_NO_USABLE_REACTION_SIGNATURE",
        )
    descriptor = asdict(descriptor_model)
    descriptor_id = str(descriptor.get("descriptor_id") or "")
    retrieval_definition_version = str(
        load_fallback_retrieval_rules()["schema_version"]
    )
    base_warnings = [
        "QUERY_TRANSFORMATION_NOT_VERIFIED",
        "UNVERIFIED_REACTION_FALLBACK_CONSIDERED",
    ]
    if unrestricted:
        base_warnings.append("UNRESTRICTED_FALLBACK_REQUESTED")
    required_source_elements = tuple(
        str(value)
        for value in descriptor.get("required_condition_source_elements") or ()
    )
    if required_source_elements:
        base_warnings.append(
            f"QUERY_PRODUCT_ATOM_SOURCE_UNVERIFIED:{','.join(required_source_elements)}"
        )
    if not bool(descriptor.get("retrieval_eligible")) and not unrestricted:
        base_warnings.extend(
            f"FALLBACK_BLOCKED:{reason}"
            for reason in descriptor.get("ineligibility_reasons") or ()
        )
        return GenericRecommendationResult(
            query_reaction_smiles=reaction_smiles,
            valid=False,
            query_fallback_descriptor_id=descriptor_id,
            recommendation_mode="abstained",
            reaction_label=_reaction_label_payload(analysis),
            transformation_class=analysis.transformation_class,
            warnings=tuple(base_warnings),
            error="QUERY_NOT_ELIGIBLE_FOR_UNVERIFIED_FALLBACK",
        )
    if unrestricted:
        base_warnings.extend(
            f"FALLBACK_GATE_OVERRIDDEN:{reason}"
            for reason in descriptor.get("ineligibility_reasons") or ()
        )
    if not index.fallback_features:
        return GenericRecommendationResult(
            query_reaction_smiles=reaction_smiles,
            valid=False,
            query_fallback_descriptor_id=descriptor_id,
            recommendation_mode="abstained",
            reaction_label=_reaction_label_payload(analysis),
            transformation_class=analysis.transformation_class,
            warnings=tuple(
                base_warnings + ["GENERIC_INDEX_REQUIRES_FALLBACK_FEATURE_REBUILD"]
            ),
            error="INDEX_HAS_NO_FALLBACK_DESCRIPTORS",
        )
    retrieval = retrieve_fallback_pool_with_trace(
        descriptor,
        index,
        minimum_pool_size=minimum_pool_size,
        unrestricted=unrestricted,
        required_source_capability_ids={
            selection.requirement_id: str(selection.capability_id)
            for selection in completion_selections
            if selection.selection_kind
            in {"compatible_source_class", "registered_substance"}
            and selection.capability_id
        },
        required_source_substance_ids={
            selection.requirement_id: str(selection.substance_id)
            for selection in completion_selections
            if selection.selection_kind == "registered_substance"
            and selection.substance_id
        },
    )
    if not retrieval.pool:
        return GenericRecommendationResult(
            query_reaction_smiles=reaction_smiles,
            valid=False,
            query_fallback_descriptor_id=descriptor_id,
            recommendation_mode="abstained",
            reaction_label=_reaction_label_payload(analysis),
            transformation_class=analysis.transformation_class,
            retrieval_definition_version=retrieval_definition_version,
            retrieval_strategy="unverified_structure_fallback",
            retrieval_level=retrieval.level,
            candidate_count=retrieval.candidate_count,
            independent_candidate_count=retrieval.independent_candidate_count,
            excluded_candidate_count=retrieval.excluded_candidate_count,
            retrieval_trace=retrieval.trace,
            warnings=tuple(base_warnings),
            error=(
                "NO_UNRESTRICTED_FALLBACK_PRECEDENT"
                if unrestricted
                else "NO_SAFE_FALLBACK_PRECEDENT"
            ),
        )
    ranking_query = {
        **descriptor,
        "spectator_groups": tuple(
            asdict(group) for group in analysis.spectator_groups
        ),
    }
    recommendations = rank_condition_recipes(
        ranking_query,
        retrieval.pool,
        retrieval_level=retrieval.level,
        top_k=top_k,
        similarity_assessor=(
            lambda query, row: assess_fallback_similarity(
                query,
                row.fallback_descriptor,
            )
        ),
        query_reaction_core=(
            asdict(analysis.reaction_core)
            if analysis.reaction_core is not None
            else None
        ),
        ranking_preferences=ranking_preferences,
    )
    fallback_cautions = (
        "Query atom correspondence and bond edits are not verified; the "
        "recommendation is an analogy, not a chemistry-confirmed match",
        "Fallback similarity uses taxonomy features, candidate hypotheses, "
        "and global structure inventories",
        *(
            (
                "Fallback eligibility, similarity, independent-support, and "
                "condition-compatibility gates were explicitly bypassed",
                "Returned recipes may be structurally weak or chemically "
                "incompatible with the query",
            )
            if unrestricted
            else (
                "Condition compatibility was checked conservatively against "
                "every observed query functional-group tag",
            )
        ),
        *(
            (
                "The query reactant structures do not contain "
                f"{', '.join(required_source_elements)}; retrieved precedents "
                "must contain a curated condition source, but direct atom "
                "origin remains unverified",
            )
            if required_source_elements
            else ()
        ),
    )
    recommendations = tuple(
        replace(
            recommendation,
            cautions=tuple(
                dict.fromkeys((*fallback_cautions, *recommendation.cautions))
            ),
        )
        for recommendation in recommendations
    )
    warnings = [
        *base_warnings,
        "UNVERIFIED_REACTION_FALLBACK_USED",
        "FALLBACK_RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW",
    ]
    if unrestricted:
        warnings.extend(
            (
                "UNRESTRICTED_FALLBACK_USED",
                "FALLBACK_SIMILARITY_SUPPORT_AND_COMPATIBILITY_GATES_BYPASSED",
            )
        )
    if retrieval.level.endswith("limited_support"):
        warnings.append("LIMITED_PRECEDENT_SUPPORT")
    if required_source_elements:
        warnings.append(
            "EXPLORATORY_PARTIAL_CORRESPONDENCE_FALLBACK_USED:"
            f"{','.join(required_source_elements)}"
        )
    if retrieval.level.startswith("source_supported_partial_transformation"):
        warnings.append("SOURCE_SUPPORTED_PARTIAL_TRANSFORMATION_USED")
    if retrieval.excluded_candidate_count:
        warnings.append(
            f"INCOMPATIBLE_PRECEDENTS_EXCLUDED:{retrieval.excluded_candidate_count}"
        )
    return GenericRecommendationResult(
        query_reaction_smiles=reaction_smiles,
        valid=True,
        query_fallback_descriptor_id=descriptor_id,
        recommendation_mode=(
            "unrestricted_unverified_structure_fallback"
            if unrestricted
            else "unverified_structure_fallback"
        ),
        reaction_label=_reaction_label_payload(analysis),
        transformation_class=analysis.transformation_class,
        retrieval_definition_version=retrieval_definition_version,
        retrieval_strategy=(
            "unrestricted_unverified_structure_fallback"
            if unrestricted
            else "unverified_structure_fallback"
        ),
        retrieval_level=retrieval.level,
        candidate_count=retrieval.candidate_count,
        independent_candidate_count=retrieval.independent_candidate_count,
        compatible_candidate_count=len(retrieval.pool),
        independent_compatible_candidate_count=(
            retrieval.independent_compatible_candidate_count
        ),
        excluded_candidate_count=retrieval.excluded_candidate_count,
        retrieval_trace=retrieval.trace,
        recommendations=recommendations,
        warnings=tuple(warnings),
    )


def recommend_generic_conditions(
    reaction_smiles: str,
    *,
    records_path: str | Path = "results/generic_conversion/records.jsonl",
    top_k: int = 5,
    minimum_pool_size: int | None = None,
    unrestricted_fallback: bool = False,
    use_rxnmapper: bool = False,
    mapping_provider: AtomMappingProvider | None = None,
    ranking_preferences: ChemistRankingPreferences | None = None,
    completion_selections: Tuple[ReactionCompletionSelection, ...] = (),
) -> GenericRecommendationResult:
    """Featurize a reaction and recommend canonical resolved recipes."""
    if use_rxnmapper and mapping_provider is None:
        mapping_provider = RxnMapperProvider()
    recommender = GenericConditionRecommender.from_path(
        records_path,
        mapping_provider=mapping_provider,
        include_review=unrestricted_fallback,
    )
    return recommender.recommend(
        reaction_smiles,
        top_k=top_k,
        minimum_pool_size=minimum_pool_size,
        unrestricted_fallback=unrestricted_fallback,
        ranking_preferences=ranking_preferences,
        completion_selections=completion_selections,
    )


__all__ = [
    "GenericConditionRecommender",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
]

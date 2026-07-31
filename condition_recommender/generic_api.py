"""Public type-agnostic condition recommendation API."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Dict, Mapping

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
    retrieve_core_shape_pool_with_trace,
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
from .models import GenericRecommendationResult
from .recipe_ranking import rank_condition_recipes


def recommend_indexed_signature(
    signature: Dict[str, Any],
    index: GenericReactionIndex,
    *,
    reaction_core: Mapping[str, Any] | None = None,
    query_reaction_smiles: str = "",
    reaction_label: str | None = None,
    reaction_label_status: str = "unavailable",
    top_k: int = 5,
    minimum_pool_size: int | None = None,
    retrieval_strategy: RetrievalStrategy = "hybrid",
    ranking_weights: Mapping[str, float] | None = None,
) -> GenericRecommendationResult:
    """Recommend from an existing signature and index without re-featurization."""
    query_context = {
        "reaction_label": reaction_label,
        "reaction_label_status": reaction_label_status,
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
    if level.startswith("reaction_core_shape"):
        warnings.append("REACTION_CORE_SHAPE_RETRIEVAL_USED")
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

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        *,
        mapping_provider: AtomMappingProvider | None = None,
    ) -> "GenericConditionRecommender":
        source = Path(path)
        return cls(
            index=load_generic_index(source),
            source_path=str(source),
            mapping_provider=mapping_provider,
        )

    def recommend(
        self,
        reaction_smiles: str,
        *,
        top_k: int = 5,
        minimum_pool_size: int | None = None,
        unrestricted_fallback: bool = False,
    ) -> GenericRecommendationResult:
        """Featurize a query and recommend without reloading the index."""
        return _recommend_with_index(
            reaction_smiles,
            self.index,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
            unrestricted_fallback=unrestricted_fallback,
            mapping_provider=self.mapping_provider,
        )


def _recommend_with_index(
    reaction_smiles: str,
    index: GenericReactionIndex,
    *,
    top_k: int,
    minimum_pool_size: int | None,
    unrestricted_fallback: bool,
    mapping_provider: AtomMappingProvider | None = None,
) -> GenericRecommendationResult:
    if top_k < 1:
        return GenericRecommendationResult(
            reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE"
        )
    base_analysis = featurize_reaction(reaction_smiles)
    assessment = (
        analyze_reaction_with_external_mapping(
            reaction_smiles,
            mapping_provider,
            base_analysis=base_analysis,
        )
        if mapping_provider is not None and base_analysis.valid
        else None
    )
    analysis = assessment.analysis if assessment is not None else base_analysis
    if not analysis.valid:
        return GenericRecommendationResult(
            reaction_smiles,
            False,
            error=analysis.error or "INVALID_REACTION",
        )
    if analysis.reaction_signature is None:
        core_attempt = None
        if analysis.edit_hypotheses:
            hypothesis_result = _recommend_hypotheses_with_index(
                analysis,
                index,
                top_k=top_k,
                minimum_pool_size=minimum_pool_size,
            )
            if (
                hypothesis_result.valid
                or hypothesis_result.error
                != "QUERY_EDIT_HYPOTHESES_NOT_RETRIEVABLE"
            ):
                return _attach_external_mapping_assessment(
                    hypothesis_result,
                    assessment,
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
            )
            if core_attempt.valid:
                return _attach_external_mapping_assessment(
                    core_attempt,
                    assessment,
                )
        result = _recommend_fallback_with_index(
            analysis,
            index,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
            unrestricted=unrestricted_fallback,
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
            reaction_label=analysis.reaction_label,
            reaction_label_status=analysis.reaction_label_status,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
        )
        result = replace(
            result,
            spectator_groups=tuple(
                asdict(group) for group in analysis.spectator_groups
            ),
        )
    return _attach_external_mapping_assessment(result, assessment)


def _recommend_core_with_index(
    analysis: Any,
    index: GenericReactionIndex,
    *,
    top_k: int,
    minimum_pool_size: int | None,
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
    retrieval = retrieve_core_shape_pool_with_trace(
        core,
        compatibility_signature,
        index,
        minimum_pool_size=minimum_pool_size,
    )
    rules = load_reaction_core_retrieval_rules()
    if not retrieval.pool:
        return GenericRecommendationResult(
            query_reaction_smiles=str(analysis.input_reaction_smiles),
            valid=False,
            query_reaction_core_id=str(core.get("core_id") or "") or None,
            recommendation_mode="reaction_core_review",
            reaction_label=(
                analysis.reaction_label
                or str(core.get("generic_label") or "")
                or None
            ),
            reaction_label_status=(
                analysis.reaction_label_status
                if analysis.reaction_label
                else "reaction_core_projection"
            ),
            transformation_class=analysis.transformation_class,
            retrieval_definition_version=str(rules["schema_version"]),
            retrieval_strategy="reaction_core_shape",
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
    )
    cautions = (
        "Query retrieval used a minimized reaction-core shape, not a verified "
        "query reaction signature",
        "Only independently admitted verified precedents supplied conditions",
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
        "REACTION_CORE_SHAPE_RETRIEVAL_USED",
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
        reaction_label=(
            analysis.reaction_label
            or str(core.get("generic_label") or "")
            or None
        ),
        reaction_label_status=(
            analysis.reaction_label_status
            if analysis.reaction_label
            else "reaction_core_projection"
        ),
        transformation_class=analysis.transformation_class,
        spectator_groups=tuple(
            asdict(group) for group in analysis.spectator_groups
        ),
        retrieval_definition_version=str(rules["schema_version"]),
        retrieval_strategy="reaction_core_shape",
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
            reaction_label=analysis.reaction_label,
            reaction_label_status=analysis.reaction_label_status,
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
    recommendations = rank_condition_recipes(
        query.to_mapping(),
        retrieval.pool,
        retrieval_level=retrieval.level,
        top_k=top_k,
        similarity_assessor=assess_hypothesis_consensus_similarity,
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
        reaction_label=analysis.reaction_label,
        reaction_label_status=analysis.reaction_label_status,
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
            reaction_label=analysis.reaction_label,
            reaction_label_status=analysis.reaction_label_status,
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
            reaction_label=analysis.reaction_label,
            reaction_label_status=analysis.reaction_label_status,
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
    )
    if not retrieval.pool:
        return GenericRecommendationResult(
            query_reaction_smiles=reaction_smiles,
            valid=False,
            query_fallback_descriptor_id=descriptor_id,
            recommendation_mode="abstained",
            reaction_label=analysis.reaction_label,
            reaction_label_status=analysis.reaction_label_status,
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
    recommendations = rank_condition_recipes(
        descriptor,
        retrieval.pool,
        retrieval_level=retrieval.level,
        top_k=top_k,
        similarity_assessor=(
            lambda query, row: assess_fallback_similarity(
                query,
                row.fallback_descriptor,
            )
        ),
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
        reaction_label=analysis.reaction_label,
        reaction_label_status=analysis.reaction_label_status,
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
) -> GenericRecommendationResult:
    """Featurize a reaction and recommend canonical resolved recipes."""
    if use_rxnmapper and mapping_provider is None:
        mapping_provider = RxnMapperProvider()
    recommender = GenericConditionRecommender.from_path(
        records_path,
        mapping_provider=mapping_provider,
    )
    return recommender.recommend(
        reaction_smiles,
        top_k=top_k,
        minimum_pool_size=minimum_pool_size,
        unrestricted_fallback=unrestricted_fallback,
    )


__all__ = [
    "GenericConditionRecommender",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
]

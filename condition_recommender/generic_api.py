"""Public type-agnostic condition recommendation API."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Mapping

from reactive_taxonomy import featurize_reaction

from .generic_indexing import (
    GenericReactionIndex,
    load_generic_index,
)
from .generic_retrieval import (
    RetrievalStrategy,
    load_generic_retrieval_rules,
    retrieve_compatible_generic_pool_with_trace,
)
from .models import GenericRecommendationResult
from .recipe_ranking import rank_condition_recipes


def recommend_indexed_signature(
    signature: Dict[str, Any],
    index: GenericReactionIndex,
    *,
    query_reaction_smiles: str = "",
    top_k: int = 5,
    minimum_pool_size: int | None = None,
    retrieval_strategy: RetrievalStrategy = "hybrid",
    ranking_weights: Mapping[str, float] | None = None,
) -> GenericRecommendationResult:
    """Recommend from an existing signature and index without re-featurization."""
    retrieval_definition_version = str(
        load_generic_retrieval_rules()["schema_version"]
    )
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
            named_family=signature.get("named_family"),
            transformation_class=signature.get("transformation_class"),
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
    if retrieval.excluded_candidate_count:
        warnings.append(
            "INCOMPATIBLE_PRECEDENTS_EXCLUDED:"
            f"{retrieval.excluded_candidate_count}"
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
        named_family=signature.get("named_family"),
        transformation_class=signature.get("transformation_class"),
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

    @classmethod
    def from_path(cls, path: str | Path) -> "GenericConditionRecommender":
        source = Path(path)
        return cls(index=load_generic_index(source), source_path=str(source))

    def recommend(
        self,
        reaction_smiles: str,
        *,
        top_k: int = 5,
        minimum_pool_size: int | None = None,
    ) -> GenericRecommendationResult:
        """Featurize a query and recommend without reloading the index."""
        return _recommend_with_index(
            reaction_smiles,
            self.index,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
        )


def _recommend_with_index(
    reaction_smiles: str,
    index: GenericReactionIndex,
    *,
    top_k: int,
    minimum_pool_size: int | None,
) -> GenericRecommendationResult:
    if top_k < 1:
        return GenericRecommendationResult(
            reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE"
        )
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid:
        return GenericRecommendationResult(
            reaction_smiles,
            False,
            error=analysis.error or "INVALID_REACTION",
        )
    if analysis.reaction_signature is None:
        return GenericRecommendationResult(
            reaction_smiles,
            False,
            error="QUERY_HAS_NO_USABLE_REACTION_SIGNATURE",
        )
    return recommend_indexed_signature(
        asdict(analysis.reaction_signature),
        index,
        query_reaction_smiles=reaction_smiles,
        top_k=top_k,
        minimum_pool_size=minimum_pool_size,
    )


def recommend_generic_conditions(
    reaction_smiles: str,
    *,
    records_path: str | Path = "results/generic_conversion/records.jsonl",
    top_k: int = 5,
    minimum_pool_size: int | None = None,
) -> GenericRecommendationResult:
    """Featurize a reaction and recommend canonical resolved recipes."""
    recommender = GenericConditionRecommender.from_path(records_path)
    return recommender.recommend(
        reaction_smiles,
        top_k=top_k,
        minimum_pool_size=minimum_pool_size,
    )


__all__ = [
    "GenericConditionRecommender",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
]

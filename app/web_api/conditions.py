"""HTTP composition for the focused, combined condition recommendation UI."""

from __future__ import annotations

import logging
import os
from typing import Any

from fastapi import APIRouter, Request
from pydantic import Field, field_validator

from condition_recommender.combined import (
    RecommendationSourceResult,
    build_automation_handoff,
    combine_recommendation_results,
    load_combined_recommendation_policy,
)

from .contracts import (
    CompletionChoiceRequest,
    RecommendationRequest,
    StrictRequest,
    envelope,
)


logger = logging.getLogger(__name__)
router = APIRouter()


class CombinedConditionsRequest(StrictRequest):
    """Chemist inputs; library selection remains server configuration."""

    reaction_smiles: str = Field(min_length=1, max_length=20_000)
    completion_choices: tuple[CompletionChoiceRequest, ...] = ()

    @field_validator("reaction_smiles")
    @classmethod
    def nonblank_reaction(cls, value: str) -> str:
        """Reject blank requests before invoking either chemistry engine."""
        if not value.strip():
            raise ValueError("A reaction is required")
        return value.strip()


@router.post("/api/v1/conditions/recommend")
def recommend_conditions(
    payload: CombinedConditionsRequest, request: Request
) -> dict[str, Any]:
    """Run canonical engines independently and preserve partial failures."""
    runtime = request.app.state.runtime
    capabilities = runtime.capabilities()
    policy = load_combined_recommendation_policy()
    library_mode = os.environ.get("CONDITION_RECOMMENDER_LIBRARY_MODE", "full")
    sources = []
    for source, mode, capability in (
        ("generic", "generic", "recommendation"),
        ("weak_label", "weak_label_fallback", "weak_label_recommendation"),
    ):
        if source == "weak_label" and payload.completion_choices:
            sources.append(
                RecommendationSourceResult(
                    source,
                    "skipped",
                    {},
                    "Screening suggestions require a complete reaction drawing; confirmed fragment sources are retained in the structural results.",
                )
            )
            continue
        if not capabilities.get(capability):
            sources.append(
                RecommendationSourceResult(
                    source,
                    "unavailable",
                    {},
                    "This recommendation library is unavailable.",
                )
            )
            continue
        try:
            query = RecommendationRequest(
                reaction_smiles=payload.reaction_smiles.strip(),
                recommendation_mode=mode,
                library_mode=library_mode,
                top_k=policy["candidate_limit_per_source"],
                use_rxnmapper=(
                    bool(capabilities.get("rxnmapper_available"))
                    and os.environ.get(
                        "CONDITION_RECOMMENDER_USE_RXNMAPPER", "true"
                    ).lower()
                    not in {"false", "0", "no"}
                ),
                completion_choices=payload.completion_choices,
            )
            result = dict(runtime.recommend(query))
            result["source_dataset_name"] = (
                capabilities.get("weak_label_dataset_name")
                if source == "weak_label"
                else f"{library_mode}/{capabilities.get('index_name', 'generic_index.sqlite')}"
            )
            status = (
                "ok"
                if result.get("valid") and result.get("recommendations")
                else "abstained"
            )
            sources.append(
                RecommendationSourceResult(
                    source,
                    status,
                    result,
                    ""
                    if status == "ok"
                    else "No supported conditions were found for this query in this source.",
                )
            )
        except Exception as exc:
            # One dataset or mapper failure must not hide the other engine's
            # valid evidence. Do not expose filesystem paths to a web client.
            logger.exception("Condition recommendation source failed: %s", source)
            rebuild_required = (
                isinstance(exc, ValueError) and "rebuild" in str(exc).lower()
            )
            sources.append(
                RecommendationSourceResult(
                    source,
                    "unavailable",
                    {},
                    (
                        "The reaction library must be rebuilt for the current chemistry schema. Contact the library administrator."
                        if rebuild_required
                        else "This source could not complete the search. Its results have not been used."
                    ),
                    "INDEX_REBUILD_REQUIRED"
                    if rebuild_required
                    else "SOURCE_UNAVAILABLE",
                )
            )
    combined = combine_recommendation_results(
        payload.reaction_smiles.strip(), tuple(sources)
    )
    result = combined.to_dict()
    result["shortlist_size"] = policy["default_shortlist_size"]
    result["automation_exports"] = {
        option.option_id: build_automation_handoff(combined, (option.option_id,))
        for option in combined.recommendations
    }
    return envelope(result)

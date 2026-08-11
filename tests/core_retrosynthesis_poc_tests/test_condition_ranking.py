"""Condition-informed retrosynthesis ranking regressions."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from core_retrosynthesis_poc.condition_ranking import (
    rank_retrosynthesis_candidates_with_conditions,
)
from core_retrosynthesis_poc.cli import _parser
from core_retrosynthesis_poc.generic_models import GenericDisconnectionCandidate


def _candidate(
    name: str,
    *,
    validation: str = "verified_signature",
    score: float = 0.5,
    level: str = "L2",
) -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles="CC",
        precursor_smiles=name,
        proposed_reaction_smiles=f"{name}>>CC",
        transformation_kind=None,
        abstraction_level=level,
        compiler_engine="test",
        template_id=f"template-{name}",
        score=score,
        context_similarity=0.5,
        product_similarity=0.5,
        precursor_similarity=0.5,
        template_specificity=0.5,
        independent_reference_support=1,
        forward_validation_status=validation,
        center_transition_key="center",
        disconnection_site_key="site",
        precedent_reaction_ids=("precedent",),
        operator_id="operator",
    )


def _result(
    reaction: str,
    *,
    valid: bool,
    fallback: bool = False,
    independent_support: int = 0,
    reference_support: int = 0,
    score: float = 0.0,
) -> SimpleNamespace:
    recommendations = (
        (
            {
                "recipe_id": f"recipe-{reaction}",
                "score": score,
                "compatibility_score": 1.0,
                "reference_support": reference_support,
                "precedent_reference_ids": ("reference",),
                "resolved_recipe": {"components": []},
            },
        )
        if valid
        else ()
    )
    return SimpleNamespace(
        query_reaction_smiles=reaction,
        valid=valid,
        recommendation_mode=("fallback_descriptor" if fallback else "verified_signature"),
        retrieval_level=("environment_neighbors" if fallback else "exact_signature"),
        candidate_count=independent_support,
        independent_candidate_count=independent_support,
        compatible_candidate_count=independent_support,
        independent_compatible_candidate_count=independent_support,
        excluded_candidate_count=0,
        recommendations=recommendations,
        warnings=(("TYPE_AGNOSTIC_FALLBACK_USED",) if fallback else ()),
        error=(None if valid else "NO_CHEMICALLY_COMPATIBLE_PRECEDENT"),
    )


class _FakeRecommender:
    def __init__(self, results: dict[str, SimpleNamespace]) -> None:
        self.results = results
        self.calls: list[tuple[str, dict]] = []

    def recommend(self, reaction_smiles: str, **kwargs: object) -> SimpleNamespace:
        self.calls.append((reaction_smiles, dict(kwargs)))
        return self.results[reaction_smiles]


def test_condition_evidence_reranks_only_verified_candidates() -> None:
    invalid = _candidate("invalid", validation="invalid")
    unsupported = _candidate("unsupported")
    direct = _candidate("direct")
    fallback = _candidate("fallback")
    recommender = _FakeRecommender(
        {
            unsupported.proposed_reaction_smiles: _result(
                unsupported.proposed_reaction_smiles,
                valid=False,
            ),
            direct.proposed_reaction_smiles: _result(
                direct.proposed_reaction_smiles,
                valid=True,
                independent_support=3,
                reference_support=2,
                score=0.7,
            ),
            fallback.proposed_reaction_smiles: _result(
                fallback.proposed_reaction_smiles,
                valid=True,
                fallback=True,
                independent_support=10,
                reference_support=8,
                score=0.9,
            ),
        }
    )

    ranked = rank_retrosynthesis_candidates_with_conditions(
        (invalid, unsupported, direct, fallback),
        recommender,
        condition_top_k=2,
    )

    assert [value.candidate.precursor_smiles for value in ranked] == [
        "direct",
        "fallback",
        "unsupported",
    ]
    assert [value.retrosynthesis_rank for value in ranked] == [3, 4, 2]
    assert [value.condition_informed_rank for value in ranked] == [1, 2, 3]
    assert [value.condition_evidence.status for value in ranked] == [
        "recommended_direct",
        "recommended_fallback",
        "insufficient_evidence",
    ]
    assert invalid.proposed_reaction_smiles not in {
        reaction for reaction, _ in recommender.calls
    }
    assert all(call[1]["top_k"] == 2 for call in recommender.calls)
    assert ranked[0].to_dict()["condition_evidence"]["recommendations"][0][
        "precedent_reference_ids"
    ] == ("reference",)


def test_condition_enrichment_can_preserve_retrosynthesis_order() -> None:
    unsupported = _candidate("unsupported")
    supported = _candidate("supported")
    recommender = _FakeRecommender(
        {
            unsupported.proposed_reaction_smiles: _result(
                unsupported.proposed_reaction_smiles,
                valid=False,
            ),
            supported.proposed_reaction_smiles: _result(
                supported.proposed_reaction_smiles,
                valid=True,
                independent_support=2,
                reference_support=2,
                score=0.8,
            ),
        }
    )

    ranked = rank_retrosynthesis_candidates_with_conditions(
        (unsupported, supported),
        recommender,
        rerank=False,
    )

    assert [value.candidate.precursor_smiles for value in ranked] == [
        "unsupported",
        "supported",
    ]
    assert [value.condition_informed_rank for value in ranked] == [1, 2]


def test_condition_support_cannot_cross_structural_score_band() -> None:
    stronger = _candidate("stronger", score=0.90)
    weaker = _candidate("weaker", score=0.80)
    recommender = _FakeRecommender(
        {
            stronger.proposed_reaction_smiles: _result(
                stronger.proposed_reaction_smiles,
                valid=False,
            ),
            weaker.proposed_reaction_smiles: _result(
                weaker.proposed_reaction_smiles,
                valid=True,
                independent_support=20,
                reference_support=10,
                score=0.95,
            ),
        }
    )

    ranked = rank_retrosynthesis_candidates_with_conditions(
        (stronger, weaker),
        recommender,
    )

    assert [value.candidate.precursor_smiles for value in ranked] == [
        "stronger",
        "weaker",
    ]
    assert [value.structural_score_band for value in ranked] == [0, 2]


def test_condition_support_reranks_within_structural_score_band() -> None:
    stronger = _candidate("stronger", score=0.90)
    supported = _candidate("supported", score=0.88)
    recommender = _FakeRecommender(
        {
            stronger.proposed_reaction_smiles: _result(
                stronger.proposed_reaction_smiles,
                valid=False,
            ),
            supported.proposed_reaction_smiles: _result(
                supported.proposed_reaction_smiles,
                valid=True,
                independent_support=3,
                reference_support=2,
                score=0.75,
            ),
        }
    )

    ranked = rank_retrosynthesis_candidates_with_conditions(
        (stronger, supported),
        recommender,
    )

    assert [value.candidate.precursor_smiles for value in ranked] == [
        "supported",
        "stronger",
    ]
    assert all(value.structural_score_band == 0 for value in ranked)
    assert all(
        value.rerank_scope
        == "same_abstraction_level_and_structural_score_band"
        for value in ranked
    )


def test_condition_support_cannot_cross_abstraction_level() -> None:
    specific = _candidate("specific", score=0.70, level="L2")
    broad = _candidate("broad", score=0.99, level="L1")
    recommender = _FakeRecommender(
        {
            specific.proposed_reaction_smiles: _result(
                specific.proposed_reaction_smiles,
                valid=False,
            ),
            broad.proposed_reaction_smiles: _result(
                broad.proposed_reaction_smiles,
                valid=True,
                independent_support=20,
                reference_support=10,
                score=0.95,
            ),
        }
    )

    ranked = rank_retrosynthesis_candidates_with_conditions(
        (specific, broad),
        recommender,
    )

    assert [value.candidate.precursor_smiles for value in ranked] == [
        "specific",
        "broad",
    ]


@pytest.mark.parametrize(
    ("condition_top_k", "minimum_pool_size"),
    ((0, None), (1, 0)),
)
def test_condition_ranking_rejects_invalid_limits(
    condition_top_k: int,
    minimum_pool_size: int | None,
) -> None:
    with pytest.raises(ValueError):
        rank_retrosynthesis_candidates_with_conditions(
            (),
            _FakeRecommender({}),
            condition_top_k=condition_top_k,
            minimum_pool_size=minimum_pool_size,
        )


def test_operator_cli_accepts_condition_ranking_options() -> None:
    arguments = _parser().parse_args(
        (
            "disconnect-operators",
            "operators.json.gz",
            "CC",
            "--condition-index",
            "generic_index.sqlite",
            "--condition-top-k",
            "4",
            "--condition-minimum-pool-size",
            "6",
            "--keep-retrosynthesis-order",
            "--no-diversity",
        )
    )

    assert arguments.condition_index == "generic_index.sqlite"
    assert arguments.condition_top_k == 4
    assert arguments.condition_minimum_pool_size == 6
    assert arguments.keep_retrosynthesis_order is True
    assert arguments.no_diversity is True

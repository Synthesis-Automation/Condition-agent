"""Hierarchical SITE1 to SYN1/REAL1 ranking regressions."""

from __future__ import annotations

from dataclasses import replace

import pytest

import core_retrosynthesis_poc.generic_search as search_module
from core_retrosynthesis_poc.generic_models import (
    GenericDisconnectionCandidate,
    GenericHandleCompletionGroup,
    GenericTemplateLibrary,
)
from core_retrosynthesis_poc.hierarchical_ranking import (
    build_completion_prior_index,
    load_hierarchical_ranking_policy,
    rank_hierarchical_candidates,
)
from core_retrosynthesis_poc.generic_search import rank_operator_site_diverse


def _group(
    group_id: str,
    *,
    operator: str,
    synthon: str,
    realization: str,
    template: str,
    support: int,
) -> GenericHandleCompletionGroup:
    return GenericHandleCompletionGroup(
        completion_group_id=group_id,
        operator_id=operator,
        completion_signature=f"handle-{group_id}",
        synthon_signatures=(synthon,),
        realization_ids=(realization,),
        template_ids=(template,),
        handle_signatures=(f"handle-{group_id}",),
        observation_support=support,
        independent_reference_support=support,
    )


def _library() -> GenericTemplateLibrary:
    return GenericTemplateLibrary(
        templates=(),
        source_row_count=15,
        accepted_observation_count=15,
        rejection_counts={},
        definition={"definition_id": "test"},
        completion_groups=(
            _group(
                "g-major",
                operator="op-a",
                synthon="syn-a",
                realization="real-major",
                template="template-major",
                support=8,
            ),
            _group(
                "g-minor",
                operator="op-a",
                synthon="syn-a",
                realization="real-minor",
                template="template-minor",
                support=2,
            ),
            _group(
                "g-other-synthon",
                operator="op-a",
                synthon="syn-b",
                realization="real-other",
                template="template-other",
                support=4,
            ),
            _group(
                "g-rare-operator",
                operator="op-b",
                synthon="syn-c",
                realization="real-rare",
                template="template-rare",
                support=1,
            ),
        ),
    )


def _candidate(
    name: str,
    score: float,
    *,
    site: str,
    operator: str,
    synthon: str,
    realization: str,
    template: str,
    support: int = 1,
    level: str = "L2",
) -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles="CC",
        precursor_smiles=name,
        proposed_reaction_smiles=f"{name}>>CC",
        transformation_kind=None,
        abstraction_level=level,
        compiler_engine="test",
        template_id=template,
        score=score,
        context_similarity=score,
        product_similarity=score,
        precursor_similarity=score,
        template_specificity=score,
        independent_reference_support=support,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key=site,
        precedent_reaction_ids=(f"precedent-{name}",),
        operator_id=operator,
        realization_id=realization,
        operator_signature=f"operator-signature-{operator}",
        synthon_signature=synthon,
    )


def test_hierarchical_policy_is_versioned_and_conservative() -> None:
    policy = load_hierarchical_ranking_policy()

    assert policy.definition_id == "hierarchical_retrosynthesis_ranking.v1"
    assert policy.schema_version == "1.0"
    assert policy.backoff_order == ("operator_synthon", "operator", "global")
    assert policy.minimum_context_independent_support == 2
    assert policy.laplace_alpha == 1.0
    assert policy.preserve_abstraction_level_order is True
    assert policy.preserve_structural_score_bands is True
    assert sum(policy.weights("realization").values()) == pytest.approx(1.0)


def test_completion_prior_uses_explicit_smoothed_backoff() -> None:
    index = build_completion_prior_index(_library())
    exact = index.assess(
        _candidate(
            "major",
            0.9,
            site="site-a",
            operator="op-a",
            synthon="syn-a",
            realization="real-major",
            template="template-major",
        )
    )
    operator = index.assess(
        _candidate(
            "other",
            0.9,
            site="site-a",
            operator="op-a",
            synthon="unseen-generated-synthon",
            realization="real-other",
            template="template-other",
        )
    )
    global_prior = index.assess(
        _candidate(
            "rare",
            0.9,
            site="site-b",
            operator="op-b",
            synthon="syn-c",
            realization="real-rare",
            template="template-rare",
        )
    )

    assert exact.backoff_level == "operator_synthon"
    assert exact.independent_support == 8
    assert exact.total_context_support == 10
    assert exact.alternative_count == 2
    assert exact.probability == 0.75
    assert operator.backoff_level == "operator"
    assert operator.total_context_support == 14
    assert operator.alternative_count == 3
    assert operator.probability == pytest.approx(5 / 17)
    assert global_prior.backoff_level == "global"
    assert global_prior.total_context_support == 15
    assert global_prior.alternative_count == 4
    assert global_prior.probability == pytest.approx(2 / 19)

    sparse_library = replace(
        _library(),
        completion_groups=(_library().completion_groups[-1],),
    )
    unsupported = build_completion_prior_index(sparse_library).assess(
        _candidate(
            "rare",
            0.9,
            site="site-b",
            operator="op-b",
            synthon="syn-c",
            realization="real-rare",
            template="template-rare",
        )
    )
    assert unsupported.completion_group_id == "g-rare-operator"
    assert unsupported.backoff_level == "unavailable"
    assert unsupported.probability is None


def test_hierarchy_ranks_sites_then_synthons_and_realizations() -> None:
    major = _candidate(
        "major",
        0.90,
        site="site-a",
        operator="op-a",
        synthon="syn-a",
        realization="real-major",
        template="template-major",
    )
    minor = _candidate(
        "minor",
        0.91,
        site="site-a",
        operator="op-a",
        synthon="syn-a",
        realization="real-minor",
        template="template-minor",
    )
    other_site = _candidate(
        "other-site",
        0.89,
        site="site-b",
        operator="op-b",
        synthon="syn-c",
        realization="real-rare",
        template="template-rare",
    )
    diverse = rank_operator_site_diverse((minor, other_site, major))

    ranked = rank_hierarchical_candidates(diverse, _library())
    reversed_input = rank_hierarchical_candidates(tuple(reversed(diverse)), _library())

    assert [candidate.precursor_smiles for candidate in ranked] == [
        "major",
        "other-site",
        "minor",
    ]
    assert [candidate.precursor_smiles for candidate in reversed_input] == [
        candidate.precursor_smiles for candidate in ranked
    ]
    assert ranked[0].completion_group_id == "g-major"
    assert ranked[0].completion_prior_backoff_level == "operator_synthon"
    assert ranked[0].completion_prior_probability == 0.75
    assert ranked[0].hierarchical_site_rank == 1
    assert ranked[0].hierarchical_synthon_rank == 1
    assert ranked[0].hierarchical_realization_rank == 1
    assert ranked[1].hierarchical_site_rank == 2
    assert ranked[2].hierarchical_realization_rank == 2
    assert [candidate.hierarchical_rank for candidate in ranked] == [1, 2, 3]
    assert all(
        candidate.hierarchical_ranking_definition_id
        == "hierarchical_retrosynthesis_ranking.v1"
        for candidate in ranked
    )


def test_hierarchy_cannot_cross_score_band_or_abstraction_level() -> None:
    strong = _candidate(
        "strong",
        0.91,
        site="site-a",
        operator="op-a",
        synthon="syn-a",
        realization="real-minor",
        template="template-minor",
    )
    weaker = _candidate(
        "weaker-with-prior",
        0.80,
        site="site-b",
        operator="op-a",
        synthon="syn-a",
        realization="real-major",
        template="template-major",
    )
    broad = _candidate(
        "broad",
        0.99,
        site="site-c",
        operator="op-a",
        synthon="syn-a",
        realization="real-major",
        template="template-major",
        level="L1",
    )

    ranked = rank_hierarchical_candidates((broad, weaker, strong), _library())

    assert [candidate.precursor_smiles for candidate in ranked] == [
        "strong",
        "weaker-with-prior",
        "broad",
    ]
    assert [candidate.hierarchical_partition_key for candidate in ranked] == [
        ("L2", 0),
        ("L2", 2),
        ("L1", 0),
    ]


def test_missing_prior_is_retained_as_uncertainty() -> None:
    candidate = _candidate(
        "unknown",
        0.9,
        site="site",
        operator="unknown-op",
        synthon="unknown-synthon",
        realization="unknown-realization",
        template="unknown-template",
    )

    ranked = rank_hierarchical_candidates((candidate,), _library())

    assert ranked[0].completion_group_id == ""
    assert ranked[0].completion_prior_probability is None
    assert ranked[0].completion_prior_backoff_level == "unavailable"
    assert ranked[0].hierarchical_realization_score > 0.0


def test_hierarchy_rejects_unvalidated_candidates() -> None:
    candidate = _candidate(
        "invalid",
        0.9,
        site="site",
        operator="op-a",
        synthon="syn-a",
        realization="real-major",
        template="template-major",
    )

    with pytest.raises(ValueError, match="forward-validated"):
        rank_hierarchical_candidates(
            (replace(candidate, forward_validation_status="unresolved"),),
            _library(),
        )


def test_operator_ladder_applies_hierarchy_after_validation(monkeypatch) -> None:
    major = _candidate(
        "major",
        0.90,
        site="site-a",
        operator="op-a",
        synthon="syn-a",
        realization="real-major",
        template="template-major",
    )
    minor = _candidate(
        "minor",
        0.91,
        site="site-a",
        operator="op-a",
        synthon="syn-a",
        realization="real-minor",
        template="template-minor",
    )
    monkeypatch.setattr(
        search_module,
        "disconnect_generic_target",
        lambda *args, **kwargs: (minor, major),
    )

    ranked = search_module.disconnect_operator_ladder(
        "CC",
        _library(),
        top_k=2,
    )

    assert [candidate.precursor_smiles for candidate in ranked] == [
        "major",
        "minor",
    ]
    assert ranked[0].hierarchical_rank == 1
    assert ranked[0].completion_prior_probability == 0.75

    legacy = search_module.disconnect_operator_ladder(
        "CC",
        _library(),
        top_k=2,
        use_hierarchical_ranking=False,
    )
    assert [candidate.precursor_smiles for candidate in legacy] == [
        "minor",
        "major",
    ]
    assert all(candidate.hierarchical_rank == 0 for candidate in legacy)

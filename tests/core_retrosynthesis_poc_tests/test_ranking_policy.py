"""General retrosynthesis diversity-policy regressions."""

from __future__ import annotations

import core_retrosynthesis_poc.generic_search as search_module
from core_retrosynthesis_poc.generic_models import GenericDisconnectionCandidate
from core_retrosynthesis_poc.generic_search import (
    disconnect_operator_ladder,
    rank_operator_site_diverse,
)
from core_retrosynthesis_poc.ranking_policy import (
    load_retrosynthesis_ranking_policy,
)
from core_retrosynthesis_poc.multistep_ranking import (
    load_multistep_ranking_policy,
)


def _candidate(
    name: str,
    score: float,
    *,
    operator: str,
    site: str,
    synthon: str,
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
        context_similarity=score,
        product_similarity=score,
        precursor_similarity=score,
        template_specificity=score,
        independent_reference_support=1,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key=site,
        precedent_reaction_ids=("precedent",),
        operator_id=operator,
        synthon_signature=synthon,
    )


def test_ranking_policy_is_versioned_and_general() -> None:
    policy = load_retrosynthesis_ranking_policy()

    assert policy.definition_id == "retrosynthesis_ranking.v1"
    assert policy.schema_version == "1.0"
    assert policy.diversity_group_fields == (
        "operator_id",
        "disconnection_site_key",
        "synthon_signature",
    )
    assert policy.candidate_pool_multiplier == 4
    assert policy.diversity_score_band_width == 0.05
    assert policy.condition_score_band_width == 0.05

    multistep = load_multistep_ranking_policy()
    assert multistep.definition_id == "multistep_ranking.v1"
    assert multistep.maximum_paths_per_state == 2
    assert multistep.minimum_candidates_per_level == 0
    assert multistep.abstraction_penalty("L0") == 0.15


def test_diversity_interleaves_groups_only_within_structural_band() -> None:
    first = _candidate("a1", 0.90, operator="A", site="site-a", synthon="syn-a")
    same_group = _candidate("a2", 0.89, operator="A", site="site-a", synthon="syn-a")
    other_group = _candidate("b", 0.87, operator="B", site="site-b", synthon="syn-b")
    weak = _candidate("c", 0.70, operator="C", site="site-c", synthon="syn-c")

    ranked = rank_operator_site_diverse((same_group, weak, other_group, first))
    reversed_input = rank_operator_site_diverse((first, other_group, weak, same_group))

    assert [value.precursor_smiles for value in ranked] == [
        "a1",
        "b",
        "a2",
        "c",
    ]
    assert [value.pre_diversity_rank for value in ranked] == [1, 3, 2, 4]
    assert [value.diversity_rank for value in ranked] == [1, 2, 3, 4]
    assert [value.structural_score_band for value in ranked] == [0, 0, 0, 4]
    assert ranked[1].diversity_group_key == ("B", "site-b", "syn-b")
    assert all(
        value.ranking_policy_definition_id == "retrosynthesis_ranking.v1"
        for value in ranked
    )
    assert [value.precursor_smiles for value in reversed_input] == [
        value.precursor_smiles for value in ranked
    ]


def test_diversity_preserves_abstraction_specificity() -> None:
    specific = _candidate(
        "specific",
        0.70,
        operator="A",
        site="site-a",
        synthon="syn-a",
        level="L2",
    )
    broad = _candidate(
        "broad",
        0.99,
        operator="B",
        site="site-b",
        synthon="syn-b",
        level="L1",
    )

    ranked = rank_operator_site_diverse((broad, specific))

    assert [value.precursor_smiles for value in ranked] == [
        "specific",
        "broad",
    ]


def test_operator_ladder_expands_pool_before_diverse_selection(
    monkeypatch,
) -> None:
    first = _candidate("a1", 0.90, operator="A", site="site-a", synthon="syn-a")
    same_group = _candidate("a2", 0.89, operator="A", site="site-a", synthon="syn-a")
    other_group = _candidate("b", 0.87, operator="B", site="site-b", synthon="syn-b")
    requested_sizes = []

    def search(*args, **kwargs):
        requested_sizes.append(kwargs["top_k"])
        return (first, same_group, other_group)[: kwargs["top_k"]]

    monkeypatch.setattr(search_module, "disconnect_generic_target", search)

    ranked = disconnect_operator_ladder("CC", object(), top_k=2)

    assert requested_sizes == [8]
    assert [value.precursor_smiles for value in ranked] == ["a1", "b"]


def test_operator_ladder_reserves_multistep_fallback_levels(monkeypatch) -> None:
    by_level = {
        "L2": (
            _candidate("specific-1", 0.9, operator="A", site="a", synthon="a"),
            _candidate("specific-2", 0.8, operator="B", site="b", synthon="b"),
        ),
        "L1": (
            _candidate("broad", 0.9, operator="C", site="c", synthon="c", level="L1"),
        ),
        "L0": (
            _candidate("generic", 0.9, operator="D", site="d", synthon="d", level="L0"),
        ),
    }
    requested_levels = []

    def search(*args, **kwargs):
        level = kwargs["levels"][0]
        requested_levels.append(level)
        return by_level[level]

    monkeypatch.setattr(search_module, "disconnect_generic_target", search)

    ranked = disconnect_operator_ladder(
        "CC",
        object(),
        top_k=4,
        minimum_candidates_per_level=1,
    )

    assert requested_levels == ["L2", "L1", "L0"]
    assert [candidate.precursor_smiles for candidate in ranked] == [
        "specific-1",
        "broad",
        "generic",
        "specific-2",
    ]

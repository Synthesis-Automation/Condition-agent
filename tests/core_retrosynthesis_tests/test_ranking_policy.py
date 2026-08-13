"""General retrosynthesis diversity-policy regressions."""

from __future__ import annotations

from dataclasses import replace

from cas_tools import PrecursorEvidence, assess_precursor_realism
import core_retrosynthesis.generic_search as search_module
from core_retrosynthesis.generic_models import (
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
)
from core_retrosynthesis.generic_search import (
    disconnect_operator_ladder,
    disconnect_operator_ladder_detailed,
    rank_operator_site_diverse,
    rank_precursor_realism,
)
from core_retrosynthesis.ranking_policy import (
    load_retrosynthesis_ranking_policy,
)
from core_retrosynthesis.multistep_ranking import (
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

    assert policy.definition_id == "retrosynthesis_ranking.v3"
    assert policy.schema_version == "3.0"
    assert policy.diversity_group_fields == (
        "operator_id",
        "disconnection_site_key",
        "synthon_signature",
    )
    assert policy.candidate_pool_multiplier == 4
    assert policy.diversity_score_band_width == 0.05
    assert policy.condition_score_band_width == 0.05
    assert policy.precursor_realism_band_penalty(0.95) == 0
    assert policy.precursor_realism_band_penalty(0.50) == 1
    assert policy.precursor_realism_band_penalty(0.25) == 2
    assert policy.precursor_realism_band_penalty(0.05) == 3
    assert policy.strategic_reserved_candidates == 2
    assert policy.strategic_maximum_band_displacement == 2

    multistep = load_multistep_ranking_policy()
    assert multistep.definition_id == "multistep_ranking.v3"
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
        value.ranking_policy_definition_id == "retrosynthesis_ranking.v3"
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


def test_realism_penalty_can_demote_across_structural_score_bands() -> None:
    unlikely = _candidate(
        "CO",
        0.90,
        operator="A",
        site="site-a",
        synthon="syn-a",
    )
    realistic = _candidate(
        "CCO",
        0.87,
        operator="B",
        site="site-b",
        synthon="syn-b",
    )
    weaker_band = _candidate(
        "CCBr",
        0.80,
        operator="C",
        site="site-c",
        synthon="syn-c",
    )
    unlikely = replace(
        unlikely,
        precursor_realism_score=assess_precursor_realism(
            "CO", PrecursorEvidence(False, False, False)
        ).score,
    )
    realistic = replace(
        realistic,
        precursor_realism_score=assess_precursor_realism(
            "CCO", PrecursorEvidence(True, True, True)
        ).score,
    )
    weaker_band = replace(
        weaker_band,
        precursor_realism_score=assess_precursor_realism(
            "CCBr", PrecursorEvidence(True, True, True)
        ).score,
    )

    ranked = rank_precursor_realism((unlikely, weaker_band, realistic))

    assert [candidate.precursor_smiles for candidate in ranked] == [
        "CCO",
        "CCBr",
        "CO",
    ]
    assert [candidate.pre_realism_rank for candidate in ranked] == [2, 3, 1]
    assert [candidate.precursor_realism_rank for candidate in ranked] == [1, 2, 3]
    assert [candidate.structural_score_band for candidate in ranked] == [0, 2, 0]
    assert [candidate.precursor_realism_band_penalty for candidate in ranked] == [
        0,
        0,
        3,
    ]
    assert [candidate.effective_structural_score_band for candidate in ranked] == [
        0,
        2,
        3,
    ]


def test_reported_free_carbamic_acid_route_is_demoted_by_realism() -> None:
    best = _candidate(
        "best",
        0.66160498,
        operator="A",
        site="site-a",
        synthon="syn-a",
    )
    carbamic_acid = _candidate(
        "CC(C)(C)I.COC(=O)[C@@H](Cc1ccc(O)cc1)NC(=O)O",
        0.60049131,
        operator="B",
        site="site-b",
        synthon="syn-b",
    )
    boc_chloride = _candidate(
        "CC(C)(C)OC(=O)Cl.COC(=O)[C@H](N)Cc1ccc(O)cc1",
        0.51744618,
        operator="C",
        site="site-c",
        synthon="syn-c",
    )
    best = replace(best, precursor_realism_score=0.975)
    carbamic_acid = replace(carbamic_acid, precursor_realism_score=0.25)
    boc_chloride = replace(boc_chloride, precursor_realism_score=0.457719)

    ranked = rank_precursor_realism((carbamic_acid, boc_chloride, best))

    assert [candidate.precursor_smiles for candidate in ranked] == [
        "best",
        "CC(C)(C)OC(=O)Cl.COC(=O)[C@H](N)Cc1ccc(O)cc1",
        "CC(C)(C)I.COC(=O)[C@@H](Cc1ccc(O)cc1)NC(=O)O",
    ]
    assert [candidate.structural_score_band for candidate in ranked] == [0, 2, 1]
    assert [candidate.precursor_realism_band_penalty for candidate in ranked] == [
        0,
        1,
        2,
    ]
    assert [candidate.effective_structural_score_band for candidate in ranked] == [
        0,
        3,
        3,
    ]


def test_realism_option_attaches_component_score_before_diversity(monkeypatch) -> None:
    unlikely = _candidate("CO", 0.90, operator="A", site="a", synthon="a")
    realistic = _candidate("CCO", 0.87, operator="B", site="b", synthon="b")

    monkeypatch.setattr(
        search_module,
        "disconnect_generic_target",
        lambda *args, **kwargs: (unlikely, realistic),
    )

    def scorer(smiles):
        evidence = PrecursorEvidence(smiles == "CCO", False, False)
        return (assess_precursor_realism(smiles, evidence),)

    ranked = disconnect_operator_ladder(
        "CC",
        object(),
        top_k=2,
        precursor_realism_scorer=scorer,
    )

    assert [candidate.precursor_smiles for candidate in ranked] == ["CCO", "CO"]
    assert ranked[0].precursor_realism_score == 0.95
    assert ranked[0].precursor_realism_assessments[0].evidence.buyable is True
    assert ranked[0].precursor_realism_aggregation is not None
    assert (
        ranked[0].precursor_realism_aggregation.known_substantial_component_bonus
        == 0.0
    )


def test_detailed_operator_ladder_applies_the_same_realism_ranking(
    monkeypatch,
) -> None:
    unlikely = _candidate("CO", 0.90, operator="A", site="a", synthon="a")
    realistic = _candidate("CCO", 0.87, operator="B", site="b", synthon="b")
    monkeypatch.setattr(
        search_module,
        "disconnect_generic_target_detailed",
        lambda *args, **kwargs: (
            (unlikely, realistic),
            GenericSearchDiagnostics(valid_candidate_count=2),
        ),
    )

    def scorer(smiles):
        evidence = PrecursorEvidence(smiles == "CCO", False, False)
        return (assess_precursor_realism(smiles, evidence),)

    ranked, diagnostics = disconnect_operator_ladder_detailed(
        "CC",
        object(),
        top_k=2,
        precursor_realism_scorer=scorer,
    )

    assert [candidate.precursor_smiles for candidate in ranked] == ["CCO", "CO"]
    assert ranked[0].precursor_realism_score == 0.95
    assert diagnostics.valid_action_count == 2


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


def test_operator_ladder_reserves_a_bounded_strategic_candidate(
    monkeypatch,
) -> None:
    routine = tuple(
        _candidate(
            f"routine-{index}",
            0.90 - index * 0.01,
            operator=f"R{index}",
            site=f"routine-site-{index}",
            synthon=f"routine-synthon-{index}",
        )
        for index in range(5)
    )
    strategic = replace(
        _candidate(
            "strategic",
            0.79,
            operator="S",
            site="strategic-site",
            synthon="strategic-synthon",
        ),
        strategic_complexity_score=0.7,
        strategic_class="scaffold_split",
        strategic_candidate=True,
    )
    monkeypatch.setattr(
        search_module,
        "disconnect_generic_target",
        lambda *args, **kwargs: (*routine, strategic),
    )

    ranked = disconnect_operator_ladder(
        "CC",
        object(),
        top_k=5,
        use_hierarchical_ranking=False,
    )

    assert len(ranked) == 5
    retained = next(
        candidate for candidate in ranked if candidate.precursor_smiles == "strategic"
    )
    assert retained.strategic_reserve_selected is True


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

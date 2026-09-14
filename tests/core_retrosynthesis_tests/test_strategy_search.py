"""Single-step STRAT1 identity and strategy-grouping regressions."""

from __future__ import annotations

from dataclasses import replace

import pytest

import core_retrosynthesis.strategy_search as strategy_search_module
from core_retrosynthesis.generic_models import (
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
)
from core_retrosynthesis.strategy_identity import build_strategy_id
from core_retrosynthesis.strategy_search import (
    disconnect_strategies,
    disconnect_strategies_detailed,
    group_strategy_candidates,
)


@pytest.mark.parametrize("strategic_score, expected", [(0.85, True), (0.1, False)])
def test_grouped_search_preserves_guarded_scaffold_reservation(
    monkeypatch, strategic_score, expected
):
    candidates = tuple(
        _candidate(f"choice-{i}", 0.90 - i * 0.01, site=f"SITE1:{i}") for i in range(5)
    )
    strategic = replace(
        _candidate("scaffold", strategic_score, site="SITE1:scaffold"),
        strategic_candidate=True,
    )
    monkeypatch.setattr(
        strategy_search_module,
        "disconnect_generic_target_detailed",
        lambda *args, **kwargs: ((*candidates, strategic), GenericSearchDiagnostics()),
    )
    result = disconnect_strategies_detailed(
        "CC", object(), top_k_strategies=5, use_hierarchical_ranking=False
    )
    assert (
        any(s.representative.strategic_reserve_selected for s in result.strategies)
        is expected
    )
    assert len({s.strategy_id for s in result.strategies}) == 5


def _candidate(
    name: str,
    score: float,
    *,
    operator: str = "OP1:operator",
    site: str = "SITE1:site",
    synthon: str = "SYN1:synthon",
    support: int = 1,
    precedent: str | None = None,
) -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles="CC",
        precursor_smiles=name,
        proposed_reaction_smiles=f"{name}>>CC",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="test",
        template_id=f"template-{name}",
        score=score,
        context_similarity=score,
        product_similarity=score,
        precursor_similarity=score,
        template_specificity=score,
        independent_reference_support=support,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key=site,
        precedent_reaction_ids=(precedent or f"precedent-{name}",),
        operator_id=operator,
        realization_id=f"REAL1:{name}",
        operator_signature=f"operator-signature-{operator}",
        synthon_signature=synthon,
    )


def test_strategy_id_is_deterministic_and_handle_independent() -> None:
    bromide = _candidate("CCBr", 0.9)
    iodide = _candidate("CCI", 0.8)
    expected = build_strategy_id("OP1:operator", "SITE1:site", "SYN1:synthon")

    assert expected.startswith("STRAT1:")
    assert bromide.strategy_id == expected
    assert iodide.strategy_id == expected
    assert (
        _candidate(
            "CCB(O)O",
            0.7,
            synthon="SYN1:other",
        ).strategy_id
        != expected
    )


def test_candidate_rejects_a_supplied_strategy_id_that_conflicts_with_graph() -> None:
    candidate = _candidate("CCBr", 0.9)

    with pytest.raises(ValueError, match="contradicts graph identity"):
        replace(candidate, strategy_id="STRAT1:not-the-graph-identity")


def test_grouping_returns_distinct_strategies_with_bounded_realizations() -> None:
    major = _candidate("major", 0.91, support=3, precedent="p-major")
    minor = _candidate("minor", 0.89, support=8, precedent="p-minor")
    truncated = _candidate("truncated", 0.87, support=2, precedent="p-third")
    other = _candidate(
        "other",
        0.88,
        operator="OP1:other",
        support=4,
        precedent="p-other",
    )

    grouped = group_strategy_candidates(
        (major, minor, other, truncated),
        top_k_strategies=2,
        max_realizations_per_strategy=2,
    )

    assert [proposal.representative.precursor_smiles for proposal in grouped] == [
        "major",
        "other",
    ]
    assert [candidate.precursor_smiles for candidate in grouped[0].realizations] == [
        "major",
        "minor",
    ]
    assert grouped[0].strategy_rank == 1
    assert grouped[0].total_realization_count == 3
    assert grouped[0].independent_reference_support == 8
    assert grouped[0].precedent_reaction_ids == ("p-major", "p-minor", "p-third")
    serialized = grouped[0].to_dict()
    assert serialized["returned_realization_count"] == 2
    assert serialized["representative"]["strategy_id"] == grouped[0].strategy_id
    assert len(serialized["alternate_realizations"]) == 1


def test_grouping_refuses_unverified_or_incomplete_candidates() -> None:
    candidate = _candidate("candidate", 0.9)

    with pytest.raises(ValueError, match="verified-signature"):
        group_strategy_candidates(
            (replace(candidate, forward_validation_status="core_only"),)
        )

    incomplete = replace(
        candidate,
        operator_id="",
        strategy_id="",
    )
    with pytest.raises(ValueError, match="complete STRAT1"):
        group_strategy_candidates((incomplete,))

    other_target = replace(candidate, target_smiles="CCC")
    with pytest.raises(ValueError, match="mix target molecules"):
        group_strategy_candidates((candidate, other_target))


def test_disconnect_strategies_uses_a_bounded_diverse_candidate_pool(
    monkeypatch,
) -> None:
    first = _candidate("first", 0.9)
    same = _candidate("same", 0.8)
    second = _candidate("second", 0.7, site="SITE1:second")
    requested: dict[str, object] = {}

    def search(*args, **kwargs):
        requested.update(kwargs)
        return (first, same, second), GenericSearchDiagnostics()

    monkeypatch.setattr(
        strategy_search_module,
        "disconnect_generic_target_detailed",
        search,
    )

    proposals = disconnect_strategies(
        "CC",
        object(),
        top_k_strategies=3,
        max_realizations_per_strategy=2,
        max_candidates_to_validate=50,
        diversify=False,
    )

    assert requested["top_k"] == 50
    assert requested["balance_operator_budget"] is True
    assert requested["max_candidates_to_validate"] == 50
    assert [proposal.representative.precursor_smiles for proposal in proposals] == [
        "first",
        "second",
    ]
    assert [candidate.precursor_smiles for candidate in proposals[0].realizations] == [
        "first",
        "same",
    ]


def test_strategy_shortfall_broadens_even_when_realizations_fill_top_k(
    monkeypatch,
) -> None:
    calls = []

    def search(*args, **kwargs):
        level = kwargs["levels"][0]
        calls.append(level)
        values = (
            (_candidate("Br", 0.9), _candidate("I", 0.8))
            if level == "L2"
            else (_candidate("other", 0.7, site="SITE1:other"),)
        )
        return values, GenericSearchDiagnostics(validation_attempt_count=len(values))

    monkeypatch.setattr(
        strategy_search_module, "disconnect_generic_target_detailed", search
    )
    result = disconnect_strategies_detailed(
        "CC", object(), top_k_strategies=2, diversify=False
    )
    assert calls == ["L2", "L1"]
    assert len(result.strategies) == 2
    assert len(result.strategies[0].realizations) == 2
    assert result.to_dict()["search_diagnostics"]["strategy_target_met"] is True
    assert result.diagnostics.validation_attempt_count == 3


def test_empty_strategy_search_reports_limits_without_claiming_library_absence(
    monkeypatch,
) -> None:
    def search(*args, **kwargs):
        return (), GenericSearchDiagnostics(
            product_query_match_count=7,
            applied_template_count=2,
            template_budget_excluded_count=5,
            generated_precursor_count=4,
            validation_attempt_count=1,
            validation_budget_excluded_count=3,
            invalid_forward_count=1,
        )

    monkeypatch.setattr(
        strategy_search_module, "disconnect_generic_target_detailed", search
    )
    result = disconnect_strategies_detailed(
        "CC", object(), include_l0=False, diversify=False
    )
    diagnostics = result.to_dict()["search_diagnostics"]
    assert diagnostics["levels_attempted"] == ["L2", "L1"]
    assert diagnostics["budget_limited"] is True
    assert diagnostics["strategy_target_met"] is False
    assert result.strategies == ()


def test_operator_scheduler_preserves_order_and_does_not_starve_other_edits() -> None:
    from core_retrosynthesis.generic_search import _interleave_operators

    values = [("A", 1), ("A", 2), ("A", 3), ("B", 1), ("C", 1), ("B", 2)]
    ordered = _interleave_operators(values, lambda value: value[0])
    assert ordered[:3] == [("A", 1), ("B", 1), ("C", 1)]
    assert [value for value in ordered if value[0] == "A"] == values[:3]


def test_validation_budget_visits_other_operators_and_sites(
    reduction_library, monkeypatch,
) -> None:
    from types import SimpleNamespace
    from core_retrosynthesis import generic_search

    template = next(
        t for t in reduction_library.templates if t.abstraction_level == "L2"
    )
    library = SimpleNamespace(
        templates=tuple(
            replace(template, template_id=key, operator_id=key, reaction_smarts=key)
            for key in ("A", "B")
        ),
        retrieval_index=None,
    )
    proposals = {
        "A": (("C", "A1"), ("CC", "A2"), ("CCC", "A3")),
        "B": (("C", "B1"),),
    }
    sites = {"A1": "site1", "A2": "site1", "A3": "site2", "B1": "site1"}
    visited = []
    monkeypatch.setattr(generic_search, "maximum_similarity", lambda *args: 0.9)
    monkeypatch.setattr(generic_search, "_apply", lambda key, target: proposals[key])
    monkeypatch.setattr(generic_search, "provisional_product_site", sites.__getitem__)

    def reject(mapped, **kwargs):
        visited.append(mapped)
        return "invalid", None, None, ""

    monkeypatch.setattr(generic_search, "_forward_analysis", reject)
    candidates, diagnostics = generic_search.disconnect_generic_target_detailed(
        "OCc1ccccc1", library, balance_operator_budget=True,
        max_candidates_to_validate=3,
    )
    assert visited == ["A1", "B1", "A3"]
    assert candidates == ()
    assert diagnostics.provisional_site_group_count == 3
    assert diagnostics.validation_budget_excluded_count == 1
    assert diagnostics.invalid_forward_count == 3


def test_incomplete_strategy_preserves_validated_proposal_for_review(
    monkeypatch,
) -> None:
    incomplete = replace(_candidate("CCBr", 0.9), synthon_signature="", strategy_id="")
    monkeypatch.setattr(
        strategy_search_module,
        "disconnect_generic_target_detailed",
        lambda *args, **kwargs: (
            (incomplete,),
            GenericSearchDiagnostics(valid_candidate_count=1),
        ),
    )
    result = disconnect_strategies_detailed("CC", object(), diversify=False)
    assert result.strategies == ()
    assert result.unresolved_candidates == (incomplete,)
    assert (
        result.to_dict()["search_diagnostics"]["incomplete_strategy_identity_count"]
        == 1
    )


@pytest.fixture(scope="module")
def reduction_library():
    from core_retrosynthesis.generic_library import build_generic_library

    return build_generic_library(
        [
            {
                "reaction_id": "amine",
                "reference_id": "reference:amine",
                "reaction_smiles": "O=Cc1ccccc1>>OCc1ccccc1",
            }
        ],
        levels=("L2", "L1", "L0"),
        admission_mode="data_driven",
    )


def test_grouped_search_reconstructs_real_mapped_chemistry_deterministically(
    reduction_library,
) -> None:
    assert reduction_library.templates
    first = disconnect_strategies_detailed(
        "OCc1ccc(F)cc1", reduction_library, top_k_strategies=2
    )
    second = disconnect_strategies_detailed(
        "Fc1ccc(CO)cc1", reduction_library, top_k_strategies=2
    )
    assert first.strategies
    assert first.to_dict() == second.to_dict()
    assert all(
        c.forward_validation_status == "verified_signature"
        for s in first.strategies
        for c in s.realizations
    )


@pytest.mark.parametrize("failure", ["invalid", "unresolved", "conflicting"])
def test_balanced_search_cannot_bypass_validation(
    reduction_library, monkeypatch, failure
) -> None:
    from core_retrosynthesis import generic_search

    if failure == "invalid":
        monkeypatch.setattr(
            generic_search,
            "_forward_analysis",
            lambda *args, **kwargs: ("invalid", None, None, ""),
        )
    elif failure == "unresolved":
        monkeypatch.setattr(
            generic_search, "analyze_generic_reaction", lambda *args: None
        )
    else:
        from types import SimpleNamespace

        monkeypatch.setattr(
            generic_search,
            "analyze_generic_reaction",
            lambda *args: SimpleNamespace(operator_signature="conflicting-edit"),
        )
    result = disconnect_strategies_detailed("OCc1ccc(F)cc1", reduction_library)
    assert not result.strategies
    assert result.diagnostics.validation_attempt_count > 0
    counter = {
        "invalid": "invalid_forward_count",
        "unresolved": "unresolved_identity_count",
        "conflicting": "operator_mismatch_count",
    }[failure]
    assert (
        sum(
            getattr(diagnostics, counter)
            for _, diagnostics in result.diagnostics.level_diagnostics
        )
        > 0
    )

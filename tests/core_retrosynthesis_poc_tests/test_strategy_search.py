"""Single-step STRAT1 identity and strategy-grouping regressions."""

from __future__ import annotations

from dataclasses import replace

import pytest

import core_retrosynthesis_poc.strategy_search as strategy_search_module
from core_retrosynthesis_poc.generic_models import GenericDisconnectionCandidate
from core_retrosynthesis_poc.strategy_identity import build_strategy_id
from core_retrosynthesis_poc.strategy_search import (
    disconnect_strategies,
    group_strategy_candidates,
)


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
    assert _candidate(
        "CCB(O)O",
        0.7,
        synthon="SYN1:other",
    ).strategy_id != expected


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
        return first, same, second

    monkeypatch.setattr(
        strategy_search_module,
        "disconnect_operator_ladder",
        search,
    )

    proposals = disconnect_strategies(
        "CC",
        object(),
        top_k_strategies=3,
        max_realizations_per_strategy=2,
        max_candidates_to_validate=50,
    )

    assert requested["top_k"] == 12
    assert requested["max_candidates_to_validate"] == 50
    assert [proposal.representative.precursor_smiles for proposal in proposals] == [
        "first",
        "second",
    ]
    assert [candidate.precursor_smiles for candidate in proposals[0].realizations] == [
        "first",
        "same",
    ]

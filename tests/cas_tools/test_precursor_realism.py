"""Deterministic precursor-realism heuristic regressions."""

from __future__ import annotations

import pytest

from cas_tools import (
    PrecursorEvidence,
    aggregate_precursor_realism,
    assess_precursor_components,
    assess_precursor_realism,
    load_precursor_realism_policy,
    validate_precursor_realism_definition,
)


def _evidence(
    *,
    buyable: bool = False,
    registry: bool = False,
    literature: bool = False,
) -> PrecursorEvidence:
    return PrecursorEvidence(
        buyable=buyable,
        in_compound_registry=registry,
        in_literature=literature,
    )


def test_all_three_sources_produce_full_realism_trace() -> None:
    assessment = assess_precursor_realism(
        "CCO",
        _evidence(buyable=True, registry=True, literature=True),
    )

    assert assessment.canonical_smiles == "CCO"
    assert assessment.evidence_tier == "buyable_registry_literature"
    assert assessment.base_score == 1.0
    assert assessment.molecular_weight_penalty == 0.0
    assert assessment.score == 1.0
    assert assessment.definition_id == "precursor_realism.v1"
    assert assessment.schema_version == "1.0"
    assert assessment.to_dict()["evidence"]["buyable"] is True


def test_buyable_evidence_prevents_small_molecule_penalty() -> None:
    assessment = assess_precursor_realism(
        "CO",
        _evidence(buyable=True),
    )

    assert assessment.evidence_tier == "buyable_only"
    assert assessment.molecular_weight_smallness == 1.0
    assert assessment.molecular_weight_penalty == 0.0
    assert assessment.score == 0.95


def test_unsupported_small_precursor_is_penalized_more_than_large_one() -> None:
    small = assess_precursor_realism("CO", _evidence())
    large = assess_precursor_realism("CCCCCCCCCCCCCCCC", _evidence())

    assert small.molecular_weight_penalty == 0.2
    assert small.score == 0.05
    assert large.molecular_weight_penalty == 0.0
    assert large.score == 0.25


def test_each_nonbuyable_evidence_tier_is_explicit() -> None:
    both = assess_precursor_realism(
        "CCBr",
        _evidence(registry=True, literature=True),
    )
    registry = assess_precursor_realism("CCBr", _evidence(registry=True))
    literature = assess_precursor_realism("CCBr", _evidence(literature=True))

    assert both.evidence_tier == "registry_literature"
    assert registry.evidence_tier == "registry_only"
    assert literature.evidence_tier == "literature_only"
    assert both.score > registry.score > literature.score


def test_candidate_aggregation_uses_weakest_required_precursor() -> None:
    strong = assess_precursor_realism("CCO", _evidence(buyable=True))
    weak = assess_precursor_realism("CO", _evidence())

    assert aggregate_precursor_realism((strong, weak)) == weak.score


def test_component_assessment_uses_source_neutral_evidence_provider() -> None:
    seen = []

    def evidence_provider(identity):
        seen.append(identity.canonical_smiles)
        return _evidence(buyable=identity.canonical_smiles == "CCO")

    assessments = assess_precursor_components("CO.CCO", evidence_provider)

    assert seen == ["CCO", "CO"]
    assert [item.evidence_tier for item in assessments] == [
        "buyable_only",
        "unsupported",
    ]
    assert aggregate_precursor_realism(assessments) == assessments[1].score


def test_invalid_structure_and_empty_aggregation_are_rejected() -> None:
    with pytest.raises(ValueError, match="valid molecular structure"):
        assess_precursor_realism("not-smiles", _evidence())
    with pytest.raises(ValueError, match="at least one"):
        aggregate_precursor_realism(())


def test_policy_and_definition_validation_are_versioned() -> None:
    policy = load_precursor_realism_policy()

    assert policy.definition_id == "precursor_realism.v1"
    assert policy.maximum_smallness_da == 50.0
    assert policy.no_smallness_da == 200.0
    with pytest.raises(ValueError, match="definition ID"):
        validate_precursor_realism_definition(
            {
                "definition_id": "wrong",
                "schema_version": "1.0",
            }
        )

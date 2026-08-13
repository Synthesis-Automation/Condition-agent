"""Planner integration for component-scoped precursor compatibility."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis_poc import (
    assess_precursor_compatibility,
    load_precursor_compatibility_policy,
)
from core_retrosynthesis_poc.generic_models import GenericDisconnectionCandidate
from core_retrosynthesis_poc.generic_search import rank_operator_site_diverse


def _candidate(name: str, score: float) -> GenericDisconnectionCandidate:
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
        independent_reference_support=1,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key=name,
        precedent_reaction_ids=("precedent",),
        operator_id=name,
        synthon_signature=name,
    )


def test_policy_demotes_strong_same_component_warning() -> None:
    policy = load_precursor_compatibility_policy()
    result = assess_precursor_compatibility("NCCC(=O)Cl")

    assert policy.definition_id == "precursor_compatibility_policy.v1"
    assert result.disposition == "demote"
    assert result.warning_strength == "strong"
    assert result.structural_band_penalty == 1
    assert result.assessments[0].scope == "same_component"


def test_separate_reactive_components_receive_no_intrinsic_penalty() -> None:
    result = assess_precursor_compatibility("CN.CC(=O)Cl")

    assert result.assessments == ()
    assert result.disposition == "pass"
    assert result.structural_band_penalty == 0
    assert result.warning_strength is None


def test_compatibility_penalty_precedes_hierarchical_partitions() -> None:
    conflicted = _candidate("conflicted", 0.90)
    compatible = _candidate("compatible", 0.87)
    conflicted = replace(
        conflicted,
        precursor_compatibility_disposition="demote",
        precursor_compatibility_warning_strength="strong",
        precursor_compatibility_band_penalty=1,
        precursor_compatibility_policy_definition_id=(
            "precursor_compatibility_policy.v1"
        ),
    )

    ranked = rank_operator_site_diverse((conflicted, compatible))

    assert [item.precursor_smiles for item in ranked] == [
        "compatible",
        "conflicted",
    ]
    assert ranked[0].effective_structural_score_band == 0
    assert ranked[1].effective_structural_score_band == 1


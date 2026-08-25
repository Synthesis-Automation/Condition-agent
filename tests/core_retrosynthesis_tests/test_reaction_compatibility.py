"""Planner integration for graph-derived reaction-regime compatibility."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis import (
    assess_candidate_reaction_compatibility,
    load_reaction_compatibility_policy,
)
from core_retrosynthesis.generic_models import GenericDisconnectionCandidate
from core_retrosynthesis.generic_search import rank_operator_site_diverse


CONFLICTING_REACTION = (
    "Brc1ccccc1.O=C(c1ccccc1)C12CC[N+](CCO)(CC1)CC2"
    ">>OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2"
)


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


def test_policy_demotes_strong_reaction_regime_conflict() -> None:
    policy = load_reaction_compatibility_policy()
    result = assess_candidate_reaction_compatibility(CONFLICTING_REACTION)

    assert policy.definition_id == "reaction_compatibility_policy.v1"
    assert result.disposition == "demote"
    assert result.warning_strength == "strong"
    assert result.structural_band_penalty == 1
    assert result.assessments[0].inferred_regime == (
        "strongly_basic_carbon_transfer"
    )


def test_reaction_compatibility_penalty_precedes_similarity_score() -> None:
    conflicted = replace(
        _candidate("conflicted", 0.90),
        reaction_compatibility_disposition="demote",
        reaction_compatibility_warning_strength="strong",
        reaction_compatibility_band_penalty=1,
        reaction_compatibility_policy_definition_id=(
            "reaction_compatibility_policy.v1"
        ),
    )
    compatible = _candidate("compatible", 0.87)

    ranked = rank_operator_site_diverse((conflicted, compatible))

    assert [item.precursor_smiles for item in ranked] == [
        "compatible",
        "conflicted",
    ]
    assert ranked[0].effective_structural_score_band == 0
    assert ranked[1].effective_structural_score_band == 1

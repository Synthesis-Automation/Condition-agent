"""Reaction-regime compatibility diagnostics from graph evidence."""

from __future__ import annotations

from copy import deepcopy

from reactive_taxonomy import (
    assess_reaction_compatibility,
    load_reaction_compatibility_definition,
    validate_reaction_compatibility_definition,
)


FREE_ALCOHOL_CARBONYL_ADDITION = (
    "Brc1ccccc1.O=C(c1ccccc1)C12CC[N+](CCO)(CC1)CC2"
    ">>OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2"
)


def test_halo_carbonyl_addition_with_free_alcohol_is_strong_conflict() -> None:
    assessments = assess_reaction_compatibility(FREE_ALCOHOL_CARBONYL_ADDITION)

    assert len(assessments) == 1
    assessment = assessments[0]
    assert assessment.rule_id == (
        "RCI001_HALO_CARBON_TRANSFER_PROTIC_CONFLICT"
    )
    assert assessment.inferred_regime == "strongly_basic_carbon_transfer"
    assert assessment.spectator_group_id == "alcohol"
    assert assessment.warning_strength == "strong"
    assert assessment.intrinsic_severity == "high"
    assert assessment.evidence_quality == "global_atom_correspondence"
    assert assessment.confidence == 0.8
    assert len(assessment.matched_edit_indices) == 3


def test_protected_alcohol_does_not_trigger_protic_conflict() -> None:
    protected = (
        "Brc1ccccc1.O=C(c1ccccc1)C12CC[N+](CCOC)(CC1)CC2"
        ">>COCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2"
    )

    assert assess_reaction_compatibility(protected) == ()


def test_free_alcohol_without_carbon_transfer_motif_does_not_trigger() -> None:
    assert assess_reaction_compatibility("O=CCCCO>>OCCCCCO") == ()


def test_reaction_compatibility_definition_rejects_unknown_spectator_tag() -> None:
    payload = deepcopy(dict(load_reaction_compatibility_definition()))
    payload["rules"] = [dict(payload["rules"][0])]
    payload["rules"][0]["spectator_tags_any"] = ["invented_tag"]

    errors = validate_reaction_compatibility_definition(payload)

    assert any("unknown_spectator_tags" in error for error in errors)

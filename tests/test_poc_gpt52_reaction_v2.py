"""Tests for general-purpose taxonomy-first reaction analysis PoC v2."""

from poc_gpt52_reaction_v2 import analyze_reaction_general


def test_v2_detects_amide_formation_without_llm() -> None:
    reaction = "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1"
    result = analyze_reaction_general(reaction, use_llm=False)

    assert result.decision.reaction_type == "Amide_formation"
    assert result.decision.source == "deterministic"
    assert result.decision.confidence >= 0.9
    assert result.validation.passed
    assert result.taxonomy_candidates


def test_v2_detects_oxidation_without_llm() -> None:
    reaction = "CCO>>CC=O"
    result = analyze_reaction_general(reaction, use_llm=False)

    assert result.decision.reaction_type == "Oxidation_primary_alcohol_to_aldehyde"
    assert result.decision.confidence >= 0.9
    assert result.validation.passed


def test_v2_returns_unknown_when_no_taxonomy_candidate() -> None:
    reaction = "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1"
    result = analyze_reaction_general(reaction, use_llm=False)

    assert result.decision.reaction_type == "unknown"
    assert result.decision.source == "deterministic"
    assert not result.taxonomy_candidates

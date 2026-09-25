"""Condition comparison preserves structural observations and source independence."""

from copy import deepcopy
from dataclasses import asdict, replace

from condition_recommender import compare_condition_evidence
from condition_recommender.generic_indexing import GenericIndexedReaction
from condition_recommender.models import PrecedentTier
from reactive_taxonomy import featurize_reaction


def precedent(recipe: dict | None = None) -> GenericIndexedReaction:
    """Small actual structural precedent without a source-corpus dependency."""
    analysis = featurize_reaction("CCBr.N>>CCN")
    return GenericIndexedReaction(
        reaction_id="r1", observation_id="o1", canonical_reaction_id="cr1",
        reaction_smiles="CCBr.N>>CCN", yield_pct=None, source_dataset="development",
        reference_id="ref1", publication_year=None, reference_condition_series_id="series1",
        scaffold_key="", scaffold_tokens=(), signature=asdict(analysis.reaction_signature),
        reaction_core={}, recipe_id="recipe1", recipe_core_id="core1",
        resolved_recipe=recipe or {"temperature_c": 0.0}, condition_uncertain=True,
        chemistry_status="verified", condition_status="resolved_partial",
        condition_stage_status="single_stage", outcome_status="missing",
        record_schema_version="10.3", converter_definition_version="generic_conversion.v10.3",
        precedent_tier=PrecedentTier.TRUSTED, core_eligibility_definition_version="core_eligibility.v1@1.0",
    )


def test_comparison_keeps_duplicates_source_counts_missingness_and_structural_differences() -> None:
    row = precedent()
    original = deepcopy(row)
    rows = (row, replace(row, observation_id="o2"), replace(row, observation_id="o3", reference_id=""))
    result = compare_condition_evidence("CCCBr.N>>CCCN", rows)
    assert result.observation_count == 3
    assert result.distinct_reference_count == 1
    assert result.missing_reference_observation_ids == ("o3",)
    entry = result.precedents[0]
    assert "temperature_c" not in entry["missing_operating_fields"]
    assert "time_h" in entry["missing_operating_fields"]
    assert entry["structural_comparison"]["formed_bond_types"]["shared"]
    assert row == original
    assert entry["observation"]["reaction_smiles"] == row.reaction_smiles


def test_unresolved_query_keeps_all_observations_and_abstains_from_comparison() -> None:
    result = compare_condition_evidence("CC>>CC", (precedent(),))
    assert len(result.precedents) == 1
    entry = result.precedents[0]
    assert entry["structural_comparison"] is None
    assert entry["compatibility"]["status"] == "unknown"
    assert not entry["compatibility"]["hard_conflicts"]

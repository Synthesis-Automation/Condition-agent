from condition_recommender.conversion.generic import convert_record
from condition_recommender.conversion.input_schema import adapt_row
from condition_recommender.generic_api import GenericConditionRecommender
from condition_recommender.generic_indexing import build_generic_index


_MAPPED_PRECEDENT = "[CH3:1][CH3:2].[NH2:3][CH3:4]>>[CH3:1][CH2:2][NH2:3]"
_CONDITION_SUPPLIED_FLUORINATION = (
    "CC(C)(I)CCOS(=O)(=O)c1ccc(F)cc1>>CC(C)(F)CCOS(=O)(=O)c1ccc(F)cc1"
)


def _precedent(index: int) -> dict:
    record = adapt_row(
        {
            "reaction_id": f"precedent-{index}",
            "reaction_smiles": _MAPPED_PRECEDENT,
            "yield_pct": "80",
            "catalyst_cas": "14221-01-3",
            "reagent_cas": "584-08-7",
            "solvent_cas": "108-88-3",
            "reference": f"Independent reference {index}",
        },
        source_dataset="fallback-test",
        source_path="fallback.csv",
        source_row_number=index,
    )
    return convert_record(record).to_dict()


def _fluorination_precedent(index: int, fluoride_cas: str) -> dict:
    record = adapt_row(
        {
            "reaction_id": f"fluorination-{index}",
            "reaction_smiles": _CONDITION_SUPPLIED_FLUORINATION,
            "yield_pct": "75",
            "reagent_cas": fluoride_cas,
            "solvent_cas": "108-88-3",
            "reference": f"Fluorination reference {index}",
        },
        source_dataset="fallback-test",
        source_path="fluorination.csv",
        source_row_number=index,
    )
    return convert_record(record).to_dict()


def test_unresolved_query_retrieves_supported_structural_analogues() -> None:
    index = build_generic_index([_precedent(1), _precedent(2)])

    result = GenericConditionRecommender(index).recommend("CC.CN>>CCN")

    assert result.valid
    assert result.query_signature_id is None
    assert result.query_fallback_descriptor_id.startswith("RFD1:")
    assert result.recommendation_mode == "unverified_structure_fallback"
    assert result.retrieval_definition_version == "1.1"
    assert result.retrieval_level == "unverified_structure_fallback"
    assert result.independent_compatible_candidate_count == 2
    assert "UNVERIFIED_REACTION_FALLBACK_USED" in result.warnings
    assert "FALLBACK_RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW" in result.warnings
    assert result.recommendations
    recommendation = result.recommendations[0]
    assert (
        recommendation.score_trace.definition_versions["fallback_retrieval.v1"] == "1.1"
    )
    assert any(
        "atom correspondence and bond edits are not verified" in caution
        for caution in recommendation.cautions
    )


def test_contradicted_query_is_blocked_before_fallback_retrieval() -> None:
    index = build_generic_index([_precedent(1), _precedent(2)])

    result = GenericConditionRecommender(index).recommend(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1"
    )

    assert not result.valid
    assert result.recommendation_mode == "abstained"
    assert result.error == "QUERY_NOT_ELIGIBLE_FOR_UNVERIFIED_FALLBACK"
    assert "FALLBACK_BLOCKED:contradicted_or_incomplete_structure" in result.warnings


def test_single_moderate_analogue_does_not_clear_support_gate() -> None:
    index = build_generic_index([_precedent(1)])

    result = GenericConditionRecommender(index).recommend("CC.CN>>CCN")

    assert not result.valid
    assert result.error == "NO_SAFE_FALLBACK_PRECEDENT"
    assert result.retrieval_level == "insufficient_safe_fallback_support"
    assert result.recommendations == ()


def test_invalid_atom_mapping_is_not_rescued_by_fallback() -> None:
    index = build_generic_index([_precedent(1), _precedent(2)])

    result = GenericConditionRecommender(index).recommend(
        "[CH3:1][Br:1].[NH2:2]C>>[CH3:1][NH:2]C"
    )

    assert not result.valid
    assert result.error == "QUERY_NOT_ELIGIBLE_FOR_UNVERIFIED_FALLBACK"
    assert "FALLBACK_BLOCKED:invalid_atom_mapping" in result.warnings


def test_condition_supplied_fluoride_partial_precedents_are_retrievable() -> None:
    records = [
        _fluorination_precedent(1, "429-41-4"),
        _fluorination_precedent(2, "13400-13-0"),
    ]

    for record in records:
        assert record["reaction_signature"] is None
        assert record["chemistry_status"] == "review"
        assert record["index_eligibility"] == "eligible"
        assert record["fallback_descriptor"]["retrieval_eligible"]
        assert record["fallback_descriptor"]["required_condition_source_elements"] == (
            "F",
        )

    index = build_generic_index(records)
    assert len(index.rows) == 2
    assert not index.bond_edits

    result = GenericConditionRecommender(index).recommend(
        _CONDITION_SUPPLIED_FLUORINATION
    )

    assert result.valid
    assert result.recommendation_mode == "unverified_structure_fallback"
    assert result.candidate_count == 2
    assert "QUERY_REQUIRES_CONDITION_SOURCE:F" in result.warnings
    assert "CONDITION_SUPPLIED_FRAGMENT_FALLBACK_USED:F" in result.warnings
    assert len(result.recommendations) == 2
    assert all(
        any(
            evidence.startswith("condition_source_supported:fluoride_source:")
            for evidence in recommendation.compatibility_evidence
        )
        for recommendation in result.recommendations
    )


def test_partial_fluorination_without_fluoride_condition_stays_ineligible() -> None:
    record = _fluorination_precedent(1, "584-08-7")

    assert record["chemistry_status"] == "rejected"
    assert record["index_eligibility"] == "ineligible"
    assert any(
        reason.startswith("condition_source_missing:fluoride_source:F")
        for reason in record["admission_reasons"]
    )
    assert not build_generic_index([record]).rows

from dataclasses import asdict

from reactive_taxonomy import featurize_reaction

from condition_recommender.conversion.generic import convert_record
from condition_recommender.conversion.input_schema import adapt_row
from condition_recommender.fallback_similarity import assess_fallback_similarity
from condition_recommender.generic_api import GenericConditionRecommender
from condition_recommender.generic_indexing import build_generic_index


_MAPPED_PRECEDENT = "[CH3:1][CH3:2].[NH2:3][CH3:4]>>[CH3:1][CH2:2][NH2:3]"
_CONDITION_SUPPLIED_FLUORINATION = (
    "CC(C)(I)CCOS(=O)(=O)c1ccc(F)cc1>>CC(C)(F)CCOS(=O)(=O)c1ccc(F)cc1"
)
_CONDITION_SUPPLIED_IODINATION = "CCCCCCCCC(C)(F)CC>>CCCCCCCCC(C)(I)CC"


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


def _iodination_precedent(index: int, iodine_cas: str) -> dict:
    record = adapt_row(
        {
            "reaction_id": f"iodination-{index}",
            "reaction_smiles": _CONDITION_SUPPLIED_IODINATION,
            "yield_pct": "90",
            "reagent_cas": iodine_cas,
            "solvent_cas": "110-54-3",
            "reference": f"Iodination reference {index}",
        },
        source_dataset="fallback-test",
        source_path="iodination.csv",
        source_row_number=index,
    )
    return convert_record(record).to_dict()


def _analogue_precedent(
    index: int,
    reaction_smiles: str,
    *,
    reference: str = "Shared analogue reference",
) -> dict:
    record = adapt_row(
        {
            "reaction_id": f"analogue-{index}",
            "reaction_smiles": reaction_smiles,
            "yield_pct": "90",
            "reagent_cas": "97-93-8",
            "solvent_cas": "110-54-3",
            "reference": reference,
        },
        source_dataset="fallback-test",
        source_path="analogues.csv",
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
    assert result.retrieval_definition_version == "1.3"
    assert result.retrieval_level == "unverified_structure_fallback"
    assert result.independent_compatible_candidate_count == 2
    assert "UNVERIFIED_REACTION_FALLBACK_USED" in result.warnings
    assert "FALLBACK_RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW" in result.warnings
    assert result.recommendations
    recommendation = result.recommendations[0]
    assert (
        recommendation.score_trace.definition_versions["fallback_retrieval.v1"] == "1.3"
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


def test_single_exact_analogue_uses_limited_support_route() -> None:
    index = build_generic_index([_precedent(1)])

    result = GenericConditionRecommender(index).recommend("CC.CN>>CCN")

    assert result.valid
    assert result.retrieval_level == "unverified_structure_fallback_limited_support"
    assert "LIMITED_PRECEDENT_SUPPORT" in result.warnings
    assert len(result.recommendations) == 1


def test_invalid_atom_mapping_is_not_rescued_by_fallback() -> None:
    index = build_generic_index([_precedent(1), _precedent(2)])

    result = GenericConditionRecommender(index).recommend(
        "[CH3:1][Br:1].[NH2:2]C>>[CH3:1][NH:2]C"
    )

    assert not result.valid
    assert result.error == "QUERY_NOT_ELIGIBLE_FOR_UNVERIFIED_FALLBACK"
    assert "FALLBACK_BLOCKED:invalid_atom_mapping" in result.warnings


def test_exploratory_fluoride_partial_precedents_are_retrievable() -> None:
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
    assert "QUERY_PRODUCT_ATOM_SOURCE_UNVERIFIED:F" in result.warnings
    assert "EXPLORATORY_PARTIAL_CORRESPONDENCE_FALLBACK_USED:F" in result.warnings
    assert len(result.recommendations) == 2


def test_partial_fluorination_keeps_unverified_source_as_a_warning() -> None:
    record = _fluorination_precedent(1, "584-08-7")

    assert record["chemistry_status"] == "review"
    assert record["index_eligibility"] == "eligible"
    assert "exploratory_partial_product_correspondence" in record["admission_reasons"]
    assert (
        "PRODUCT_ATOM_SOURCE_UNVERIFIED:F" in record["fallback_descriptor"]["warnings"]
    )
    assert len(build_generic_index([record]).rows) == 1


def test_exploratory_iodine_partial_precedents_are_retrievable() -> None:
    records = [
        _iodination_precedent(1, "97-93-8"),
        _iodination_precedent(2, "584-08-7"),
    ]

    for record in records:
        assert record["reaction_signature"] is None
        assert record["chemistry_status"] == "review"
        assert record["index_eligibility"] == "eligible"
        assert record["fallback_descriptor"]["retrieval_eligible"]
        assert record["fallback_descriptor"]["required_condition_source_elements"] == (
            "I",
        )

    result = GenericConditionRecommender(build_generic_index(records)).recommend(
        _CONDITION_SUPPLIED_IODINATION
    )

    assert result.valid
    assert result.recommendation_mode == "unverified_structure_fallback"
    assert result.candidate_count == 2
    assert result.retrieval_definition_version == "1.3"
    assert "QUERY_PRODUCT_ATOM_SOURCE_UNVERIFIED:I" in result.warnings
    assert "EXPLORATORY_PARTIAL_CORRESPONDENCE_FALLBACK_USED:I" in result.warnings
    assert len(result.recommendations) == 2
    assert all(
        recommendation.similarity_score == 1.0
        for recommendation in result.recommendations
    )


def test_partial_iodination_keeps_unverified_source_as_a_warning() -> None:
    record = _iodination_precedent(1, "97-93-8")

    assert record["chemistry_status"] == "review"
    assert record["index_eligibility"] == "eligible"
    assert "exploratory_partial_product_correspondence" in record["admission_reasons"]
    assert (
        "PRODUCT_ATOM_SOURCE_UNVERIFIED:I" in record["fallback_descriptor"]["warnings"]
    )
    assert len(build_generic_index([record]).rows) == 1


def test_local_center_signature_ranks_nonidentical_same_edit_analogues() -> None:
    query = asdict(
        featurize_reaction(_CONDITION_SUPPLIED_IODINATION).fallback_descriptor
    )
    secondary = asdict(
        featurize_reaction("CCCCCCCCCCC(C)F>>CCCCCCCCCCC(C)I").fallback_descriptor
    )
    bridged_tertiary = asdict(
        featurize_reaction(
            "FC12CC3CC(CC(C3)C1)C2>>IC12CC3CC(CC(C3)C1)C2"
        ).fallback_descriptor
    )

    assert "center:substitution_class:tertiary" in query["reaction_center_core_tokens"]
    assert (
        "center:substitution_class:secondary"
        in secondary["reaction_center_core_tokens"]
    )
    assert (
        "center:substitution_class:tertiary"
        in bridged_tertiary["reaction_center_core_tokens"]
    )
    assert query["reaction_center_radius_1_tokens"]
    assert query["reaction_center_radius_2_tokens"]
    assert query["reaction_center_radius_3_tokens"]

    secondary_score = assess_fallback_similarity(query, secondary).score
    bridged_score = assess_fallback_similarity(query, bridged_tertiary).score
    assert secondary_score > 0.65
    assert bridged_score > secondary_score


def test_leave_one_reaction_out_returns_nonidentical_edit_analogues() -> None:
    records = [
        _analogue_precedent(
            1,
            "CCCCCCCCCCC(C)F>>CCCCCCCCCCC(C)I",
        ),
        _analogue_precedent(
            2,
            "FC12CC3CC(CC(C3)C1)C2>>IC12CC3CC(CC(C3)C1)C2",
        ),
    ]

    result = GenericConditionRecommender(build_generic_index(records)).recommend(
        _CONDITION_SUPPLIED_IODINATION
    )

    assert result.valid
    assert result.candidate_count == 2
    assert result.independent_compatible_candidate_count == 1
    assert (
        result.retrieval_level
        == "unverified_structure_fallback_exploratory_limited_support"
    )
    assert "LIMITED_PRECEDENT_SUPPORT" in result.warnings
    assert result.recommendations
    recommendation = result.recommendations[0]
    assert set(recommendation.precedent_reaction_ids) == {
        "analogue-1",
        "analogue-2",
    }
    assert (
        recommendation.score_trace.similarity_components["reaction_center_core"] > 0.0
    )
    assert any(
        note.startswith("Reaction-center similarity:")
        for note in recommendation.explanation
    )
    assert all(
        context["reaction_center_core"]
        for context in recommendation.precedent_reaction_contexts
    )

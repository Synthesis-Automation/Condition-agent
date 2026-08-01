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
_ACYL_FLUORIDE_QUERY = (
    "N#Cc1ccc(C(=O)O)cc1>>N#Cc1ccc(C(=O)F)cc1"
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


def _acyl_fluoride_precedent(index: int, reaction_smiles: str) -> dict:
    record = adapt_row(
        {
            "reaction_id": f"acyl-fluoride-{index}",
            "reaction_smiles": reaction_smiles,
            "yield_pct": "90",
            "reagent_cas": "63517-29-3",
            "solvent_cas": "141-78-6",
            "reference": f"Acyl fluoride reference {index}",
        },
        source_dataset="fallback-test",
        source_path="acyl-fluoride.csv",
        source_row_number=index,
    )
    return convert_record(record).to_dict()


def _analogue_precedent(
    index: int,
    reaction_smiles: str,
    *,
    reference: str = "Shared analogue reference",
    reagent_cas: str = "97-93-8",
) -> dict:
    record = adapt_row(
        {
            "reaction_id": f"analogue-{index}",
            "reaction_smiles": reaction_smiles,
            "yield_pct": "90",
            "reagent_cas": reagent_cas,
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
    assert result.query_fallback_descriptor_id.startswith("RFD2:")
    assert result.recommendation_mode == "unverified_structure_fallback"
    assert result.retrieval_definition_version == "2.0"
    assert result.retrieval_level == "unverified_structure_fallback"
    assert result.independent_compatible_candidate_count == 2
    assert "UNVERIFIED_REACTION_FALLBACK_USED" in result.warnings
    assert "FALLBACK_RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW" in result.warnings
    assert result.recommendations
    recommendation = result.recommendations[0]
    assert (
        recommendation.score_trace.definition_versions["fallback_retrieval.v1"] == "2.0"
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


def test_unrestricted_fallback_bypasses_all_chemistry_gates() -> None:
    records = [
        _analogue_precedent(
            1,
            "CC(=O)O.CN>>CC(=O)NC",
            reference="Independent amide reference 1",
        ),
        _analogue_precedent(
            2,
            "CCC(=O)O.CN>>CCC(=O)NC",
            reference="Independent amide reference 2",
        ),
    ]
    index = build_generic_index(records)
    query = "c1ccc(c(N)c1)O.CC(=O)O>>Cc1nc2ccccc2o1"

    gated = GenericConditionRecommender(index).recommend(query)
    unrestricted = GenericConditionRecommender(index).recommend(
        query,
        unrestricted_fallback=True,
    )

    assert not gated.valid
    assert gated.error == "QUERY_NOT_ELIGIBLE_FOR_UNVERIFIED_FALLBACK"
    assert unrestricted.valid
    assert (
        unrestricted.recommendation_mode
        == "unrestricted_unverified_structure_fallback"
    )
    assert (
        unrestricted.retrieval_level
        == "unrestricted_unverified_structure_fallback"
    )
    assert unrestricted.recommendations
    assert "UNRESTRICTED_FALLBACK_REQUESTED" in unrestricted.warnings
    assert "UNRESTRICTED_FALLBACK_USED" in unrestricted.warnings
    assert (
        "FALLBACK_GATE_OVERRIDDEN:contradicted_or_incomplete_structure"
        in unrestricted.warnings
    )
    assert (
        "FALLBACK_SIMILARITY_SUPPORT_AND_COMPATIBILITY_GATES_BYPASSED"
        in unrestricted.warnings
    )
    assert any(
        "condition-compatibility gates were explicitly bypassed" in caution
        for recommendation in unrestricted.recommendations
        for caution in recommendation.cautions
    )


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
    assert result.retrieval_level == "source_supported_partial_transformation"
    assert result.candidate_count == 2
    assert "QUERY_PRODUCT_ATOM_SOURCE_UNVERIFIED:F" in result.warnings
    assert "EXPLORATORY_PARTIAL_CORRESPONDENCE_FALLBACK_USED:F" in result.warnings
    assert len(result.recommendations) == 2


def test_partial_fluorination_without_capable_source_is_not_indexed() -> None:
    record = _fluorination_precedent(1, "584-08-7")

    assert record["chemistry_status"] == "rejected"
    assert record["index_eligibility"] == "ineligible"
    assert "missing_condition_fragment_source" in record["admission_reasons"]
    assert record["fragment_source_support"][0]["status"] == "unsupported"
    assert len(build_generic_index([record]).rows) == 0


def test_acyl_fluoride_query_prefers_exact_partial_transformation_to_amidation() -> None:
    acyl_fluorides = [
        _acyl_fluoride_precedent(
            1,
            "COc1cccc(C(=O)O)c1>>COc1cccc(C(=O)F)c1",
        ),
        _acyl_fluoride_precedent(
            2,
            "O=C(O)c1ccccc1>>O=C(F)c1ccccc1",
        ),
    ]
    amide = _analogue_precedent(
        3,
        "COC(=O)c1ccc(C(=O)O)cc1.CN>>"
        "COC(=O)c1ccc(C(=O)NC)cc1",
        reference="Amidation reference",
    )

    result = GenericConditionRecommender(
        build_generic_index([*acyl_fluorides, amide])
    ).recommend(_ACYL_FLUORIDE_QUERY)

    assert result.valid
    assert result.retrieval_level == "source_supported_partial_transformation"
    assert result.candidate_count == 2
    assert "SOURCE_SUPPORTED_PARTIAL_TRANSFORMATION_USED" in result.warnings
    hit_ids = {
        reaction_id
        for recommendation in result.recommendations
        for reaction_id in recommendation.precedent_reaction_ids
    }
    assert hit_ids == {"acyl-fluoride-1", "acyl-fluoride-2"}
    assert all(
        recommendation.precedent_reaction_contexts[0][
            "fragment_source_support"
        ][0]["status"]
        == "supported"
        for recommendation in result.recommendations
    )


def test_exploratory_iodine_partial_precedents_are_retrievable() -> None:
    records = [
        _iodination_precedent(1, "7553-56-2"),
        _iodination_precedent(2, "7681-11-0"),
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
    assert result.retrieval_definition_version == "2.0"
    assert result.retrieval_level == "source_supported_partial_transformation"
    assert "QUERY_PRODUCT_ATOM_SOURCE_UNVERIFIED:I" in result.warnings
    assert "EXPLORATORY_PARTIAL_CORRESPONDENCE_FALLBACK_USED:I" in result.warnings
    assert len(result.recommendations) == 2
    assert all(
        recommendation.similarity_score == 1.0
        for recommendation in result.recommendations
    )


def test_partial_iodination_without_capable_source_is_not_indexed() -> None:
    record = _iodination_precedent(1, "97-93-8")

    assert record["chemistry_status"] == "rejected"
    assert record["index_eligibility"] == "ineligible"
    assert "missing_condition_fragment_source" in record["admission_reasons"]
    assert len(build_generic_index([record]).rows) == 0


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
            reagent_cas="7553-56-2",
        ),
        _analogue_precedent(
            2,
            "FC12CC3CC(CC(C3)C1)C2>>IC12CC3CC(CC(C3)C1)C2",
            reagent_cas="7553-56-2",
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
        == "source_supported_partial_transformation_limited_support"
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

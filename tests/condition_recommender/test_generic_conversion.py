import csv
import json
from pathlib import Path

import condition_recommender.conversion.generic as generic_module
import condition_recommender.conversion.sharded as sharded_module
import pytest
from condition_recommender.conversion.concise_review import (
    CONCISE_REACTION_REVIEW_FIELDS,
    ConciseReviewConversionCancelled,
    concise_reaction_review_row,
    convert_dataset_folder_to_concise_review_csv,
    export_concise_reaction_review_csv,
)
from condition_recommender.conversion.engine import convert_datasets
from condition_recommender.conversion.flatten import flatten_generic_record
from condition_recommender.conversion.generic import (
    GenericConversionCache,
    convert_record,
)
from condition_recommender.conversion.input_schema import (
    adapt_row,
    discover_conversion_datasets,
    discover_csv_datasets,
)
from condition_recommender.conversion.sharded import (
    ShardedConversionCancelled,
    convert_datasets_sharded,
    validate_sharded_conversion,
)
from condition_recommender.generic_indexing import (
    build_generic_index,
    load_generic_index,
)
from condition_recommender.models import (
    AdmissionTier,
    ChemistryStatus,
    ConditionStageStatus,
    ConditionStatus,
    IndexEligibility,
    OutcomeStatus,
)


def _raw(
    reaction: str,
    *,
    reaction_id: str = "record-1",
    yield_pct: str = "80",
    catalyst_cas: str = "14221-01-3",
    reagent_cas: str = "584-08-7",
    solvent_cas: str = "108-88-3",
    stages: str = "1",
):
    row = {
        "reaction_id": reaction_id,
        "reaction_type": "untrusted-source-label",
        "reaction_smiles": reaction,
        "yield_pct": yield_pct,
        "catalyst_cas": catalyst_cas,
        "reagent_cas": reagent_cas,
        "solvent_cas": solvent_cas,
        "reference": "reference",
        "stages": stages,
    }
    return adapt_row(
        row,
        source_dataset="mixed",
        source_path="mixed.csv",
        source_row_number=2,
    )


def _csv_row(
    reaction_id: str,
    reaction: str,
    *,
    reaction_type: str,
    solvent_cas: str = "108-88-3",
) -> dict[str, str]:
    return {
        "reaction_id": reaction_id,
        "reaction_type": reaction_type,
        "reaction_smiles": reaction,
        "yield_pct": "80",
        "temperature_c": "",
        "time_h": "",
        "reference": "reference",
        "reactant_cas": "",
        "product_cas": "",
        "reagent_cas": "584-08-7",
        "catalyst_cas": "",
        "solvent_cas": solvent_cas,
        "experimental_procedure": "",
        "stages": "1",
        "steps": "1",
        "notes": "",
    }


def test_exact_signature_is_verified_without_trusting_source_family() -> None:
    record = convert_record(_raw("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"))

    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.source_declared_family == "untrusted-source-label"
    assert record.named_family == "suzuki_miyaura"
    assert record.reaction_signature is not None
    assert record.canonical_reaction_id.startswith("CRX1:")
    assert record.observation_id.startswith("OBS1:")
    assert record.raw_recipe_id.startswith("RAWCOND1:")
    assert record.resolved_recipe_id.startswith("RCR1:")
    assert record.resolved_recipe["catalysts"][0]["primary_role"] == ("metal_catalyst")
    assert record.resolved_recipe["bases"][0]["primary_role"] == "base"
    assert record.condition_resolution["component_count"] == 3
    assert record.schema_version == "9.0"
    assert record.converter_definition_version == "generic_conversion.v9.0"
    assert record.reaction_signature["schema_version"] == "3.4"
    assert record.reaction_observation is not None
    assert record.reaction_interpretation is not None
    assert record.reaction_interpretation["primary_pattern_id"] == (
        "organoboron_c_c_coupling_like"
    )
    assert record.reaction_signature["topology"]["reaction_scope"] == ("intermolecular")
    assert record.reference_id.startswith("REF1:")
    assert record.reference_identity["resolution_status"] == "bibliographic_text"
    assert record.chemistry_status == ChemistryStatus.VERIFIED
    assert record.condition_status == ConditionStatus.RESOLVED_COMPLETE
    assert record.outcome_status == OutcomeStatus.USABLE
    assert record.index_eligibility == IndexEligibility.ELIGIBLE
    assert record.resolved_recipe_core_id.startswith("RCORE1:")
    assert record.reference_condition_series_id.startswith("RCS1:")
    payload = record.to_dict()
    assert payload["reaction_label"]["text"]
    assert payload["reaction_label"]["basis"]
    assert set(payload) & {
        "reaction_label_status",
        "reaction_display_label",
    } == set()


def test_mapped_unknown_reaction_serializes_reaction_core_for_review() -> None:
    reaction = (
        "[CH3:1][OH:2].O[CH3:5]."
        "[CH:3](=[O:4])[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
        ">>[CH3:1][O:2][CH:3]([O:4][CH3:5])"
        "[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
    )

    record = convert_record(_raw(reaction))

    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.reaction_core is not None
    assert record.reaction_core["schema_version"] == "3.0"
    assert record.reaction_core["algorithm_version"] == (
        "reaction_core_projection.v11"
    )
    assert record.reaction_core["shape_core_key"].startswith("RSH2:")
    assert "generic_label" not in record.reaction_core
    assert "presentation" not in record.reaction_core


def test_cycloaddition_ring_observation_is_chemist_readable_in_review() -> None:
    record = convert_record(_raw("CC#C.CN=[N+]=[N-]>>Cc1nnn(C)c1"))

    assert record.reaction_signature is not None
    topology = record.reaction_signature["topology"]
    assert topology["formed_ring_sizes"] == (5,)
    assert topology["ring_count_delta"] == 1
    assert len(topology["ring_changes"]) == 1
    change = topology["ring_changes"][0]
    assert change["element_sequence"] == ("C", "C", "N", "N", "N")
    assert change["formed_bond_types"] == ("C-N", "C-N")
    review = concise_reaction_review_row(record.to_dict())
    assert record.reaction_label is not None
    assert record.reaction_label["basis"] == "ring_topology"
    assert review["reaction_label"].endswith(
        "aromatic 5-membered C₂N₃ ring"
    )
    assert review["ring_change_count"] == "1"
    assert review["formed_ring_sizes"] == "5"
    assert review["ring_count_delta"] == "1"
    assert review["ring_change_summaries"] == (
        "aromatic 5-membered C2N3 ring; formed=C-N,C-N"
    )


def test_ambiguous_alcohol_displacement_retains_structural_hypotheses() -> (
    None
):
    reaction = "CNCC[C@H](O)c1ccccc1.Cc1ccccc1O>>CNCC[C@@H](Oc1ccccc1C)c1ccccc1"

    record = convert_record(_raw(reaction))

    assert record.chemistry_status == ChemistryStatus.REVIEW
    assert record.index_eligibility == IndexEligibility.REVIEW_ONLY
    assert record.reaction_signature is None
    assert record.named_family is None
    assert record.transformation_class is None
    assert len(record.reaction_edit_hypotheses) == 2
    assert all(
        hypothesis["stereo_changes"][0]["change_type"]
        == "descriptor_changed"
        for hypothesis in record.reaction_edit_hypotheses
    )
    review = concise_reaction_review_row(record.to_dict())
    assert review["edit_hypothesis_count"] == "2"


def test_conversion_cache_reuses_deterministic_reaction_analysis(
    monkeypatch,
) -> None:
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    original = generic_module.featurize_reaction
    calls = []

    def counted(value):
        calls.append(value)
        return original(value)

    monkeypatch.setattr(generic_module, "featurize_reaction", counted)
    cache = GenericConversionCache()
    convert_record(_raw(reaction, reaction_id="first"), cache=cache)
    convert_record(_raw(reaction, reaction_id="second"), cache=cache)

    assert calls == [reaction]


def test_mapped_unknown_family_signature_is_verified() -> None:
    record = convert_record(
        _raw(
            "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]",
            catalyst_cas="",
        )
    )

    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.named_family is None
    assert record.evidence_quality == "validated_atom_mapping"
    assert record.reaction_label["text"] == "CH₂=CH₂ → CH₃–CH₃"
    assert record.reaction_label["status"] == "structural_equation"
    assert record.reaction_label["basis"] == "local_context"
    assert record.reaction_signature["order_changes"] == ("C-C:DOUBLE>SINGLE",)
    assert record.resolved_recipe_id.startswith("RCR1:")


def test_exact_multi_event_signature_is_verified() -> None:
    record = convert_record(_raw("CO.CS.Fc1ccc(F)cc1>>COc1ccc(SC)cc1"))

    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.evidence_quality == "global_atom_correspondence"
    assert record.named_family is None
    assert record.reaction_signature is not None
    assert record.reaction_signature["event_count"] == 2
    assert record.reaction_signature["event_scope"] == "multi_event"
    assert record.reaction_label is not None
    assert record.reaction_label["event_count"] == 2
    assert "event_labels" not in record.reaction_label
    assert " + " in record.reaction_label["text"]


def test_global_correspondence_is_verified_structural_evidence() -> None:
    record = convert_record(_raw("CC=O.CN>>CC=NC"))

    assert record.evidence_quality == "global_atom_correspondence"
    assert record.reaction_signature is not None
    assert record.reaction_completeness is not None
    assert record.reaction_completeness["status"] == "verified"
    assert record.chemistry_status == ChemistryStatus.VERIFIED
    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.index_eligibility == IndexEligibility.ELIGIBLE
    assert record.admission_reasons == (
        "verified_reaction_signature_with_conditions",
    )
    index = build_generic_index([record.to_dict()])
    assert len(index.rows) == 1
    assert index.rows[0].chemistry_status == "verified"


def test_optional_patterns_do_not_create_stereochemical_conflicts() -> None:
    reaction = (
        "COc1ccc(C)cc1B(O)O."
        "O=C1c2ccccc2C(=O)N1C/C=C(/Br)c1ccccc1>>"
        "COc1ccc(C)cc1/C(=C\\CN1C(=O)c2ccccc2C1=O)c1ccccc1"
    )
    record = convert_record(_raw(reaction))

    assert record.reaction_signature is not None
    assert record.evidence_quality == "global_atom_correspondence"
    assert record.chemistry_status == ChemistryStatus.VERIFIED
    assert record.index_eligibility == IndexEligibility.ELIGIBLE
    assert "conflicting_stereochemical_evidence" not in record.admission_reasons
    assert len(build_generic_index([record.to_dict()]).rows) == 1


def test_unresolved_record_is_rejected_even_with_fallback_descriptor() -> None:
    record = convert_record(_raw("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1"))

    assert record.admission_tier == AdmissionTier.REJECTED
    assert record.admission_reasons == (
        "no_usable_transformation_evidence",
        "reaction_completeness_unresolved",
    )
    assert record.fallback_descriptor is not None
    assert record.fallback_descriptor["schema_version"] == "3.0"
    assert record.fallback_descriptor["retrieval_eligible"]


def test_ambiguous_edit_hypotheses_are_serialized_but_not_indexed() -> None:
    record = convert_record(
        _raw(
            "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
            ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
        )
    )

    assert record.reaction_signature is None
    assert len(record.reaction_edit_hypotheses) == 2
    assert {
        value["provider"] for value in record.reaction_evidence_candidates
    } == {"global_correspondence"}
    assert all(
        value["hypothesis_id"].startswith("REH1:")
        for value in record.reaction_edit_hypotheses
    )
    assert record.index_eligibility == IndexEligibility.REVIEW_ONLY
    assert not build_generic_index([record.to_dict()]).rows

    review = concise_reaction_review_row(record.to_dict())
    assert review["edit_hypothesis_count"] == "2"
    assert "edit_hypothesis_ids" not in review


def test_ambiguous_record_serializes_ineligible_fallback_descriptor() -> None:
    record = convert_record(_raw("CC.CN>>CCN"))

    assert record.reaction_signature is None
    assert record.fallback_descriptor is not None
    assert record.fallback_descriptor["descriptor_id"].startswith("RFD3:")
    assert record.fallback_descriptor["evidence_mode"] == ("structure_inventory_only")
    assert not record.fallback_descriptor["retrieval_eligible"]
    assert "ambiguous_edit_hypotheses" in (
        record.fallback_descriptor["ineligibility_reasons"]
    )
    assert record.fallback_descriptor["reactant_component_tokens"] == (
        "CC",
        "CN",
    )


def test_unaccounted_product_atoms_are_ineligible_and_serialized() -> None:
    record = convert_record(_raw("[CH2:1]=[CH2:2]>>[CH3:1][CH2:2]C"))

    assert record.admission_tier == AdmissionTier.REJECTED
    assert record.chemistry_status == ChemistryStatus.REJECTED
    assert record.index_eligibility == IndexEligibility.INELIGIBLE
    assert record.admission_reasons == (
        "suspected_missing_reactant",
        "unaccounted_product_heavy_atoms",
    )
    assert record.reaction_signature is None
    assert record.reaction_completeness is not None
    assert record.reaction_completeness["status"] == "incomplete"
    assert record.reaction_completeness["product_element_excess"] == {"C": 1}


def test_partial_product_observation_is_serialized_but_not_indexed() -> None:
    record = convert_record(_raw("O=C(O)c1cccc(I)c1>>O=C(Cl)c1cccc(I)c1"))

    assert record.admission_tier == AdmissionTier.REJECTED
    assert record.chemistry_status == ChemistryStatus.REJECTED
    assert record.index_eligibility == IndexEligibility.INELIGIBLE
    assert record.reaction_signature is None
    assert record.reaction_completeness["status"] == "incomplete"
    assert record.reaction_label["text"] == "R–C(=O)–O → R–C(=O)–Cl"
    assert record.reaction_label["status"] == "partial_product_correspondence"
    assert record.transformation_class == "acyl_heteroatom_substitution"
    assert record.transformation_confidence == 0.8
    assert record.partial_product_transformation is not None
    assert record.partial_product_transformation["missing_product_atom_elements"] == (
        "Cl",
    )

    review = concise_reaction_review_row(record.to_dict())
    assert review["reaction_label_status"] == "partial_product_correspondence"
    assert review["reaction_label_basis"] == "partial_product_correspondence"
    assert review["partial_transformation_class"] == "acyl_heteroatom_substitution"
    assert review["missing_product_atom_elements"] == "Cl"
    assert "PRODUCT_ATOM_SOURCE_UNRESOLVED:Cl" in review["warnings"]


def test_multi_atom_product_origin_gap_is_serialized_for_review() -> None:
    record = convert_record(_raw("Brc1ccccc1>>N#Cc1ccccc1"))

    assert record.chemistry_status == ChemistryStatus.REJECTED
    assert record.index_eligibility == IndexEligibility.INELIGIBLE
    assert record.reaction_signature is None
    partial = record.partial_product_transformation
    assert partial is not None
    assert partial["schema_version"] == "2.0"
    assert partial["transformation_type"] == ("attachment_fragment_replacement")
    fragment = partial["installed_fragment"]
    assert fragment["canonical_fragment_smiles"] == "C#N"
    assert fragment["rooted_fragment_smiles"] == "*C#N"
    assert fragment["fragment_key"].startswith("PFG1:")
    assert fragment["source_status"] == "unresolved"
    assert fragment["source_candidates"] == ()
    provenance = partial["product_atom_provenance"]
    assert len(provenance) == 8
    assert (
        sum(item["source_kind"] == "reactant_correspondence" for item in provenance)
        == 6
    )
    assert sum(item["source_kind"] == "unresolved" for item in provenance) == 2

    review = concise_reaction_review_row(record.to_dict())
    assert review["removed_fragment"] == "Br"
    assert review["installed_fragment"] == "C#N"
    assert "installed_fragment_key" not in review
    assert review["fragment_source_status"] == "unresolved"
    assert review["fragment_source_candidates"] == ""
    assert review["missing_product_atom_elements"] == "C; N"


def test_inconsistent_product_mapping_is_review_only() -> None:
    record = convert_record(_raw("[CH3:1][OH:2]>>[CH2:1]=[O:3]"))

    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.chemistry_status == ChemistryStatus.REVIEW
    assert record.index_eligibility == IndexEligibility.REVIEW_ONLY
    assert "inconsistent_product_atom_mapping" in record.admission_reasons
    assert record.reaction_completeness is not None
    assert (
        "PRODUCT_MAPS_MISSING_FROM_REACTANTS"
        in record.reaction_completeness["warnings"]
    )


def test_unresolved_transformation_and_missing_conditions_are_rejected() -> None:
    unresolved = convert_record(_raw("CC>>CCC"))
    no_conditions = convert_record(
        _raw(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            catalyst_cas="",
            reagent_cas="",
            solvent_cas="",
        )
    )

    assert unresolved.admission_tier == AdmissionTier.REJECTED
    assert unresolved.admission_reasons == (
        "suspected_missing_reactant",
        "unaccounted_product_heavy_atoms",
    )
    assert no_conditions.admission_tier == AdmissionTier.REJECTED
    assert no_conditions.admission_reasons == ("no_condition_identifiers",)


def test_unresolved_condition_identifier_is_retained_for_review() -> None:
    record = convert_record(
        _raw(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            catalyst_cas="not-a-cas",
            reagent_cas="",
            solvent_cas="",
        )
    )

    assert record.admission_tier == AdmissionTier.REVIEW
    assert "condition_identifier_uncertainty" in record.admission_reasons
    assert "unresolved_condition_identifiers" in record.admission_reasons
    component = record.condition_resolution["components"][0]
    assert component["raw_identifier"] == "not-a-cas"
    assert component["status"] == "invalid_identifier"
    assert record.chemistry_status == ChemistryStatus.VERIFIED
    assert record.condition_status == ConditionStatus.RESOLVED_PARTIAL
    assert record.index_eligibility == IndexEligibility.REVIEW_ONLY


def test_missing_yield_does_not_discard_a_usable_condition_precedent() -> None:
    record = convert_record(
        _raw(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            yield_pct="",
        )
    )

    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.chemistry_status == ChemistryStatus.VERIFIED
    assert record.condition_status == ConditionStatus.RESOLVED_COMPLETE
    assert record.outcome_status == OutcomeStatus.MISSING
    assert record.index_eligibility == IndexEligibility.ELIGIBLE
    index = build_generic_index([record.to_dict()])
    assert len(index.rows) == 1
    assert index.rows[0].yield_pct is None


def test_valid_unknown_condition_identity_is_retained_for_retrieval() -> None:
    record = convert_record(
        _raw(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            catalyst_cas="999999-99-4",
            reagent_cas="",
            solvent_cas="",
        )
    )

    assert record.condition_status == ConditionStatus.UNRESOLVED_RETAINED
    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.index_eligibility == IndexEligibility.ELIGIBLE
    assert len(build_generic_index([record.to_dict()]).rows) == 1


def test_resolved_unassigned_multistage_conditions_are_indexable_with_caution() -> None:
    record = convert_record(
        _raw(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            stages="2",
        )
    )

    assert record.condition_status == ConditionStatus.RESOLVED_COMPLETE
    assert record.condition_stage_status == ConditionStageStatus.UNASSIGNED_MULTISTAGE
    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.index_eligibility == IndexEligibility.ELIGIBLE
    index = build_generic_index([record.to_dict()])
    assert len(index.rows) == 1
    assert index.rows[0].condition_stage_status == "unassigned_multistage"
    assert index.rows[0].condition_uncertain


def test_unresolved_unassigned_multistage_conditions_remain_review_only() -> None:
    record = convert_record(
        _raw(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            catalyst_cas="999999-99-4",
            reagent_cas="",
            solvent_cas="",
            stages="2",
        )
    )

    assert record.condition_status == ConditionStatus.UNRESOLVED_RETAINED
    assert record.condition_stage_status == ConditionStageStatus.UNASSIGNED_MULTISTAGE
    assert record.index_eligibility == IndexEligibility.REVIEW_ONLY
    assert len(build_generic_index([record.to_dict()]).rows) == 0


def test_mixed_engine_writes_canonical_jsonl_and_review_views(tmp_path) -> None:
    dataset = tmp_path / "mixed.csv"
    rows = [
        _csv_row(
            "exact",
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            reaction_type="Suzuki",
        ),
        _csv_row(
            "mapped",
            "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]",
            reaction_type="Unknown",
        ),
        _csv_row(
            "review",
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1",
            reaction_type="Suzuki",
        ),
        _csv_row("rejected", "CC>>CCC", reaction_type="Unknown"),
    ]
    with dataset.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    output = tmp_path / "converted"
    report = convert_datasets(dataset, output)

    assert report["tier_counts"] == {
        "verified": 2,
        "review": 0,
        "rejected": 2,
    }
    assert report["signature_count"] == 2
    assert report["resolved_recipe_count"] > 0
    assert sum(report["role_confidence_counts"].values()) > 0
    records = [
        json.loads(line)
        for line in (output / "records.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    assert len(records) == 4
    assert records[1]["named_family"] is None
    assert records[1]["reaction_signature"]["order_changes"] == ["C-C:DOUBLE>SINGLE"]
    assert records[1]["resolved_recipe_id"].startswith("RCR1:")
    with (output / "verified.csv").open(encoding="utf-8-sig", newline="") as handle:
        verified = list(csv.DictReader(handle))
    assert len(verified) == 2
    assert all(row["reaction_signature_id"].startswith("RS3:") for row in verified)
    assert {row["reaction_event_count"] for row in verified} == {"1"}
    assert {row["reaction_event_scope"] for row in verified} == {"single_event"}
    assert verified[0]["reaction_scope"] == "intermolecular"
    assert {row["reaction_scope"] for row in verified} == {
        "intermolecular",
        "unimolecular",
    }
    assert json.loads((output / "conversion_report.json").read_text()) == report
    assert report["schema_version"] == "2.1"
    assert report["reaction_signature_schema_version"] == "3.4"
    assert report["reaction_scope_counts"] == {
        "intermolecular": 1,
        "unimolecular": 1,
    }
    assert report["reaction_completeness_status_counts"] == {
        "incomplete": 1,
        "unresolved": 1,
        "verified": 2,
    }
    assert (output / "conversion_report.md").exists()


def test_concise_reaction_review_export_has_only_requested_columns(
    tmp_path: Path,
) -> None:
    dataset = tmp_path / "mixed.csv"
    rows = [
        _csv_row(
            "exact",
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            reaction_type="Original Suzuki Label",
        )
    ]
    with dataset.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    converted = tmp_path / "converted"
    convert_datasets(dataset, converted)
    output = tmp_path / "concise_review.csv"

    report = export_concise_reaction_review_csv(
        converted / "records.jsonl",
        output,
    )

    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        review_rows = list(csv.DictReader(handle))
    assert report["schema_version"] == "11.0"
    assert report["row_count"] == 1
    assert tuple(review_rows[0]) == CONCISE_REACTION_REVIEW_FIELDS
    label_index = CONCISE_REACTION_REVIEW_FIELDS.index(
        "reaction_label"
    )
    assert CONCISE_REACTION_REVIEW_FIELDS[:2] == (
        "canonical_reaction_smiles",
        "reaction_label",
    )
    assert CONCISE_REACTION_REVIEW_FIELDS[2:10] == (
        "primary_reaction_pattern",
        "primary_reaction_pattern_count",
        "reaction_pattern_matches",
        "co_occurring_reaction_patterns",
        "identified_reaction_type",
        "compatible_reaction_types",
        "reaction_pattern_confidence",
        "reaction_pattern_requires_condition_evidence",
    )
    assert CONCISE_REACTION_REVIEW_FIELDS[label_index : label_index + 1] == (
        "reaction_label",
    )
    assert review_rows[0]["canonical_reaction_smiles"]
    assert review_rows[0]["reaction_label"]
    assert review_rows[0]["original_reaction_type"] == "Original Suzuki Label"
    assert review_rows[0]["primary_reaction_pattern"] == (
        "organoboron_c_c_coupling_like"
    )
    assert review_rows[0]["primary_reaction_pattern_count"] == "1"
    assert review_rows[0]["reaction_pattern_matches"].startswith(
        "organoboron_c_c_coupling_like"
    )
    assert review_rows[0]["identified_reaction_type"] == "suzuki_miyaura"
    assert review_rows[0]["compatible_reaction_types"] == "suzuki_miyaura"
    assert review_rows[0][
        "reaction_pattern_requires_condition_evidence"
    ] == "False"
    assert review_rows[0]["reaction_label_status"] == "structural_equation"
    assert review_rows[0]["transformation_class"] == (
        "generic_graph_transformation"
    )
    assert review_rows[0]["fallback_retrieval_eligible"] == "True"
    assert not {
        "signature_id",
        "reaction_core_id",
        "reaction_core_shape_key",
        "reaction_core_mapping_equivalence_key",
        "fallback_descriptor_id",
        "external_mapping_matched_hypotheses",
        "edit_hypothesis_ids",
        "installed_fragment_key",
        "partial_transformation_key",
    }.intersection(CONCISE_REACTION_REVIEW_FIELDS)
    assert review_rows[0]["evidence_quality"] == "global_atom_correspondence"
    assert review_rows[0]["reaction_completeness_status"] == "verified"
    assert review_rows[0]["reaction_core_status"] == "available_inferred"
    assert review_rows[0]["reaction_label"] == "Ar–Br + Ar–B(OH)₂ → Ar–Ar"
    assert review_rows[0]["reaction_core_unavailability_reasons"] == ""
    assert review_rows[0]["chemistry_status"] == "verified"
    assert review_rows[0]["condition_stage_status"] == "single_stage"
    assert review_rows[0]["index_eligibility"] == "eligible"
    assert review_rows[0]["reactivity_profile"] == ""


def test_review_exports_ambiguous_structural_pattern_candidates() -> None:
    record = convert_record(_raw("Brc1ccccc1.CN>>CNc1ccccc1"))

    row = concise_reaction_review_row(record.to_dict())

    assert row["primary_reaction_pattern"] == "sp2_c_n_substitution_like"
    assert row["identified_reaction_type"] == ""
    assert row["compatible_reaction_types"] == (
        "buchwald_hartwig_c_n; snar_amination; ullmann_c_n"
    )
    assert row["reaction_pattern_requires_condition_evidence"] == "True"


def test_review_exports_multi_site_pattern_count() -> None:
    record = convert_record(
        _raw(
            "Brc1ccc(Br)cc1.CC(N)=O"
            ">>CC(=O)Nc1ccc(NC(C)=O)cc1"
        )
    )

    row = concise_reaction_review_row(record.to_dict())

    assert row["primary_reaction_pattern"] == "sp2_c_n_substitution_like"
    assert row["primary_reaction_pattern_count"] == "2"
    assert row["reaction_label"] == (
        "Ar–Br + R–C(O)–NH₂ + R–C(O)–NH₂ "
        "→ Ar–NH–R–C(O) + Ar–NH–R–C(O)"
    )


def test_review_exports_co_occurring_tandem_patterns() -> None:
    record = convert_record(
        _raw(
            "CC(C)(C)OC(=O)NCCN."
            "O=c1oc2cc(Br)ccc2cc1-c1ccccc1"
            ">>NCCNc1ccc2cc(-c3ccccc3)c(=O)oc2c1"
        )
    )

    row = concise_reaction_review_row(record.to_dict())

    assert row["co_occurring_reaction_patterns"] == (
        "sp2_c_n_substitution_like; amine_deprotection_like"
    )
    assert row["reaction_label"] == "HetAr–Br + RO–C(O)–NHR → HetAr–NH–Alk"


def test_review_exports_carbonyl_alpha_c_h_coupling_not_heck() -> None:
    record = convert_record(
        _raw(
            "CC(=O)c1ccc(Br)cc1.O=C1CCCC1"
            ">>CC(=O)c1ccc(C2CCCC2=O)cc1"
        )
    )

    row = concise_reaction_review_row(record.to_dict())

    assert row["primary_reaction_pattern"] == (
        "carbonyl_alpha_c_h_sp2_c_c_coupling_like"
    )
    assert row["identified_reaction_type"] == (
        "carbonyl_alpha_c_h_coupling"
    )
    assert "heck_coupling_like" not in row["reaction_pattern_matches"]


def test_concise_review_formats_spectators_and_partner_environment() -> None:
    row = concise_reaction_review_row(
        {
            "chemistry_status": ChemistryStatus.REVIEW,
            "condition_status": ConditionStatus.RESOLVED_COMPLETE,
            "condition_stage_status": ConditionStageStatus.UNASSIGNED_MULTISTAGE,
            "index_eligibility": IndexEligibility.ELIGIBLE,
            "reaction_signature": {
                "stereo_changes": [
                    {
                        "stereo_type": "atom",
                        "old_descriptor": "S",
                        "new_descriptor": "R",
                        "change_type": "inverted",
                    }
                ],
                "spectator_groups": [
                    {
                        "group_id": "ether",
                        "chemist_label": "R–O–R",
                        "graph_distance": 1,
                    },
                    {
                        "group_id": "ether",
                        "chemist_label": "R–O–R",
                        "graph_distance": 3,
                    },
                ],
                "partners": [
                    {
                        "component_index": 0,
                        "role": "nitrogen_partner",
                        "reactivity_profile": {
                            "context_kind": "heteroatom",
                            "context": {
                                "context_kind": "heteroatom",
                                "element": "N",
                                "substitution_class": "primary",
                                "resonance_class": "aryl_delocalized",
                            },
                            "steric": {
                                "accessibility_class": "moderate",
                                "accessibility_score": 0.4,
                            },
                            "electronic": {
                                "activation_axis": "lone_pair_availability",
                                "activation_class": "low",
                                "activation_score": 0.55,
                            },
                            "reactive_center": {
                                "element": "N",
                                "substitution_class": "primary",
                                "hydrogen_count": 1,
                                "lone_pair_class": "aryl_delocalized",
                                "lone_pair_availability": "low",
                            },
                            "modifiers": [],
                        },
                    }
                ],
            },
        }
    )

    assert row["stereochemical_changes"] == "atom: S→R (inverted)"
    assert row["spectators"] == "2× R–O–R [ether] (d=1/3)"
    assert row["reactivity_profile"] == (
        "nitrogen partner: primary N, aryl delocalized; access moderate; "
        "lone-pair availability low; primary N-H; "
        "lone pair aryl delocalized/low"
    )
    assert row["chemistry_status"] == "review"
    assert row["condition_status"] == "resolved_complete"
    assert row["condition_stage_status"] == "unassigned_multistage"
    assert row["index_eligibility"] == "eligible"


def test_concise_review_uses_the_single_reaction_label() -> None:
    row = concise_reaction_review_row(
        {
            "reaction_label": {
                "text": "R–X + R′ → R–R′",
                "status": "structural_equation",
                "basis": "reaction_sites",
            },
            "reaction_core": {
                "evidence_status": "verified",
            },
            "reaction_signature": {
                "spectator_groups": [],
                "partners": [],
            },
        }
    )

    assert row["reaction_label"] == "R–X + R′ → R–R′"
    assert row["reaction_core_status"] == "available_verified"
    assert row["reaction_core_unavailability_reasons"] == ""


def test_review_csv_rows_leave_unavailable_reaction_label_blank() -> None:
    record = convert_record(_raw("CC.O>>N#N"))

    concise_row = concise_reaction_review_row(record.to_dict())
    tier_row = flatten_generic_record(record)

    assert record.reaction_label is not None
    assert record.reaction_label["status"] == "unavailable"
    assert record.reaction_label["text"] == "Unavailable"
    assert concise_row["reaction_label"] == ""
    assert concise_row["reaction_label_status"] == "unavailable"
    assert tier_row["reaction_label"] == ""
    assert tier_row["reaction_label_status"] == "unavailable"


def test_recursive_dataset_folder_converts_to_one_concise_review_csv(
    tmp_path: Path,
) -> None:
    source = tmp_path / "datasets"
    nested = source / "nested"
    nested.mkdir(parents=True)
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    for path, reaction_id in (
        (source / "root.csv", "root"),
        (nested / "child.csv", "child"),
    ):
        row = _csv_row(reaction_id, reaction, reaction_type="Suzuki source")
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
    progress = []
    output = source / "concise_review.csv"

    report = convert_dataset_folder_to_concise_review_csv(
        source,
        output,
        progress_callback=progress.append,
        progress_interval=1,
    )

    assert [
        path.relative_to(source).as_posix()
        for path in discover_csv_datasets(source)
        if path != output
    ] == ["nested/child.csv", "root.csv"]
    assert report["source_file_count"] == 2
    assert report["row_count"] == 2
    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        assert len(list(csv.DictReader(handle))) == 2
    assert progress[0].phase == "discovered"
    assert progress[-1].phase == "completed"


def test_conversion_input_discovery_accepts_files_and_folders_without_duplicates(
    tmp_path: Path,
) -> None:
    folder = tmp_path / "folder"
    folder.mkdir()
    first = folder / "first.csv"
    second = folder / "second.observations.jsonl.gz"
    ignored = folder / "notes.txt"
    first.write_text("reaction_smiles\n", encoding="utf-8")
    second.write_bytes(b"")
    ignored.write_text("not a dataset", encoding="utf-8")

    discovered = discover_conversion_datasets((first, folder, ignored))

    assert discovered == (first, second)


def test_recursive_concise_review_cancellation_removes_temporary_file(
    tmp_path: Path,
) -> None:
    dataset = tmp_path / "dataset.csv"
    row = _csv_row("cancel", "CC>>CC", reaction_type="Unknown")
    with dataset.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    output = tmp_path / "review.csv"

    with pytest.raises(ConciseReviewConversionCancelled):
        convert_dataset_folder_to_concise_review_csv(
            tmp_path,
            output,
            cancel_check=lambda: True,
        )

    assert not output.exists()
    assert not output.with_suffix(".csv.tmp").exists()


def test_sharded_conversion_is_restartable_and_integrity_checked(
    tmp_path,
    monkeypatch,
) -> None:
    dataset = tmp_path / "mixed.csv"
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    rows = [
        _csv_row(
            f"reaction-{index}",
            reaction,
            reaction_type="untrusted",
        )
        for index in range(4)
    ]
    with dataset.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    output = tmp_path / "sharded"
    progress = []

    first = convert_datasets_sharded(
        dataset,
        output,
        shard_size=2,
        mode="test",
        workers=2,
        progress_callback=progress.append,
    )
    first_catalog = (output / "recipe_catalog.jsonl.gz").read_bytes()
    second = convert_datasets_sharded(dataset, output, shard_size=2, mode="test")

    assert first["shard_count"] == 2
    assert any(
        update.phase == "source_processing"
        and update.message == "Processing 1/1 (total): mixed.csv"
        for update in progress
    )
    assert first["output_row_count"] == 4
    assert first["index_eligibility_counts"] == {"eligible": 4}
    assert first["transformation_class_counts"] == {
        "generic_graph_transformation": 4
    }
    assert first["named_family_counts"] == {"suzuki_miyaura": 4}
    assert first["integrity"]["valid"]
    assert second["reused_shard_count"] == 2
    assert (output / "recipe_catalog.jsonl.gz").read_bytes() == first_catalog
    assert len(load_generic_index(output / "records.jsonl.gz").rows) == 4
    assert len(load_generic_index(output / "shard_manifest.json").rows) == 4

    monkeypatch.setattr(
        sharded_module,
        "reaction_label_definition_versions",
        lambda: {"chemist_notation.v1.json": "changed"},
    )
    invalidated = convert_datasets_sharded(
        dataset,
        output,
        shard_size=2,
        mode="test",
    )
    assert invalidated["reused_shard_count"] == 0

    manifest = json.loads((output / "shard_manifest.json").read_text())
    first_shard = output / manifest["shards"][0]["output_path"]
    first_shard.write_bytes(first_shard.read_bytes() + b"tamper")
    integrity = validate_sharded_conversion(
        output / "shard_manifest.json",
        verify_rows=False,
    )
    assert not integrity["valid"]
    assert any(
        issue.startswith("output_checksum_mismatch") for issue in integrity["issues"]
    )


def test_cancelled_sharded_conversion_checkpoints_and_resumes(
    tmp_path: Path,
) -> None:
    dataset = tmp_path / "mixed.csv"
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    rows = [
        _csv_row(
            f"reaction-{index}",
            reaction,
            reaction_type="untrusted",
        )
        for index in range(3)
    ]
    with dataset.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    output = tmp_path / "sharded"
    cancel = {"requested": False}

    def on_progress(progress) -> None:
        if progress.phase == "shard_completed":
            cancel["requested"] = True

    with pytest.raises(ShardedConversionCancelled):
        convert_datasets_sharded(
            dataset,
            output,
            shard_size=1,
            checkpoint_interval=1,
            progress_callback=on_progress,
            cancel_check=lambda: cancel["requested"],
        )

    partial_manifest = json.loads(
        (output / "shard_manifest.json").read_text(encoding="utf-8")
    )
    assert len(partial_manifest["shards"]) == 1
    resumed = convert_datasets_sharded(dataset, output, shard_size=1)
    assert resumed["reused_shard_count"] == 1
    assert resumed["output_row_count"] == 3

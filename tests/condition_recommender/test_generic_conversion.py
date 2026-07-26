import csv
import json
from pathlib import Path

import condition_recommender.conversion.generic as generic_module
import pytest
from condition_recommender.conversion.concise_review import (
    CONCISE_REACTION_REVIEW_FIELDS,
    ConciseReviewConversionCancelled,
    convert_dataset_folder_to_concise_review_csv,
    export_concise_reaction_review_csv,
)
from condition_recommender.conversion.engine import convert_datasets
from condition_recommender.conversion.generic import (
    GenericConversionCache,
    convert_record,
)
from condition_recommender.conversion.input_schema import (
    adapt_row,
    discover_csv_datasets,
)
from condition_recommender.conversion.sharded import (
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
    assert record.schema_version == "1.9"
    assert record.converter_definition_version == "generic_conversion.v1.6"
    assert record.reaction_signature["schema_version"] == "1.3"
    assert record.reaction_signature["topology"]["reaction_scope"] == ("intermolecular")
    assert record.reference_id.startswith("REF1:")
    assert record.reference_identity["resolution_status"] == "bibliographic_text"
    assert record.chemistry_status == ChemistryStatus.VERIFIED
    assert record.condition_status == ConditionStatus.RESOLVED_COMPLETE
    assert record.outcome_status == OutcomeStatus.USABLE
    assert record.index_eligibility == IndexEligibility.ELIGIBLE
    assert record.resolved_recipe_core_id.startswith("RCORE1:")
    assert record.reference_condition_series_id.startswith("RCS1:")


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
    assert record.reaction_label == "H2C=CH2 → H3C–CH3"
    assert record.reaction_label_status == "mapped_generic_pattern"
    assert record.reaction_display_label is not None
    assert record.reaction_display_label["status"] == "generic_pattern"
    assert record.reaction_display_label["pattern_id"] == "hydrogenation"
    assert record.reaction_display_label["transformation_label"] == (
        "C=C hydrogenation"
    )
    assert record.reaction_display_label["product_context_label"] == "H3C–CH3"
    assert record.reaction_display_label["structural_label"] == (
        "C=C → C–C; 2 × H gain at C"
    )
    assert len(record.reaction_display_label["clauses"]) == 3
    assert record.reaction_signature["order_changes"] == ("C-C:DOUBLE>SINGLE",)
    assert record.resolved_recipe_id.startswith("RCR1:")


def test_exact_multi_event_signature_is_verified() -> None:
    record = convert_record(
        _raw("CO.CS.Fc1ccc(F)cc1>>COc1ccc(SC)cc1")
    )

    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.evidence_quality == "exact_multi_event_reconstruction"
    assert record.named_family is None
    assert record.reaction_signature is not None
    assert record.reaction_signature["event_count"] == 2
    assert record.reaction_signature["event_scope"] == "multi_event"
    assert record.reaction_label is not None
    assert record.reaction_label.count("substitution") == 2
    assert " + " in record.reaction_label


def test_grammar_only_record_is_review_not_rejected() -> None:
    record = convert_record(_raw("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1"))

    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.admission_reasons == ("missing_verified_reaction_signature",)


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
    assert unresolved.admission_reasons == ("no_usable_transformation_evidence",)
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


def test_unstructured_multistage_conditions_are_review_only() -> None:
    record = convert_record(
        _raw(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            stages="2",
        )
    )

    assert record.condition_status == ConditionStatus.MULTISTAGE_AMBIGUOUS
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
        "review": 1,
        "rejected": 1,
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
    assert all(row["reaction_signature_id"].startswith("RS1:") for row in verified)
    assert {row["reaction_event_count"] for row in verified} == {"1"}
    assert {row["reaction_event_scope"] for row in verified} == {"single_event"}
    assert verified[0]["reaction_scope"] == "intermolecular"
    assert {row["reaction_scope"] for row in verified} == {
        "intermolecular",
        "unimolecular",
    }
    assert json.loads((output / "conversion_report.json").read_text()) == report
    assert report["schema_version"] == "1.2"
    assert report["reaction_signature_schema_version"] == "1.3"
    assert report["reaction_scope_counts"] == {
        "intermolecular": 1,
        "unimolecular": 1,
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
    assert report["row_count"] == 1
    assert tuple(review_rows[0]) == CONCISE_REACTION_REVIEW_FIELDS
    assert review_rows[0]["canonical_reaction_smiles"]
    assert review_rows[0]["reaction_display_label_detailed"]
    assert review_rows[0]["original_reaction_type"] == "Original Suzuki Label"
    assert review_rows[0]["detected_reaction_family"] == "suzuki_miyaura"
    assert review_rows[0]["detection_status"] == "family_overlay"


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

    first = convert_datasets_sharded(
        dataset,
        output,
        shard_size=2,
        mode="test",
        workers=2,
    )
    first_catalog = (output / "recipe_catalog.jsonl.gz").read_bytes()
    second = convert_datasets_sharded(dataset, output, shard_size=2, mode="test")

    assert first["shard_count"] == 2
    assert first["output_row_count"] == 4
    assert first["index_eligibility_counts"] == {"eligible": 4}
    assert first["transformation_class_counts"] == {
        "c_c_transfer_coupling": 4
    }
    assert first["named_family_counts"] == {"suzuki_miyaura": 4}
    assert first["integrity"]["valid"]
    assert second["reused_shard_count"] == 2
    assert (output / "recipe_catalog.jsonl.gz").read_bytes() == first_catalog
    assert len(load_generic_index(output / "records.jsonl.gz").rows) == 4
    assert len(load_generic_index(output / "shard_manifest.json").rows) == 4

    manifest = json.loads((output / "shard_manifest.json").read_text())
    first_shard = output / manifest["shards"][0]["output_path"]
    first_shard.write_bytes(first_shard.read_bytes() + b"tamper")
    integrity = validate_sharded_conversion(
        output / "shard_manifest.json",
        verify_rows=False,
    )
    assert not integrity["valid"]
    assert any(
        issue.startswith("output_checksum_mismatch")
        for issue in integrity["issues"]
    )

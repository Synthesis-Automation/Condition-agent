import csv
import json

from condition_recommender.conversion.engine import convert_datasets
from condition_recommender.conversion.generic import convert_record
from condition_recommender.conversion.input_schema import adapt_row
from condition_recommender.models import AdmissionTier


def _raw(
    reaction: str,
    *,
    reaction_id: str = "record-1",
    yield_pct: str = "80",
    catalyst_cas: str = "14221-01-3",
    reagent_cas: str = "584-08-7",
    solvent_cas: str = "108-88-3",
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
    record = convert_record(
        _raw("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
    )

    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.source_declared_family == "untrusted-source-label"
    assert record.named_family == "suzuki_miyaura"
    assert record.reaction_signature is not None
    assert record.canonical_reaction_id.startswith("CRX1:")
    assert record.observation_id.startswith("OBS1:")
    assert record.raw_recipe_id.startswith("RAWCOND1:")
    assert record.resolved_recipe_id.startswith("RCR1:")
    assert record.resolved_recipe["catalysts"][0]["primary_role"] == (
        "metal_catalyst"
    )
    assert record.resolved_recipe["bases"][0]["primary_role"] == "base"
    assert record.condition_resolution["component_count"] == 3


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
    assert record.reaction_label == "C=C → C–C; 2 × H gain at C"
    assert record.reaction_label_status == "mapped_edit_summary"
    assert record.reaction_display_label is not None
    assert record.reaction_display_label["status"] == "observed_edits"
    assert len(record.reaction_display_label["clauses"]) == 3
    assert record.reaction_signature["order_changes"] == (
        "C-C:DOUBLE>SINGLE",
    )
    assert record.resolved_recipe_id.startswith("RCR1:")


def test_grammar_only_record_is_review_not_rejected() -> None:
    record = convert_record(
        _raw("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1")
    )

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
        for line in (output / "records.jsonl").read_text().splitlines()
    ]
    assert len(records) == 4
    assert records[1]["named_family"] is None
    assert records[1]["reaction_signature"]["order_changes"] == [
        "C-C:DOUBLE>SINGLE"
    ]
    assert records[1]["resolved_recipe_id"].startswith("RCR1:")
    with (output / "verified.csv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        verified = list(csv.DictReader(handle))
    assert len(verified) == 2
    assert all(row["reaction_signature_id"].startswith("RS1:") for row in verified)
    assert json.loads((output / "conversion_report.json").read_text()) == report
    assert (output / "conversion_report.md").exists()

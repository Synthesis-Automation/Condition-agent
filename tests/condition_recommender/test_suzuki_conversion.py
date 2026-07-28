import csv
import json

from condition_recommender.condition_normalization import (
    is_valid_cas,
    normalize_cas_list,
)
from condition_recommender.conversion.suzuki import (
    convert_file,
    convert_row,
    flatten_record,
)
from condition_recommender.models import AdmissionTier


def _row(reaction_smiles: str, **overrides: str) -> dict[str, str]:
    row = {
        "reaction_id": "test-1",
        "reaction_type": "Suzuki_Miyaura",
        "reaction_smiles": reaction_smiles,
        "yield_pct": "82",
        "temperature_c": "80",
        "time_h": "12",
        "catalyst_cas": "14221-01-3",
        "reagent_cas": "584-08-7",
        "solvent_cas": "108-88-3",
        "reference": "test",
    }
    row.update(overrides)
    return row


def test_condition_identifiers_are_deduplicated_and_sorted() -> None:
    assert normalize_cas_list("584-08-7, 108-88-3;584-08-7") == ("108-88-3", "584-08-7")
    assert is_valid_cas("7732-18-5")
    assert not is_valid_cas("7732-18-4")


def test_exact_suzuki_with_complete_conditions_is_verified() -> None:
    record = convert_row(_row("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"), 2)
    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.family_environment is not None
    assert record.product_connection is not None
    assert record.product_connection["concise_label"] == "Ar1–Ar2"
    assert record.conditions.complete
    assert record.reaction_signature is not None
    assert record.reaction_signature["signature_id"].startswith("RS2:")
    assert record.transformation_class == "c_c_transfer_coupling"
    assert record.transformation_confidence == 1.0
    flat = flatten_record(record)
    assert flat["signature_l0_exact"].startswith("L0:")
    assert (
        json.loads(flat["reaction_signature_json"])["signature_id"]
        == (record.reaction_signature["signature_id"])
    )


def test_unverified_product_is_sent_to_review() -> None:
    record = convert_row(_row("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1"), 2)
    assert record.admission_tier == AdmissionTier.REVIEW
    assert "reaction_not_exactly_verified" in record.admission_reasons
    assert record.reaction_label == "Ar1–Br + Ar2–B(OH)2 →"
    assert record.reaction_label_status == "source_family_reactant_only"


def test_source_family_removes_non_suzuki_partial_label_alternatives() -> None:
    record = convert_row(
        _row("Nc1ccccc1N.OCc1ccc(Br)cc1.OB(O)c1cccc2ccccc12>>c1ccccc1"), 2
    )
    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.reaction_label == "Ar1–Br + Ar2–B(OH)2 →"
    assert record.reaction_label_status == "source_family_reactant_only"


def test_missing_yield_is_rejected() -> None:
    record = convert_row(
        _row(
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            yield_pct="",
        ),
        2,
    )
    assert record.admission_tier == AdmissionTier.REJECTED
    assert record.admission_reasons == ("missing_or_invalid_yield",)


def test_convert_file_writes_three_tiers_and_report(tmp_path) -> None:
    input_path = tmp_path / "input.csv"
    rows = [
        _row("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"),
        _row("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1", reaction_id="test-2"),
        _row("not-a-reaction", reaction_id="test-3"),
    ]
    with input_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = convert_file(input_path, tmp_path / "out")
    assert report["tier_counts"] == {"verified": 1, "review": 1, "rejected": 1}
    assert (
        json.loads((tmp_path / "out" / "conversion_report.json").read_text()) == report
    )
    assert (tmp_path / "out" / "conversion_report.md").exists()
    for tier in ("verified", "review", "rejected"):
        assert (tmp_path / "out" / f"{tier}.csv").exists()
    with (tmp_path / "out" / "verified.csv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        verified = next(csv.DictReader(handle))
    assert verified["reaction_signature_id"].startswith("RS2:")
    assert verified["signature_l3_bond_edit"].startswith("L3:")

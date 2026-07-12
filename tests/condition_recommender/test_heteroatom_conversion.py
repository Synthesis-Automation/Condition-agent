from condition_recommender.conversion.heteroatom import convert_row
from condition_recommender.models import AdmissionTier


def _row(reaction: str) -> dict[str, str]:
    return {
        "reaction_id": "test", "reaction_smiles": reaction, "yield_pct": "80",
        "catalyst_cas": "14221-01-3", "reagent_cas": "584-08-7", "solvent_cas": "108-88-3",
    }


def test_exact_co_record_is_verified() -> None:
    record = convert_row(_row("Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1"), 2, element="O")
    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.product_connection["connection_type"] == "C_O"
    assert record.family_environment["family_id"] == "c_o_coupling"


def test_exact_cs_record_is_verified() -> None:
    record = convert_row(_row("Brc1ccccc1.CCS>>CCSc1ccccc1"), 2, element="S")
    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.product_connection["connection_type"] == "C_S"
    assert record.family_environment["family_id"] == "c_s_coupling"


def test_unverified_heteroatom_product_is_reviewed_with_blank_product() -> None:
    record = convert_row(_row("Brc1ccccc1.CCO>>c1ccccc1"), 2, element="O")
    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.reaction_label.endswith("→")
    assert record.product_connection is None

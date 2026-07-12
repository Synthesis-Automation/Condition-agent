from condition_recommender.conversion.cn import convert_row
from condition_recommender.models import AdmissionTier


def _row(reaction: str, **updates: str) -> dict[str, str]:
    row = {
        "reaction_id": "cn-1", "reaction_type": "C_N_Coupling",
        "reaction_smiles": reaction, "yield_pct": "80",
        "catalyst_cas": "14221-01-3", "reagent_cas": "584-08-7",
        "solvent_cas": "108-88-3",
    }
    row.update(updates); return row


def test_exact_cn_record_is_verified() -> None:
    record = convert_row(_row("Brc1ccccc1.CNC>>CN(C)c1ccccc1"), 2)
    assert record.admission_tier == AdmissionTier.VERIFIED
    assert record.product_connection["connection_type"] == "C_N"
    assert record.product_connection["concise_label"] == "Ar–NR1R2"
    assert record.family_environment["family_id"] == "c_n_coupling"


def test_partial_cn_record_keeps_product_empty() -> None:
    record = convert_row(_row("Brc1ccccc1.CNC>>c1ccccc1"), 2)
    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.reaction_label.endswith("→")
    assert record.product_connection is None


def test_multiple_product_cn_record_requires_review() -> None:
    reaction = "Brc1ccccc1.CNC>>CN(C)c1ccccc1.C"
    record = convert_row(_row(reaction), 2)
    assert record.admission_tier == AdmissionTier.REVIEW
    assert "multiple_products" in record.admission_reasons

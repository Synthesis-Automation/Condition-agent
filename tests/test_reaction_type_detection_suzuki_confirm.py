import pytest

from chemtools.featurizers import reaction_detection as rd
from chemtools.util.rdkit_helpers import rdkit_available


def test_confirm_suzuki_product_by_attachment_positive() -> None:
    if not rdkit_available():
        pytest.skip("rdkit not available")

    reactants = ["Brc1ccccc1", "c1ccc(B(O)O)cc1"]
    products = ["c1ccc(-c2ccccc2)cc1"]
    ok, reason = rd.confirm_suzuki_product_by_attachment(reactants, products)
    assert ok
    assert reason == "substructure_match"


def test_confirm_suzuki_product_by_attachment_negative() -> None:
    if not rdkit_available():
        pytest.skip("rdkit not available")

    reactants = ["Brc1ccccc1", "c1ccc(B(O)O)cc1"]
    products = ["c1ccc(Br)cc1"]
    ok, reason = rd.confirm_suzuki_product_by_attachment(reactants, products)
    assert not ok
    assert reason in {"no_substructure_match", "ambiguous_sites", "candidate_build_failed"}


def test_confirm_cn_coupling_by_attachment_positive() -> None:
    if not rdkit_available():
        pytest.skip("rdkit not available")

    reactants = ["Brc1ccccc1", "Nc1ccccc1"]
    products = ["c1ccc(Nc2ccccc2)cc1"]
    ok, reason = rd.confirm_coupling_product_by_attachment(
        reactants,
        products,
        "C_N_Coupling",
    )
    assert ok
    assert reason == "substructure_match"

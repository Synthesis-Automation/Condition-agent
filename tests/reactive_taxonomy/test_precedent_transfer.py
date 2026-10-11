"""Structural transfer review must not confuse product identity with a precedent."""

from reactive_taxonomy.precedent_transfer import review_precedent_transfer

BORON = "OB(O)c1ccccc1.Brc1ccccc1"
TIN = "C[Sn](C)(C)c1ccccc1.Brc1ccccc1"
PRODUCT = "c1ccc(-c2ccccc2)cc1"


def test_boron_tin_change_with_identical_product_requires_transfer_review() -> None:
    review = review_precedent_transfer(BORON, PRODUCT, TIN, PRODUCT)
    assert review.structure_scope == "same_product_different_precursors"
    assert "TRANSFER_GROUP_INVENTORY_DIFFERS_CHECK_DONOR_CHEMISTRY" in review.warnings
    assert any("SnR3" in item.signature for item in review.precursor_site_differences)
    assert review.proposed_only_components and review.recorded_only_components
    assert review.experimental_transfer == "not_established"


def test_exact_structures_ignore_order_and_mapping_without_claiming_experimental_transfer() -> None:
    review = review_precedent_transfer("CCBr.N", "CCN", "N.Br[CH2:1][CH3:2]", "NCC")
    assert review.structure_scope == "exact_structures"
    assert not review.precursor_site_differences
    assert review.experimental_transfer == "not_established"


def test_stereo_conflict_is_not_exact_even_when_site_inventories_agree() -> None:
    review = review_precedent_transfer("CCBr.N", "N[C@H](C)CO", "N.CCBr", "N[C@@H](C)CO")
    assert review.structure_scope == "different_product_requires_review"
    assert review.same_product is False
    assert not review.precursor_site_differences


def test_component_multiplicity_and_missing_structures_are_not_silently_ignored() -> None:
    review = review_precedent_transfer("N.N.CCBr", "CCN", "N.CCBr", "CCN")
    assert review.proposed_only_components == ("N",)
    gap = review_precedent_transfer("CCBr.N", "CCN", "", "CCN")
    assert gap.structure_scope == "unresolved"
    assert gap.same_product is gap.same_precursors is None


def test_no_detected_handle_does_not_prove_transferability() -> None:
    review = review_precedent_transfer("CC", "CCC", "CCCC", "CCCCC")
    assert review.experimental_transfer == "not_established"
    assert review.structure_scope == "different_product_requires_review"

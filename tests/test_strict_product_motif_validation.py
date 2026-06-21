from __future__ import annotations

import pytest

from chemtools.reaction import featurize_reaction
from chemtools.core.rdkit import rdkit_available


RXN_MISSING_ALKYNE_PRODUCT = (
    "C#Cc1ccccc1.COc1ccc(NC(=N)c2ccccc2Cl)cc1OC"
    ">>COc1ccc(-n2c(-c3ccccc3)cc3ccccc3c2=O)cc1OC"
)


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_crk_validation_requires_observed_product_motif_by_default() -> None:
    result = featurize_reaction(
        RXN_MISSING_ALKYNE_PRODUCT,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )

    assert result.get("reaction_type") != "Sonogashira"
    detection = result.get("detection") or {}
    validation = detection.get("validation") or {}
    assert validation.get("validated_detection") != "Sonogashira"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_crk_validation_can_opt_in_to_inferred_product_motifs() -> None:
    result = featurize_reaction(
        RXN_MISSING_ALKYNE_PRODUCT,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
            "strict_product_motif_validation": False,
        },
    )

    assert result.get("reaction_type") == "Sonogashira"
    detection = result.get("detection") or {}
    validation = detection.get("validation") or {}
    assert validation.get("validated_detection") == "Sonogashira"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_include_reaction_type_does_not_keep_initial_cn_when_crk_has_no_taxonomy_match() -> None:
    result = featurize_reaction(
        RXN_MISSING_ALKYNE_PRODUCT,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )

    assert result.get("reaction_type") != "C_N_Coupling"
    detection = result.get("detection") or {}
    validation = detection.get("validation") or {}
    assert validation.get("validated_detection") == "Unknown"

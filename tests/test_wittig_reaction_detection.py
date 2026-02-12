import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_wittig_reaction_detected_for_phosphonium_olefination_example() -> None:
    rxn = (
        "O=Cc1ccccc1."
        "[CH2-][P+](c1ccccc1)(c1ccccc1)c1ccccc1"
        ">>C=Cc1ccccc1"
    )
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Wittig_reaction"
    assert "Ar-PR3+" in reaction_key
    assert "Ar-Alkenyl_terminal" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_wittig_reaction_detected_for_standard_ylide_case_without_explicit_pr3_motif() -> None:
    rxn = (
        "COC(=O)C=P(c1ccccc1)(c1ccccc1)c1ccccc1."
        "COc1cc(/C(C)=C/C=O)ccc1C"
        ">>COC(=O)/C=C/C=C(\\C)c1ccc(C)c(OC)c1"
    )
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Wittig_reaction"
    assert "bond_broken: C-P" in reaction_key

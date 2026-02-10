from __future__ import annotations

import pytest

from chemtools.featurizers.formatters.detection_validation import (
    validate_detection_with_crk_key,
)
from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


def test_validate_detection_with_crk_key_prefers_quinazolinone_annulation() -> None:
    crk_raw = (
        "|Ar-NH2|Ar-NHCOR|CH3-OH -> Pyrimidine "
        "| bond_formed: C-N; C-N | bond_broken: C-O"
    )
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "Quinazolinone_annulation"
    assert result["validation_method"] == "crk_pattern"


def test_validate_detection_with_crk_key_keeps_sandmeyer_for_simple_case() -> None:
    crk_raw = "|Ar-NH2 -> Ar-Cl"
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "Sandmeyer"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_featurize_reaction_quinazolinone_annulation_key_is_centered() -> None:
    rxn = "CO.Nc1cc(Cl)ccc1C(=O)Nc1ccccc1>>O=c1c2ccc(Cl)cc2ncn1-c1ccccc1"
    result = featurize_reaction(
        rxn,
        options={"detailed": True, "confirm_coupling_products": True},
    )
    assert result.get("reaction_type") == "Quinazolinone_annulation"
    reaction_key = result.get("reaction_key") or ""
    assert "-> Pyrimidine" in reaction_key
    assert "-> HeteroAr-Cl" not in reaction_key
    detection = result.get("detection") or {}
    evidence = detection.get("evidence") or {}
    selected = evidence.get("selected") or {}
    assert selected.get("reaction_id") == "Quinazolinone_annulation"

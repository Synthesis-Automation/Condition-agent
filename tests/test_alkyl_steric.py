import pytest

from chemtools.featurizers.structural import featurize_molecule
from chemtools.util.rdkit_helpers import rdkit_available


def _analysis_for_prefix(result, prefix: str):
    for analysis in result["analyses"]:
        if analysis["compound_id"].startswith(prefix):
            return analysis
    return None


def _scaffold_score(smiles: str, prefix: str):
    result = featurize_molecule(smiles)
    analysis = _analysis_for_prefix(result, prefix)
    assert analysis is not None
    return analysis["steric"]["score_0_10"], analysis


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_alkyl_scaffold_trends():
    methyl, _ = _scaffold_score("CCl", "R-")
    ethyl, _ = _scaffold_score("CCBr", "R-")
    isopropyl, _ = _scaffold_score("CC(C)Br", "R-")
    tert_butyl, _ = _scaffold_score("CC(C)(C)Cl", "R-")

    assert tert_butyl > isopropyl > ethyl > methyl

    neopentyl, _ = _scaffold_score("CC(C)(C)CBr", "R-")
    assert neopentyl > ethyl


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_benzyl_ring_penalty():
    result = featurize_molecule("c1ccccc1CBr", options={"include_steric_details": True})
    analysis = _analysis_for_prefix(result, "Bn-") or _analysis_for_prefix(result, "R-")
    assert analysis is not None
    substituents = analysis["steric"]["details"]["substituents"]
    assert any(sub["has_ring"] for sub in substituents)


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_allyl_flag():
    result = featurize_molecule("C=CCBr", options={"include_steric_details": True})
    analysis = _analysis_for_prefix(result, "Allyl-") or _analysis_for_prefix(result, "R-")
    assert analysis is not None
    assert analysis["steric"]["details"]["alpha"]["allylic"] is True

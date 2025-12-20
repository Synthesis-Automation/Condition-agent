import pytest

from chemtools.taxonomy.v2.analyze import analyze_smiles
from chemtools.util.rdkit_helpers import rdkit_available


def _analysis_for_prefix(result, prefix: str):
    for analysis in result["analyses"]:
        if analysis["compound_id"].startswith(prefix):
            return analysis
    return None


def _scaffold_score(smiles: str, prefix: str):
    result = analyze_smiles(smiles)
    analysis = _analysis_for_prefix(result, prefix)
    assert analysis is not None
    return analysis["steric"]["scaffold_score_0_10"], analysis


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
    result = analyze_smiles("c1ccccc1CBr")
    analysis = _analysis_for_prefix(result, "Bn-") or _analysis_for_prefix(result, "R-")
    assert analysis is not None
    assert any(sub["has_ring"] for sub in analysis["steric"]["substituents"])


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_allyl_flag():
    result = analyze_smiles("C=CCBr")
    analysis = _analysis_for_prefix(result, "Allyl-") or _analysis_for_prefix(result, "R-")
    assert analysis is not None
    assert analysis["steric"]["alpha"]["allylic"] is True

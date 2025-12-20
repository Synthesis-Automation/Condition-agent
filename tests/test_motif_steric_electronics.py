from pathlib import Path

import pytest

from chemtools.taxonomy.v2.analyze import analyze_smiles
from chemtools.util.rdkit_helpers import rdkit_available

ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATHS = {
    "groups": ROOT / "chemtools" / "taxonomy" / "v2" / "organic_groups.v1.2.json",
    "compounds": ROOT / "chemtools" / "taxonomy" / "v2" / "organic_compounds.v1.2.json",
    "templates": ROOT / "chemtools" / "taxonomy" / "v2" / "smarts_templates.v1.json",
}


def _analysis_for(result, compound_id: str):
    for analysis in result["analyses"]:
        if analysis["compound_id"] == compound_id:
            return analysis
    return None


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_bromobenzene_steric_electronics():
    result = analyze_smiles("Brc1ccccc1", registry_paths=REGISTRY_PATHS)
    assert "Ar-Br" in {hit["compound_id"] for hit in result["motifs"]}
    analysis = _analysis_for(result, "Ar-Br")
    assert analysis is not None
    assert analysis["steric"]["score_0_10"] <= 1.0
    score = analysis["electronic"]["score_0_10"]
    assert 4.5 <= score <= 5.5


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_benzaldehyde_electronics_including_ipso():
    result = analyze_smiles("O=CC1=CC=CC=C1", registry_paths=REGISTRY_PATHS)
    assert "Ar-CHO" in {hit["compound_id"] for hit in result["motifs"]}
    analysis = _analysis_for(result, "Ar-CHO")
    assert analysis is not None
    assert analysis["steric"]["score_0_10"] <= 1.0
    assert analysis["electronic"]["score_0_10"] > 5.0


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_ortho_alkyl_bromobenzene_steric():
    result = analyze_smiles("Cc1c(Br)cccc1", registry_paths=REGISTRY_PATHS)
    analysis = _analysis_for(result, "Ar-Br")
    assert analysis is not None
    assert analysis["steric"]["score_0_10"] > 0.0


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_para_nitro_bromobenzene_electronics():
    result = analyze_smiles("O=[N+]([O-])c1ccc(Br)cc1", registry_paths=REGISTRY_PATHS)
    analysis = _analysis_for(result, "Ar-Br")
    assert analysis is not None
    assert analysis["electronic"]["score_0_10"] > 6.0

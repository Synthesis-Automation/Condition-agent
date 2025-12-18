import pytest

from chemtools.util.rdkit_helpers import rdkit_available


pytestmark = pytest.mark.skipif(not rdkit_available(), reason="RDKit required for electronics POC v1 tests")


def test_electronics_v1_bromobenzene_detects_ar_br_and_score_in_range():
    from chemtools.taxonomy.v1.ar_context_electronics_v1 import analyze_smiles_reactivity_electronics_v1

    res = analyze_smiles_reactivity_electronics_v1("Brc1ccccc1")
    assert res.get("error") is None
    assert res["aryl_anchor_site_count"] >= 1
    assert any(s.get("label") == "Ar-Br" for s in res.get("sites", []))
    assert 0 <= res["aryl_electron_poor_score_0_10_max"] <= 10


def test_electronics_v1_anisyl_bromide_more_electron_rich_than_bromobenzene():
    from chemtools.taxonomy.v1.ar_context_electronics_v1 import analyze_smiles_reactivity_electronics_v1

    bromo = analyze_smiles_reactivity_electronics_v1("Brc1ccccc1")
    anisyl = analyze_smiles_reactivity_electronics_v1("COc1ccc(Br)cc1")
    assert bromo.get("error") is None
    assert anisyl.get("error") is None
    assert anisyl["aryl_electron_poor_score_0_10_max"] < bromo["aryl_electron_poor_score_0_10_max"]


def test_electronics_v1_nitro_bromobenzene_more_electron_poor_than_bromobenzene():
    from chemtools.taxonomy.v1.ar_context_electronics_v1 import analyze_smiles_reactivity_electronics_v1

    bromo = analyze_smiles_reactivity_electronics_v1("Brc1ccccc1")
    nitro = analyze_smiles_reactivity_electronics_v1("O=[N+]([O-])c1ccc(Br)cc1")
    assert bromo.get("error") is None
    assert nitro.get("error") is None
    assert nitro["aryl_anchor_site_count"] >= 1
    assert nitro["aryl_electron_poor_score_0_10_max"] > bromo["aryl_electron_poor_score_0_10_max"]


def test_electronics_v1_benzaldehyde_detects_ar_cho():
    from chemtools.taxonomy.v1.ar_context_electronics_v1 import analyze_smiles_reactivity_electronics_v1

    res = analyze_smiles_reactivity_electronics_v1("O=Cc1ccccc1")
    assert res.get("error") is None
    assert res["aryl_anchor_site_count"] >= 1
    assert any(s.get("label") == "Ar-CHO" for s in res.get("sites", []))


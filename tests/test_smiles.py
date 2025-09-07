from chemtools.smiles import normalize
from chemtools.util.rdkit_helpers import rdkit_available

def test_normalize_basic():
    res = normalize("Brc1ccc(F)cc1.O")
    assert res["largest_smiles"]
    assert res["smiles_norm"]

def test_normalize_salt_hcl():
    res = normalize("Nc1ccccc1.Cl")
    assert res["smiles_norm"].startswith("Nc1ccccc1")

def test_normalize_carboxylate_to_acid():
    res = normalize("O=C([O-])c1ccccc1")
    assert "[O-]" not in (res["smiles_norm"] or "")

def test_normalize_quat_ammonium_kept():
    res = normalize("C[N+](C)(C)C")
    # Permanent cation should remain charged when RDKit available
    if rdkit_available():
        assert '+' in (res["smiles_norm"] or '')
    else:
        assert res["smiles_norm"]

def test_invalid_smiles():
    res = normalize("notasmiles")
    if rdkit_available():
        assert res.get("error") == "INVALID_SMILES"
    else:
        # Fallback mode: heuristics only; no strict validation possible
        assert res.get("smiles_norm") == "notasmiles"

def test_normalize_hydrate_and_fragments():
    res = normalize("Nc1ccccc1.O")
    assert res["smiles_norm"].startswith("Nc1ccccc1")

def test_normalize_tautomers_equal():
    # Keto-enol pair for acetylacetone-like pattern
    a = normalize("CC(=O)CH2C(=O)C")
    b = normalize("CC(=O)C=C(O)C")
    if rdkit_available():
        assert a["smiles_norm"] == b["smiles_norm"]

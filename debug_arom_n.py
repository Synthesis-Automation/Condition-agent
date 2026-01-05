from rdkit import Chem
from chemtools.featurizers.motif_detect import detect_motifs
from chemtools.featurizers.motif_registry import build_compound_registry, _default_registry_paths

# Load registry
registry_paths = _default_registry_paths()
registry = build_compound_registry(registry_paths)
compiled_compounds = registry["compiled_compounds"]

def test_mol(smiles, label):
    mol = Chem.MolFromSmiles(smiles)
    motifs = detect_motifs(mol, compiled_compounds)
    ids = [m['compound_id'] for m in motifs]
    print(f"{label}: {smiles}")
    print(f"Motifs: {ids}")
    return ids

# Case 1: Imidazole reactants
test_mol("Cc1ncc[nH]1", "Imidazole Reactant 1")
test_mol("O=C(Nc1cccc(Br)n1)c1ccccc1", "Imidazole Reactant 2")
test_mol("Cc1nccn1-c1cccc(NC(=O)c2ccccc2)n1", "Imidazole Product")

# Case 2: Carbazole product
test_mol("c1ccc2c(c1)c1ccccc1n2-c1ccc(-c2nncnn2)cc1", "Carbazole Product")

# Case 3: N-phenylpyrrole (should work)
test_mol("c1ccc(N2C=CC=C2)cc1", "N-phenylpyrrole")


from chemtools.featurizers.motif_registry import build_compound_detect_registry, _default_registry_paths
from chemtools.util.rdkit_helpers import parse_smiles

paths = _default_registry_paths()
registry = build_compound_detect_registry(paths)
compiled = registry.get("compiled_compounds", {})

smiles_list = [
    "OCC1CCCCC1", # Alcohol
    "O=CC1CCCCC1", # Carbonyl
    "O=[N+]([O-])C1CCCCC1", # Nitro
    "c1ccccc1" # Benzene
]

motifs_to_check = ["Any-OH", "Any-CHO", "Any-NO2", "Ar-H"]

for smiles in smiles_list:
    mol = parse_smiles(smiles)
    print(f"SMILES: {smiles}")
    for mid in motifs_to_check:
        if mid in compiled:
            match = any(mol.HasSubstructMatch(q) for q in compiled[mid])
            print(f"  {mid}: {match}")
        else:
            print(f"  {mid}: NOT IN REGISTRY")

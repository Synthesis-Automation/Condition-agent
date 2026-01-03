from chemtools.featurizers.molecule import featurize_molecule
from chemtools.featurizers.motif_registry import build_compound_registry, _default_registry_paths
from rdkit import Chem

smiles = "NCc1ccccc1"
mol = Chem.MolFromSmiles(smiles)
options = {"discovery_mode": True}

paths = _default_registry_paths()
registry = build_compound_registry(paths)
compiled_compounds = registry["compiled_compounds"]

print(f"Total compiled compounds: {len(compiled_compounds)}")

bn_nh2 = next((c for c in compiled_compounds if c.compound_id == "Bn-NH2"), None)
if bn_nh2:
    print(f"Bn-NH2 SMARTS: {bn_nh2.smarts}")
    matches = mol.GetSubstructMatches(bn_nh2.query)
    print(f"Bn-NH2 matches: {matches}")
else:
    print("Bn-NH2 not found in registry")

payload = featurize_molecule(smiles, options=options)

print("\nDetected Motifs:")
for m in payload.get("motifs", []):
    print(f"ID: {m['compound_id']}, Priority: {m['priority']}, Undocumented: {m.get('undocumented', False)}, Atoms: {m['atoms']}")

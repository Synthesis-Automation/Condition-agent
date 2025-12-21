
from chemtools.featurizers.motif_detect import detect_motifs
from chemtools.featurizers.reaction_detection import _load_compound_registry
from chemtools.util.rdkit_helpers import parse_smiles
import json

registry = _load_compound_registry()
smiles_list = ["Brc1ccccc1", "CCNCC"]

for s in smiles_list:
    mol = parse_smiles(s)
    print(f"\nMotifs for {s}:")
    hits = detect_motifs(mol, registry["compiled_compounds"])
    for hit in hits:
        print(f"  - {hit['compound_id']}")


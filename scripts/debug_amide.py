import sys
from pathlib import Path
from rdkit import Chem

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.featurizers.motif_detect import detect_motifs
from chemtools.featurizers.motif_registry import build_compound_registry
from chemtools import detect_reaction

def debug_amide():
    reaction = "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1"
    print(f"Testing reaction: {reaction}")
    
    # 1. Test full detection
    res = detect_reaction(reaction)
    print(f"\nDetection result: {res['family']} (Confidence: {res['confidence']})")
    if 'details' in res:
        print(f"Slot evidence: {res['details'].get('slot_evidence')}")
        print(f"Matched slots: {res['details'].get('matched_slots')}/{res['details'].get('required_slots')}")

    # 2. Inspect motifs
    base = PROJECT_ROOT / "chemtools" / "taxonomy" / "data"
    registry_paths = {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
    }
    registry = build_compound_registry(registry_paths)
    compiled = registry["compiled_compounds"]
    
    parts = reaction.split(">>")
    reactants = parts[0].split(".")
    products = parts[1].split(".") if len(parts) > 1 else []
    
    print("\n--- Reactant Motifs ---")
    for s in reactants:
        mol = Chem.MolFromSmiles(s)
        if mol:
            mol = Chem.AddHs(mol)
            hits = detect_motifs(mol, compiled)
            print(f"\n{s}:")
            for h in hits:
                print(f"  - {h['compound_id']}")
                
    print("\n--- Product Motifs ---")
    for s in products:
        mol = Chem.MolFromSmiles(s)
        if mol:
            mol = Chem.AddHs(mol)
            hits = detect_motifs(mol, compiled)
            print(f"\n{s}:")
            for h in hits:
                print(f"  - {h['compound_id']}")

if __name__ == "__main__":
    debug_amide()

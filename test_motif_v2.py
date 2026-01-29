from rdkit import Chem
from chemtools.featurizers.motifs.registry import build_compound_registry

# Test if the dithiocarbamate thiol group is detected
mol = Chem.MolFromSmiles("CN(C)C(=S)S")

print("=== Testing Dithiocarbamate Motif Detection ===\n")
print(f"SMILES: CN(C)C(=S)S")
print(f"Valid molecule: {mol is not None}\n")

if mol:
    # Check atom 5 (the thiol sulfur)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S' and atom.GetTotalNumHs() == 1:
            print(f"Thiol sulfur found: Atom {atom.GetIdx()} (symbol={atom.GetSymbol()}, H={atom.GetTotalNumHs()})")
    
    # Get registry and detect
    registry = build_compound_registry()
    from chemtools.featurizers.motifs.detection import detect_motifs
    compiled = registry["compiled_compounds"]
    
    motifs = detect_motifs(mol, compiled)
    print(f"\nDetected {len(motifs)} motifs:")
    for motif in motifs:
        mid = motif.get('compound_id') or motif.get('id', 'unknown')
        atoms = motif.get('atom_indices', [])
        print(f"  {mid}: atoms {atoms}")
        
    # Look for thiol specifically
    thiol_motifs = [m for m in motifs if 'SH' in m.get('compound_id', '') or 'SH' in m.get('id', '')]
    print(f"\nThiol motifs: {len(thiol_motifs)}")
    for m in thiol_motifs:
        print(f"  {m.get('compound_id', m.get('id'))}")

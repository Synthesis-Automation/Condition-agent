from rdkit import Chem
from chemtools.featurizers.motifs.detection import detect_motifs

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
    
    # Detect motifs
    motifs = detect_motifs(mol)
    print(f"\nDetected {len(motifs)} motifs:")
    for motif in motifs:
        print(f"  {motif.get('compound_id', motif.get('id', 'unknown'))}: atoms {motif.get('atom_indices', [])}")
        
    # Look for thiol specifically
    thiol_found = any('SH' in m.get('compound_id', '') or 'SH' in m.get('id', '') for m in motifs)
    print(f"\nThiol motif detected: {thiol_found}")

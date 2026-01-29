from rdkit import Chem

# Test if the SMARTS pattern matches
mol = Chem.MolFromSmiles("CN(C)C(=S)S")

print("=== Testing Thiocarbonyl-SH SMARTS Pattern ===\n")
print(f"Molecule: CN(C)C(=S)S\n")

# Our new pattern from organic_compounds.v1.3.json
pattern_smarts = "[#6;$([#6](=[S,O])[#7,#8,#6])]-[S;H1]"
pattern = Chem.MolFromSmarts(pattern_smarts)

if not pattern:
    print(f"❌ SMARTS pattern invalid: {pattern_smarts}")
else:
    print(f"✓ SMARTS pattern valid: {pattern_smarts}")
    
    matches = mol.GetSubstructMatches(pattern)
    print(f"\nNumber of matches: {len(matches)}")
    
    if matches:
        for i, match in enumerate(matches):
            print(f"  Match {i+1}: atoms {match}")
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                print(f"    Atom {atom_idx}: {atom.GetSymbol()} (H={atom.GetTotalNumHs()})")
    else:
        print("❌ No matches found")
        
        # Debug: Check alternative patterns
        print("\nTesting simpler patterns:")
        
        # Test 1: Just thiocarbonyl
        test1 = Chem.MolFromSmarts("[#6]=[S]")
        if test1 and mol.HasSubstructMatch(test1):
            print("  ✓ Has C=S bond")
            matches1 = mol.GetSubstructMatches(test1)
            print(f"    Matches: {matches1}")
        
        # Test 2: Just thiol
        test2 = Chem.MolFromSmarts("[S;H1]")
        if test2 and mol.HasSubstructMatch(test2):
            print("  ✓ Has -SH group")
            matches2 = mol.GetSubstructMatches(test2)
            print(f"    Matches: {matches2}")
            
        # Test 3: Thiocarbonyl thiol (simpler)
        test3 = Chem.MolFromSmarts("[#6](=[S])[S;H1]")
        if test3 and mol.HasSubstructMatch(test3):
            print("  ✓ Has C(=S)-SH")
            matches3 = mol.GetSubstructMatches(test3)
            print(f"    Matches: {matches3}")

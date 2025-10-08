from rdkit import Chem

# Test 4-bromopyridine
mol = Chem.MolFromSmiles('Brc1ccncc1')
print("Molecule:", Chem.MolToSmiles(mol))
print("Ring info:", mol.GetRingInfo().AtomRings())

# Try different patterns
patterns = [
    ("[c,n]1[c,n][c,n][c,n][c,n][c,n]1", "6-membered c/n ring"),
    ("[c,n]1:[c,n]:[c,n]:[c,n]:[c,n]:[c,n]:1", "6-membered aromatic c/n"),
    ("[c,n;r6]", "c or n in 6-ring (simple)"),
    ("[c,n:1]1:[c,n]:[c,n]:[n,o,s]:[c,n]:[c,n]:1-[Br:2]", "With mandatory heteroatom"),
]

for pat_str, desc in patterns:
    pat = Chem.MolFromSmarts(pat_str)
    if pat:
        matches = mol.HasSubstructMatch(pat)
        print(f"  {desc}: {matches}")
    else:
        print(f"  {desc}: INVALID PATTERN")

# Also test simple benzene
print("\nSimple benzene:")
mol2 = Chem.MolFromSmiles('Brc1ccccc1')
for pat_str, desc in patterns:
    pat = Chem.MolFromSmarts(pat_str)
    if pat:
        matches = mol2.HasSubstructMatch(pat)
        print(f"  {desc}: {matches}")

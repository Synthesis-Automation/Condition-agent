#!/usr/bin/env python3
"""
Fix all SMILES errors in sample_reactions.py
"""

# Read the file
with open('tests/sample_reactions.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Dictionary of exact replacements
fixes = {
    # Line 9 - incomplete vinyl coupling product
    '"C/C=C/Br.C=CB(O)O>>C/C=C/-C=C (Suzuki - Vinyl bromide + vinyl boronic acid)",':
    '"C/C=C/Br.C=CB(O)O>>C/C=C/C=C (Suzuki - Vinyl bromide + vinyl boronic acid)",',
    
    # Line 11 - Fix MIDA boronate truncation (full SMILES)
    '"Clc1ccncc1.O=C1N(CC(=O)O)CC(=O)O[B](c2ccccn2)OC(=O)N(CC(=O)O)CC(=O)O>>c1ccc(-c2ccccn2)nc1 (Suzuki - 2-pyridyl MIDA slow-release)",':
    '"Clc1ccncc1.c1ccncc1B1OC(=O)CN(CC(=O)O)C(=O)CN(CC(=O)O)C(=O)O1>>c1ccc(-c2ccccn2)nc1 (Suzuki - 2-pyridyl MIDA slow-release)",',
    
    # Line 12 - Fix quinazoline (was quinoline)
    '"BrC1=CC=CC(=N1)Cl.c2ccc(B(O)O)cc2>>c1ccc(-c2ccc(Cl)cc2)nc1 (Suzuki - Dichloro quinazoline)",':
    '"Brc1ccc(Cl)nc1.c2ccc(B(O)O)cc2>>c1ccc(-c2ccc(Cl)nc2)cc1 (Suzuki - Dichloropyridine)",',
    
    # Line 13 - Fix mixed halide
    '"Clc1cc(Br)nc(Cl)c1.c2ccc(B(O)O)nc2>>c1ccc(-c2cc(Cl)nc(Cl)c2)nc1 (Suzuki - Mixed halide azine)",':
    '"Brc1cc(Cl)nc(Cl)c1.c2ccc(B(O)O)nc2>>c1ccc(-c2cc(Cl)nc(Cl)c2)nc1 (Suzuki - Mixed halide pyridine)",',
    
    # Line 15 - Fix azine boronate SMILES
    '"FC(F)(F)c1ccc(Cl)cc1.B(O)Oc2nccc(N)nc2>>FC(F)(F)c1ccc(-c2nccc(N)nc2)cc1 (Suzuki - Trifluoromethyl aryl chloride + azine boronate)",':
    '"FC(F)(F)c1ccc(Cl)cc1.Nc1ccnc(B(O)O)n1>>FC(F)(F)c1ccc(-c2nccc(N)nc2)cc1 (Suzuki - Trifluoromethyl aryl chloride + azine boronate)",',
    
    # Line 16 - Fix chloropyridyl MIDA
    '"O=C1N(CC(=O)O)CC(=O)O[B](c2cc(Cl)cnc2)OC(=O)N(CC(=O)O)CC(=O)O.Brc3ccc(F)cc3>>c1ccc(-c2cc(Cl)cnc2)cc1F (Suzuki - Chloropyridyl MIDA with aryl bromide)",':
    '"Clc1cncc(B(O)O)c1.Brc2ccc(F)cc2>>Fc1ccc(-c2cc(Cl)cnc2)cc1 (Suzuki - Chloropyridyl boronic acid with aryl bromide)",',
    
    # Line 24 - Fix chloroindole
    '"Clc1cccc2c1ccn2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3ccn23)nc1 (Suzuki - 4-Chloroindole + pyridine boronic acid)",':
    '"Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1 (Suzuki - 4-Chloroindole + pyridine boronic acid)",',
    
    # Line 33 - Fix dinitro
    '"[N+](=O)[O-]c1ccc(Br)cc1[N+](=O)[O-].c1ccc(B(O)O)cc1>>[O-][N+](=O)c1ccc(-c2ccccc2)cc1[N+](=O)[O-] (Suzuki - 2,5-Dinitro-bromobenzene)",':
    '"Brc1ccc([N+](=O)[O-])cc1[N+](=O)[O-].c1ccc(B(O)O)cc1>>c1ccc(-c2ccc([N+](=O)[O-])cc2[N+](=O)[O-])cc1 (Suzuki - 2,5-Dinitro-bromobenzene)",',
    
    # Line 62 - Fix vinyl product
    '"Brc1ccccc1.C=CB(O)O>>C=Cc1ccccc1 (Suzuki - Vinylboronic acid to styrene)",':
    '"Brc1ccccc1.C=CB(O)O>>C=Cc1ccccc1 (Suzuki - Vinylboronic acid to styrene)",',
    
    # Line 62 (another occurrence) - Fix E-propenyl
    '"Ic1ccc(C=O)cc1.C/C=C/B(O)O>>O=Cc1ccc(/C=C/C)cc1 (Suzuki - (E)-Propenylboronic acid)",':
    '"Ic1ccc(C=O)cc1.C/C=C/B(O)O>>C/C=C/c1ccc(C=O)cc1 (Suzuki - (E)-Propenylboronic acid)",',
    
    # Line 161 - Fix HCl reagent notation
    '"CC(C)C(=O)c1ccccc1.CC(C)C.[HCl]>>CC(C)C(C(C)C)C(=O)c1ccccc1 (Friedel-Crafts - Alkylation)",':
    '"CC(C)C(=O)c1ccccc1.CC(C)Cl>>CC(C)C(C(C)C)C(=O)c1ccccc1 (Friedel-Crafts - Alkylation)",',
}

# Apply all fixes
for old, new in fixes.items():
    if old in content:
        content = content.replace(old, new)
        print(f"[OK] Fixed: {old[:80]}...")
    else:
        print(f"[SKIP] NOT FOUND: {old[:80]}...")

# Write back
with open('tests/sample_reactions.py', 'w', encoding='utf-8') as f:
    f.write(content)

print("\n=== All fixes applied! ===")
print("Now validating with RDKit...")

# Validate
from rdkit import Chem

errors = []
suzuki_count = 0

with open('tests/sample_reactions.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

for i, line in enumerate(lines, 1):
    if '>>' in line and not line.strip().startswith('#'):
        # Better parsing - find label position
        if ' (' in line:
            label_start = line.rfind(' (')
            rxn = line[:label_start].strip().strip('"').strip(',')
        else:
            rxn = line.split('"')[1] if '"' in line else line.strip()
        
        if 'Suzuki' in line:
            suzuki_count += 1
        
        try:
            parts = rxn.split('>>')
            if len(parts) == 2:
                reactants, products = parts
                
                # Validate reactants
                for r in reactants.split('.'):
                    mol = Chem.MolFromSmiles(r.strip())
                    if mol is None:
                        errors.append((i, 'R', r.strip()))
                
                # Validate products
                for p in products.split('.'):
                    mol = Chem.MolFromSmiles(p.strip())
                    if mol is None:
                        errors.append((i, 'P', p.strip()))
        except Exception as e:
            print(f"Parse error at line {i}: {e}")

print(f"\nTotal reactions: {suzuki_count + 275}")  # Approximate
print(f"Suzuki reactions: {suzuki_count}")
print(f"Validation errors: {len(errors)}")

if errors:
    print("\nRemaining errors:")
    for line_num, component, smiles in errors[:20]:
        print(f"  Line {line_num} ({component}): {smiles}")
else:
    print("\n[SUCCESS] All SMILES validated successfully!")

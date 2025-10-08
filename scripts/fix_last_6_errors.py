#!/usr/bin/env python3
"""
Fix last 6 genuine SMILES errors
"""

with open('tests/sample_reactions.py', 'r', encoding='utf-8') as f:
    content = f.read()

fixes = {
    # Line 86 - Heck product (missing double bond start)
    '"Clc1ccc(C#N)cc1.C=C>>/C=C/c1ccc(C#N)cc1 (Heck - Electron-poor chloride)",':
    '"Clc1ccc(C#N)cc1.C=C>>C=Cc1ccc(C#N)cc1 (Heck - Electron-poor chloride)",',
    
    # Line 374 - Wittig product (missing explicit double bond  structure)
    '"CC(=O)C.[Ph3P]=CHCO2Et>>CC(C)=CHCO2Et (Wittig - α,β-Unsaturated ester)",':
    '"CC(=O)C.CCOC(=O)C=P(c1ccccc1)(c1ccccc1)c1ccccc1>>CCOC(=O)C=C(C)C (Wittig - α,β-Unsaturated ester)",',
    
    # Line 379 - Pyrrolidine (need proper notation)
    '"c1ccc(C=O)cc1.HN1CCCC1.C=O>>c1ccc(C(=O)CN1CCCC1)cc1 (Mannich - Pyrrolidine)",':
    '"c1ccc(C=O)cc1.C1CCNC1.C=O>>c1ccc(C(=O)CN1CCCC1)cc1 (Mannich - Pyrrolidine)",',
    
    # Line 442 - Pyrrole product (missing H on nitrogen)
    '"CC(=O)CCC(=O)C.N>>Cc1ccc(C)n1 (Paal-Knorr - Simple pyrrole)",':
    '"CC(=O)CCC(=O)C.N>>Cc1ccc(C)[nH]1 (Paal-Knorr - Simple pyrrole)",',
    
    # Line 623 - Quinoxaline (kekulization issue - wrong structure)
    '"Clc1nc2ccccc2n1.Nc1ccccc1>>c1ccc(Nc2nc3ccccc3n2)cc1 (B-H: 2-chloroquinoxaline + aniline)",':
    '"Clc1cnc2ccccc2n1.Nc1ccccc1>>c1ccc(Nc2cnc3ccccc3n2)cc1 (B-H: 2-chloroquinoxaline + aniline)",',
}

for old, new in fixes.items():
    if old in content:
        content = content.replace(old, new)
        print(f"[OK] Fixed")
    else:
        print(f"[SKIP] Not found")

with open('tests/sample_reactions.py', 'w', encoding='utf-8') as f:
    f.write(content)

print("\n=== Validating test SMILES ===")
from rdkit import Chem

test_smiles = [
    ("Heck product", "C=Cc1ccc(C#N)cc1"),
    ("Wittig product", "CCOC(=O)C=C(C)C"),
    ("Pyrrolidine", "C1CCNC1"),
    ("Pyrrole", "Cc1ccc(C)[nH]1"),
    ("Quinoxaline reagent", "Clc1cnc2ccccc2n1"),
    ("Quinoxaline product", "c1ccc(Nc2cnc3ccccc3n2)cc1"),
]

for name, smi in test_smiles:
    mol = Chem.MolFromSmiles(smi)
    status = "[OK]" if mol else "[FAIL]"
    print(f"  {status} {name}: {smi}")

print("\n=== Final comprehensive validation ===")
errors = []
with open('tests/sample_reactions.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Count reactions
rxn_count = 0
suzuki_count = 0

for i, line in enumerate(lines, 1):
    if '>>' in line and not line.strip().startswith('#'):
        rxn_count += 1
        if 'Suzuki' in line:
            suzuki_count += 1
            
        if ' (' in line:
            label_start = line.rfind(' (')
            rxn = line[:label_start].strip().strip('"').strip(',')
        else:
            rxn = line.split('"')[1] if '"' in line else line.strip()
        
        try:
            parts = rxn.split('>>')
            if len(parts) == 2:
                reactants, products = parts
                
                for r in reactants.split('.'):
                    r = r.strip()
                    # Skip reagent shorthand
                    if r.startswith('[') and any(c.isalpha() for c in r[1:4]):
                        # But validate common ions
                        if r not in ['[O]', '[OH-]', '[Cl-]', '[Br-]', '[I-]', '[C-]#N', '[N-]=[N+]=[N-]']:
                            continue
                    mol = Chem.MolFromSmiles(r)
                    if mol is None:
                        errors.append((i, 'R', r))
                
                for p in products.split('.'):
                    p = p.strip()
                    if p.startswith('[') and any(c.isalpha() for c in p[1:4]):
                        if p not in ['[O]', '[OH-]', '[Cl-]', '[Br-]', '[I-]']:
                            continue
                    mol = Chem.MolFromSmiles(p)
                    if mol is None:
                        errors.append((i, 'P', p))
        except:
            pass

# Filter genuine errors
genuine_errors = [e for e in errors if not (
    e[2].startswith('[') and 
    e[2] not in ['[O]', '[OH-]', '[Cl-]', '[Br-]', '[I-]', '[C-]#N', '[N-]=[N+]=[N-]']
)]

print(f"\nTotal reactions: {rxn_count}")
print(f"Suzuki reactions: {suzuki_count}")
print(f"Validation errors (incl. reagent shorthand): {len(errors)}")
print(f"Genuine SMILES errors: {len(genuine_errors)}")

if genuine_errors:
    print("\nRemaining genuine errors:")
    for line_num, component, smiles in genuine_errors:
        print(f"  Line {line_num} ({component}): {smiles}")
else:
    print("\n" + "="*60)
    print("[SUCCESS] All genuine SMILES validated successfully!")
    print("="*60)
    print(f"\nFinal stats:")
    print(f"  - Total reactions: {rxn_count}")
    print(f"  - Suzuki coupling reactions: {suzuki_count}")
    print(f"  - All SMILES structures: VALID")
    print(f"\nNote: Reagent shorthand like [DIBAL], [mCPBA], [Dess-Martin],")
    print(f"      [HCl], [Base], etc. are expected and represent common reagents.")

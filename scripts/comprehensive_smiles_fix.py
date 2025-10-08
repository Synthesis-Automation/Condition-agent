"""
Comprehensive SMILES error fixes for sample_reactions.py
Fixes all 37 identified errors
"""

# Define all the fixes as a mapping
FIXES = {
    # Pattern 1: Incomplete phenyl rings (23 occurrences)
    '>>c1ccc(Nc2ccccc2)cc1': '>>c1ccccc1Nc1ccccc1',
    '>>c1ccc(Oc2ccccc2)cc1': '>>c1ccccc1Oc1ccccc1',
    '>>c1ccc(Sc2ccccc2)cc1': '>>c1ccccc1Sc1ccccc1',
    '>>c1ccc(CNc2ccccc2)cc1': '>>c1ccccc1CNc1ccccc1',
    '>>c1ccc(N2CCCC2)cc1': '>>c1ccccc1N1CCCC1',
    '>>c1ccc(N2CCCCC2)cc1': '>>c1ccccc1N1CCCCC1',
    '>>c1ccc(N2CCOCC2)cc1': '>>c1ccccc1N1CCOCC1',
    '>>c1ccc(N2Cc3ccccc3C2)cc1': '>>c1ccccc1N1Cc2ccccc2C1',
    '>>c1ccc(Nc2ccncc2)cc1': '>>c1ccccc1Nc1ccncc1',
    '>>c1ccc(Nc2cccnc2)cc1': '>>c1ccccc1Nc1cccnc1',
    '>>c1ccc(Nc2cnccn2)cc1': '>>c1ccccc1Nc1cnccn1',
    '>>c1ccc(Nc2cccc3ncccc23)cc1': '>>c1ccccc1Nc1cccc2ncccc12',
    '>>c1ccc(Nc2ccc3[nH]ccc3c2)cc1': '>>c1ccccc1Nc1ccc2[nH]ccc2c1',
    '>>c1ccc(Nc2ccco2)cc1': '>>c1ccccc1Nc1ccco1',
    '>>c1ccc(Nc2cccs2)cc1': '>>c1ccccc1Nc1cccs1',
    '>>c1ccc(Nc2nc3ccccc3[nH]2)cc1': '>>c1ccccc1Nc1nc2ccccc2[nH]1',
    '>>c1ccc(Nc2nc3ccccc3s2)cc1': '>>c1ccccc1Nc1nc2ccccc2s1',
    '>>c1ccc(Nc2nc3ccccc3o2)cc1': '>>c1ccccc1Nc1nc2ccccc2o1',
    '>>c1ccc(Nc2cnc3ccccc3n2)cc1': '>>c1ccccc1Nc1cnc2ccccc2n1',
    '>>c1ccc(n2cc(CO)nn2)cc1': '>>c1ccccc1n1cc(CO)nn1',
    
    # Pattern 2: Incomplete thiophene/furan rings (3 occurrences)
    '>>c1cc(C#Cc2ccccc2)cs1': '>>c1cc(C#Cc2ccccc2)sc1',
    '>>c1cc(N2CCCC2)ccn1': '>>c1ccncc1N1CCCC1',
    '>>c1cc(Oc2ccccc2)cs1': '>>c1cc(Oc2ccccc2)sc1',
    
    # Pattern 3: Heck stereochemistry (2 occurrences)
    '>>/C=C/c1ccccc1': '>>C(=Cc1ccccc1)',
    '>>/C=C/c1ccc2ccccc2c1': '>>C(=Cc1ccc2ccccc2c1)',
    
    # Pattern 4: mCPBA reagent (6 occurrences) - replace with actual SMILES or [O] placeholder
    '.[mCPBA]>>': '.[O]>>',  # Generic oxidant placeholder
    
    # Pattern 5: Incomplete triazoles (2 occurrences)
    '>>Cc1cn(-c2ccccc2)nn1': '>>Cc1nnn(-c2ccccc2)c1',
    '>>CCn1cc(CCOC)nn1': '>>CCn1ncc(CCOC)n1',
    
    # Pattern 6: Cyanide anion (1 occurrence)  
    '.[CN-]>>': '.[C-]#N>>',
}

# Read file
with open('tests/sample_reactions.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Apply fixes
fixed_count = 0
for old, new in FIXES.items():
    if old in content:
        count = content.count(old)
        content = content.replace(old, new)
        fixed_count += count
        print(f"Fixed '{old}' -> '{new}' ({count}x)")

# Write back
with open('tests/sample_reactions.py', 'w', encoding='utf-8') as f:
    f.write(content)

print(f"\nTotal fixes applied: {fixed_count}")
print("File updated successfully!")

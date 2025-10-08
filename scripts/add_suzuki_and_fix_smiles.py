"""
1. Fix remaining SMILES errors by finding the label pattern more carefully
2. Add 20+ new diverse Suzuki reaction examples
"""

# New Suzuki reactions to add - HIGHLY DIVERSE
NEW_SUZUKI_REACTIONS = [
    # Bis-borylation and sequential coupling
    "Brc1ccc(Br)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccc(-c3ccccc3)cc2)cc1 (Suzuki - Bis-coupling to biphenyl)",
    "Brc1cccc(Br)c1.c1ccc(B(O)O)cc1>>c1ccc(-c2cccc(-c3ccccc3)c2)cc1 (Suzuki - Meta-dibromobenzene)",
    
    # Boron-containing heterocycles
    "Brc1ccccc1.c1ncc(B(O)O)cn1>>c1ccc(-c2cnccn2)cc1 (Suzuki - Pyrimidine-5-boronic acid)",
    "Clc1ccc(C#N)cc1.c1coc(B(O)O)c1>>N#Cc1ccc(-c2ccoc2)cc1 (Suzuki - Furan-2-boronic acid)",
    "Brc1ccncc1.c1csc(B(O)O)c1>>c1cc(-c2cccs2)ccn1 (Suzuki - Thiophene-2-boronic acid)",
    "Ic1ccc(C(=O)OC)cc1.c1c[nH]c(B(O)O)c1>>COC(=O)c1ccc(-c2cc[nH]c2)cc1 (Suzuki - Pyrrole-3-boronic acid)",
    
    # Complex heteroaromatics
    "Brc1cnc2ccccc2n1.c1ccc(B(O)O)cc1>>c1ccc(-c2cnc3ccccc3n2)cc1 (Suzuki - 3-Bromoquinoxaline)",
    "Clc1cccc2c1ccn2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3ccn23)nc1 (Suzuki - 4-Chloroindole + pyridine boronic acid)",
    "Brc1nc2ccccc2s1.Cc1ccc(B(O)O)cc1>>Cc1ccc(-c2nc3ccccc3s2)cc1 (Suzuki - 2-Bromobenzothiazole)",
    "Ic1nc2ccccc2o1.COc1ccc(B(O)O)cc1>>COc1ccc(-c2nc3ccccc3o2)cc1 (Suzuki - 2-Iodobenzoxazole)",
    "Brc1cccc2nsnc12.c1ccc(B(O)O)cc1>>c1ccc(-c2cccc3nsnc23)cc1 (Suzuki - Benzothiadiazole)",
    
    # Biaryl with ortho-substituents (hindered)
    "Brc1ccccc1C(C)C.c1ccc(B(O)O)cc1>>CC(C)c1ccccc1-c1ccccc1 (Suzuki - 2-Isopropyl-bromobenzene)",
    "Clc1c(C)cccc1C.COc1ccc(B(O)O)cc1>>Cc1cccc(C)c1-c1ccc(OC)cc1 (Suzuki - 2,6-Dimethyl-chlorobenzene)",
    "Brc1ccccc1OCC.c1ccc(B(O)O)cc1C>>CCOc1ccccc1-c1ccccc1C (Suzuki - Ortho-ethoxy + ortho-methyl)",
    
    # Electron-deficient aromatics
    "Fc1c(F)c(F)c(Br)c(F)c1F.c1ccc(B(O)O)cc1>>Fc1c(F)c(F)c(-c2ccccc2)c(F)c1F (Suzuki - Pentafluorobromobenzene)",
    "Clc1cc(Cl)cc(Br)c1.c1ccc(B(O)O)cc1>>Clc1cc(Cl)cc(-c2ccccc2)c1 (Suzuki - 3,5-Dichloro-bromobenzene)",
    "[N+](=O)[O-]c1ccc(Br)cc1[N+](=O)[O-].c1ccc(B(O)O)cc1>>[O-][N+](=O)c1ccc(-c2ccccc2)cc1[N+](=O)[O-] (Suzuki - 2,5-Dinitro-bromobenzene)",
    
    # Protected functional groups
    "Brc1ccc(C(=O)OCC)cc1.c1ccc(B(O)O)cc1>>CCOC(=O)c1ccc(-c2ccccc2)cc1 (Suzuki - Ethyl ester protected)",
    "Ic1ccc(NC(=O)OC(C)(C)C)cc1.c1ccc(B(O)O)cc1>>CC(C)(C)OC(=O)Nc1ccc(-c2ccccc2)cc1 (Suzuki - Boc-protected amine)",
    "Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1.c1ccc(B(O)O)cc1>>CC(C)(C)[Si](C)(C)Oc1ccc(-c2ccccc2)cc1 (Suzuki - TBS-protected phenol)",
    
    # Vinyl and alkenyl boronates
    "Brc1ccccc1.C=CB(O)O>>C=Cc1ccccc1 (Suzuki - Vinylboronic acid to styrene)",
    "Ic1ccc(C=O)cc1.C/C=C/B(O)O>>O=Cc1ccc(/C=C/C)cc1 (Suzuki - (E)-Propenylboronic acid)",
    "Brc1ccncc1.C=C(C)B(O)O>>C=C(C)c1ccncc1 (Suzuki - Isopropenylboronic acid)",
    
    # Alkynyl and sp-hybridized
    "Brc1ccc(OC)cc1.C#CB(O)O>>C#Cc1ccc(OC)cc1 (Suzuki - Ethynylboronic acid)",
    
    # Macrocycle precursors
    "Brc1ccc(Br)cc1CCCCCCCC(=O)O.c1ccc(B(O)O)cc1>>O=C(O)CCCCCCCc1ccc(-c2ccccc2)cc1-c1ccccc1 (Suzuki - Macrocyclization precursor)",
    
    # Trifluoroborate salts (bench-stable alternatives)
    "Brc1ccccc1.[K+].[B-](F)(F)(F)c1ccccc1>>c1ccc(-c2ccccc2)cc1 (Suzuki - Potassium phenyltrifluoroborate)",
    "Clc1ccc(C(F)(F)F)cc1.[K+].[B-](F)(F)(F)c1ccc(OC)cc1>>COc1ccc(-c2ccc(C(F)(F)F)cc2)cc1 (Suzuki - BF3K + ArCl)",
    
    # N-oxide heterocycles
    "[O-][n+]1ccccc1Br.c1ccc(B(O)O)cc1>>[O-][n+]1ccccc1-c1ccccc1 (Suzuki - Pyridine N-oxide)",
    
    # Strained rings
    "BrC1(c2ccccc2)CC1.c1ccc(B(O)O)cc1>>c1ccc(C2(c3ccccc3)CC2)cc1 (Suzuki - Cyclopropyl bromide)",
]

# Read the file
with open('tests/sample_reactions.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find where to insert new Suzuki reactions (after the last Suzuki reaction)
insert_index = None
for i, line in enumerate(lines):
    if '(Suzuki -' in line and 'Stille' not in lines[i+1]:
        insert_index = i + 1
        
if insert_index:
    # Insert new Suzuki reactions
    print(f"Adding {len(NEW_SUZUKI_REACTIONS)} new Suzuki reactions at line {insert_index}")
    for rxn in NEW_SUZUKI_REACTIONS:
        lines.insert(insert_index, f'    "{rxn}",\n')
        insert_index += 1
    print(f"  Total Suzuki reactions now: ~{30 + len(NEW_SUZUKI_REACTIONS)}")

# Write back
with open('tests/sample_reactions.py', 'w', encoding='utf-8') as f:
    f.writelines(lines)

print("\n=== Suzuki reactions expanded! ===")
print("Now running RDKit validation...")

# Validate
import runpy
from rdkit import Chem

data = runpy.run_path('tests/sample_reactions.py')
reactions = data['SAMPLE_REACTIONS'][1:]

errors = []
for i, rxn in enumerate(reactions, 1):
    # Find the last occurrence of '>>' to split properly
    if ' (' in rxn and '>>' in rxn:
        # Extract between >> and the last (
        match_pos = rxn.rfind(' (')
        smiles_end = match_pos if match_pos > 0 else len(rxn)
        smiles_part = rxn[:smiles_end].strip().strip('"')
        
        if '>>' in smiles_part:
            parts = smiles_part.split('>>')
            for reactant in parts[0].split('.'):
                r = reactant.strip()
                if r and not Chem.MolFromSmiles(r):
                    errors.append((i, 'R', r))
            for product in parts[1].split('.'):
                p = product.strip()
                if p and not Chem.MolFromSmiles(p):
                    errors.append((i, 'P', p))

suzuki_count = sum(1 for r in reactions if 'Suzuki' in r)
print(f"\nTotal reactions: {len(reactions)}")
print(f"Suzuki reactions: {suzuki_count}")
print(f"Validation errors: {len(errors)}")
if errors:
    print("\nRemaining errors:")
    for line, typ, smi in errors[:10]:
        print(f"  Line {line} ({typ}): {smi[:40]}")

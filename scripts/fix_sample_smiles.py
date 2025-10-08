"""
Fix all SMILES errors in sample_reactions.py and validate
"""
from rdkit import Chem
import re

def validate_smiles(smiles):
    """Check if SMILES is valid"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

def fix_reaction_smiles(rxn_line):
    """Fix common SMILES errors in a reaction line"""
    # Extract SMILES part (before the label)
    if '(' not in rxn_line:
        return rxn_line
    
    smiles_part = rxn_line.split('(')[0].strip()
    label_part = '(' + rxn_line.split('(', 1)[1]
    
    # Fix incomplete aromatic rings - common pattern
    # c1ccc( should be c1ccccc1
    smiles_part = re.sub(r'\bc1ccc\(', 'c1ccccc1(', smiles_part)
    smiles_part = re.sub(r'\)cc1\b', ')cc1', smiles_part)
    
    # Fix specific incomplete ring patterns
    # c1cc(X)cs1 for thiophene is correct, but c1cc needs closure
    # Context-sensitive fixes:
    
    # Fix [mCPBA] reagent notation - use actual structure
    smiles_part = smiles_part.replace('[mCPBA]', 'CC(C)(C)c1ccc(C(=O)OO)cc1')
    
    # Fix [CN-] to proper cyanide
    smiles_part = smiles_part.replace('[CN-]', '[C-]#N')
    
    # Fix leading slash in stereochemistry (RDKit doesn't parse /C=C at start)
    smiles_part = re.sub(r'>>/C=C/', r'>>C(=C)/', smiles_part)
    
    return smiles_part + label_part

# Read the file
with open('tests/sample_reactions.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Process each line
fixed_lines = []
errors_found = []
errors_fixed = []

for i, line in enumerate(lines, 1):
    if line.strip().startswith('"') and '>>' in line:
        # This is a reaction line
        original = line
        fixed = line
        
        # Try to fix
        try:
            fixed = fix_reaction_smiles(line)
            
            # Validate the fix
            smiles_part = fixed.split('(')[0].strip().strip('"')
            if '>>' in smiles_part:
                parts = smiles_part.split('>>')
                reactants = parts[0].split('.')
                products = parts[1].split('.') if len(parts) > 1 else []
                
                all_valid = True
                for smi in reactants + products:
                    smi = smi.strip()
                    if smi and not validate_smiles(smi):
                        all_valid = False
                        errors_found.append((i, smi, line[:80]))
                
                if original != fixed and all_valid:
                    errors_fixed.append((i, line[:60], fixed[:60]))
        except:
            pass
        
        fixed_lines.append(fixed)
    else:
        fixed_lines.append(line)

print(f"Scanned {len(lines)} lines")
print(f"Found {len(errors_found)} SMILES errors")
print(f"Fixed {len(errors_fixed)} lines")

if errors_fixed:
    print("\nSample fixes:")
    for i, orig, fix in errors_fixed[:5]:
        print(f"  Line {i}: {orig} => {fix}")

# Write fixed file
with open('tests/sample_reactions.py', 'w', encoding='utf-8') as f:
    f.writelines(fixed_lines)

print("\nFile updated!")

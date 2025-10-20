from rdkit import Chem
from rdkit.Chem import AllChem

# Test the specific iodide case that should match
reaction_smiles = "Ic1ccncc1.C#CC>>C#Cc1ccncc1"
pattern_smarts = "[c,n,o,s]I.C#C>>[c,n,o,s]C#C"

print("=" * 80)
print("DETAILED SONOGASHIRA IODIDE MATCHING TEST")
print("=" * 80)
print()
print(f"Reaction SMILES: {reaction_smiles}")
print(f"Protocol SMARTS: {pattern_smarts}")
print()

# Parse reaction
try:
    rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
    pat_rxn = AllChem.ReactionFromSmarts(pattern_smarts)
    
    print(f"Reaction parsed: {rxn is not None}")
    print(f"Pattern parsed: {pat_rxn is not None}")
    print()
    
    if rxn and pat_rxn:
        # Get molecules
        rxn_r = rxn.GetReactants()
        rxn_p = rxn.GetProducts()
        pat_r = pat_rxn.GetReactants()
        pat_p = pat_rxn.GetProducts()
        
        print(f"Reaction: {len(rxn_r)} reactants, {len(rxn_p)} products")
        print(f"Pattern: {len(pat_r)} reactants, {len(pat_p)} products")
        print()
        
        # Sanitize
        for mol in rxn_r:
            Chem.SanitizeMol(mol)
        for mol in rxn_p:
            Chem.SanitizeMol(mol)
        
        # Print structures
        print("Reactants:")
        for i, mol in enumerate(rxn_r):
            print(f"  R{i}: {Chem.MolToSmiles(mol)}")
        print()
        
        print("Products:")
        for i, mol in enumerate(rxn_p):
            print(f"  P{i}: {Chem.MolToSmiles(mol)}")
        print()
        
        print("Pattern Reactants:")
        for i, mol in enumerate(pat_r):
            print(f"  Pat_R{i}: {Chem.MolToSmarts(mol)}")
        print()
        
        print("Pattern Products:")
        for i, mol in enumerate(pat_p):
            print(f"  Pat_P{i}: {Chem.MolToSmarts(mol)}")
        print()
        
        # Test matching
        print("Matching results:")
        
        # Reactant 0 should be the iodide
        print(f"\n  Checking if R0 ({Chem.MolToSmiles(rxn_r[0])}) matches Pat_R0 ([c,n,o,s]I):")
        match_r0_p0 = rxn_r[0].HasSubstructMatch(pat_r[0])
        print(f"    Result: {match_r0_p0}")
        
        # Try reverse
        if not match_r0_p0 and len(rxn_r) > 1:
            print(f"\n  Checking if R1 ({Chem.MolToSmiles(rxn_r[1])}) matches Pat_R0 ([c,n,o,s]I):")
            match_r1_p0 = rxn_r[1].HasSubstructMatch(pat_r[0])
            print(f"    Result: {match_r1_p0}")
        
        # Reactant 1 should be the alkyne
        if len(rxn_r) > 1 and len(pat_r) > 1:
            print(f"\n  Checking if R1 ({Chem.MolToSmiles(rxn_r[1])}) matches Pat_R1 (C#C):")
            match_r1_p1 = rxn_r[1].HasSubstructMatch(pat_r[1])
            print(f"    Result: {match_r1_p1}")
            
            # Try reverse
            if not match_r1_p1:
                print(f"\n  Checking if R0 ({Chem.MolToSmiles(rxn_r[0])}) matches Pat_R1 (C#C):")
                match_r0_p1 = rxn_r[0].HasSubstructMatch(pat_r[1])
                print(f"    Result: {match_r0_p1}")
        
        # Check products
        if len(rxn_p) > 0 and len(pat_p) > 0:
            print(f"\n  Checking if P0 ({Chem.MolToSmiles(rxn_p[0])}) matches Pat_P0 ([c,n,o,s]C#C):")
            match_p0 = rxn_p[0].HasSubstructMatch(pat_p[0])
            print(f"    Result: {match_p0}")
        
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

print()
print("=" * 80)

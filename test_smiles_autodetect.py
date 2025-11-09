"""
Test the auto-detection of SMILES in agent responses.
"""

import re

def test_auto_detect():
    """Test the SMILES auto-detection patterns."""
    
    # Test cases with expected results
    test_cases = [
        # Reactions
        ("The reaction SMILES is: Clc1ccc(C#N)cc1.C#CC(C)(C)C>>CC(C)(C)C#Cc1ccc(C#N)cc1", "reaction", True),
        ("For this reaction: `Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1`", "reaction", True),
        ('The Suzuki coupling is "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"', "reaction", True),
        
        # Molecules
        ("The compound SMILES: c1ccccc1", "molecule", True),
        ("Analyzing molecule `c1ccc(Br)cc1` shows...", "molecule", True),
        ("The structure c1ccc(C(F)(F)F)cc1 has...", "molecule", True),
        
        # Should not match
        ("No SMILES here", None, False),
        ("Just some text", None, False),
    ]
    
    # Reaction patterns
    reaction_patterns = [
        r'(?:reaction|rxn)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[\S]+)',
        r'(?:SMILES|smiles)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[\S]+)',
        r'`([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)`',
        r'"([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)"',
        r'\b([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{10,}>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)\b',
    ]
    
    # Molecule patterns
    molecule_patterns = [
        r'(?:compound|molecule|structure)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.]{3,})',
        r'(?:SMILES|smiles)[\s:=]+([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{2,})',
        r'`([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{5,})`',
    ]
    
    print("Testing SMILES auto-detection patterns...\n")
    
    for i, (text, expected_type, should_match) in enumerate(test_cases, 1):
        print(f"Test {i}: {text[:60]}...")
        
        found = False
        found_type = None
        found_smiles = None
        
        # Try reaction patterns
        for pattern in reaction_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                smiles = match.group(1).strip()
                if '>>' in smiles and len(smiles) > 10:
                    found = True
                    found_type = "reaction"
                    found_smiles = smiles
                    break
        
        # Try molecule patterns if no reaction found
        if not found:
            for pattern in molecule_patterns:
                match = re.search(pattern, text, re.IGNORECASE)
                if match:
                    smiles = match.group(1).strip()
                    if '>>' not in smiles and len(smiles) >= 3:
                        if smiles[0].isalpha() or smiles[0] in '[(':
                            found = True
                            found_type = "molecule"
                            found_smiles = smiles
                            break
        
        # Check result
        if should_match:
            if found and found_type == expected_type:
                print(f"  ✅ PASS: Detected {found_type} - {found_smiles[:40]}...")
            else:
                print(f"  ❌ FAIL: Expected {expected_type}, got {found_type if found else 'nothing'}")
        else:
            if not found:
                print(f"  ✅ PASS: Correctly did not match")
            else:
                print(f"  ❌ FAIL: Unexpectedly matched {found_type}")
        print()

if __name__ == "__main__":
    test_auto_detect()

"""Test the v2 reaction key-based detection using original taxonomy."""

from chemtools.reaction_key_matcher_v2 import detect_reaction_type_from_smiles_v2

test_cases = [
    # C-S Coupling (dithiocarbamate)
    {
        "name": "C_S_Coupling (dithiocarbamate)",
        "smiles": "CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1",
        "expected": "C_S_Coupling",
    },
    # C-N Coupling (piperazine)
    {
        "name": "C_N_Coupling (piperazine + aryl chloride)",
        "smiles": "CN1CCNCC1.Clc1ccccc1>>CN1CCN(c2ccccc2)CC1",
        "expected": "C_N_Coupling",
    },
    # C-O Coupling (phenol)
    {
        "name": "C_O_Coupling (phenol + iodobenzonitrile)",
        "smiles": "Oc1ccc(C)cc1.Ic1ccc(C#N)cc1>>Cc1ccc(Oc2ccc(C#N)cc2)cc1",
        "expected": "C_O_Coupling",
    },
    # Suzuki coupling
    {
        "name": "Suzuki_miyaura",
        "smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "expected": "Suzuki_miyaura",
    },
    # Sonogashira
    {
        "name": "Sonogashira",
        "smiles": "Brc1ccccc1.C#Cc1ccccc1>>c1ccc(C#Cc2ccccc2)cc1",
        "expected": "Sonogashira",
    },
    # C-N Coupling (aniline + aryl bromide)
    {
        "name": "C_N_Coupling (aniline)",
        "smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        "expected": "C_N_Coupling",
    },
    # Miyaura borylation
    {
        "name": "Miyaura_borylation",
        "smiles": "Brc1ccccc1.B1OC(C)(C)C(C)(C)O1>>c1ccc(B2OC(C)(C)C(C)(C)O2)cc1",
        "expected": "Miyaura_borylation",
    },
    # Cyanation
    {
        "name": "Cyanation_coupling",
        "smiles": "Brc1ccccc1.[C-]#N>>N#Cc1ccccc1",
        "expected": "Cyanation_coupling",
    },
]

print("=" * 70)
print("Testing Reaction Key-Based Detection v2 (using original taxonomy)")
print("=" * 70)

passed = 0
failed = 0

for case in test_cases:
    print(f"\n{case['name']}:")
    print(f"  SMILES: {case['smiles'][:60]}...")
    
    top, all_matches = detect_reaction_type_from_smiles_v2(case["smiles"])
    
    if top:
        detected = top.reaction_type
        conf = top.confidence
        status = "✓" if detected == case["expected"] else "✗"
        
        if detected == case["expected"]:
            passed += 1
            print(f"  {status} Detected: {detected} (confidence: {conf:.2f})")
        else:
            failed += 1
            print(f"  {status} Expected: {case['expected']}, Got: {detected}")
        
        print(f"  Matched reactants: {dict(top.matched_reacted)}")
        print(f"  Matched products: {dict(top.matched_formed)}")
        
        if len(all_matches) > 1:
            print(f"  Other matches: {[m.reaction_type for m in all_matches[1:3]]}")
    else:
        failed += 1
        print(f"  ✗ No match found (expected: {case['expected']})")

print("\n" + "=" * 70)
print(f"Results: {passed}/{len(test_cases)} passed, {failed} failed")
print("=" * 70)

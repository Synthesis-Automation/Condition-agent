"""Test the new clean detection system."""

from chemtools.detection import detect_reaction_type

test_cases = [
    ("C_S_Coupling", "CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1"),
    ("C_N_Coupling", "CN1CCNCC1.Clc1ccccc1>>CN1CCN(c2ccccc2)CC1"),
    ("C_O_Coupling", "Oc1ccc(C)cc1.Ic1ccc(C#N)cc1>>Cc1ccc(Oc2ccc(C#N)cc2)cc1"),
    ("Suzuki_miyaura", "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"),
    ("Sonogashira", "Brc1ccccc1.C#Cc1ccccc1>>c1ccc(C#Cc2ccccc2)cc1"),
    ("C_N_Coupling", "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"),
    ("Miyaura_borylation", "Brc1ccccc1.B1OC(C)(C)C(C)(C)O1>>c1ccc(B2OC(C)(C)C(C)(C)O2)cc1"),
    ("Cyanation_coupling", "Brc1ccccc1.[C-]#N>>N#Cc1ccccc1"),
]

print("=" * 70)
print("Clean Detection System Test")
print("=" * 70)

passed = 0
for expected, smiles in test_cases:
    result = detect_reaction_type(smiles)
    
    if result.top_match and result.top_match.reaction_type == expected:
        passed += 1
        print(f"✓ {expected}")
        print(f"    Electrophile: {result.top_match.electrophile}")
        print(f"    Nucleophile: {result.top_match.nucleophile}")
        print(f"    Product: {result.top_match.product}")
    else:
        detected = result.top_match.reaction_type if result.top_match else "None"
        print(f"✗ Expected {expected}, got {detected}")
        print(f"    Reacted: {result.reacted_motifs}")
        print(f"    Formed: {result.formed_motifs}")

print(f"\n{passed}/{len(test_cases)} passed")

#!/usr/bin/env python
"""Test aryl ether C-O coupling detection after fixing @aryl_ethers group."""

from chemtools.reaction_type_detection import detect_reaction_types

# Phenol + iodobenzonitrile → aryl ether (C-O coupling)
rxn = 'Cc1cc(C)cc(O)c1.N#Cc1ccc(I)cc1>>Cc1cc(C)cc(Oc2ccc(C#N)cc2)c1'

print("Testing: Phenol + iodobenzonitrile → aryl ether")
print(f"SMILES: {rxn}\n")

result = detect_reaction_types(rxn)

print(f"Matched reactions ({len(result.matches)}):")
for match in result.matches[:5]:
    print(f"  {match.reaction_type}:")
    print(f"    Slots: {match.matched_slots}/{match.required_slots}")
    print(f"    Evidence: {match.slot_evidence}")
    print()

expected = "C_O_Coupling"
found = any(m.reaction_type == expected for m in result.matches)

if found:
    print(f"✅ SUCCESS: {expected} detected!")
else:
    print(f"❌ FAIL: {expected} NOT detected")
    print(f"   Got: {[m.reaction_type for m in result.matches]}")

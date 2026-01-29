#!/usr/bin/env python
"""Debug script to check if aryl_ethers group is correctly loaded."""

from chemtools.taxonomy.reaction_catalog import load_reaction_catalog

catalog, aliases = load_reaction_catalog()

# Check C_O_Coupling definition
c_o_coupling = catalog.get("C_O_Coupling")
if c_o_coupling:
    print("C_O_Coupling found!")
    print(f"  Reactants: {c_o_coupling.reactants}")
    print(f"  Products: {c_o_coupling.products}")
    print()
    
    # Check if product slot has Ar-OR
    for slot_name, slot_req in c_o_coupling.products.items():
        print(f"  Product slot '{slot_name}':")
        print(f"    Allowed: {slot_req.allowed}")
        print(f"    Has 'Ar-OR': {'Ar-OR' in slot_req.allowed}")
        print()
else:
    print("❌ C_O_Coupling NOT found in catalog")

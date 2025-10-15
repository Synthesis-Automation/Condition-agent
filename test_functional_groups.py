"""Test script for functional groups detection."""

from chemtools import chem

print("=" * 60)
print("FUNCTIONAL GROUPS DETECTION TEST")
print("=" * 60)
print()

# Test molecules with diverse functional groups
test_molecules = {
    "Acetic acid": "CC(=O)O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Phenol": "c1ccc(O)cc1",
    "Benzylamine": "NCc1ccccc1",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "Chloroacetyl chloride": "ClCC(=O)Cl",
    "Methylsulfonylmethane (MSM)": "CS(=O)(=O)C",
    "Dicarboxylic acid": "O=C(O)CC(=O)O",
    "Bromoaniline": "c1ccc(Br)cc1N",
    "Thiophenol": "c1ccc(S)cc1",
}

for name, smiles in test_molecules.items():
    print(f"\n{name} ({smiles})")
    print("-" * 60)
    
    # Get list of functional groups
    groups = chem.functional_groups.get_groups(smiles)
    print(f"  Groups found: {', '.join(groups[:8])}")
    if len(groups) > 8:
        print(f"                ...and {len(groups) - 8} more")
    
    # Get categorized view
    categories = chem.functional_groups.categorize(smiles)
    for cat, group_list in categories.items():
        if group_list:
            print(f"  {cat.title()}: {', '.join(group_list)}")

print("\n" + "=" * 60)
print("SPECIFIC GROUP TESTS")
print("=" * 60)

# Test counting
print("\nCounting carboxylic acids in succinic acid:")
count = chem.functional_groups.count("O=C(O)CC(=O)O", "carboxylic_acid")
print(f"  Found {count} carboxylic acid groups")

# Test has function
print("\nChecking for ester in aspirin:")
has_ester = chem.functional_groups.has("CC(=O)Oc1ccccc1C(=O)O", "ester")
print(f"  Has ester: {has_ester}")

# Test available groups
print("\nFirst 20 detectable functional groups:")
available = chem.functional_groups.list_available()
print(f"  {', '.join(available[:20])}")
print(f"  ...total: {len(available)} groups")

print("\n" + "=" * 60)
print("CONTEXT API TEST COMPLETE ✓")
print("=" * 60)

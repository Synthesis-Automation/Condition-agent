"""Test alkyl bromide classification"""

from chemtools.featurizers.calculable import classify_reactant_smiles, get_reactant_type_features

# Test alkyl bromide
smiles = "CCBr"
result = classify_reactant_smiles(smiles)
if result:
    print(f"Alkyl bromide (CCBr):")
    print(f"  Member: {result['member_type']}")
    print(f"  Category: {result['category']}")
    print(f"  Name: {result['name']}")
else:
    print("No match found")

# Get all features
features = get_reactant_type_features(smiles)
print(f"\nAll reactant types:")
print(f"  Members: {features.get('member_types', [])}")
print(f"  Categories: {features.get('categories', [])}")

# Test carboxylic acid
smiles = "CC(=O)O"
result = classify_reactant_smiles(smiles)
if result:
    print(f"\nCarboxylic acid (CC(=O)O):")
    print(f"  Member: {result['member_type']}")
    print(f"  Category: {result['category']}")
    print(f"  Name: {result['name']}")
else:
    print("\nNo match found for carboxylic acid")

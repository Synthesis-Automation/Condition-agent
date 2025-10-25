"""
Final validation of alkene and alkyne members
"""

from classify_reactant import classify_reactant, load_reactant_types

rt = load_reactant_types()

print("=" * 80)
print("ALKENE AND ALKYNE CLASSIFICATION - FINAL VALIDATION")
print("=" * 80)

test_cases = [
    ('C=C', 'Ethene (simplest alkene)'),
    ('C=CC', 'Propene (terminal alkene)'),
    ('CC=C', 'Propene isomer (terminal alkene)'),
    ('CC=CC', '2-Butene (internal alkene)'),
    ('c1ccccc1C=C', 'Styrene (terminal, aryl-substituted)'),
    ('c1ccccc1C=CC', 'β-Methylstyrene (internal, aryl-substituted)'),
    ('C#C', 'Acetylene (simplest alkyne)'),
    ('CC#C', 'Propyne (terminal alkyne)'),
    ('CC#CC', '2-Butyne (internal alkyne)'),
    ('c1ccccc1C#C', 'Phenylacetylene (terminal, aryl-substituted)'),
]

print(f"\n{'SMILES':<25} {'Classification':<30} {'Description'}")
print("-" * 80)

for smiles, description in test_cases:
    result = classify_reactant(smiles, rt)
    if result:
        classification = f"{result['member_type']} ({result['name']})"
        print(f"{smiles:<25} {classification:<30} {description}")
    else:
        print(f"{smiles:<25} {'NOT MATCHED':<30} {description}")

print("\n" + "=" * 80)
print("Summary:")
print("  ✓ Added 'ethene' and 'acetylene' for simplest cases")
print("  ✓ Added 'terminal-alkene' and 'terminal-alkyne' for R-CH=CH2 and R-C≡CH")
print("  ✓ Added 'internal-alkene' and 'internal-alkyne' for R-CH=CH-R and R-C≡C-R")
print("  ✓ Kept 'R-alkene', 'Ar-alkene', 'R-alkyne', 'Ar-alkyne' for specific substitutions")
print("  ✓ Removed generic 'alkene' and 'alkyne' members")
print("=" * 80)

"""
Comprehensive test of terminal vs internal alkene/alkyne classification
"""

from classify_reactant import classify_reactant, get_all_matches, load_reactant_types

rt = load_reactant_types()

print("=" * 80)
print("TERMINAL vs INTERNAL ALKENE/ALKYNE TEST")
print("=" * 80)

# Test cases: (SMILES, expected_member_type, description)
test_cases = [
    # Terminal alkenes
    ('C=C', 'terminal-alkene', 'ethene'),
    ('C=CC', 'terminal-alkene', 'propene'),
    ('C=CCCC', 'terminal-alkene', '1-pentene'),
    ('c1ccccc1C=C', 'terminal-alkene', 'styrene'),
    ('C=CC(C)C', 'terminal-alkene', '3-methylbutene'),
    
    # Internal alkenes
    ('CC=CC', 'internal-alkene', '2-butene'),
    ('CCC=CCC', 'internal-alkene', '3-hexene'),
    ('c1ccccc1C=Cc1ccccc1', 'internal-alkene', 'stilbene'),
    
    # Alkyl-substituted (should still match)
    ('CCCC=CCCC', 'R-alkene', 'dialkyl internal alkene'),
    
    # Aryl-substituted
    ('c1ccccc1C=CC', 'Ar-alkene', 'trans-β-methylstyrene'),
    
    # Terminal alkynes
    ('C#C', 'terminal-alkyne', 'acetylene'),
    ('CC#C', 'terminal-alkyne', 'propyne'),
    ('CCC#C', 'terminal-alkyne', '1-butyne'),
    ('c1ccccc1C#C', 'terminal-alkyne', 'phenylacetylene'),
    
    # Internal alkynes
    ('CC#CC', 'internal-alkyne', '2-butyne'),
    ('CCC#CCC', 'internal-alkyne', '3-hexyne'),
    ('c1ccccc1C#Cc1ccccc1', 'internal-alkyne', 'diphenylacetylene'),
]

passed = 0
failed = 0

for smiles, expected, description in test_cases:
    result = classify_reactant(smiles, rt)
    
    if result:
        actual = result['member_type']
        status = '✓' if actual == expected else '✗'
        
        if actual == expected:
            passed += 1
            print(f"{status} {smiles:25s} → {actual:25s} ({description})")
        else:
            failed += 1
            print(f"{status} {smiles:25s} → {actual:25s} (expected: {expected}, desc: {description})")
            # Show all matches for failed cases
            all_matches = get_all_matches(smiles, rt)
            if len(all_matches) > 1:
                print(f"   All matches: {', '.join([m['member_type'] for m in all_matches[:3]])}")
    else:
        failed += 1
        print(f"✗ {smiles:25s} → NOT MATCHED (expected: {expected})")

print("\n" + "=" * 80)
print(f"Results: {passed} passed, {failed} failed out of {len(test_cases)} tests")
print("=" * 80)

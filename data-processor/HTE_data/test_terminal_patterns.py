"""
Test terminal alkene and alkyne patterns
"""

from classify_reactant import classify_reactant, load_reactant_types

rt = load_reactant_types()

test_cases = [
    ('C=C', 'ethene (terminal)'),
    ('C=CC', 'propene (terminal)'),
    ('CC=CC', 'butene (internal)'),
    ('c1ccccc1C=C', 'styrene (terminal)'),
    ('C#C', 'acetylene (terminal)'),
    ('CC#C', 'propyne (terminal)'),
    ('CC#CC', 'butyne (internal)'),
    ('c1ccccc1C#C', 'phenylacetylene (terminal)'),
]

print('Testing Alkene and Alkyne patterns:')
print('=' * 70)

for smiles, description in test_cases:
    result = classify_reactant(smiles, rt)
    if result:
        print(f'{smiles:20s} → {result["member_type"]:25s} ({result["name"]})')
    else:
        print(f'{smiles:20s} → NOT MATCHED')

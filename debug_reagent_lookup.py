"""
Debug script to understand why reagent filtering is too aggressive.
"""

from chemtools import reagent_lookup

# Test some common reagents from the precedents
test_reagents = [
    ('Pd/Tri-tert-butylphosphinetetrafluoroborate', 'metal_precursor'),
    ('Pd/diphenylphosphino', 'metal_precursor'),
    ('Pd/RuPhos', 'metal_precursor'),
    ('534-17-8', 'base'),  # Carbonate base
    ('1310-58-3', 'base'),  # Sodium hydroxide
    ('7778-53-2', 'base'),  # Potassium phosphate
    ('108-88-3', 'solvent'),  # Toluene
    ('68-12-2', 'solvent'),  # DMF
    ('phen', 'ligand'),  # Phenanthroline
    ('Cu', 'metal_precursor'),  # Copper
]

print("Testing reagent lookup for common reagents:")
print("=" * 80)

for reagent, role in test_reagents:
    result = reagent_lookup.enrich_reagent_info(reagent, role)
    found = result.get('found', False)
    name = result.get('name', 'N/A')
    print(f"{reagent:50s} ({role:15s}) -> Found: {str(found):5s} | Name: {name}")

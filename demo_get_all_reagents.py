"""
Quick demo: How to get all reagents by type
"""

from chemtools.reagent_lookup import (
    get_all_reagents_by_type,
    count_reagents_by_type,
    get_all_reagent_types
)

print("="*70)
print("  Quick Demo: Get All Reagents by Type")
print("="*70)

# Show available types
print("\n1. Available Reagent Types:")
types = get_all_reagent_types()
for t in types:
    count = count_reagents_by_type(t)
    print(f"   {t:30s} ({count:3d} reagents)")

# Get all solvents
print("\n2. Get All Solvents:")
solvents = get_all_reagents_by_type('solvent')
print(f"   Total: {len(solvents)} solvents")
print(f"   Examples: {', '.join([s['name'] for s in solvents[:5]])}")

# Get all ligands
print("\n3. Get All Ligands:")
ligands = get_all_reagents_by_type('ligand')
print(f"   Total: {len(ligands)} ligands")
abbr_examples = []
for lig in ligands[:5]:
    abbr = lig.get('abbreviation', [])
    if abbr:
        abbr_examples.append(abbr[0])
print(f"   Examples: {', '.join(abbr_examples)}")

# Filter example
print("\n4. Filter Example - Find All Phosphine Ligands:")
phosphines = [
    lig for lig in ligands 
    if 'PPh' in lig.get('name', '') or 
       'phosphine' in lig.get('name', '').lower()
]
print(f"   Found {len(phosphines)} phosphine ligands")
for p in phosphines[:8]:
    abbr = p.get('abbreviation', [])
    abbr_str = f" ({abbr[0]})" if abbr else ""
    print(f"   • {p['name']}{abbr_str}")

print("\n" + "="*70)
print("  ✅ Demo Complete!")
print("="*70)
print("\nUsage:")
print("  from chemtools.reagent_lookup import get_all_reagents_by_type")
print("  solvents = get_all_reagents_by_type('solvent')")
print("  print(f'Total: {len(solvents)} solvents')")

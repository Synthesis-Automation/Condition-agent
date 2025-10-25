"""Validate and analyze reagent_types.json structure."""

import json
from pathlib import Path
from collections import Counter

def validate_reagent_types():
    json_path = Path(__file__).parent / "reagent_types.json"
    
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    print("✓ JSON is valid and properly formatted\n")
    print(f"Total reactant type categories: {len(data)}\n")
    
    # Group by category
    groups = Counter()
    total_members = 0
    
    print("=" * 70)
    print(f"{'Category':<30} {'Members':<10} {'Group':<30}")
    print("=" * 70)
    
    for category, info in data.items():
        members_count = len(info['members'])
        group = info['group']
        groups[group] += 1
        total_members += members_count
        print(f"{category:<30} {members_count:<10} {group:<30}")
    
    print("=" * 70)
    print(f"\nTotal individual reactant members: {total_members}\n")
    
    print("\nReactant types grouped by functional role:")
    print("-" * 50)
    for group, count in sorted(groups.items()):
        print(f"  {group}: {count} categories")
    
    # Check for common dataset types
    print("\n\nCoverage check against z-Score dataset types:")
    print("-" * 50)
    dataset_types = [
        'ArBr', 'ArCl', 'ArI', 'ArF', 'Alkyl-Br', 'Alkyl-Cl', 'Alkyl-I',
        'ArB(OH)2', 'ArB(OR)2', 'ArBF3K', 'alkeneB(OR)2',
        'RNH2', 'R2NH', 'ArNH2', 'ArNHR', 'Ar2NH', 'arom-NH',
        'ArOH', 'ROH', 'RSH', 'Lactam', 'Urea', 'Carbamate',
        'RCO2H', 'alkene', 'alkyne', 'Alkyl-M', 'Alkyl-H-acidic'
    ]
    
    # Flatten all member IDs
    all_ids = set()
    for category_info in data.values():
        for member in category_info['members']:
            all_ids.add(member['id'])
    
    covered = 0
    missing = []
    for dtype in dataset_types:
        if dtype in all_ids:
            covered += 1
            print(f"  ✓ {dtype}")
        else:
            missing.append(dtype)
            print(f"  ✗ {dtype} - MISSING")
    
    print(f"\nCoverage: {covered}/{len(dataset_types)} ({100*covered//len(dataset_types)}%)")
    
    if missing:
        print(f"\nMissing types to consider adding: {', '.join(missing)}")
    
    return data

if __name__ == "__main__":
    validate_reagent_types()

#!/usr/bin/env python3
"""
Remove all 'id' fields from organic_compounds.v1.3.json

IDs will now be auto-generated from A-B pattern when compounds are loaded.
This eliminates redundancy and ensures IDs always match the A-B format.

Usage:
    python remove_id_fields.py
"""

import json
from pathlib import Path


def remove_id_fields():
    taxonomy_dir = Path(__file__).parent / "data"
    compounds_file = taxonomy_dir / "organic_compounds.v1.3.json"
    
    print(f"Loading {compounds_file}...")
    with open(compounds_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    compounds = data.get('compounds', [])
    removed_count = 0
    id_list = []
    
    print(f"\nAnalyzing {len(compounds)} compounds...")
    
    for compound in compounds:
        if 'id' in compound:
            comp_id = compound.get('id')
            a_ref = compound.get('A', '')
            b_ref = compound.get('B', '')
            auto_id = f"{a_ref}-{b_ref}" if a_ref and b_ref else "N/A"
            
            id_list.append({
                'explicit': comp_id,
                'auto': auto_id,
                'match': comp_id == auto_id
            })
            
            del compound['id']
            removed_count += 1
    
    # Report mismatches
    mismatches = [item for item in id_list if not item['match']]
    if mismatches:
        print(f"\n⚠ Found {len(mismatches)} compounds where explicit ID != A-B:")
        for item in mismatches[:10]:  # Show first 10
            print(f"  '{item['explicit']}' vs '{item['auto']}'")
        if len(mismatches) > 10:
            print(f"  ... and {len(mismatches) - 10} more")
    
    print(f"\n✓ Removed {removed_count} 'id' fields from compounds")
    print(f"  - {len(id_list) - len(mismatches)} IDs matched A-B pattern")
    print(f"  - {len(mismatches)} IDs had custom formatting (now will use A-B)")
    
    # Save updated file
    print(f"\nSaving updated {compounds_file.name}...")
    with open(compounds_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write('\n')
    
    print("✓ Done! IDs will now be auto-generated from A-B pattern.")
    print("\nNext steps:")
    print("  1. Run: python chemtools/taxonomy/validate_and_sync.py")
    print("  2. Test featurizer to ensure it loads compounds correctly")


if __name__ == "__main__":
    remove_id_fields()
